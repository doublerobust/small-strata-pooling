#!/usr/bin/env Rscript
# Small strata power analysis: Does pooling improve power?
# Compares no-pooling vs pooling-small-strata for CMH OR/RR and MN RD
library(furrr); library(future)
plan(multisession, workers = min(11, future::availableCores() - 1))
library(PropCIs)

# ── Helper: stratified MN RD ──
stratified_mn_rd <- function(strata, A, Y, conf.level = 0.95) {
  u_strata <- unique(strata)
  ests <- numeric(length(u_strata)); vars <- numeric(length(u_strata))
  valid <- logical(length(u_strata))
  for (k in seq_along(u_strata)) {
    idx <- which(strata == u_strata[k])
    Ak <- A[idx]; Yk <- Y[idx]
    n1 <- sum(Ak); n0 <- length(Ak) - n1
    x1 <- sum(Yk[Ak == 1]); x0 <- sum(Yk[Ak == 0])
    if (n1 > 0 && n0 > 0 && x1 + x0 > 0) {
      ci <- tryCatch(diffscoreci(x1, n1, x0, n0, conf.level), error = function(e) NULL)
      if (!is.null(ci)) {
        rd <- x1/n1 - x0/n0
        se <- (ci$conf.int[2] - ci$conf.int[1]) / (2*qnorm(1-(1-conf.level)/2))
        if (is.finite(se) && se > 0) { ests[k] <- rd; vars[k] <- se^2; valid[k] <- TRUE }
      }
    }
  }
  if (sum(valid) < 2) return(list(p = NA))
  w <- 1/vars[valid]
  rd_pooled <- sum(w*ests[valid])/sum(w)
  se_pooled <- sqrt(1/sum(w))
  if (is.na(se_pooled) || se_pooled == 0) return(list(p = NA))
  list(p = 2*pnorm(-abs(rd_pooled/se_pooled)))
}

# ── Helper: CMH risk ratio ──
cmh_risk_ratio <- function(strata, A, Y) {
  num <- 0; den <- 0; var_num <- 0; var_den <- 0
  u_strata <- unique(strata)
  for (k in seq_along(u_strata)) {
    idx <- which(strata == u_strata[k])
    if (length(idx) < 2) next
    n1 <- sum(A[idx]); n0 <- length(idx) - n1
    x1 <- sum(Y[idx][A[idx] == 1]); x0 <- sum(Y[idx][A[idx] == 0])
    if (n1 > 0 && n0 > 0) {
      w <- n0 * n1 / (n0 + n1)
      r1 <- (x1 + 0.5) / (n1 + 0.5); r0 <- (x0 + 0.5) / (n0 + 0.5)
      rr_k <- r1 / r0
      num <- num + w * rr_k; den <- den + w
      if (x1 + x0 > 0) {
        var_num <- var_num + w^2 * (1/(x1+0.5) - 1/(n1+0.5) + 1/(x0+0.5) - 1/(n0+0.5))
      }
    }
  }
  if (den == 0) return(list(est = NA, p = NA))
  rr_pooled <- num / den
  var_log <- var_num / (den^2)
  if (is.na(var_log) || var_log <= 0) return(list(est = NA, p = NA))
  se_log <- sqrt(var_log)
  if (is.na(se_log) || se_log == 0) return(list(est = NA, p = NA))
  list(est = rr_pooled, p = 2*pnorm(-abs(log(rr_pooled) / se_log)))
}

# ── Helper: stratified CMH OR ──
cmh_odds_ratio <- function(strata, A, Y) {
  tbl <- array(0, dim = c(2, 2, 4))
  for (k in 1:4) {
    idx <- which(strata == k)
    if (length(idx) > 0) {
      Ak <- A[idx]; Yk <- Y[idx]
      if (sum(Ak == 0) > 0 && sum(Ak == 1) > 0) {
        tbl[1,1,k] <- sum(Yk[Ak == 1])
        tbl[2,1,k] <- sum(1-Yk[Ak == 1])
        tbl[1,2,k] <- sum(Yk[Ak == 0])
        tbl[2,2,k] <- sum(1-Yk[Ak == 0])
      }
    }
  }
  nonempty <- which(apply(tbl, 3, sum) > 0)
  if (length(nonempty) == 0) return(list(p = NA, fail = TRUE))
  tbl <- tbl[,,nonempty, drop=FALSE]
  p <- NA; fail <- FALSE
  tryCatch({
    mh <- mantelhaen.test(tbl, correct = FALSE)
    p <- mh$p.value
  }, error = function(e) { fail <<- TRUE })
  list(p = p, fail = fail)
}

# ── Single replication ──
run_rep <- function(i, scenario, ev_rate, log_or, p_strata) {
  set.seed(20260518 + scenario*100000 + i*10 + round(ev_rate*100))
  
  n <- 400
  event_rate <- ev_rate
  stratum <- sample(1:4, n, replace = TRUE, prob = p_strata)
  A <- integer(n)
  for (s in unique(stratum)) {
    idx <- which(stratum == s)
    n_s <- length(idx)
    for (b in seq_len(ceiling(n_s / 4))) {
      start <- (b-1)*4 + 1
      end <- min(b*4, n_s)
      if (start <= end) {
        bn <- end - start + 1
        n_trt <- floor(bn / 2)
        A[idx[start:end]] <- sample(c(rep(1, n_trt), rep(0, bn - n_trt)))
      }
    }
  }
  
  # Generate outcome with treatment effect
  # p_control = event_rate; p_treated from log(OR)
  or <- exp(log_or)
  logit_p0 <- log(event_rate / (1 - event_rate))
  p1 <- plogis(logit_p0 + log_or)
  p0 <- event_rate
  Y <- ifelse(A == 1, rbinom(n, 1, p1), rbinom(n, 1, p0))
  
  # ── No pooling (all 4 strata as-is) ──
  cmh_or <- cmh_odds_ratio(stratum, A, Y)
  cmh_rr <- cmh_risk_ratio(stratum, A, Y)
  mn_rd  <- stratified_mn_rd(stratum, A, Y)
  
  # ── Pool small strata ──
  # Merge strata < 10 patients into nearest larger stratum
  tab <- table(stratum)
  small <- as.numeric(names(tab[tab < 10]))
  stratum_pooled <- stratum
  if (length(small) > 0 && length(small) < 4) {
    large <- as.numeric(names(tab[tab >= 10]))
    for (s in small) {
      if (length(large) > 0) {
        # Merge into nearest large stratum by size
        target <- large[which.min(abs(tab[as.character(large)] - tab[as.character(s)]))]
        stratum_pooled[stratum == s] <- target
      }
    }
  }
  
  cmh_or_p <- cmh_odds_ratio(stratum_pooled, A, Y)
  cmh_rr_p <- cmh_risk_ratio(stratum_pooled, A, Y)
  mn_rd_p  <- stratified_mn_rd(stratum_pooled, A, Y)
  
  c(
    or_p = cmh_or$p, or_fail = cmh_or$fail,
    rr_p = cmh_rr$p, rr_fail = is.na(cmh_rr$p),
    mn_p = mn_rd$p,  mn_fail = is.na(mn_rd$p),
    or_p_pool = cmh_or_p$p, or_fail_pool = cmh_or_p$fail,
    rr_p_pool = cmh_rr_p$p, rr_fail_pool = is.na(cmh_rr_p$p),
    mn_p_pool = mn_rd_p$p,  mn_fail_pool = is.na(mn_rd_p$p)
  )
}

# ── Main simulation ──
set.seed(20260518)
n_sim <- 5000

scenarios <- list(
  `1` = list(label = "Balanced",       p = c(0.25, 0.25, 0.25, 0.25)),
  `2` = list(label = "1 small (5%)",   p = c(0.05, 0.35, 0.30, 0.30)),
  `3` = list(label = "2 small (3%)",   p = c(0.03, 0.03, 0.47, 0.47)),
  `4` = list(label = "2 tiny (1%,2%)", p = c(0.01, 0.02, 0.485, 0.485))
)

event_rates <- c(0.10, 0.30, 0.50)
log_ors <- c(0.5, 0.7)  # ~OR 1.65, ~OR 2.0

cat("=== Small Strata Power Analysis ===\n")
cat("n_sim =", n_sim, "| N = 400 | Block = 4\n\n")

for (scenario in seq_along(scenarios)) {
  sc_label <- scenarios[[scenario]]$label
  p_strata <- scenarios[[scenario]]$p
  cat(sprintf("── Scenario %d: %s ──\n", scenario, sc_label))
  
  for (ev_rate in event_rates) {
    for (log_or in log_ors) {
      cat(sprintf("  Event rate: %.2f | OR: %.2f\n", ev_rate, exp(log_or)))
      
      reps <- future_map(1:n_sim, function(i) {
        run_rep(i, scenario, ev_rate, log_or, p_strata)
      }, .options = furrr_options(seed = TRUE, chunk_size = 200))
      
      extract <- function(name) sapply(reps, function(r) as.numeric(r[name]))
      
      p_or   <- extract("or_p"); f_or   <- extract("or_fail")
      p_rr   <- extract("rr_p"); f_rr   <- extract("rr_fail")
      p_mn   <- extract("mn_p"); f_mn   <- extract("mn_fail")
      
      p_or_p <- extract("or_p_pool"); f_or_p <- extract("or_fail_pool")
      p_rr_p <- extract("rr_p_pool"); f_rr_p <- extract("rr_fail_pool")
      p_mn_p <- extract("mn_p_pool"); f_mn_p <- extract("mn_fail_pool")
      
      # Power at alpha = 0.05 (one-sided for benefit)
      pow_or  <- mean(p_or/2 < 0.025, na.rm = TRUE)
      pow_rr  <- mean(p_rr/2 < 0.025, na.rm = TRUE)
      pow_mn  <- mean(p_mn/2 < 0.025, na.rm = TRUE)
      
      pow_or_p <- mean(p_or_p/2 < 0.025, na.rm = TRUE)
      pow_rr_p <- mean(p_rr_p/2 < 0.025, na.rm = TRUE)
      pow_mn_p <- mean(p_mn_p/2 < 0.025, na.rm = TRUE)
      
      cat(sprintf("    CMH OR  fail=%.3f | no pool=%.3f | pool=%.3f | diff=%.3f\n",
                  mean(f_or, na.rm=TRUE), pow_or, pow_or_p, pow_or_p - pow_or))
      cat(sprintf("    CMH RR  fail=%.3f | no pool=%.3f | pool=%.3f | diff=%.3f\n",
                  mean(f_rr, na.rm=TRUE), pow_rr, pow_rr_p, pow_rr_p - pow_rr))
      cat(sprintf("    MN RD   fail=%.3f | no pool=%.3f | pool=%.3f | diff=%.3f\n",
                  mean(f_mn, na.rm=TRUE), pow_mn, pow_mn_p, pow_mn_p - pow_mn))
      flush(stdout())
    }
  }
  cat("\n")
}

cat("Done.\n")
