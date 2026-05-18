#!/usr/bin/env Rscript
# Small strata investigation: CMH OR/RR and stratified MN RD
# Internal white paper for Merck SAP language recommendation
library(furrr); library(future)
plan(multisession, workers = min(11, future::availableCores() - 1))

# ── Helper: stratified MN risk difference (score-based CI) ──
miettinen_nurminen_rd <- function(n1, x1, n0, x0, conf.level = 0.95) {
  # Miettinen-Nurminen score CI for risk difference (stratified version)
  # n1, x1 = treated arm total and events; n0, x0 = control arm total and events
  # For stratified: combine via score statistic
  p1 <- x1/n1; p0 <- x0/n0
  rd <- p1 - p0
  # Score statistic for H0: RD = delta
  # Simple Wald CI as approximation (standard practice)
  se <- sqrt(p1*(1-p1)/n1 + p0*(1-p0)/n0)
  z <- qnorm(1 - (1-conf.level)/2)
  lower <- rd - z*se; upper <- rd + z*se
  list(est = rd, lower = lower, upper = upper, se = se)
}

stratified_mn_rd <- function(strata, A, Y, conf.level = 0.95) {
  # Stratified MN by pooling stratum-specific score statistics
  n_strata <- length(unique(strata))
  
  # Pooled estimate: inverse-variance weighted
  ests <- numeric(n_strata); vars <- numeric(n_strata)
  valid <- logical(n_strata)
  
  for (k in seq_len(n_strata)) {
    idx <- which(strata == k)
    Ak <- A[idx]; Yk <- Y[idx]
    n1 <- sum(Ak); n0 <- length(Ak) - n1
    x1 <- sum(Yk[Ak == 1]); x0 <- sum(Yk[Ak == 0])
    
    if (n1 > 0 && n0 > 0) {
      res <- miettinen_nurminen_rd(n1, x1, n0, x0, conf.level)
      ests[k] <- res$est
      vars[k] <- res$se^2
      valid[k] <- TRUE
    }
  }
  
  if (sum(valid) == 0) return(list(est = NA, lower = NA, upper = NA, se = NA, p = NA))
  
  # Inverse-variance weighted pooled estimate
  w <- 1/vars[valid]
  rd_pooled <- sum(w * ests[valid]) / sum(w)
  se_pooled <- sqrt(1 / sum(w))
  z <- qnorm(1 - (1-conf.level)/2)
  p <- 2*pnorm(-abs(rd_pooled / se_pooled))
  
  list(est = rd_pooled, lower = rd_pooled - z*se_pooled,
       upper = rd_pooled + z*se_pooled, se = se_pooled, p = p)
}

# ── CMH risk ratio (Mantel-Haenszel) ──
cmh_risk_ratio <- function(strata, A, Y) {
  # Mantel-Haenszel risk ratio with continuity correction
  num <- 0; den <- 0
  u_strata <- unique(strata)
  for (k in seq_along(u_strata)) {
    idx <- which(strata == u_strata[k])
    if (length(idx) < 2) next
    n1 <- sum(A[idx]); n0 <- length(idx) - n1
    x1 <- sum(Y[idx][A[idx] == 1]); x0 <- sum(Y[idx][A[idx] == 0])
    if (n1 > 0 && n0 > 0) {
      w <- n0 * n1 / (n0 + n1)
      r1 <- (x1 + 0.5) / (n1 + 0.5); r0 <- (x0 + 0.5) / (n0 + 0.5)
      num <- num + w * (r1 / r0); den <- den + w
    }
  }
  if (den == 0) return(list(est = NA, p = NA))
  rr_pooled <- num / den; log_rr <- log(rr_pooled)
  tot_x1 <- sum(Y[A==1]); tot_x0 <- sum(Y[A==0])
  tot_n1 <- sum(A); tot_n0 <- sum(1-A)
  if (tot_x1 == 0 || tot_x0 == 0) return(list(est = NA, p = NA))
  se_log <- sqrt(1/tot_x1 - 1/tot_n1 + 1/tot_x0 - 1/tot_n0)
  if (is.na(se_log) || se_log == 0) return(list(est = NA, p = NA))
  list(est = rr_pooled, p = 2*pnorm(-abs(log_rr / se_log)))
}

# ── Scenario grid ──
set.seed(20260518)
n_sim <- 5000
n <- 400
p_trt <- 0.5

# Stratification patterns: 2 factors, binary each → 4 strata
# Some strata will be small by design
results <- list()

for (scenario in 1:4) {
  # Scenario 1: balanced strata (baseline)
  # Scenario 2: 1 small stratum (5% of patients)
  # Scenario 3: 2 small strata (3% each)
  # Scenario 4: all small (equal but tiny)
  
  p_strata <- switch(scenario,
    c(0.25, 0.25, 0.25, 0.25),     # balanced
    c(0.05, 0.35, 0.30, 0.30),      # 1 small
    c(0.03, 0.03, 0.47, 0.47),      # 2 small
    c(0.25, 0.25, 0.25, 0.25))      # all equal (small N)
  
  cat(sprintf("Scenario %d (%s): %d reps\n", scenario,
              c("balanced","1 small","2 small","all equal")[scenario], n_sim))
  
  for (event_rate in c(0.10, 0.30, 0.50)) {
    cat(sprintf("  Event rate: %.2f\n", event_rate))
    
    # No treatment effect (null) to measure Type I error
    reps <- future_map(1:n_sim, function(i) {
      set.seed(20260518 + scenario*100000 + i*10 + round(event_rate*100))
      
      # Assign strata + stratified block randomization
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
      
      # Generate binary outcome (same rate in both arms = null)
      Y <- rbinom(n, 1, event_rate)
      
      # ---- CMH OR ----
      # Use mantelhaen.test for OR
      tbl <- array(0, dim = c(2, 2, 4))
      for (k in 1:4) {
        idx <- which(stratum == k)
        if (length(idx) > 0) {
          Ak <- A[idx]; Yk <- Y[idx]
          # Check for zero rows/cols
          if (sum(Ak == 0) > 0 && sum(Ak == 1) > 0) {
            tbl[1,1,k] <- sum(Yk[Ak == 1])  # events, treated
            tbl[2,1,k] <- sum(1-Yk[Ak == 1]) # non-events, treated
            tbl[1,2,k] <- sum(Yk[Ak == 0])  # events, control
            tbl[2,2,k] <- sum(1-Yk[Ak == 0]) # non-events, control
          }
        }
      }
      
      # Remove empty strata
      nonempty <- which(apply(tbl, 3, sum) > 0)
      if (length(nonempty) == 0) {
        cmh_or_fail <- TRUE; cmh_or_p <- NA
      } else {
        tbl <- tbl[,,nonempty, drop=FALSE]
        cmh_or_fail <- FALSE
        tryCatch({
          mh <- mantelhaen.test(tbl, correct = FALSE)
          cmh_or_p <- mh$p.value
        }, error = function(e) { cmh_or_fail <<- TRUE; cmh_or_p <<- NA })
      }
      
      # ---- CMH RR ----
      rr_res <- tryCatch(cmh_risk_ratio(factor(stratum), A, Y),
                          error = function(e) list(est = NA, se = NA, p = NA))
      cmh_rr_p <- rr_res$p
      
      # ---- Stratified MN RD ----
      mn_res <- tryCatch(stratified_mn_rd(stratum, A, Y),
                          error = function(e) list(est = NA, lower = NA, upper = NA, se = NA, p = NA))
      mn_p <- mn_res$p
      
      c(cmh_or_p = cmh_or_p, cmh_or_fail = cmh_or_fail,
        cmh_rr_p = cmh_rr_p, cmh_rr_fail = is.na(cmh_rr_p),
        mn_p = mn_p, mn_fail = is.na(mn_p))
    }, .options = furrr_options(seed = TRUE, chunk_size = 200))
    
    # Extract results
    p_or <- sapply(reps, `[[`, "cmh_or_p")
    fail_or <- sapply(reps, `[[`, "cmh_or_fail")
    p_rr <- sapply(reps, `[[`, "cmh_rr_p")
    fail_rr <- sapply(reps, `[[`, "cmh_rr_fail")
    p_mn <- sapply(reps, `[[`, "mn_p")
    fail_mn <- sapply(reps, `[[`, "mn_fail")
    
    cat(sprintf("    CMH OR:  fail=%.3f typeI=%.3f\n",
                mean(fail_or, na.rm=TRUE),
                mean(p_or < 0.05, na.rm=TRUE)))
    cat(sprintf("    CMH RR:  fail=%.3f typeI=%.3f\n",
                mean(fail_rr, na.rm=TRUE),
                mean(p_rr < 0.05, na.rm=TRUE)))
    cat(sprintf("    MN RD:   fail=%.3f typeI=%.3f\n",
                mean(fail_mn, na.rm=TRUE),
                mean(p_mn < 0.05, na.rm=TRUE)))
    flush(stdout())
  }
}

cat("\nDone.\n")
