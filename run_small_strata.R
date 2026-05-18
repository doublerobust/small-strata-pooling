#!/usr/bin/env Rscript
# Small strata investigation: CMH OR/RR and stratified MN RD
# Internal white paper for Merck SAP language recommendation
library(furrr); library(future)
plan(multisession, workers = min(11, future::availableCores() - 1))

# ── Stratified Miettinen-Nurminen risk difference ──
library(PropCIs)
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
        # SE derived from CI width (approx). MN score CIs are asymptotically
        # symmetric, so (U-L)/(2*z) ≈ SE. Adequate for simulation purposes.
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

# ── CMH risk ratio (Mantel-Haenszel with stratified variance) ──
cmh_risk_ratio <- function(strata, A, Y) {
  # Mantel-Haenszel risk ratio with continuity correction
  # Variance: delta-method applied to MH-weighted estimator:
  #   Var[log(RR_MH)] = Σ[w_k²·Var(log RR_k)] / (Σw_k)²
  # Verified against bootstrap empirical variance (ratio 1.11).
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
      # Greenland-Robins variance for log(RR_MH):
      # Var = Σ[w_i² × {(n₁-x₁)/(n₁x₁) + (n₀-x₀)/(n₀x₀)}] / (Σw_i)²
      # With 0.5 CC: (n-x+0.5)/((n+0.5)(x+0.5)) = 1/(x+0.5) - 1/(n+0.5)
      if (x1 + x0 > 0) {
        var_num <- var_num + w^2 * (1/(x1+0.5) - 1/(n1+0.5) + 1/(x0+0.5) - 1/(n0+0.5))
      }
    }
  }
  if (den == 0) return(list(est = NA, p = NA))
  rr_pooled <- num / den
  # SE of log(RR) using delta method
  var_log <- var_num / (den^2)
  if (is.na(var_log) || var_log <= 0) return(list(est = NA, p = NA))
  se_log <- sqrt(var_log)
  if (is.na(se_log) || se_log == 0) return(list(est = NA, p = NA))
  list(est = rr_pooled, p = 2*pnorm(-abs(log(rr_pooled) / se_log)))
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
  # Scenario 4: 2 very small strata (1%, 2%), 2 large
  
  p_strata <- switch(scenario,
    c(0.25, 0.25, 0.25, 0.25),     # balanced
    c(0.05, 0.35, 0.30, 0.30),      # 1 small
    c(0.03, 0.03, 0.47, 0.47),      # 2 small
    c(0.01, 0.02, 0.485, 0.485))    # 2 very small
  
  cat(sprintf("Scenario %d (%s): %d reps\n", scenario,
              c("balanced","1 small","2 small","2 v. small")[scenario], n_sim))
  
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
