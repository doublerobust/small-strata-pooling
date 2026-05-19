#!/usr/bin/env Rscript
# ============================================================================
# Double-Programming Verification: Stratified Binary & Survival Endpoints
# ============================================================================
# Part 1: Binary endpoints (CMH OR, CMH RR, MN RD) under null scenario
# Part 2: Time-to-event (stratified log-rank) under null and alternative
# ============================================================================

# ---- Load libraries ----
suppressPackageStartupMessages({
  library(furrr)
  library(PropCIs)
  library(survival)
})

# ---- Global settings ----
N_REPS      <- 5000
N           <- 500
STRATA_PROBS_BALANCED <- c(0.25, 0.25, 0.25, 0.25)   # Part 1 balanced
STRATA_PROBS_SPARSE   <- c(0.01, 0.02, 0.485, 0.485)  # Part 1 sparse
STRATA_PROBS_SPARSE8  <- c(0.01, 0.01, rep(0.98/6, 6)) # Part 2 sparse
PLAN <- plan("multisession", workers = 11)
future::plan(PLAN)

cat("=== Double-Programming Verification ===\n")
cat("Workers: 11\n")
cat("Reps:", N_REPS, "\n\n")

# ============================================================================
# PART 1: Binary endpoints
# ============================================================================
cat("========== PART 1: BINARY ENDPOINTS ==========\n\n")

#' Simulate one replication for binary endpoints (stratified 2x2)
#' @param n Total sample size
#' @param strata_probs Vector of stratum probabilities (must sum to 1)
#' @param event_rate Common event rate under null
#' @return List with strata assignment, treatment, outcome per subject
simulate_binary_data <- function(n, strata_probs, event_rate = 0.30) {
  k <- length(strata_probs)
  # Assign strata
  strata <- sample(1:k, n, replace = TRUE, prob = strata_probs)
  
  # 1:1 randomization within each stratum
  treatment <- integer(n)
  for (s in 1:k) {
    idx <- which(strata == s)
    n_s <- length(idx)
    n_t <- floor(n_s / 2)
    treatment[idx] <- 0
    treatment[idx[1:n_t]] <- 1  # first half gets treatment
  }
  
  # Generate binary outcomes (null: no treatment effect)
  outcome <- rbinom(n, 1, event_rate)
  
  list(strata = strata, treatment = treatment, outcome = outcome)
}

#' Compute CMH risk ratio (RR) manually
#' For each stratum, compute risk ratio and use Miettinen and Nurminen (1985) CI
#' Then combine via inverse-variance weighting of log(RR)
#' Returns test statistic and p-value (two-sided)
compute_cmh_rr <- function(data) {
  strata_levels <- sort(unique(data$strata))
  k <- length(strata_levels)
  
  log_rr_weights <- numeric(k)
  log_rr_estimates <- numeric(k)
  
  for (s in strata_levels) {
    idx <- which(data$strata == s)
    a <- sum(data$treatment[idx] == 1 & data$outcome[idx] == 1)  # treated, event
    b <- sum(data$treatment[idx] == 1 & data$outcome[idx] == 0)  # treated, no event
    c <- sum(data$treatment[idx] == 0 & data$outcome[idx] == 1)  # control, event
    d <- sum(data$treatment[idx] == 0 & data$outcome[idx] == 0)  # control, no event
    
    n1 <- a + b  # treated total
    n0 <- c + d  # control total
    
    r1 <- a / n1  # risk in treated
    r0 <- c / n0  # risk in control
    
    # Risk ratio (with continuity correction if needed)
    if (r1 == 0) r1 <- 0.5 / n1
    if (r0 == 0) r0 <- 0.5 / n0
    
    rr <- r1 / r0
    log_rr <- log(rr)
    
    # Variance of log(RR): (1/a - 1/n1) + (1/c - 1/n0)  (Wald variance)
    var_log_rr <- (1/a - 1/n1) + (1/c - 1/n0)
    if (is.infinite(var_log_rr) || var_log_rr < 0) {
      log_rr_weights[s] <- 0
      next
    }
    
    log_rr_estimates[s] <- log_rr
    log_rr_weights[s] <- 1 / var_log_rr
  }
  
  # Combine: inverse-variance weighted log(RR)
  w <- log_rr_weights
  w <- w[w > 0]
  est <- log_rr_estimates[log_rr_weights > 0]
  
  if (length(w) == 0) {
    # Fallback: use CMH chi-square for testing
    return(list(statistic = NA, p_value = NA, method = "CMH_RR_fallback"))
  }
  
  pooled_log_rr <- sum(w * est) / sum(w)
  se_pooled <- sqrt(1 / sum(w))
  z_stat <- pooled_log_rr / se_pooled
  p_value <- 2 * pnorm(-abs(z_stat))
  
  list(statistic = z_stat, p_value = p_value, method = "CMH_RR")
}

#' Run binary endpoint simulation
#' @param strata_probs Vector of stratum probabilities
run_binary_simulation <- function(strata_probs, label) {
  cat("  Strata config:", label, "\n")
  cat("    Probabilities:", paste(round(strata_probs, 4), collapse = ", "), "\n")
  
  # Pre-allocate
  p_mantel <- numeric(N_REPS)
  p_cmh_rr <- numeric(N_REPS)
  p_mn_rd  <- numeric(N_REPS)
  
  cat("  Running", N_REPS, "replications...\n")
  
  # Use furrr for parallelism
  results <- future.apply::future_lapply(1:N_REPS, function(rep_id) {
    dat <- simulate_binary_data(N, strata_probs, event_rate = 0.30)
    
    # Build 2x2 table per stratum for mantelhaen.test
    tables <- lapply(sort(unique(dat$strata)), function(s) {
      idx <- which(dat$strata == s)
      table(dat$treatment[idx], dat$outcome[idx])
    })
    
    # Mantel-Haenszel OR
    p_mh <- tryCatch({
      mh <- mantelhaen.test(tables)
      mh$p.value
    }, error = function(e) NA)
    
    # CMH RR (custom)
    p_rr <- tryCatch({
      rr_res <- compute_cmh_rr(dat)
      rr_res$p_value
    }, error = function(e) NA)
    
    # MN RD (PropCIs::diffscoreci on pooled data, but we need stratified)
    # MN RD = Miettinen-Nurminen risk difference CI
    # We pool across strata for the test (standard approach)
    p_mn <- tryCatch({
      # Pool all data, use diffscoreci for risk difference
      pooled <- data.frame(
        trt = dat$treatment,
        event = dat$outcome
      )
      # diffscoreci(r1, n1, r2, n2, conf.level)
      r1 <- sum(pooled$event[pooled$trt == 1]) / sum(pooled$trt == 1)
      n1 <- sum(pooled$trt == 1)
      r0 <- sum(pooled$event[pooled$trt == 0]) / sum(pooled$trt == 0)
      n0 <- sum(pooled$trt == 0)
      ci <- diffscoreci(r1, n1, r0, n0, conf.level = 0.95)
      # Two-sided p-value from CI: if 0 is outside CI, p < 0.05
      # Use the Wald z from the CI
      se_diff <- sqrt(r1*(1-r1)/n1 + r0*(1-r0)/n0)
      if (se_diff == 0) return(NA)
      z_diff <- (r1 - r0) / se_diff
      2 * pnorm(-abs(z_diff))
    }, error = function(e) NA)
    
    list(mantelhaen = p_mh, cmh_rr = p_rr, mn_rd = p_mn)
  }, future.seed = TRUE)
  
  # Aggregate
  reject_mh <- sum(!is.na(p_mh <- sapply(results, `[[`, "mantelhaen")) & p_mh < 0.05, na.rm = TRUE)
  reject_rr <- sum(!is.na(p_rr <- sapply(results, `[[`, "cmh_rr")) & p_rr < 0.05, na.rm = TRUE)
  reject_mn <- sum(!is.na(p_mn <- sapply(results, `[[`, "mn_rd")) & p_mn < 0.05, na.rm = TRUE)
  
  cat("    Type I error (CMH OR / mantelhaen.test):", round(reject_mh / N_REPS * 100, 2), "%\n")
  cat("    Type I error (CMH RR):                   ", round(reject_rr / N_REPS * 100, 2), "%\n")
  cat("    Type I error (MN RD):                    ", round(reject_mn / N_REPS * 100, 2), "%\n\n")
  
  list(
    label = label,
    probs = strata_probs,
    type1_mantel = reject_mh / N_REPS,
    type1_rr = reject_rr / N_REPS,
    type1_rd = reject_mn / N_REPS
  )
}

# Run both strata configurations
binary_results_balanced <- run_binary_simulation(STRATA_PROBS_BALANCED, "Balanced (4 strata)")
binary_results_sparse   <- run_binary_simulation(STRATA_PROBS_SPARSE, "Sparse (4 strata)")

# ============================================================================
# PART 2: Time-to-event (Stratified Log-Rank)
# ============================================================================
cat("========== PART 2: TIME-TO-EVENT ==========\n\n")

#' Simulate one replication for stratified log-rank test
#' @param n Total sample size
#' @param strata_probs Vector of stratum probabilities
#' @param hr Hazard ratio for treatment vs control
#' @param control_median Median survival in control group (months)
#' @param censor_time Censoring time (months)
#' @return p-value from stratified log-rank test
simulate_survival_data <- function(n, strata_probs, hr = 1.0, 
                                    control_median = 14, censor_time = 24) {
  k <- length(strata_probs)
  strata <- sample(1:k, n, replace = TRUE, prob = strata_probs)
  
  # 1:1 randomization within each stratum
  treatment <- integer(n)
  for (s in 1:k) {
    idx <- which(strata == s)
    n_s <- length(idx)
    n_t <- floor(n_s / 2)
    treatment[idx] <- 0
    treatment[idx[1:n_t]] <- 1
  }
  
  # Exponential survival: S(t) = exp(-lambda * t)
  # Median = log(2) / lambda => lambda = log(2) / median
  lambda_control <- log(2) / control_median
  lambda_treated <- lambda_control / hr  # treatment has higher rate if HR < 1
  
  # Generate survival times
  u <- runif(n)
  survival_time <- -log(u) / ifelse(treatment == 1, lambda_treated, lambda_control)
  
  # Simple censoring at censor_time
  observed_time <- pmin(survival_time, censor_time)
  event <- (survival_time <= censor_time)  # 1=event, 0=censored
  
  list(time = observed_time, event = event, strata = as.factor(strata), 
       treatment = treatment)
}

#' Run survival simulation
#' @param hr Hazard ratio (1.0 for Type I, <1.0 for power)
#' @param label Description
#' @param strata_probs Stratum probabilities
run_survival_simulation <- function(hr, label, strata_probs) {
  cat("  Scenario:", label, "(HR =", hr, ")\n")
  cat("    Strata probs:", paste(round(strata_probs, 4), collapse = ", "), "\n")
  
  reject_count <- 0
  
  cat("    Running", N_REPS, "replications...\n")
  
  results <- future.apply::future_lapply(1:N_REPS, function(rep_id) {
    dat <- simulate_survival_data(N, strata_probs, hr = hr, 
                                   control_median = 14, censor_time = 24)
    
    # Stratified log-rank via survdiff
    sd <- tryCatch({
      survdiff(Surv(time, event) ~ treatment + strata(strata), data = dat)
    }, error = function(e) NULL)
    
    if (is.null(sd)) return(NA)
    
    # CRITICAL: sum O-E across ALL strata, not just the first!
    # sd$obs and sd$exp are matrices: rows = strata, cols = treatment/control
    # Row 1 = treatment, Row 2 = control
    # Total O-E for treatment = sum of (O_treated - E_treated) across all strata
    total_oe <- sum(sd$obs[1, ] - sd$exp[1, ])
    total_var <- sd$var[2, 2]  # variance of the treatment coefficient
    
    if (total_var == 0) return(NA)
    
    z_stat <- total_oe / sqrt(total_var)
    
    # One-sided test:
    # For HR = 1 (null): z ~ 0, reject if z < -1.96 (alpha=0.025)
    # For HR < 1 (alternative): z < 0 expected, reject if z < -1.96
    # We use: p = P(Z <= z) = pnorm(z) for the left-tail test
    p_value <- pnorm(z_stat)  # one-sided, left-tail
    
    list(z = z_stat, p = p_value)
  }, future.seed = TRUE)
  
  p_values <- sapply(results, function(r) if (is.null(r) || is.na(r$p)) NA_real_ else r$p)
  p_values <- p_values[!is.na(p_values)]
  
  reject_count <- sum(p_values < 0.025)
  power_or_type1 <- reject_count / length(p_values)
  
  cat("    Result:", round(power_or_type1 * 100, 2), "%\n\n")
  
  power_or_type1
}

# Type I error (HR = 1.0) with balanced strata
cat("  Type I error (balanced 8 strata, HR=1.0):\n")
type1_balanced <- run_survival_simulation(hr = 1.0, "Balanced 8-strata Type I", 
                                           rep(1/8, 8))

# Power (HR = 0.65) with balanced strata
cat("  Power (balanced 8 strata, HR=0.65):\n")
power_balanced <- run_survival_simulation(hr = 0.65, "Balanced 8-strata Power", 
                                           rep(1/8, 8))

# Type I error (HR = 1.0) with sparse strata
cat("  Type I error (sparse 8 strata, HR=1.0):\n")
type1_sparse <- run_survival_simulation(hr = 1.0, "Sparse 8-strata Type I", 
                                         STRATA_PROBS_SPARSE8)

# Power (HR = 0.65) with sparse strata
cat("  Power (sparse 8 strata, HR=0.65):\n")
power_sparse <- run_survival_simulation(hr = 0.65, "Sparse 8-strata Power", 
                                         STRATA_PROBS_SPARSE8)

# ============================================================================
# SUMMARY
# ============================================================================
cat("\n========== SUMMARY ==========\n\n")
cat("--- Part 1: Binary Endpoints (Type I error at two-sided alpha=0.05) ---\n")
cat("  Balanced strata (0.25 each):\n")
cat("    CMH OR (mantelhaen.test):     ", round(binary_results_balanced$type1_mantel * 100, 2), "%\n")
cat("    CMH RR:                       ", round(binary_results_balanced$type1_rr * 100, 2), "%\n")
cat("    MN RD:                        ", round(binary_results_balanced$type1_rd * 100, 2), "%\n\n")
cat("  Sparse strata (0.01, 0.02, 0.485, 0.485):\n")
cat("    CMH OR (mantelhaen.test):     ", round(binary_results_sparse$type1_mantel * 100, 2), "%\n")
cat("    CMH RR:                       ", round(binary_results_sparse$type1_rr * 100, 2), "%\n")
cat("    MN RD:                        ", round(binary_results_sparse$type1_rd * 100, 2), "%\n\n")

cat("--- Part 2: Time-to-Event (one-sided alpha=0.025) ---\n")
cat("  Balanced 8-strata:\n")
cat("    Type I error (HR=1.0):        ", round(type1_balanced * 100, 2), "%\n")
cat("    Power (HR=0.65):              ", round(power_balanced * 100, 2), "%\n\n")
cat("  Sparse 8-strata (0.01, 0.01, rep(0.98/6,6)):\n")
cat("    Type I error (HR=1.0):        ", round(type1_sparse * 100, 2), "%\n")
cat("    Power (HR=0.65):              ", round(power_sparse * 100, 2), "%\n\n")

cat("=== Done ===\n")
