# ============================================================
# Independent simulation: Does sparse stratification reduce
# power of the stratified log-rank test?
# ============================================================
#
# Key insight: With no censoring and HR=0.65, power is ~100%
# in both scenarios. Need to add realistic censoring and test
# a wider range of effect sizes.

library(survival)
library(furrr)

N            <- 500
n_strata     <- 4
n_reps       <- 5000
alpha        <- 0.025

# ---- Helper ----
run_one <- function(seed, probs, HR, censor_rate) {
  set.seed(seed)
  rate_ctrl    <- log(2) / 14
  rate_treat   <- rate_ctrl * HR
  stratum      <- sample(1:n_strata, N, replace = TRUE, prob = probs)
  trt          <- rbinom(N, 1, 0.5)
  rate         <- ifelse(trt == 1, rate_treat, rate_ctrl)
  time         <- rexp(N, rate = rate)
  
  # Add uniform censoring over [0, 24] months
  censor_time  <- runif(N, 0, 24)
  time         <- pmin(time, censor_time)
  event        <- as.integer(time < censor_time)
  
  # Check that stratified log-rank can run (need events in each stratum for each group)
  if (any(table(stratum, event) == 0)) {
    return(list(power = 0, events = sum(event), min_n = min(table(stratum)),
                zero_cells = TRUE))
  }
  
  fit <- tryCatch(
    survdiff(Surv(time, event) ~ trt + strata(stratum)),
    error = function(e) NULL
  )
  if (is.null(fit)) {
    return(list(power = 0, events = sum(event), min_n = min(table(stratum)),
                zero_cells = TRUE))
  }
  
  O1  <- fit$n[1]; E1  <- sum(fit$exp[1,]); V11 <- fit$var[1, 1]
  if (V11 <= 0) return(list(power = 0, events = sum(event), min_n = min(table(stratum)),
                             zero_cells = TRUE))
  z   <- (O1 - E1) / sqrt(V11)
  pval <- pnorm(z, lower.tail = FALSE)
  
  list(power   = as.numeric(pval < alpha),
       events  = sum(event),
       min_n   = min(table(stratum)),
       zero_cells = FALSE)
}

plan(multisession, workers = 11)

# ---- Primary analysis: HR=0.65 with realistic censoring ----
cat("=== PRIMARY ANALYSIS: HR=0.65 with uniform censoring [0,24] months ===\n\n")

HR <- 0.65
censor_rate <- 1/24
seeds <- sample.int(1e8, n_reps)

cat("Balanced...\n")
res_bal <- future_map(seeds, ~ run_one(.x, rep(0.25, n_strata), HR, censor_rate),
                      .options = furrr_options(seed = NULL))
cat("Sparse...\n")
res_spa <- future_map(seeds, ~ run_one(.x, c(0.01, 0.02, 0.485, 0.485), HR, censor_rate),
                      .options = furrr_options(seed = NULL))

power_bal  <- mean(sapply(res_bal, `[[`, "power"))
power_spa  <- mean(sapply(res_spa, `[[`, "power"))
events_bal <- mean(sapply(res_bal, `[[`, "events"))
events_spa <- mean(sapply(res_spa, `[[`, "events"))
min_bal    <- median(sapply(res_bal, `[[`, "min_n"))
min_spa    <- median(sapply(res_spa, `[[`, "min_n"))
fail_bal   <- sum(sapply(res_bal, `[[`, "zero_cells"))
fail_spa   <- sum(sapply(res_spa, `[[`, "zero_cells"))

cat(sprintf("  Balanced power:   %.4f\n", power_bal))
cat(sprintf("  Sparse power:     %.4f\n", power_spa))
cat(sprintf("  Power diff:       %.4f\n", power_bal - power_spa))
cat(sprintf("  Mean events:      %.1f / %.1f\n", events_bal, events_spa))
cat(sprintf("  Min-stratum (med): %.0f / %.0f\n", min_bal, min_spa))
cat(sprintf("  Failed (zero cells): %d / %d\n\n", fail_bal, fail_spa))

# ---- Power across HRs with censoring ----
cat("=== POWER ACROSS HR VALUES (with censoring) ===\n\n")

HR_values <- c(0.60, 0.65, 0.70, 0.75, 0.80, 0.85, 0.90, 0.95)

for (HR in HR_values) {
  seeds <- sample.int(1e8, n_reps)
  
  res_bal <- future_map(seeds, ~ run_one(.x, rep(0.25, n_strata), HR, censor_rate),
                        .options = furrr_options(seed = NULL))
  res_spa <- future_map(seeds, ~ run_one(.x, c(0.01, 0.02, 0.485, 0.485), HR, censor_rate),
                        .options = furrr_options(seed = NULL))
  
  power_bal  <- mean(sapply(res_bal, `[[`, "power"))
  power_spa  <- mean(sapply(res_spa, `[[`, "power"))
  events_bal <- mean(sapply(res_bal, `[[`, "events"))
  events_spa <- mean(sapply(res_spa, `[[`, "events"))
  min_bal    <- median(sapply(res_bal, `[[`, "min_n"))
  min_spa    <- median(sapply(res_spa, `[[`, "min_n"))
  fail_bal   <- sum(sapply(res_bal, `[[`, "zero_cells"))
  fail_spa   <- sum(sapply(res_spa, `[[`, "zero_cells"))
  
  cat(sprintf("HR=%.2f:  bal=%.4f  spa=%.4f  diff=%.4f  ev=%.0f/%.0f  min=%.0f/%.0f  fail=%d/%d\n",
              HR, power_bal, power_spa, power_bal-power_spa,
              events_bal, events_spa, min_bal, min_spa, fail_bal, fail_spa))
}

cat("\n=== END ===\n")
