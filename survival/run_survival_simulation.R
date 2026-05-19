#!/usr/bin/env Rscript
# Survival: Small Strata Pooling — Power & Type I Error with Group Sequential Design
# Stratified Cox PH + Stratified Log-rank, oncology trial setup
library(survival); library(gsDesign); library(furrr); library(future)
plan(multisession, workers = min(11, future::availableCores() - 1))

# ── O'Brien-Fleming boundaries ──
gs <- gsDesign(k = 3, test.type = 2, alpha = 0.05, sfu = "OF")
obf_z <- gs$upper$bound  # z-value boundaries
obf_p <- 2 * pnorm(-obf_z)  # two-sided nominal alpha

# ── Simulation parameters ──
N <- 500
accrual_months <- 18
cutoff_months <- 36
block_size <- 4
hr_alt <- 0.65
n_sim <- list(type1 = 10000, power = 5000)

# ── Strata scenarios ──
scenarios <- list(
  list(label = "Balanced",       p = c(0.25, 0.25, 0.25, 0.25)),
  list(label = "1 small (5%)",   p = c(0.05, 0.35, 0.30, 0.30)),
  list(label = "2 small (3%)",   p = c(0.03, 0.03, 0.47, 0.47)),
  list(label = "2 tiny (1%,2%)", p = c(0.01, 0.02, 0.485, 0.485))
)

# ── Generate one trial ──
gen_trial <- function(seed, p_strata, hr, n = N) {
  set.seed(seed)
  
  # Stratum membership
  stratum <- sample(1:4, n, replace = TRUE, prob = p_strata)
  
  # Stratified block randomization
  trt <- integer(n)
  for (s in unique(stratum)) {
    idx <- which(stratum == s)
    ns <- length(idx)
    for (b in seq_len(ceiling(ns / block_size))) {
      start <- (b - 1) * block_size + 1
      end <- min(b * block_size, ns)
      bn <- end - start + 1
      n_trt <- floor(bn / 2)
      if (n_trt > 0 && bn - n_trt > 0) {
        trt[idx[start:end]] <- sample(c(rep(1, n_trt), rep(0, bn - n_trt)))
      } else {
        trt[idx[start:end]] <- 0  # all to control if odd block
      }
    }
  }
  
  # Survival times (Weibull, shape=1 = exponential)
  scale_c <- 14 / log(2)  # ~20.20 → median 14 months
  scale_t <- scale_c / hr
  surv_time <- rweibull(n, shape = 1, scale = ifelse(trt == 1, scale_t, scale_c))
  
  # Accrual time
  accrual_time <- runif(n, 0, accrual_months)
  
  # Calendar time of event if uncensored
  cal_time <- accrual_time + surv_time
  
  # Random censoring (5% annual dropout)
  cens_rand <- rexp(n, rate = -log(0.95) / 12)
  
  # Observed time and event
  obs_time <- pmin(surv_time, cens_rand, cutoff_months - accrual_time)
  obs_event <- as.integer(obs_time == surv_time & obs_time < cutoff_months - accrual_time)
  
  data.frame(stratum = factor(stratum), trt = trt,
             time = pmax(obs_time, 0.01), event = obs_event,
             stringsAsFactors = FALSE)
}

# ── Pool small strata ──
pool_strata <- function(d, threshold = 10) {
  tab <- table(d$stratum)
  small <- as.numeric(names(tab[tab < threshold]))
  if (length(small) == 0 || length(small) >= length(tab)) return(d$stratum)
  large <- as.numeric(names(tab[tab >= threshold]))
  new_stratum <- as.character(d$stratum)
  for (s in small) {
    if (length(large) > 0) {
      # Merge into largest large stratum (by count)
      target <- large[which.max(tab[as.character(large)])]
      new_stratum[new_stratum == as.character(s)] <- as.character(target)
    }
  }
  factor(new_stratum)
}

# ── Analyze at one look ──
analyze_look <- function(d) {
  # Stratified Cox
  cox_fit <- tryCatch(
    coxph(Surv(time, event) ~ trt + strata(stratum), data = d),
    error = function(e) NULL, warning = function(w) NULL
  )
  if (is.null(cox_fit) || is.null(cox_fit$coefficients) || any(is.na(coef(cox_fit)))) {
    cox_z <- NA; cox_hr <- NA; cox_se <- NA; cox_conv <- FALSE
  } else {
    cox_z <- coef(cox_fit) / sqrt(diag(vcov(cox_fit)))
    cox_hr <- exp(coef(cox_fit))
    cox_se <- sqrt(diag(vcov(cox_fit)))
    cox_conv <- TRUE
  }
  
  # Stratified log-rank
  lr_fit <- tryCatch(
    survdiff(Surv(time, event) ~ trt + strata(stratum), data = d),
    error = function(e) NULL, warning = function(w) NULL
  )
  if (is.null(lr_fit)) {
    lr_p <- NA; lr_conv <- FALSE
  } else {
    # One-sided p from chi-square statistic
    lr_chisq <- lr_fit$chisq
    lr_p <- if (is.null(lr_chisq) || is.na(lr_chisq) || lr_chisq < 0) NA else
      pchisq(lr_chisq, df = 1, lower.tail = FALSE) / 2
    # Check direction: observed events in treatment vs expected
    # survdiff returns obs and exp for each arm; check treatment arm (second element)
    if (!is.na(lr_p) && !is.null(lr_fit$obs) && !is.null(lr_fit$exp) && length(lr_fit$obs) >= 2) {
      lr_oe <- lr_fit$obs[2] - lr_fit$exp[2]  # treatment: O - E (negative = benefit)
      if (!is.na(lr_oe) && lr_oe > 0) lr_p <- 1 - lr_p  # harmful direction
    }
    lr_conv <- TRUE
  }
  
  c(cox_z = unname(cox_z), cox_hr = unname(cox_hr), cox_se = unname(cox_se),
    cox_conv = cox_conv, cox_p = if (is.na(cox_z)) NA else 2 * pnorm(-abs(cox_z)),
    lr_p = lr_p, lr_conv = lr_conv)
}

# ── Analyze trial at all looks ──
analyze_trial <- function(d) {
  nevents <- sum(d$event)
  # Target: 33%, 66%, 100% of events
  targets <- round(c(0.33, 0.66, 1.00) * nevents)
  # Ensure milestones don't exceed total events
  targets <- pmin(targets, nevents)
  # Remove duplicates and zeros
  targets <- unique(targets[targets > 0])
  
  # Sort by event time
  d <- d[order(d$time), ]
  d$cum_events <- cumsum(d$event)
  
  results <- list()
  for (look in seq_along(targets)) {
    tgt <- targets[look]
    look_idx <- which(d$cum_events >= tgt)[1]
    if (is.na(look_idx)) break
    d_look <- d[1:look_idx, ]
    
    # Unpooled
    res_no <- analyze_look(d_look)
    
    # Pooled
    d_look$stratum_pooled <- pool_strata(d_look)
    res_pool <- analyze_look(transform(d_look, stratum = stratum_pooled))
    
    results[[look]] <- list(
      look = look, events = look_idx,
      no_pool = res_no, pool = res_pool
    )
  }
  results
}

# ── Single replication ──
run_rep <- function(i, scenario_idx, hr, p_strata) {
  seed <- 20260519 + scenario_idx * 1e6 + ifelse(hr == 1, 0, 1) * 1e5 + i
  d <- gen_trial(seed, p_strata, hr)
  
  res <- analyze_trial(d)
  nlooks <- length(res)
  n_events <- sum(d$event)
  
  # Initialize output
  out <- list(
    seed = seed, scenario = scenario_idx, hr = hr, n_events = n_events,
    ever_reject_cox_no = 0, ever_reject_lr_no = 0,
    ever_reject_cox_pool = 0, ever_reject_lr_pool = 0,
    cox_conv_no = 0, cox_conv_pool = 0,
    lr_conv_no = 0, lr_conv_pool = 0,
    cox_hr_last = NA, cox_se_last = NA,
    pool_hr_last = NA, pool_se_last = NA
  )
  
  for (lk in seq_len(min(nlooks, 3))) {
    r <- res[[lk]]
    z_no <- r$no_pool["cox_z"]
    p_lr_no <- r$no_pool["lr_p"]
    z_pool <- r$pool["cox_z"]
    p_lr_pool <- r$pool["lr_p"]
    
    # Check OBF boundaries
    if (!is.na(z_no) && abs(z_no) >= obf_z[lk]) {
      # Check direction for power (HR < 1)
      if (hr < 1 && z_no < 0) out$ever_reject_cox_no <- 1
      if (hr == 1) out$ever_reject_cox_no <- 1  # Type I error = any rejection
    }
    if (!is.na(p_lr_no) && p_lr_no < 0.025) {
      out$ever_reject_lr_no <- 1
    }
    if (!is.na(z_pool) && abs(z_pool) >= obf_z[lk]) {
      if (hr < 1 && z_pool < 0) out$ever_reject_cox_pool <- 1
      if (hr == 1) out$ever_reject_cox_pool <- 1
    }
    if (!is.na(p_lr_pool) && p_lr_pool < 0.025) {
      out$ever_reject_lr_pool <- 1
    }
    
    # Convergence (last look only)
    if (lk == nlooks) {
      out$cox_conv_no <- r$no_pool["cox_conv"]
      out$cox_conv_pool <- r$pool["cox_conv"]
      out$lr_conv_no <- r$no_pool["lr_conv"]
      out$lr_conv_pool <- r$pool["lr_conv"]
      out$cox_hr_last <- r$no_pool["cox_hr"]
      out$cox_se_last <- r$no_pool["cox_se"]
      out$pool_hr_last <- r$pool["cox_hr"]
      out$pool_se_last <- r$pool["cox_se"]
    }
  }
  
  out
}

# ── Run simulation ──
run_sim <- function(hr, n_reps, sim_label) {
  cat(sprintf("\n═══ %s (HR=%.2f, %d reps) ═══\n", sim_label, hr, n_reps))
  
  for (sc in seq_along(scenarios)) {
    sc_label <- scenarios[[sc]]$label
    p_strata <- scenarios[[sc]]$p
    cat(sprintf("\n── Scenario %d: %s ──\n", sc, sc_label))
    
    reps <- future_map(1:n_reps, function(i) {
      run_rep(i, sc, hr, p_strata)
    }, .options = furrr_options(seed = TRUE, chunk_size = 200))
    
    # Aggregate
    extract <- function(nm) sapply(reps, function(r) as.numeric(r[[nm]]))
    
    pow_c_no <- mean(extract("ever_reject_cox_no"), na.rm = TRUE)
    pow_c_pool <- mean(extract("ever_reject_cox_pool"), na.rm = TRUE)
    pow_l_no <- mean(extract("ever_reject_lr_no"), na.rm = TRUE)
    pow_l_pool <- mean(extract("ever_reject_lr_pool"), na.rm = TRUE)
    
    conv_c_no <- mean(extract("cox_conv_no"), na.rm = TRUE)
    conv_c_pool <- mean(extract("cox_conv_pool"), na.rm = TRUE)
    conv_l_no <- mean(extract("lr_conv_no"), na.rm = TRUE)
    conv_l_pool <- mean(extract("lr_conv_pool"), na.rm = TRUE)
    
    hr_no <- mean(extract("cox_hr_last"), na.rm = TRUE)
    hr_pool <- mean(extract("pool_hr_last"), na.rm = TRUE)
    se_no <- mean(extract("cox_se_last"), na.rm = TRUE)
    se_pool <- mean(extract("pool_se_last"), na.rm = TRUE)
    
    n_ev <- mean(extract("n_events"), na.rm = TRUE)
    
    cat(sprintf("  Events: %.0f\n", n_ev))
    cat(sprintf("  Cox:    no pool=%.4f | pool=%.4f | diff=%.4f | conv=%.3f,%.3f\n",
                pow_c_no, pow_c_pool, pow_c_pool - pow_c_no, conv_c_no, conv_c_pool))
    cat(sprintf("  LR:     no pool=%.4f | pool=%.4f | diff=%.4f | conv=%.3f,%.3f\n",
                pow_l_no, pow_l_pool, pow_l_pool - pow_l_no, conv_l_no, conv_l_pool))
    cat(sprintf("  HR:     no pool=%.4f | pool=%.4f | SE: %.4f,%.4f\n",
                hr_no, hr_pool, se_no, se_pool))
    flush(stdout())
  }
}

# ── Phase 1: Code test (2 reps) ──
cat("═══ Phase 1: Code Test (2 reps) ═══\n")
for (sc in 1:2) {
  d <- gen_trial(20260519 + sc * 1e6 + 1, scenarios[[sc]]$p, 1.0)
  cat(sprintf("  Scenario %d: n=%d events=%.0f strata=%s\n",
              sc, nrow(d), sum(d$event),
              paste(table(d$stratum), collapse = ",")))
}
cat("  OK: Data generation works\n")

# Test analysis on one trial
d_test <- gen_trial(20260519 + 1e6 + 1, scenarios[[1]]$p, 0.65)
res <- analyze_trial(d_test)
cat(sprintf("  Analysis: %d looks\n", length(res)))
for (lk in seq_along(res)) {
  cat(sprintf("    Look %d: Cox z=%.3f LR p=%.3f\n", lk,
              res[[lk]]$no_pool["cox_z"], res[[lk]]$no_pool["lr_p"]))
}
cat("  OK: Analysis pipeline works\n\n")

# ── Phase 2: Proof of concept (100 reps, first scenario only) ──
cat("═══ Phase 2: Proof of Concept (100 reps) ═══\n")
for (hr in c(1.0, 0.65)) {
  cat(sprintf("\n── HR=%.2f ──\n", hr))
  for (sc in 1:2) {
    p_strata <- scenarios[[sc]]$p
    pc_reps <- future_map(1:100, function(i) {
      run_rep(i, sc, hr, p_strata)
    }, .options = furrr_options(seed = TRUE, chunk_size = 25))
    
    ever_c <- mean(sapply(pc_reps, function(r) r$ever_reject_cox_no), na.rm = TRUE)
    conv_c <- mean(sapply(pc_reps, function(r) r$cox_conv_no), na.rm = TRUE)
    n_ev <- mean(sapply(pc_reps, function(r) r$n_events), na.rm = TRUE)
    cat(sprintf("  Sc%d: events=%.0f power=%.3f conv=%.3f\n", sc, n_ev, ever_c, conv_c))
    flush(stdout())
  }
}
cat("\nProof of concept complete — proceeding to full runs.\n")
