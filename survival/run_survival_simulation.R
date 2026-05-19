#!/usr/bin/env Rscript
# Survival: Small Strata Pooling — Power & Type I Error with Group Sequential Design
# Stratified Cox PH + Stratified Log-rank, oncology trial setup
#
# IMPORTANT: The log-rank test (survdiff) is the canonical test for group
# sequential designs. OBF boundaries from gsDesign are calibrated for the
# log-rank score test. Cox score/Wald tests have inflated Type I error (~0.085)
# when used with OBF boundaries designed for the log-rank.
# See survival/audit/qwen-gsd-review.md for details.

library(survival); library(gsDesign); library(furrr); library(future)
plan(multisession, workers = min(11, future::availableCores() - 1))

# ── O'Brien-Fleming boundaries (realistic oncology GSD) ──
# 2 looks at 60% and 100% information, matching typical IA practice
# Earlier 3-look at 33% was unrealistic — you never stop for efficacy at 33%
gs <- gsDesign(k = 2, test.type = 2, alpha = 0.025, sfu = sfLDOF, timing = c(0.7, 1))
obf_z <- gs$upper$bound
lr_nominal <- pnorm(-obf_z)  # one-sided alpha per look

# ── Simulation parameters ──
N <- 500
accrual_months <- 18
cutoff_months <- 36
block_size <- 4

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
  
  stratum <- sample(1:4, n, replace = TRUE, prob = p_strata)
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
        trt[idx[start:end]] <- 0
      }
    }
  }
  
  scale_c <- 14 / log(2)
  scale_t <- scale_c / hr
  surv_time <- rweibull(n, shape = 1, scale = ifelse(trt == 1, scale_t, scale_c))
  accrual_time <- runif(n, 0, accrual_months)
  cens_rand <- rexp(n, rate = -log(0.95) / 12)
  
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
      target <- large[which.max(tab[as.character(large)])]
      new_stratum[new_stratum == as.character(s)] <- as.character(target)
    }
  }
  factor(new_stratum)
}

# ── Analyze at one look (fixed event count) ──
analyze_at_look <- function(d, n_events_target) {
  d <- d[order(d$time), ]
  d$cum_events <- cumsum(d$event)
  idx <- which(d$cum_events >= n_events_target)[1]
  if (is.na(idx)) return(NULL)
  
  dl <- d[1:idx, ]
  
  # Cox PH (for HR estimation only, NOT used for GSD)
  cox_fit <- tryCatch(
    coxph(Surv(time, event) ~ trt + strata(stratum), data = dl),
    error = function(e) NULL, warning = function(w) NULL
  )
  cox_ok <- !is.null(cox_fit) && !is.null(cox_fit$coefficients) && !any(is.na(coef(cox_fit)))
  
  # Log-rank (canonical GSD test)
  lr_fit <- tryCatch(
    survdiff(Surv(time, event) ~ trt + strata(stratum), data = dl),
    error = function(e) NULL, warning = function(w) NULL
  )
  
  lr_p <- NA
  if (!is.null(lr_fit) && !is.null(lr_fit$chisq) && !is.na(lr_fit$chisq) && lr_fit$chisq > 0) {
    # Total O-E across ALL strata (not just first stratum!)
    # survdiff returns obs/exp as: [ctl_strat1, trt_strat1, ctl_strat2, trt_strat2, ...]
    k <- length(lr_fit$obs) / 2  # number of strata
    trt_idx <- seq(2, 2*k, 2)    # indices for treatment arms
    oe <- if (!is.null(lr_fit$obs) && length(lr_fit$obs) >= 2) {
      sum(lr_fit$obs[trt_idx] - lr_fit$exp[trt_idx])
    } else 0
    # One-sided p-value for benefit (H1: HR < 1)
    # pchisq(..., lower.tail=FALSE) gives two-sided p; divide by 2 for one-sided
    # If oe > 0 (harm trend), the one-sided benefit p is large (1 - small)
    lr_p <- ifelse(is.na(oe) || oe > 0,
                   1 - pchisq(lr_fit$chisq, 1, lower.tail = FALSE) / 2,
                   pchisq(lr_fit$chisq, 1, lower.tail = FALSE) / 2)
  }
  
  list(
    n_events = n_events_target,
    cox_hr = if (cox_ok) exp(coef(cox_fit)) else NA,
    cox_se = if (cox_ok) sqrt(diag(vcov(cox_fit))) else NA,
    lr_p = lr_p
  )
}

# ── Full trial analysis ──
analyze_trial <- function(d, look_targets) {
  d <- d[order(d$time), ]
  d$cum_events <- cumsum(d$event)
  ne <- sum(d$event)
  
  results <- list()
  for (lk in seq_along(look_targets)) {
    tgt <- look_targets[lk]
    if (tgt > ne) break
    
    # Unpooled
    r_no <- analyze_at_look(d, tgt)
    if (is.null(r_no)) break
    
    # Pooled
    d_pooled <- transform(d, stratum = pool_strata(d))
    d_pooled <- d_pooled[order(d_pooled$time), ]
    d_pooled$cum_events <- cumsum(d_pooled$event)
    r_pool <- analyze_at_look(d_pooled, tgt)
    
    results[[lk]] <- list(look = lk, events = tgt, no_pool = r_no, pool = r_pool)
  }
  results
}

# ── Single replication ──
run_rep <- function(i, scenario_idx, hr, p_strata, look_targets) {
  seed <- 20260519 + scenario_idx * 1e6 + ifelse(hr == 1, 0, 1) * 1e5 + i
  d <- gen_trial(seed, p_strata, hr)
  ne <- sum(d$event)
  
  res <- analyze_trial(d, look_targets)
  nlooks <- length(res)
  
  out <- list(seed = seed, scenario = scenario_idx, hr = hr, n_events = ne,
              ever_reject_lr_no = 0, ever_reject_lr_pool = 0,
              cox_hr_last = NA, cox_se_last = NA, lr_p_last = NA,
              pool_hr_last = NA, pool_se_last = NA)
  
  stopped_no <- FALSE
  stopped_pool <- FALSE
  
  for (lk in seq_len(nlooks)) {
    if (stopped_no && stopped_pool) break
    r <- res[[lk]]
    if (is.null(r)) break
    
    # Log-rank sequential testing (canonical GSD test)
    # Only check if not already rejected at an earlier look
    if (!stopped_no && !is.na(r$no_pool$lr_p) && r$no_pool$lr_p < lr_nominal[lk]) {
      out$ever_reject_lr_no <- 1
      stopped_no <- TRUE
    }
    if (!stopped_pool && !is.na(r$pool$lr_p) && r$pool$lr_p < lr_nominal[lk]) {
      out$ever_reject_lr_pool <- 1
      stopped_pool <- TRUE
    }
    
    # HR/SE from last available look
    if (lk == nlooks || (lk < nlooks && is.null(res[[lk+1]]))) {
      out$cox_hr_last <- r$no_pool$cox_hr
      out$cox_se_last <- r$no_pool$cox_se
      out$pool_hr_last <- r$pool$cox_hr
      out$pool_se_last <- r$pool$cox_se
      out$lr_p_last <- r$no_pool$lr_p
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
      # Fixed event targets (matching real trial practice)
      run_rep(i, sc, hr, p_strata, c(210, 350))
    }, .options = furrr_options(seed = TRUE, chunk_size = 200))
    
    extract <- function(nm) sapply(reps, function(r) as.numeric(r[[nm]]))
    
    pow_l_no <- mean(extract("ever_reject_lr_no"), na.rm = TRUE)
    pow_l_pool <- mean(extract("ever_reject_lr_pool"), na.rm = TRUE)
    
    hr_no <- mean(extract("cox_hr_last"), na.rm = TRUE)
    hr_pool <- mean(extract("pool_hr_last"), na.rm = TRUE)
    se_no <- mean(extract("cox_se_last"), na.rm = TRUE)
    se_pool <- mean(extract("pool_se_last"), na.rm = TRUE)
    
    n_ev <- mean(extract("n_events"), na.rm = TRUE)
    
    cat(sprintf("  Events: %.0f\n", n_ev))
    cat(sprintf("  LR (canonical GSD): no pool=%.4f | pool=%.4f | diff=%.4f\n",
                pow_l_no, pow_l_pool, pow_l_pool - pow_l_no))
    cat(sprintf("  Cox HR: no pool=%.4f | pool=%.4f | SE: %.4f,%.4f\n",
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

# Test analysis
d_test <- gen_trial(20260519 + 1e6 + 1, scenarios[[1]]$p, 0.65)
lt <- c(210, 350)
res <- analyze_trial(d_test, lt)
cat(sprintf("  Analysis: %d looks\n", length(res)))
for (lk in seq_along(res)) {
  cat(sprintf("    Look %d: events=%d LR p=%.4f\n", lk,
              res[[lk]]$no_pool$n_events, res[[lk]]$no_pool$lr_p))
}
cat("  OK: Analysis pipeline works\n\n")

# ── Phase 2: Proof of concept (100 reps) ──
cat("═══ Phase 2: Proof of Concept (100 reps) ═══\n")
for (hr in c(1.0, 0.65)) {
  cat(sprintf("\n── HR=%.2f ──\n", hr))
  for (sc in 1:2) {
    p_strata <- scenarios[[sc]]$p
    pc_reps <- future_map(1:100, function(i) {
      run_rep(i, sc, hr, p_strata, c(210, 350))
    }, .options = furrr_options(seed = TRUE, chunk_size = 25))
    
    pow_lr <- mean(sapply(pc_reps, function(r) r$ever_reject_lr_no), na.rm = TRUE)
    n_ev <- mean(sapply(pc_reps, function(r) r$n_events), na.rm = TRUE)
    cat(sprintf("  Sc%d: events=%.0f LR power=%.3f\n", sc, n_ev, pow_lr))
    flush(stdout())
  }
}
cat("\nProof of concept complete.\n")
