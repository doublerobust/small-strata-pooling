#!/usr/bin/env Rscript
# Survival simulation with fixed event targets c(203, 290)
# Sources the original file (which auto-runs Phase 1 & 2 with old c(210,350) targets),
# then runs the actual simulation with c(203, 290).

source("/home/yue-shentu/.openclaw/workspace/merck-internal/small-strata/survival/run_survival_simulation.R")

# ── Fixed look targets ──
LOOK_TARGETS <- c(203, 290)

cat("\n\n")
cat("═══════════════════════════════════════════════════════════\n")
cat("  MAIN SIMULATION — Fixed targets: c(203, 290)\n")
cat("═══════════════════════════════════════════════════════════\n\n")

# ── Type I Error (HR=1.0, 10000 reps) ──
cat("️ Type I Error (HR=1.00, 10000 reps) ️\n")
for (sc in seq_along(scenarios)) {
  sc_label <- scenarios[[sc]]$label
  p_strata <- scenarios[[sc]]$p
  cat(sprintf("\n── Scenario %d: %s ──\n", sc, sc_label))
  
  reps <- future_map(1:10000, function(i) {
    run_rep(i, sc, 1.0, p_strata, LOOK_TARGETS)
  }, .options = furrr_options(seed = TRUE, chunk_size = 200))
  
  extract <- function(nm) sapply(reps, function(r) as.numeric(r[[nm]]))
  
  pow_l_no    <- mean(extract("ever_reject_lr_no"), na.rm = TRUE)
  pow_l_pool  <- mean(extract("ever_reject_lr_pool"), na.rm = TRUE)
  n_ev        <- mean(extract("n_events"), na.rm = TRUE)
  
  cat(sprintf("  Scenario %d: events=%.0f LR: no pool=%.4f | pool=%.4f | diff=%.4f\n",
              sc, n_ev, pow_l_no, pow_l_pool, pow_l_pool - pow_l_no))
  flush(stdout())
}

# ── Power (HR=0.65, 5000 reps) ──
cat("\n\n️ Power (HR=0.65, 5000 reps) ️\n")
for (sc in seq_along(scenarios)) {
  sc_label <- scenarios[[sc]]$label
  p_strata <- scenarios[[sc]]$p
  cat(sprintf("\n── Scenario %d: %s ──\n", sc, sc_label))
  
  reps <- future_map(1:5000, function(i) {
    run_rep(i, sc, 0.65, p_strata, LOOK_TARGETS)
  }, .options = furrr_options(seed = TRUE, chunk_size = 200))
  
  extract <- function(nm) sapply(reps, function(r) as.numeric(r[[nm]]))
  
  pow_l_no    <- mean(extract("ever_reject_lr_no"), na.rm = TRUE)
  pow_l_pool  <- mean(extract("ever_reject_lr_pool"), na.rm = TRUE)
  n_ev        <- mean(extract("n_events"), na.rm = TRUE)
  
  cat(sprintf("  Scenario %d: events=%.0f LR: no pool=%.4f | pool=%.4f | diff=%.4f\n",
              sc, n_ev, pow_l_no, pow_l_pool, pow_l_pool - pow_l_no))
  flush(stdout())
}

cat("\n═══════════════════════════════════════════════════════════\n")
cat("  SIMULATION COMPLETE\n")
cat("═══════════════════════════════════════════════════════════\n")
