#!/usr/bin/env Rscript
# Survival: Small Strata Pooling — Full Simulation Run
# 2-look OBF (60%, 100%) — realistic oncology GSD design
# Log-rank canonical test for GSD decisions
source("run_survival_simulation.R")

# ── Full Type I error run (HR = 1.0, 10,000 reps) ──
run_sim(hr = 1.0, n_reps = 10000, sim_label = "Type I Error")

# ── Full Power run (HR = 0.65, 5,000 reps) ──
run_sim(hr = 0.65, n_reps = 5000, sim_label = "Power")

cat("\n═══════════════════════════════════════════════\n")
cat("FULL SIMULATION COMPLETE\n")
cat(Sys.time(), "\n")
cat("═══════════════════════════════════════════════\n")
