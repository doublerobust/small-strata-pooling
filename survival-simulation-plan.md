# Time-to-Event: Small Strata Pooling — Simulation Plan (v2)

## Objective
Investigate whether pooling small strata affects power, Type I error, or convergence for stratified Cox PH and stratified log-rank tests in oncology trials.

## Trial Design

| Parameter | Value |
|-----------|-------|
| N | 500 |
| Randomization | 1:1, stratified block (block size 4) |
| Stratification factors | 2 binary factors → 4 strata |
| Accrual period | 18 months, uniform enrollment (~28/month) |
| Follow-up | Administrative cutoff at 36 months from start |
| Random censoring | 5% annual loss to follow-up (exponential) |
| Control median survival | 14 months |
| Treatment effect | HR = 0.65 (realistic oncology effect) |
| Target events | ~350 at cutoff (~70% event rate) |

## Strata Proportions (4 scenarios, same as binary study)

| Scenario | Label | S1 | S2 | S3 | S4 |
|:--------:|:------|:--:|:--:|:--:|:--:|
| 1 | Balanced | 0.25 | 0.25 | 0.25 | 0.25 |
| 2 | 1 small (5%) | 0.05 | 0.35 | 0.30 | 0.30 |
| 3 | 2 small (3% ea) | 0.03 | 0.03 | 0.47 | 0.47 |
| 4 | 2 tiny (1%, 2%) | 0.01 | 0.02 | 0.485 | 0.485 |

## Survival Generation (Weibull — corrected)

**Shape parameter:** γ = 1 (exponential, constant hazard — simplest, adequate for this comparison)

**Scale parameter computation:**
- Control median = 14 months → scale_c = 14 / ln(2) = 20.20
- Treatment HR = 0.65 → scale_t = scale_c / 0.65 = 31.08
- Both have shape = 1 (exponential)

```r
scale_c <- 14 / log(2)  # ≈ 20.20
scale_t <- scale_c / HR  # ≈ 31.08 for HR = 0.65
T <- rweibull(n, shape = 1, scale = ifelse(trt == 1, scale_t, scale_c))
```

**Proportional hazards assumption:** Common HR across all strata (the standard PH assumption). The simulation tests precision in small strata, not PH violations.

## Accrual and Censoring

1. **Accrual time:** A ~ Uniform(0, 18) months
2. **Calendar time:** C = A + T (date of event if uncensored)
3. **Administrative censoring at cutoff (36 months):** event if C < 36
4. **Random censoring:** C_rand ~ Exp(rate = -ln(0.95)/12) → ~5% annual dropout
5. **Observed time:** min(T, C_rand, 36 - A)
6. **Observed event:** 1 if T is the minimum (and < cutoff)

## Group Sequential Design
- **Looks:** 3
- **Boundaries:** O'Brien-Fleming via `gsDesign::gsDesign(k=3, test.type=2, alpha=0.05, sfu="OF")`
- **Information fractions at looks:** 33%, 66%, 100% of target events
- **Target events for information fraction:** ~350 (approximate; adaptive to actual events observed)
- **Look triggers:** Looks fire when cumulative observed events ≥ milestone fractions of final events (not calendar-based)

## Methods

### Primary Analysis (No Pooling)
1. **Stratified Cox PH:** `coxph(Surv(time, event) ~ trt + strata(stratum))`
2. **Stratified log-rank:** `survdiff(Surv(time, event) ~ trt + strata(stratum))`

### Comparison (Pooling Small Strata)
- Pool strata with < 10 patients into the largest stratum by patient count
- **Tie-breaking:** If multiple large strata tie for size, merge into the one with the smallest stratum index (1→2→3→4)
- Re-run both analyses with pooled stratum assignment
- This rule matches the binary analysis for consistency

## Output Metrics
- **Cumulative power (ever-reject):** Proportion of reps where p < stopping boundary at ANY look for correct direction — primary metric
- **Power by look:** Proportion rejecting at each individual look
- **Type I error:** Same under null (HR = 1.0)
- **Convergence failure (Cox):** `tryCatch` wrapping `coxph` — if error/warning on convergence, flag as failure (`coxph.convergence != 0`, singular Hessian, or `iter.max` exceeded)
- **Log-rank failure:** Stratum with zero events in both arms → dropped silently by `survdiff`; track proportion of reps affected
- **HR bias:** mean(log(HR_hat) − true_log(HR)) — assess whether pooling biases effect estimate
- **Empirical SE of log(HR):** SD of log(HR_hat) across reps — assess precision impact
- **CI coverage:** proportion of reps where true HR is within coxph 95% CI
- **Power difference:** (pooled − unpooled) with paired delta-method SE
- **Stratum balance check:** mean absolute imbalance by stratum (treatment n − control n)

## Simulation Pipeline

| Phase | Reps | Purpose |
|:-----:|:----:|:--------|
| 1 | 2 | Code test — syntax, seed, data generation |
| 2 | 100 | Proof of concept — bug fixes, timing estimate |
| 3 | 10,000 | Type I error (HR = 1.0) |
| 4 | 5,000 | Power (HR = 0.65) — 5K sufficient for binary outcome precision |

## Reproducibility: Seed Scheme
```r
base_seed <- 20260519  # study date
# Per rep: seed = base_seed + scenario*1e6 + hr_id*1e5 + rep
# hr_id: 0 = null (HR=1), 1 = alt (HR=0.65)
# This ensures unique, non-overlapping seeds across all conditions
```

## R Packages Required
- `survival` — Cox PH, log-rank
- `gsDesign` — O'Brien-Fleming boundaries
- `furrr` — parallel execution

## Edge Cases
- **All strata small:** If pooling merges everything to 1 stratum, `strata()` produces a single-level factor. `coxph` with `+ strata()` on a single stratum = unstratified analysis. Log-rank similarly degrades. This is conservative (pooling effect exaggerated).
- **Cox with tiny strata:** No convergence issues expected with `coxph.fit`; tiny strata just contribute noise. If convergence fails, flag as failure.
- **Log-rank with zero-event strata:** `survdiff` handles these silently. Stratum with zero events in both arms contributes nothing.

## Data Generation Diagram
```
For each replication i:
  For each patient j in 1:n:
    stratum[j] ← sample(1:4, prob=p_strata)
    trt[j] ← block_randomization(stratum[j], block_size=4)
    T[j] ← rweibull(scale_c or scale_t, shape=1)
    A[j] ← runif(0, 18)  # accrual time
    C_adm[j] ← 36        # admin cutoff
    C_rand[j] ← rexp(rate=-ln(0.95)/12)
    obs_time[j] ← min(T[j], C_rand[j], 36 - A[j])
    obs_event[j] ← (T[j] == obs_time[j])
  
  For each look in 1:3:
    info_frac ← events_observed / 350
    if info_frac >= milestone[look]:
      Analyze with stratified Cox + log-rank (pooled & unpooled)
      Check stopping boundary
```

## References
- O'Brien & Fleming (1979). *A multiple testing procedure for clinical trials.* Biometrics.
- Greenland & Robins (1985). *Estimation of a common effect parameter from sparse follow-up data.* Biometrics.
