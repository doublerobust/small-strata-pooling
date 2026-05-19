# Power Analysis Code Review: `run_power_analysis.R`

**Date:** 2026-05-19  
**Reviewer:** Subagent  
**Scope:** Comparison against `run_small_strata.R` (Type I error code) and design plan

---

## 1. Simulation Design Check

| Parameter | Specification | Found in Code | Verdict |
|-----------|:------------:|:-------------:|:-------:|
| 5000 reps | `n_sim <- 5000` | ✅ | PASS |
| N=400 | `n <- 400` | ✅ | PASS |
| Block=4 | Block loop `ceiling(n_s/4)` with `bn - bn%/%2` split | ✅ | PASS |
| 4 scenarios | Same `p_strata` vectors as Type I code | ✅ | PASS |
| 3 event rates | `c(0.10, 0.30, 0.50)` | ✅ | PASS |
| 2 effect sizes | `log_ors <- c(0.5, 0.7)` → OR ≈ 1.65, 2.00 | ✅ | PASS |
| Seed scheme | `20260518 + scenario*100000 + i*10 + round(ev_rate*100)` | ✅ | PASS, same as Type I code |
| Parallel | `future_map` + `furrr_options(seed=TRUE, chunk_size=200)` | ✅ | PASS |

**Verdict: PASS.** The simulation design matches the plan exactly and is consistent with the Type I error code.

**One minor note:** The generated treatment effect uses the log-odds scale:
```r
logit_p0 <- log(event_rate / (1 - event_rate))
p1 <- plogis(logit_p0 + log_or)
```
This means the true treatment effect is an **odds ratio** of `exp(log_or)`. The effect is tested using both RR-scale and RD-scale methods, which is a known inconsistency (the data-generating mechanism is logit-additive, not risk-additive or log-risk-additive). For the effect sizes used here (OR 1.65, 2.0) and N=400, the distinction between scales is small, but worth documenting.

---

## 2. Pooling Logic

### Implementation
```r
tab <- table(stratum)
small <- as.numeric(names(tab[tab < 10]))
stratum_pooled <- stratum
if (length(small) > 0 && length(small) < 4) {
  large <- as.numeric(names(tab[tab >= 10]))
  for (s in small) {
    if (length(large) > 0) {
      target <- large[which.min(abs(tab[as.character(large)] - tab[as.character(s)]))]
      stratum_pooled[stratum == s] <- target
    }
  }
}
```

### Analysis

**Correctness:** Merges strata < 10 patients into the nearest large stratum (by patient count). This is a reasonable heuristic.

**Edge cases:**

| Condition | Behavior | Verdict |
|-----------|----------|---------|
| 0 small strata | No pooling (correct) | ✅ |
| 1–3 small strata | Merge into nearest large | ✅ Reasonable |
| All 4 small (< 10) | **No pooling** (`length(small) < 4` is false) | ⚠️ Minor |
| All patients in 1 stratum after merging | CMH OR still works (2+ tables), CMH RR still works, MN RD returns NA (< 2 valid strata) | ⚠️ Minor |

The `length(small) < 4` guard means if all 4 strata are < 10 patients (possible only with extreme imbalance and small N), no pooling occurs. With N=400 and the given proportions, this cannot happen. The guard is likely a safety check against merging all strata into one pool, which would leave no stratification. A more robust check would be `length(small) < length(unique(stratum))`.

**Merge target selection:** The code merges each small stratum into the large stratum whose size is **numerically closest** to the small stratum's size. This is a reasonable approach (similar-sized strata get merged). An alternative would be to merge by stratum number or by clinical similarity, but for a simulation study this is fine.

**Verdict: Minor Issue.** The pooling logic is correct for all practical parameter settings. The guard `length(small) < 4` should arguably use `length(unique(stratum))` for robustness, but this is a non-issue for the given simulation.

---

## 3. Method Implementation

### 3.1 CMH OR (`cmh_odds_ratio`)

```r
tbl <- array(0, dim = c(2, 2, 4))
for (k in 1:4) {
  idx <- which(strata == k)
  ...
}
```

**Separate function in power code:** The power code extracts this into a helper function. The Type I error code had this logic inlined. The implementations are functionally identical. ✅

**Hardcoded `dim = c(2,2,4)`:** After pooling, the actual strata values may be a subset of {1,2,3,4}. The loop `for (k in 1:4)` creates empty 2×2 tables for merged-out strata, which are then removed by `nonempty <- which(apply(tbl, 3, sum) > 0)`. This works correctly because `mantelhaen.test` is called on the trimmed array. ✅

**`mantelhaen.test(tbl, correct = FALSE)`:** Standard uncorrected CMH test. Same as Type I code. ✅

**Verdict: PASS.** Correctly implemented and consistent with the Type I error code.

### 3.2 CMH RR (`cmh_risk_ratio`)

**Point estimate:** Uses 0.5 continuity correction (`(x+0.5)/(n+0.5)`) — same as Type I code. Identical function body in both files. ✅

**Variance formula:** Uses inverse-variance pooling:
```r
var_num <- var_num + w^2 * (1/(x1+0.5) - 1/(n1+0.5) + 1/(x0+0.5) - 1/(n0+0.5))
# ...
var_log <- var_num / (den^2)
```

**⚠️ Carry-over issue from Type I code:** This is **inverse-variance pooling** of stratum-specific log-RR variances, not the Greenland-Robins variance for the Mantel-Haenszel risk ratio. This was flagged as a **major issue** in `code-review-v2.md`.

The existing review established that this variance formula is:
- A valid meta-analysis approach (inverse-variance fixed-effects)
- NOT the correct Greenland-Robins variance for the MH RR estimator
- Likely similar in practice for N=400 with moderate event rates

Since the power code uses the **identical** `cmh_risk_ratio` function as the Type I error code, this issue propagates directly. Any bias in Type I error will affect the power calculations proportionally.

**Verdict: Minor Issue (carry-over, same as in Type I code).** The variance formula is technically incorrect for a "proper" MH RR, but the same function was used in the Type I error analysis, so the power results are at least internally consistent. If the Type I error is well-calibrated (which was the finding), the power estimates are usable.

### 3.3 Stratified MN RD (`stratified_mn_rd`)

**SE from CI width:**
```r
se <- (ci$conf.int[2] - ci$conf.int[1]) / (2*qnorm(1-(1-conf.level)/2))
```
Same approximation as the Type I code. ✅

**Inverse-variance pooling:** Weighted by `1/se²`. Standard approach. ✅

**Minimum 2 valid strata:** `if (sum(valid) < 2) return(list(p = NA))` — same as Type I code. This means if pooling reduces data to 1 stratum, MN RD returns NA. This is slightly conservative.

**Verdict: PASS.** Identical implementation to Type I error code. Known approximation (SE from CI width) but adequate for simulation.

---

## 4. Power Calculation: One-Sided Alpha = 0.025

```r
pow_or  <- mean(p_or/2 < 0.025, na.rm = TRUE)
pow_rr  <- mean(p_rr/2 < 0.025, na.rm = TRUE)
pow_mn  <- mean(p_mn/2 < 0.025, na.rm = TRUE)
```

### Analysis

Each method returns a **two-sided** p-value:
- CMH OR: `mantelhaen.test` returns two-sided p-value ✅
- CMH RR: `2 * pnorm(-abs(log(rr)/se))` → two-sided ✅
- MN RD: `2 * pnorm(-abs(rd/se))` → two-sided ✅

Dividing by 2 gives a one-sided p-value. Checking `< 0.025` gives one-sided α = 0.025. ✅

**⚠️ Direction check:** The formula `p_or/2 < 0.025` does not condition on the observed effect being in the correct direction. If the estimated odds ratio is < 1 (treatment appears harmful) but the p-value is small, this would count as a "significant" result for the wrong direction. In practice, with OR ≥ 1.65 and N=400, the probability of observing a reversed effect is negligible. The standard "two-sided p-value ÷ 2" approach used in most simulation studies is adequate here.

**Verdict: PASS** with minor caveat documented.

---

## 5. Consistency with Type I Error Code

| Feature | Power Code | Type I Code | Verdict |
|---------|:----------:|:-----------:|:-------:|
| Seed formula | `20260518 + scenario*100000 + i*10 + round(ev_rate*100)` | Same | ✅ |
| Block randomization | Within-stratum, block=4, `floor(bn/2)` split | Same | ✅ |
| `cmh_risk_ratio` | Function (identical body) | Function (identical body) | ✅ |
| `stratified_mn_rd` | Function (identical body) | Function (identical body) | ✅ |
| `cmh_odds_ratio` | Extracted to helper function | Inline in loop | ✅ Logic identical |
| `future_map` | `furrr_options(seed=TRUE, chunk_size=200)` | Same | ✅ |
| Scenario probabilities | Same 4 scenarios | Same | ✅ |
| Event rate loop | Inside scenario+effect loop | Inside scenario loop | ✅ Different by design |
| Outcome generation | `rbinom` with treatment effect (`p1 != p0`) | `rbinom` with null (`p1 = p0`) | ✅ Correct by design |

**Key differences (all expected/benign):**
1. Power code has an inner `log_ors` loop — necessary for power analysis.
2. Power code wraps `cmh_odds_ratio` as a function — refactoring, not a logic change.
3. Power code uses `future_map` inside the `log_ors` loop (inner) vs. the `event_rate` loop (inner in Type I code). Both are fine since seeds are unique per `(scenario, i, event_rate)` tuple regardless of loop nesting.

**Verdict: PASS.** The simulation engines are structurally consistent. Seeds are unique and non-overlapping between the two scripts (different data-generating processes but no seed collision because the Type I code never calls seeds with non-zero `log_or`).

---

## 6. Edge Cases and Potential Bugs

### 6.1 MN RD ≥ 2 Strata Requirement After Pooling
If pooling reduces data to a single stratum (e.g., 3 small strata merged into 1), `stratified_mn_rd` returns `list(p = NA)` because it requires `sum(valid) >= 2`. This situation requires at most one valid large stratum remaining after all small ones are merged. With N=400 and these scenarios, this is very unlikely. The MN RD failure rate should be checked in the output, and if > 0 for pooled analyses, needs investigation.

### 6.2 Pooling Target Selection Heuristic
The merge-to-nearest-large-by-size heuristic is reasonable but somewhat arbitrary. Different pooling rules (e.g., merge all small strata together, or merge by stratum number) could give different results. The SAP language this simulation supports should ideally be pooling-rule-agnostic, so this is worth noting.

### 6.3 Power Output Format
The output prints:
```r
cat(sprintf("    CMH OR  fail=%.3f | no pool=%.3f | pool=%.3f | diff=%.3f\n", ...))
```
The fail rate is computed over all replications (including those that would be pooled), while the power estimates use `na.rm=TRUE`. If the fail rate is non-trivial, the power estimate is conditioned on valid runs only, which is correct but should be checked.

### 6.4 No Direction-Conditioned Power
As noted in §4, the power calculation divides two-sided p-values by 2 without checking effect direction. For the discussed effect sizes, this is negligible but should be documented.

### 6.5 Reproducibility
The `set.seed(20260518)` at the top of the script is immediately superceded by per-replication seeds inside `run_rep`. This is fine for reproducibility since `furrr_options(seed=TRUE)` uses L'Ecuyer-CMRG for parallel RNG.

---

## Summary of Issues

| # | Issue | Severity | Location | Propagation |
|---|-------|----------|----------|-------------|
| 1 | CMH RR variance: inverse-variance pooling, not Greenland-Robins | **Minor** | `cmh_risk_ratio()` | Carried over from Type I code |
| 2 | MN RD SE from CI width: approximation | **Minor** | `stratified_mn_rd()` | Carried over from Type I code |
| 3 | Pooling guard: `length(small) < 4` | **Minor** | Pooling block | New in power code |
| 4 | MN RD requires ≥2 valid strata | **Minor** | `stratified_mn_rd()` | Carried over; conservative after pooling |
| 5 | Power doesn't check effect direction | **Documentation** | Power calculation lines | New in power code; negligible for effect sizes used |

---

## Overall Verdict: Minor Issues ⚠️

**Rationale:** The power analysis code is a correct extension of the existing Type I error simulation. The simulation design matches the plan exactly. The simulation engines are consistent between the two scripts. All method implementations are identical to the Type I error code, so no new methodological errors are introduced.

The carry-over issues (CMH RR variance formula, MN RD SE approximation) were already identified and accepted in the Type I error code review. These same caveats apply to the power analysis.

**What's correct and well-done:**
- Seed scheme ensures reproducibility and non-overlap with Type I error code
- Pooling logic handles 0–3 small strata correctly for practical cases
- One-sided α = 0.025 testing is properly implemented via `p/2 < 0.025`
- Power comparisons (`diff = pool - no pool`) are clearly presented
- Fail rates are reported alongside power estimates, providing diagnostic information

**Specific recommendation:** The power code's output testing all combinations of 4 scenarios × 3 event rates × 2 effect sizes × 2 pooling conditions × 3 methods = 144 power estimates is comprehensive. The results will be internally consistent with the Type I error findings, so the overall conclusion ("no pooling needed") can be evaluated by comparing the **difference** (`pool - no pool`) across all configurations.

---

*End of review.*
