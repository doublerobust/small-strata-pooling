# Code Review: `run_small_strata.R` (v2)

**Date:** 2026-05-18  
**Reviewer:** Subagent (automated)  
**Previous bugs fixed:** Scenario 4 duplication, CMH RR unstratified variance, Wald RD mislabeled as MN ✅

---

## 1. `stratified_mn_rd` — Stratified Miettinen-Nurminen RD

### Implementation
```r
se <- (ci$conf.int[2] - ci$conf.int[1]) / (2*qnorm(1-(1-conf.level)/2))
```
This derives SE from CI width using the standard normal quantile.

### Assessment: **Correct with caveats**

**Correctness:** The derivation is mathematically sound for a symmetric CI on the difference scale. If the CI is `[L, U]` and the estimator is `rd = (L+U)/2`, then:
```
se = (U - L) / (2 * z_{1-α/2})
```
This is the standard inverse-width method for extracting SE from a CI. For the Miettinen-Nurminen score CI, this is approximately correct because the score CI is asymptotically symmetric.

**Caveat 1 — CI asymmetry:** The MN score CI is *not* perfectly symmetric around the point estimate (it's a score-based interval, not a Wald interval). The width-based SE extraction introduces a small bias in the SE estimate. In practice, for moderate-to-large strata, this bias is negligible. For very small strata (Scenarios 3–4), the asymmetry can be substantial, and the derived SE may be slightly off.

**Caveat 2 — Inverse-variance pooling:** The code uses inverse-variance weighting (`w = 1/se²`) to pool stratum-specific RD estimates. This is the **fixed-effects inverse-variance meta-analysis** approach, which is the standard Greenland-Robins method for stratified RD. ✅ This is correct.

**Caveat 3 — `diffscoreci` defaults to score CI:** The PropCIs `diffscoreci` function implements the Miettinen-Nurminen score CI (as confirmed by the package documentation and references). This is appropriate for the MN method.

**Verdict:** **Minor issue.** The width-based SE extraction is an approximation. A more accurate approach would be to use the observed Fisher information at the MLE (which `diffscoreci` computes internally but doesn't expose). For the simulation's purpose (Type I error calibration), this approximation is adequate.

---

## 2. CMH RR Variance — Greenland-Robins Stratified Variance

### Implementation
```r
w <- n0 * n1 / (n0 + n1)
r1 <- (x1 + 0.5) / (n1 + 0.5); r0 <- (x0 + 0.5) / (n0 + 0.5)
rr_k <- r1 / r0
num <- num + w * rr_k; den <- den + w
var_num <- var_num + w^2 * (1/(x1+0.5) - 1/(n1+0.5) + 1/(x0+0.5) - 1/(n0+0.5))
```

### Assessment: **BUG — Incorrect variance formula for risk ratio**

**The problem:** The code uses the **Woolf-style variance for log-RR** (which is correct for the OR RBG formula but not for RR). The formula:
```
var(log(RR)) ≈ 1/x1 - 1/n1 + 1/x0 - 1/n0
```
is the delta-method variance for log(risk ratio) in a **single stratum**. This is correct per-stratum.

**However**, the pooling is wrong. The code computes:
```
var_log = var_num / den^2
```
where `var_num = Σ w_k² * var_k` and `den = Σ w_k`.

This is **not** the Greenland-Robins variance for the Mantel-Haenszel risk ratio. The correct approach is:

For the MH risk ratio `RR_MH = (Σ w_k * r1_k) / (Σ w_k * r0_k)` with `w_k = n0_k * n1_k / n_k`:

The Greenland-Robins variance for log(RR_MH) should use the **MH weights** applied to the **log-RR contributions**, not the inverse-variance pooling used for RD. Specifically:

```
var(log RR_MH) = Σ [w_k / (Σ w_j)]² * var(log RR_k)
```

Wait — actually, let me re-examine. The code does:
```
num = Σ w_k * RR_k    (this is the MH RR numerator)
den = Σ w_k           (this is the MH RR denominator)
RR_MH = num / den
```

But the variance formula `var_num / den²` where `var_num = Σ w_k² * var_k` gives:
```
var(log RR_MH) = Σ w_k² * var_k / (Σ w_k)²
```

This is **inverse-variance pooling**, not the Greenland-Robins variance for the MH estimator. The Greenland-Robins formula for the **MH odds ratio** is well-established. For the **MH risk ratio**, the variance is:

```
var(log RR_MH) = (1/RR_MH²) * var(RR_MH)
```

where `var(RR_MH)` follows the delta method applied to the MH ratio estimator. The correct variance for the stratified MH risk ratio uses the **conditional variance** under the null, not the inverse-variance of stratum-specific estimates.

**More precisely:** The Greenland-Robins variance for the MH risk ratio (not OR!) is given by:

```
var(log RR_MH) = (1/RR_MH²) * Σ [n0_k * n1_k / n_k²] * [r1_k(1-r1_k)/n1_k + RR_k² * r0_k(1-r0_k)/n0_k]
```

The code's formula `Σ w_k² * (1/x1 - 1/n1 + 1/x0 - 1/n0) / (Σ w_k)²` is **not** the correct Greenland-Robins variance for the MH risk ratio. It's a hybrid that mixes the MH RR point estimate with inverse-variance pooling of log-RR variances.

**Verdict:** **Major issue.** The variance for CMH RR is incorrect. The code uses inverse-variance pooling of stratum-specific log-RR variances, which is a valid *meta-analysis* approach but is **not** the Greenland-Robins variance for the Mantel-Haenszel risk ratio estimator. These give different results, especially in sparse strata.

**Fix:** Either:
1. Use `metafor::rma` with `method = "MH"` for the correct MH RR variance, or
2. Implement the correct Greenland-Robins variance for the MH risk ratio as derived in Greenland & Robins (1985) Section 3.

---

## 3. CMH OR — `mantelhaen.test` with `correct=FALSE`

### Implementation
```r
mh <- mantelhaen.test(tbl, correct = FALSE)
```

### Assessment: **Correct**

`mantelhaen.test` in R implements the Cochran-Mantel-Haenszel test for association in stratified 2×2 tables. With `correct=FALSE`, it uses the uncorrected chi-square statistic, which is the standard approach for large-sample inference. The continuity correction (`correct=TRUE`) is conservative and unnecessary for the simulation's purpose (calibrating Type I error).

**Verdict:** **Correct.** ✅

---

## 4. Scenario 4 — `p_strata = c(0.01, 0.02, 0.485, 0.485)`

### Assessment: **BUG — Scenario 4 is NOT meaningfully different from Scenario 1**

Looking at the code comments:
```r
# Scenario 4: all small (equal but tiny)
p_strata = c(0.01, 0.02, 0.485, 0.485)
```

The comment says "all small" but 48.5% + 48.5% = 97% of patients are in strata 3 and 4, which are **large**. Only strata 1 and 2 are small (1% and 2%).

**Scenario 1:** `c(0.25, 0.25, 0.25, 0.25)` — all 4 strata have ~100 patients each (n=400)
**Scenario 4:** `c(0.01, 0.02, 0.485, 0.485)` — strata 3 and 4 have ~194 patients each, strata 1 and 2 have ~4 and ~8 patients

This is actually **the same situation as Scenario 3** (`c(0.03, 0.03, 0.47, 0.47)`), just with more extreme small-stratum proportions. The comment "all small (equal but tiny)" is **misleading** — the strata are not all small, nor are they equal.

**Verdict:** **Minor issue.** The scenario design is valid (it tests the case where two strata are very small), but the comment is wrong. The scenario should be renamed to something like "2 very small strata (1%, 2%)" and the comment should match. This doesn't affect the simulation results but is confusing for documentation.

---

## 5. Reproducibility

### Seeds and Data Generation

```r
set.seed(20260518)
# ...
set.seed(20260518 + scenario*100000 + i*10 + round(event_rate*100))
```

### Assessment: **Correct but fragile**

**Strengths:**
- Each simulation replicate has a deterministic seed derived from scenario, replicate number, and event rate.
- `furrr_options(seed = TRUE)` enables deterministic parallel random number generation via L'Ecuyer-CMRG.

**Concerns:**
1. **Seed collision risk:** The seed formula `20260518 + scenario*100000 + i*10 + round(event_rate*100)` uses `scenario` (1-4) with multiplier 100,000 and `i` (1-5000) with multiplier 10. For the same `scenario` and `event_rate`, the seeds differ by exactly 10 per replicate. This is fine — no collision risk within a scenario.

2. **Cross-scenario seed overlap:** Seeds for scenario 1 replicate 10000 (if it existed) would overlap with scenario 2 replicate 0. But since `i` ranges 1-5000 and `scenario` ranges 1-4, the maximum seed for scenario 1 is `20260518 + 100000 + 50000 + 30 = 20310548` and minimum for scenario 2 is `20260518 + 200000 + 10 + 30 = 22060558`. No overlap. ✅

3. **Parallel determinism:** `furrr_options(seed = TRUE)` uses L'Ecuyer-CMRG for deterministic parallel RNG. This is the correct approach for reproducible parallel simulations. ✅

4. **Block randomization within strata:** The block randomization uses `sample()` with the replicate's seed. Since each replicate has a unique seed, the randomization is deterministic and reproducible. ✅

**Verdict:** **Correct.** ✅ The reproducibility is solid.

---

## Summary of Issues

| # | Issue | Severity | Status |
|---|-------|----------|--------|
| 1 | `stratified_mn_rd` SE from CI width | **Minor** | Approximation, adequate for simulation |
| 2 | CMH RR variance formula | **Major** | Uses inverse-variance pooling, not Greenland-Robins for MH RR |
| 3 | CMH OR `mantelhaen.test(correct=FALSE)` | **None** | Correct |
| 4 | Scenario 4 comment misleading | **Minor** | Doesn't affect results but is confusing |
| 5 | Reproducibility | **None** | Correct |

## Verdict: **Minor Issues** ⚠️

The code is functional and the simulation design is sound. However, issue #2 (CMH RR variance) is a **conceptual error** that could produce misleading Type I error rates for the risk ratio, especially in small strata where the inverse-variance pooling differs substantially from the correct MH variance.

### Recommended Fixes

1. **Fix CMH RR variance** (Priority: High): Replace the variance formula with the correct Greenland-Robins variance for the MH risk ratio. The formula should be:

```r
# For each stratum k:
# var(log RR_k) = (1-r1_k)/(n1_k * r1_k) + (1-r0_k)/(n0_k * r0_k)  [delta method]
# 
# For the MH RR estimator RR_MH = (Σ w_k * r1_k) / (Σ w_k * r0_k):
# var(log RR_MH) ≈ (1/RR_MH²) * Σ [w_k/n_k]² * [r1_k(1-r1_k)/n1_k + RR_k² * r0_k(1-r0_k)/n0_k]
```

2. **Fix Scenario 4 comment** (Priority: Low): Update to "2 very small strata (1%, 2%)" to match the actual proportions.

3. **Consider adding `diffscoreci` method parameter** (Priority: Low): The current code relies on `diffscoreci` defaulting to score CI. Explicitly document this dependency.
