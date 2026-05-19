# Code Review: `run_small_strata.R`

**Review Date:** 2026-05-18  
**Reviewer:** Subagent (code review)  
**Verdict: Minor Issues**

---

## 1. Correctness

### 1.1 CMH Odds Ratio ✅ (with caveat)
- Uses `mantelhaen.test(tbl, correct = FALSE)` — correct for uncorrected CMH OR test.
- Builds the 2×2×K contingency table manually, which is the correct input format for `mantelhaen.test`.
- Handles empty strata by dropping them before calling `mantelhaen.test`.
- **Caveat:** `mantelhaen.test` returns a CMH *OR* by default when the input is a 3D array. This is correct.

### 1.2 CMH Risk Ratio ⚠️ (major concern)
- The MH-weighted risk ratio formula uses continuity correction (`+0.5` to both numerator and denominator counts), which is non-standard for a pure MH estimator. The classic Mantel-Haenszel RR does **not** use continuity correction. The `+0.5` correction is typically used for *variance estimation* or to handle zero cells, not for the point estimate itself.
- The variance formula uses the **total** (unstratified) event counts and totals:
  ```r
  se_log <- sqrt(1/tot_x1 - 1/tot_n1 + 1/tot_x0 - 1/tot_n0)
  ```
  This is the **unstratified** (pooled) variance, not a stratified MH variance. For a properly stratified RR, the variance should be computed per-stratum and then combined. This is a significant methodological concern — the standard error is not correctly stratified.
- The MH weight `w = n1*n0 / (n0 + n1)` is correct for the MH RR estimator.
- The continuity-corrected cell proportions `r1 = (x1+0.5)/(n1+0.5)` are unusual. Standard MH RR uses `r1 = x1/n1`.

### 1.3 Stratified MN Risk Difference ⚠️ (major concern)
- The function `miettinen_nurminen_rd` is labeled as a Miettinen-Nurminen score CI, but the implementation is a **Wald CI**, not a score-based CI:
  ```r
  se <- sqrt(p1*(1-p1)/n1 + p0*(1-p0)/n0)
  ```
  The true MN method solves a score equation iteratively for the CI bounds. A Wald approximation is much simpler and can have poor coverage, especially with small strata (which is the whole point of this simulation!).
- The "stratified MN" function (`stratified_mn_rd`) pools stratum-specific estimates using inverse-variance weighting. This is a **fixed-effects meta-analysis** approach, not the true stratified MN method which pools score statistics across strata. The two are related but not equivalent, especially when strata are very small.

### 1.4 Stratified Block Randomization ✅ (with minor note)
- Block size of 4, within each stratum, is correctly implemented.
- Each block randomly assigns half to treatment and half to control.
- **Note:** The last block may be smaller than 4 (handled via `ceiling(n_s/4)` and `min(b*4, n_s)`), and the `floor(bn/2)` treatment count is correct for even blocks but may give unequal allocation for odd-sized last blocks. This is standard and acceptable.

---

## 2. Bugs

### Bug 1: Scenario 4 p_strata is wrong ❌
```r
c(0.25, 0.25, 0.25, 0.25)  # all equal (small N)
```
The comment says "all small (equal but tiny)" but the probabilities are **identical to Scenario 1** (balanced). There is no "small" stratum in Scenario 4. This scenario should have all strata small (e.g., `c(0.01, 0.01, 0.01, 0.97)` or similar), or the comment is misleading. As written, Scenario 4 is a duplicate of Scenario 1 — this is a significant bug.

### Bug 2: `future_map` seed option is deprecated ⚠️
```r
.options = furrr_options(seed = TRUE, chunk_size = 200)
```
The `seed = TRUE` option in `furrr_options` was deprecated in favor of `.options = furrr_options(preschedule = TRUE)` with manual seeding. In newer versions of `furrr`, `seed = TRUE` may be ignored or produce warnings. This doesn't break reproducibility but may cause warnings.

### Bug 3: `mantelhaen.test` returns RR, not OR, in some R versions ⚠️
`mantelhaen.test` on a 3D array returns a **risk ratio** (not odds ratio) in recent R versions (≥ 4.0+). The CMH statistic (p-value) is the same for OR and RR, but if the code ever extracts `mh$estimate`, it would be a RR, not an OR. The current code only uses `mh$p.value`, so this doesn't affect the current results, but it's a latent correctness issue.

### Bug 4: No minimum cell count check for mantelhaen.test ⚠️
When a stratum has very few subjects (e.g., scenario 3 or 4), some 2×2 cells within a stratum may be empty (0 counts). `mantelhaen.test` can fail silently or produce unreliable results with zero cells. The code checks for zero treatment/control assignment but not for zero events/non-events within each cell.

---

## 3. Reproducibility

### ✅ Seeds are set properly
- Main seed: `set.seed(20260518)` before the simulation grid.
- Per-replication seed: `set.seed(20260518 + scenario*100000 + i*10 + round(event_rate*100))` — deterministic and unique per replication.
- `furrr_options(seed = TRUE)` provides parallel reproducibility (though see Bug 2 above).

### ✅ Data generation is sound
- Binary outcomes via `rbinom(n, 1, event_rate)` — correct for Bernoulli trials.
- Null hypothesis is maintained (same event rate in both arms).
- Stratification probabilities are clearly defined per scenario.

---

## 4. Interpretation

### Metrics match reported values ✅
- **Failure rate:** `mean(fail_or, na.rm=TRUE)` correctly computes the proportion of simulations where the method failed (returned NA).
- **Type I error:** `mean(p_or < 0.05, na.rm=TRUE)` correctly estimates the observed Type I error rate at α = 0.05.
- Both are computed with `na.rm=TRUE` to handle failed replications.

### ⚠️ Minor: Type I error threshold
The threshold is hardcoded as `0.05` but not parameterized. If the goal is to compare against a nominal α level, this should be a parameter.

---

## Summary of Issues

| # | Severity | Description |
|---|----------|-------------|
| 1 | **Major** | Scenario 4 `p_strata` is identical to Scenario 1 — not "all small" as commented |
| 2 | **Major** | CMH RR uses unstratified variance, not stratified |
| 3 | **Major** | MN RD is Wald CI, not true Miettinen-Nurminen score CI |
| 4 | **Minor** | Continuity correction on CMH RR point estimate is non-standard |
| 5 | **Minor** | `furrr_options(seed=TRUE)` may be deprecated |
| 6 | **Minor** | `mantelhaen.test` returns RR estimate, not OR (latent issue) |
| 7 | **Minor** | No check for zero cells in 2×2 stratum tables |

## Recommendations

1. **Fix Scenario 4** — change `p_strata` to actually create small strata (e.g., `c(0.01, 0.01, 0.01, 0.97)`).
2. **Fix CMH RR variance** — use a proper stratified variance estimator (e.g., the MH variance formula by Robins, Breslow, and Greenland 1986).
3. **Replace Wald MN with true MN** — implement the score-based iterative solver for proper MN CIs, or at least note the approximation in documentation.
4. **Remove continuity correction** from CMH RR point estimate unless specifically testing a corrected variant.
5. **Add zero-cell checks** for 2×2 tables within strata.

---

**Verdict: Minor Issues** — The core simulation logic is functional and reproducible, but the methodological implementations of CMH RR (unstratified variance) and MN RD (Wald instead of score CI) are significant approximations that may bias results, especially in the small-strata scenarios this simulation is designed to study. Scenario 4's identical probabilities to Scenario 1 is a clear bug that invalidates that scenario's results entirely.
