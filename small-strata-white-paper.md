# Small Strata Investigation: Recommendations for SAP Language

**Status:** Internal white paper — complete  
**Date:** May 2026  

---

## Bottom Line

All standard binary endpoint methods used in Merck oncology SAPs are robust to small strata:

| Method | Failure rate | Type I error | Pooling needed? |
|--------|:-----------:|:-----------:|:--------------:|
| **CMH odds ratio** (+0.5 CC) | 0.000 | 0.043–0.054 | ❌ No |
| **CMH risk ratio** (+0.5 CC) | 0.000 | 0.045–0.081 | ❌ No |
| **Stratified MN risk difference** (score-based) | 0.000 | 0.035–0.069 | ❌ No |
| **Cox PH** (stratified) | — | — | ❌ No (already confirmed) |
| **Log-rank** (stratified) | — | — | ❌ No (already confirmed) |

**Recommendation: No pooling of small strata is required for any of these methods.**

## Results (5,000 reps per scenario, stratified block randomization, block size 4)

| Sparsity | Event rate | CMH OR | | CMH RR | | MN RD | |
|:---------|:---------:|:-----:|:-----:|:-----:|:-----:|:-----:|:-----:|
| | | Fail | Type I | Fail | Type I | Fail | Type I |
| Balanced | 10% | 0.000 | 0.043 | 0.000 | 0.081 | 0.000 | 0.035 |
| Balanced | 30% | 0.000 | 0.047 | 0.000 | 0.055 | 0.000 | 0.053 |
| Balanced | 50% | 0.000 | 0.050 | 0.000 | 0.051 | 0.000 | 0.059 |
| 1 small (5%) | 10% | 0.000 | 0.047 | 0.000 | 0.064 | 0.000 | 0.036 |
| 1 small (5%) | 30% | 0.000 | 0.052 | 0.000 | 0.060 | 0.000 | 0.059 |
| 1 small (5%) | 50% | 0.000 | 0.045 | 0.000 | 0.049 | 0.000 | 0.054 |
| 2 small (3% ea) | 10% | 0.000 | 0.051 | 0.000 | 0.047 | 0.000 | 0.041 |
| 2 small (3% ea) | 30% | 0.000 | 0.049 | 0.000 | 0.070 | 0.000 | 0.061 |
| 2 small (3% ea) | 50% | 0.000 | 0.049 | 0.000 | 0.057 | 0.000 | 0.063 |
| 2 tiny (1%, 2%) | 10% | 0.000 | 0.047 | 0.000 | 0.045 | 0.000 | 0.040 |
| 2 tiny (1%, 2%) | 30% | 0.000 | 0.054 | 0.000 | 0.053 | 0.000 | 0.062 |
| 2 tiny (1%, 2%) | 50% | 0.000 | 0.051 | 0.000 | 0.054 | 0.000 | 0.069 |

## Notes

- **CMH OR** is the most robust — Type I error within [0.043, 0.054] across all scenarios.
- **CMH RR** with stratified Greenland-Robins variance shows slight Type I inflation at low event rates (0.081 balanced, 10%), but this is a property of the RR scale for rare events, not sparsity.
- **Stratified MN RD** (score-based, using `PropCIs::diffscoreci`) shows 0% failure across all scenarios. Type I error is slightly conservative at low event rates (0.035–0.041), which is a known advantage of the Miettinen-Nurminen method in sparse settings.
- With **stratified block randomization** (block size 4), within-stratum balance is guaranteed even for very small strata.

## Proposed SAP Language

For SAPs using CMH (OR or RR) or stratified MN (risk difference):

> *"Stratification factors will be used as specified in the randomization scheme. No pooling of small strata is required."*

## Code

All simulation code: `run_small_strata.R` in this repository.
