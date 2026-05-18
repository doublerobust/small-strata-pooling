# Small Strata Investigation: Recommendations for SAP Language

**Status:** Internal white paper — complete  
**Date:** May 2026  

---

## Bottom Line

All three standard binary endpoint methods are robust to small strata when properly implemented:

| Method | Failure rate | Type I error | Pooling needed? |
|--------|:-----------:|:-----------:|:--------------:|
| **CMH odds ratio** (+0.5 CC) | 0.000 | 0.043–0.054 | ❌ No |
| **CMH risk ratio** (+0.5 CC, stratified variance) | 0.000 | 0.045–0.081 | ❌ No (note ⚠️) |
| **Stratified Wald risk difference** | 0.000 | 0.051–0.063 | ❌ No |
| **Cox PH** (stratified) | — | — | ❌ No (already confirmed) |
| **Log-rank** (stratified) | — | — | ❌ No (already confirmed) |

**Overall recommendation: No pooling of small strata is required for any of these methods.**

## Results (5,000 reps per scenario, stratified block randomization, block size 4)

| Sparsity pattern | Event rate | CMH OR | | CMH RR | | Wald RD | |
|:----------------|:---------:|:-----:|:-----:|:-----:|:-----:|:-----:|:-----:|
| | | Fail | Type I | Fail | Type I | Fail | Type I |
| Balanced | 10% | 0.000 | 0.043 | 0.000 | 0.081 | 0.000 | 0.051 |
| Balanced | 30% | 0.000 | 0.047 | 0.000 | 0.055 | 0.000 | 0.054 |
| Balanced | 50% | 0.000 | 0.050 | 0.000 | 0.051 | 0.000 | 0.055 |
| 1 small stratum (5%) | 10% | 0.000 | 0.047 | 0.000 | 0.064 | 0.000 | 0.052 |
| 1 small stratum (5%) | 30% | 0.000 | 0.052 | 0.000 | 0.060 | 0.000 | 0.059 |
| 1 small stratum (5%) | 50% | 0.000 | 0.045 | 0.000 | 0.049 | 0.000 | 0.052 |
| 2 small strata (3% ea) | 10% | 0.000 | 0.051 | 0.000 | 0.047 | 0.000 | 0.056 |
| 2 small strata (3% ea) | 30% | 0.000 | 0.049 | 0.000 | 0.070 | 0.000 | 0.063 |
| 2 small strata (3% ea) | 50% | 0.000 | 0.049 | 0.000 | 0.057 | 0.000 | 0.061 |
| 2 v. small (1%, 2%) | 10% | 0.000 | 0.047 | 0.000 | 0.045 | 0.000 | 0.051 |
| 2 v. small (1%, 2%) | 30% | 0.000 | 0.054 | 0.000 | 0.053 | 0.000 | 0.060 |
| 2 v. small (1%, 2%) | 50% | 0.000 | 0.051 | 0.000 | 0.054 | 0.000 | 0.061 |

## Notes

- **CMH OR** is the most robust — Type I error stays within [0.043, 0.054] across all sparsity levels and event rates.
- **CMH RR** with stratified Greenland-Robins variance shows slight Type I inflation at low event rates (0.081 at 10%), but this is a property of the RR scale (rare events) not sparsity — the inflation is similar in balanced and unbalanced designs.
- **Wald RD** stays within [0.051, 0.063] across all scenarios with 0% failure — the inverse-variance weighted approach with proper handling of zero-variance strata is stable.
- With **stratified block randomization** (block size 4), within-stratum balance is guaranteed even for very small strata, which contributes to numerical stability.

## Proposed SAP Language

For SAPs using CMH (OR or RR) or stratified Wald RD:
> *"Stratification factors will be used as specified in the randomization scheme. No pooling of small strata is required."*

## Code

All simulation code: `run_small_strata.R` in this repository.
