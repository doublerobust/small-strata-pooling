# Small Strata Investigation: Recommendations for SAP Language

**Status:** Internal white paper — complete  
**Date:** May 2026  

---

## Bottom Line

| Method | Pooling needed? | Notes |
|--------|:--------------:|-------|
| **Cox PH** (stratified) | ❌ No | Already confirmed — no issues |
| **Log-rank** (stratified) | ❌ No | Already confirmed — no issues |
| **CMH odds ratio** (+0.5 CC) | ❌ No | 0% failure, Type I error ~0.05 across all sparsity levels |
| **CMH risk ratio** (+0.5 CC) | ❌ No (with note) | 0% failure, but Type I error up to 0.097 at 10% event rate |
| **MN risk difference** | ✅ Yes | Fails 14-51% with small strata at low event rates |

## Findings (5,000 reps per scenario)

### CMH Odds Ratio — Fully Robust
| Sparsity | Event rate | Failure rate | Type I error |
|:---------|:----------:|:------------:|:-----------:|
| Balanced | 10% | 0.000 | 0.052 |
| Balanced | 30% | 0.000 | 0.052 |
| Balanced | 50% | 0.000 | 0.053 |
| 1 small stratum | 10% | 0.000 | 0.047 |
| 1 small stratum | 30% | 0.000 | 0.048 |
| 1 small stratum | 50% | 0.000 | 0.047 |
| 2 small strata | 10% | 0.000 | 0.056 |
| 2 small strata | 30% | 0.000 | 0.054 |
| 2 small strata | 50% | 0.000 | 0.053 |

**Recommendation:** No pooling required for CMH OR. The +0.5 continuity correction and Mantel-Haenszel weighting handle sparsity without issues. Standard SAP language is adequate.

### CMH Risk Ratio — Stable but Slight Type I Inflation
| Sparsity | Event rate | Failure rate | Type I error |
|:---------|:----------:|:------------:|:-----------:|
| Balanced | 10% | 0.000 | 0.097 |
| Balanced | 30% | 0.000 | 0.066 |
| Balanced | 50% | 0.000 | 0.055 |
| 1 small stratum | 10% | 0.000 | 0.070 |
| 1 small stratum | 30% | 0.000 | 0.063 |
| 1 small stratum | 50% | 0.000 | 0.052 |
| 2 small strata | 10% | 0.000 | 0.052 |
| 2 small strata | 30% | 0.000 | 0.064 |
| 2 small strata | 50% | 0.000 | 0.063 |

**Recommendation:** No pooling required for numerical stability (0% failure). The slight Type I inflation at low event rates is a property of the variance estimator, not sparsity. A more conservative variance estimator could be used if concern exists.

### Stratified MN Risk Difference — Fails with Small Strata + Low Events
| Sparsity | Event rate | Failure rate | Type I error |
|:---------|:----------:|:------------:|:-----------:|
| Balanced | 10% | 0.000 | 0.062 |
| Balanced | 30% | 0.000 | 0.059 |
| Balanced | 50% | 0.000 | 0.062 |
| 1 small stratum | 10% | **0.141** | 0.055 |
| 1 small stratum | 30% | 0.003 | 0.059 |
| 1 small stratum | 50% | 0.000 | 0.059 |
| 2 small strata | 10% | **0.514** | 0.063 |
| 2 small strata | 30% | **0.059** | 0.077 |
| 2 small strata | 50% | 0.020 | 0.079 |

**Recommendation:** Pooling IS needed for MN RD. When a stratum has very few patients (≤5 per arm) and event rates are low (≤30%), the inverse-variance weighted estimator fails 14-51% of the time. SAPs should specify pooling for any stratum with fewer than 5 patients per arm when MN RD is the primary method.

## Proposed SAP Language

For SAPs using **CMH** (OR or RR):
> *"Stratification factors will be used as specified in the randomization scheme. No pooling of small strata is required."*

For SAPs using **stratified MN** (risk difference):
> *"If any stratum has fewer than 5 patients in either treatment arm, that stratum will be pooled with the adjacent stratum according to the hierarchy specified in Table X."*

## Code

**Simulation details:** 5,000 reps per scenario using stratified block randomization (block size 4) to match real trial practice. 

All simulation code: `simulation/run_small_strata.R` in the ridge-cal repository.
