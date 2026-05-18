# Small Strata Pooling: Internal Investigation

**Author:** Yue Shentu  
**Status:** Internal white paper — complete  
**Date:** May 2026  
**Repository:** github.com/doublerobust/small-strata-pooling

---

## Bottom Line

**No pooling of small strata is required** for any of the standard methods used in Merck oncology SAPs:

| Method | Failure Rate | Type I Error | Pooling Needed? |
|--------|:-----------:|:-----------:|:---------------:|
| CMH Odds Ratio (+0.5 CC) | 0.000 | 0.043–0.054 | ❌ No |
| CMH Risk Ratio (Greenland-Robins) | 0.000 | 0.045–0.081 | ❌ No |
| Stratified Miettinen-Nurminen Risk Diff | 0.000 | 0.035–0.069 | ❌ No |
| Cox PH (stratified) | — | — | ❌ No (confirmed) |
| Log-rank (stratified) | — | — | ❌ No (confirmed) |

## Motivation

In stratified randomized trials, SAPs routinely require pre-specified pooling rules for small strata. This creates operational complexity: statisticians review blinded data pre-interim to identify small strata and determine pooling. Internal investigation found this is unnecessary for Cox/log-rank. This investigation extends the same question to binary endpoint methods.

## Methods Investigated

- **CMH Odds Ratio** — Mantel-Haenszel estimator with +0.5 continuity correction
- **CMH Risk Ratio** — MH-weighted with Greenland-Robins stratified variance
- **Stratified Miettinen-Nurminen Risk Difference** — Score-based CIs via `PropCIs::diffscoreci`, inverse-variance pooled

## Simulation Design

- 5,000 reps per scenario
- Stratified block randomization (block size 4)
- N = 400, 2–4 stratification factors
- Event rates: 10%, 30%, 50%
- Sparsity: from balanced to extreme (1% stratum size)

## Key Finding

Type I error departures from nominal are **inherent to the methods**, not caused by small strata. Pooling would not address them.

## Repository Contents

| File | Description |
|------|-------------|
| `small-strata-white-paper.md` | Final white paper with SAP language |
| `small-strata-proposal.md` | Original research proposal |
| `run_small_strata.R` | R simulation code |
| `code-review.md` | Independent code review v1 |
| `code-review-v2.md` | Independent code review v2 |
| `final-review.md` | Final white paper + code review |

## Reproducing Results

Run in R with the `PropCIs` package installed:

```r
install.packages("PropCIs")
source("run_small_strata.R")
```

The simulation takes ~2-3 minutes with 11 parallel workers.

## Proposed SAP Language

For SAPs using CMH (OR or RR) or stratified MN (risk difference):

> *"Stratification factors will be used as specified in the randomization scheme. No pooling of small strata is required."*
