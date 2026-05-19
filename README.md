# Small Strata Pooling Investigation

**Author:** Yue Shentu  
**Repository:** github.com/doublerobust/small-strata-pooling

Comprehensive investigation of whether pooling small strata is necessary for stratified analyses in Merck oncology SAPs.

## Key Findings

| Method | Pooling Needed? |
|--------|:--------------:|
| CMH OR, CMH RR, MN RD (binary) | ❌ No |
| Stratified Cox PH (time-to-event) | ❌ No |
| Stratified Log-rank (time-to-event) | ⚠️ Pool strata < 10 patients |

## Documents

| File | Description |
|:-----|:------------|
| `small-strata-white-paper.md` | Comprehensive white paper with all findings and recommendations |
| `elstic-guidance-update.md` | Proposed ELSTIC guidance update |
| `logrank-pooling-explanation.md` | Technical explanation of why log-rank is affected |
| `binary/` | Binary endpoint simulations + audit trail |
| `survival/` | Time-to-event simulations + audit trail |
| `README.md` | This file |

## Reproducing Results

### Binary
```r
source("binary/run_small_strata.R")      # Type I error (~2 min)
source("binary/run_power_analysis.R")     # Power analysis (~3 min)
```

### Survival
```r
source("survival/run_survival_full.R")    # Full 10K+5K reps (~15 min with 11 workers)
```
