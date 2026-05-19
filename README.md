# Small Strata Pooling Investigation

**Author:** Yue Shentu  
**Repository:** github.com/doublerobust/small-strata-pooling

Comprehensive investigation of whether pooling small strata is necessary for stratified analyses in Merck oncology SAPs.

## Key Findings

**No pooling required for any standard analysis method.** Confirmed by double-programming with independent implementation.

| Method | Pooling Needed? |
|--------|:--------------:|
| CMH OR, CMH RR, MN RD (binary) | ❌ No |
| Stratified Cox PH (time-to-event) | ❌ No |
| Stratified Log-rank (time-to-event) | ❌ No |

## Repository Structure

```
small-strata/
├── README.md
├── small-strata-white-paper.md          # Comprehensive white paper
├── binary/                               # Binary endpoint analysis
│   ├── run_small_strata.R               # Type I error simulation
│   ├── run_power_analysis.R             # Power simulation
│   └── audit/                           # Independent reviews
└── survival/                            # Time-to-event analysis
    ├── run_survival_simulation.R        # Simulation code
    ├── audit/                           # Independent reviews
    └── survival_final_corrected.txt     # Final results
```

## Reproducing Results

### Binary
```r
source("binary/run_small_strata.R")      # Type I error (~2 min)
source("binary/run_power_analysis.R")     # Power analysis (~3 min)
```

### Survival
```r
source("survival/run_survival_full.R")    # Full 10K+5K reps (~15 min)
```
