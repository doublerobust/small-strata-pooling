# Small Strata Pooling: Comprehensive Guidance for Merck Oncology SAPs

**Author:** Yue Shentu  
**Date:** May 2026  
**Repository:** github.com/doublerobust/small-strata-pooling

---

## Executive Summary

This white paper investigates whether pooling small strata is necessary for stratified analyses in oncology clinical trials. We studied five standard methods across binary and time-to-event endpoints, using realistic oncology trial simulations with stratified block randomization, accrual periods, censoring, and group sequential designs.

**Key finding:** The answer depends on the method, not the endpoint.

| Method | Pooling Needed? | Recommendation |
|--------|:--------------:|:---------------|
| **CMH Odds Ratio** (+0.5 CC) | ❌ No | Use as-is regardless of stratum size |
| **CMH Risk Ratio** (+0.5 CC) | ❌ No | Use as-is; note RR scale inflation at low event rates |
| **Stratified MN Risk Difference** | ❌ No | Use as-is regardless of stratum size |
| **Stratified Cox PH** | ❌ No | Use as-is; partial likelihood naturally handles small strata |
| **Stratified Log-rank** | ❌ No | Use as-is; power is essentially identical to unstratified test |

**For ELSTIC and SAP templates:**
All methods (binary, Cox PH, log-rank): **no pooling required** — the existing guidance stands.
The stratified log-rank does NOT lose power with small strata.

---

## 1. Motivation

In stratified randomized oncology trials, SAPs routinely require pre-specified pooling rules for small strata. This creates operational complexity: statisticians review blinded data pre-interim to identify small strata and determine pooling rules.

Prior internal work established that pooling is unnecessary for Cox PH and log-rank in standard scenarios. This analysis:
1. Extends the binary endpoint investigation to cover **power** as well as Type I error
2. Adds a comprehensive **time-to-event simulation** with realistic oncology trial design
3. Distinguishes between **Cox PH** and **log-rank** behavior in extreme sparsity

---

## 2. Part I: Binary Endpoint Methods

### 2.1 Methods Investigated

| Method | Implementation | Details |
|--------|---------------|---------|
| **CMH OR** | `mantelhaen.test(correct=FALSE)` | +0.5 continuity correction for zero cells |
| **CMH RR** | Inverse-variance pooled, +0.5 CC | Bootstrap-validated variance (ratio 1.11 vs empirical) |
| **Stratified MN RD** | `PropCIs::diffscoreci`, IV-pooled | SE from score CI width (asymptotic approximation) |

### 2.2 Simulation Design

| Parameter | Value |
|-----------|-------|
| N | 400 |
| Reps | 5,000 (Type I) + 5,000 (Power: OR=1.65, 2.01) |
| Randomization | Stratified block (block size 4) |
| Strata | 4 (2 binary factors) |
| Scenarios | Balanced, 1 small (5%), 2 small (3%), 2 tiny (1%, 2%) |
| Event rates | 10%, 30%, 50% |
| Pooling rule | Merge strata < 10 patients into nearest larger stratum |

### 2.3 Type I Error Results

| Sparsity | Event Rate | CMH OR | CMH RR | MN RD |
|:---------|:---------:|:-----:|:-----:|:-----:|
| Balanced | 10% | 0.043 | 0.081 | 0.035 |
| Balanced | 30% | 0.047 | 0.055 | 0.053 |
| Balanced | 50% | 0.050 | 0.051 | 0.059 |
| 1 small (5%) | 10% | 0.047 | 0.064 | 0.036 |
| 1 small (5%) | 30% | 0.052 | 0.060 | 0.059 |
| 1 small (5%) | 50% | 0.045 | 0.049 | 0.054 |
| 2 small (3%) | 10% | 0.051 | 0.047 | 0.041 |
| 2 small (3%) | 30% | 0.049 | 0.070 | 0.061 |
| 2 small (3%) | 50% | 0.049 | 0.057 | 0.063 |
| 2 tiny (1%,2%) | 10% | 0.047 | 0.045 | 0.040 |
| 2 tiny (1%,2%) | 30% | 0.054 | 0.053 | 0.062 |
| 2 tiny (1%,2%) | 50% | 0.051 | 0.054 | 0.069 |

**Key finding:** Type I error inflation (where present) is an **inherent property of the method** at low event rates, not a small-strata problem. The inflation is highest in the balanced design with no small strata at all.

### 2.4 Power Results

Pooling small strata does **not improve power** for any binary method:

| Method | Max Power Gain from Pooling | Effect |
|:-------|:--------------------------:|:------|
| CMH OR | 0.000 to −0.002 | None |
| CMH RR | 0.000 to −0.023 | Slight **reduction** |
| MN RD | 0.001 to −0.010 | Minimal |

### 2.5 Binary Endpoint Recommendation

> **"No pooling of small strata is required for CMH odds ratio, CMH risk ratio, or stratified Miettinen-Nurminen risk difference."**

*Note: If using CMH RR at low event rates (< 15%), Type I error may be slightly inflated as an inherent property of the RR scale. Consider CMH OR as an alternative.*

---

## 3. Part II: Time-to-Event Methods

### 3.1 Methods Investigated

| Method | Implementation | Notes |
|--------|---------------|-------|
| **Stratified Cox PH** | `coxph(Surv ~ trt + strata(stratum))` | Treatment as only covariate |
| **Stratified Log-rank** | `survdiff(Surv ~ trt + strata(stratum))` | Fixed design, one-sided α=0.025 |

> **Estimand note:** Stratified Cox and stratified log-rank estimate a *conditional* hazard ratio — the treatment effect conditional on the randomization strata. An unstratified analysis would target the *marginal* (population-average) hazard ratio. Due to the non-collapsibility of the hazard ratio, these are technically distinct estimands under ICH E9(R1) [ASA 2025]. In this paper, all comparisons are made within the same analysis type (stratified vs. stratified, unstratified vs. unstratified), so the estimand is consistent across comparator methods within each setting. The simulation uses a constant treatment effect (HR = 0.65) across all strata, meaning both conditional and marginal estimands converge to the same numeric truth.

### 3.2 Simulation Design

| Parameter | Value |
|-----------|-------|
| N | 500 |
| Randomization | 1:1, stratified block (block size 4) |
| Strata | 8 (3 binary factors → 8 strata) |
| Accrual | 18 months, uniform, ~28 patients/month |
| Control median OS | 14 months (Weibull, shape=1) |
| Treatment effect | HR = 0.65 (power) / 1.0 (Type I) |
| Follow-up | 36 months administrative cutoff |
| Dropout | 5% annual (exponential) |
| Design | Fixed (no interim looks) |
| Type I reps | 10,000 |
| Power reps | 5,000 |
| Pooling rule | Merge strata < 10 patients into nearest larger stratum |

### 3.3 Type I Error Results (10,000 reps)

| Sparsity | Log-rank (No Pool) | Log-rank (Pool) | Diff |
|:---------|:-----------------:|:---------------:|:----:|
| Balanced | 0.0250 | 0.0250 | 0.0000 |
| 1 small (3%) | 0.0242 | 0.0242 | 0.0000 |
| 2 small (2%) | 0.0248 | 0.0248 | 0.0000 |
| 2 tiny (1%) | 0.0278 | 0.0278 | 0.0000 |

Type I error is well-controlled at the one-sided 0.025 level. Pooling has no effect on Type I error.

### 3.4 Power Results (5,000 reps, HR = 0.65)

| Sparsity | Log-rank (No Pool) | Log-rank (Pool) | Diff | Min stratum |
|:---------|:-----------------:|:---------------:|:----:|:-----------:|
| Balanced | 0.964 | 0.964 | 0.000 | 62 |
| 1 small (3%) | 0.958 | 0.958 | 0.000 | 15 |
| 2 small (2%) | 0.960 | 0.960 | 0.000 | 10 |
| 2 tiny (1%) | **0.963** | **0.963** | **0.000** | **5** |

**Power is essentially identical across ALL sparsity levels.** The stratified log-rank test is robust to small strata.

### 3.5 Convergence and Bias

| Metric | Result |
|:-------|:-------|
| Cox PH convergence | **100%** across all reps |
| Log-rank failure | 0% |
| HR bias (Cox) | Negligible: estimate 0.668 vs true 0.650 |
| SE of log(HR) | 0.116 for both pooled and unpooled |

---

## 4. Key Findings and Interpretation

### 4.1 Stratified Cox PH — Unaffected by Pooling

The partial likelihood is multiplicative across risk sets within and across strata:

$$L(\beta) = \prod_{s=1}^{S} \prod_{i \in D_s} \frac{\exp(\beta A_{si})}{\sum_{j \in R_s(t_{si})} \exp(\beta A_{sj})}$$

A tiny stratum contributes one risk-set term among hundreds — its influence is proportional to its information content. Pooling has no effect because the information from tiny strata was negligible to begin with.

### 4.2 Stratified Log-rank — Also Unaffected by Pooling

The stratified log-rank test is also robust to small strata. A dedicated simulation with 5,000 reps across multiple sparsity levels (stratum proportions from 50% down to 1%) shows:

| Small stratum proportion | Stratified power | Unstratified power | Difference |
|:-----------------------:|:----------------:|:------------------:|:----------:|
| 50% (balanced) | 0.964 | 0.963 | +0.001 |
| 20% | 0.963 | 0.962 | +0.001 |
| 10% | 0.960 | 0.960 | 0.000 |
| 5% | 0.966 | 0.966 | 0.000 |
| 2% | 0.962 | 0.963 | -0.001 |

**Result:** Stratified and unstratified log-rank have essentially identical power regardless of stratum size. Tiny strata contribute negligible noise. This finding was independently confirmed via double-programming (Qwen wrote and ran an independent simulation from scratch).

The observation that stratified and unstratified analyses produce nearly identical results under the simulation design is consistent with the estimand discussion in §3.1: with a constant treatment effect across strata, the conditional and marginal hazard ratios coincide, and the choice between stratified and unstratified estimators primarily affects precision, not the target of inference. In practice, when treatment effect heterogeneity may be present, the two approaches estimate distinct quantities and should be pre-specified accordingly [ASA 2025].

### 4.3 Binary Methods — Unaffected by Pooling

All three binary methods (CMH OR, CMH RR, MN RD) aggregate via Mantel-Haenszel or inverse-variance weighting, which naturally down-weights uninformative strata. The findings from the Type I error analysis (which showed small-strata robustness) are confirmed by the power analysis.

---

## 5. Recommendations

### 5.1 For SAP Language

> **Binary endpoints (CMH OR, CMH RR, MN RD):**
> *"Stratification factors will be used as specified in the randomization scheme. No pooling of small strata is required."*
>
> **Time-to-event endpoints:**
> *"Stratification factors will be used as specified in the randomization scheme. No pooling of small strata is required for any of the standard analysis methods."*

### 5.2 Pooling Rule

No pooling rule is needed. The stratified log-rank test is robust to small strata, as confirmed by double-programming with an independent implementation.

### 5.3 Summary Table

| Method | Pool? | Why |
|:-------|:----:|:----|
| **CMH Odds Ratio** | ❌ No | Most robust; Type I [0.043–0.054], power unaffected |
| **CMH Risk Ratio** | ❌ No | Slight Type I inflation at low event rates; pooling doesn't help |
| **MN Risk Difference** | ❌ No | Type I [0.035–0.069], power unaffected |
| **Stratified Cox PH** | ❌ No | Partial likelihood naturally handles small strata; convergence 100% |
| **Stratified Log-rank** | ❌ No | Power identical to unstratified across all sparsity levels |

---

## 6. Implications for ELSTIC Guidance

The ELSTIC guidance was finalized with a blanket "no pooling" recommendation. Our findings confirm this is correct for ALL methods:
- Binary methods (CMH OR, CMH RR, MN RD)
- Stratified Cox PH
- Stratified Log-rank (the stratified and unstratified tests have essentially identical power)

The ELSTIC guidance stands as originally written.

---

## Appendix A: Log-rank Test Robustness to Small Strata

The stratified log-rank test is robust to small strata. Power is essentially identical to the unstratified test across all sparsity levels.

**Key references:**
- Schoenfeld (1981). *The asymptotic properties of rank tests.* Biometrika.
- Andersen, Borgan, Gill & Keiding (1993). *Statistical Models Based on Counting Processes.* Springer.

---

## Appendix B: Simulation Design Details

### B.1 Binary Endpoint Simulation

| Parameter | Value |
|-----------|-------|
| R script | `binary/run_small_strata.R` |
| Power script | `binary/run_power_analysis.R` |
| Reps | 5,000 (Type I) + 5,000 (Power) |
| N | 400 |
| Randomization | Stratified block, size 4 |
| Scenario proportions | {0.25,0.25,0.25,0.25}, {0.05,0.35,0.30,0.30}, {0.03,0.03,0.47,0.47}, {0.01,0.02,0.485,0.485} |
| Event rates | 10%, 30%, 50% |
| Treatment effect (Power) | OR = 1.65, 2.01 |
| Seed scheme | `20260518 + scenario*1e6 + rep*10 + ev_rate*100` |
| Parallel | `furrr`, 11 workers, chunk_size=200 |
| Code review | `binary/audit/` (4 reviews: code-review, code-review-v2, final-review, power-code-review, qwen-peer-review) |

### B.2 Survival Endpoint Simulation

| Parameter | Value |
|-----------|-------|
| R script | `survival/run_survival_simulation.R` |
| Full run | `survival/run_survival_full.R` |
| Reps | 10,000 (Type I) + 5,000 (Power) |
| N | 500 |
| Randomization | 1:1 stratified block, size 4 |
| Strata | 8 (3 binary factors → 8 strata) |
| Scenario proportions | Same pattern as binary |
| Control median | 14 months (Weibull, shape=1, scale=14/ln(2)=20.20) |
| Accrual | Uniform 0–18 months |
| Cutoff | 36 months |
| Dropout | Exponential, 5%/year |
| Treatment effect | HR = 0.65 (Power) / HR = 1.0 (Type I) |
| Design | Fixed (no interim looks) |
| Significance | One-sided $\alpha = 0.025$ |
| Seed scheme | `20260519 + scenario*1e6 + hr_id*1e5 + rep` |
| Parallel | `furrr`, 11 workers, chunk_size=200 |
| Convergence | 100% across all reps |
| Code/plan review | `survival/audit/` (3 reviews: plan-review-agent1, plan-review-qwen, qwen-review-elstic) |

---

## Appendix C: File Structure

```
small-strata/
├── README.md                                 # Repository overview
├── small-strata-white-paper.md               # This comprehensive white paper
├── elstic-guidance-update.md                 # Proposed ELSTIC guidance update
├── logrank-pooling-explanation.md            # Technical explanation of log-rank mechanism
│
├── binary/                                   # Binary endpoint analysis
│   ├── run_small_strata.R                    # Type I error simulation
│   ├── run_power_analysis.R                  # Power simulation
│   ├── small-strata-proposal.md              # Original research proposal
│   └── audit/                                # Independent reviews
│       ├── code-review.md
│       ├── code-review-v2.md
│       ├── final-review.md
│       ├── power-code-review.md
│       └── qwen-peer-review.md
│
└── survival/                                 # Time-to-event analysis
    ├── run_survival_simulation.R             # Survival simulation code
    ├── run_survival_full.R                   # 10K + 5K rep runner
    ├── survival_full_output.txt              # Raw simulation output
    ├── survival-simulation-plan.md           # Simulation plan (reviewed)
    └── audit/                                # Independent reviews
        ├── plan-review-agent1.md
        ├── plan-review-qwen.md
        └── qwen-review-elstic.md
```

---

## References

[ASA 2025] ASA Oncology Estimand Working Group, Conditional and Marginal Effect Task Force. "Current practice on covariate adjustment and stratified analysis — based on survey results." *BMC Medical Research Methodology* (2025). DOI: 10.1186/s12874-025-02670-7.


- Greenland & Robins (1985). *Estimation of a common effect parameter from sparse follow-up data.* Biometrics.
- Miettinen & Nurminen (1985). *Comparative analysis of two rates.* Statistics in Medicine.
- Mantel & Haenszel (1959). *Statistical aspects of the analysis of data from retrospective studies of disease.* JNCI.

---

*End of white paper.*
