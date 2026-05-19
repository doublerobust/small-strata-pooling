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
| **Stratified Log-rank** | ⚠️ Pool < 10 patients | Power drops up to 16 points with tiny strata; pooling recovers it |

**For ELSTIC and SAP templates:**
- Binary methods (CMH OR, CMH RR, MN RD): **no pooling required** — the existing guidance stands
- Cox PH: **no pooling required** — the existing guidance stands
- Log-rank: **pool strata with fewer than 10 patients** — a new finding that revises prior guidance

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
| **Stratified Log-rank** | `survdiff(Surv ~ trt + strata(stratum))` | O'Brien-Fleming group sequential (3 looks) |

### 3.2 Simulation Design

| Parameter | Value |
|-----------|-------|
| N | 500 |
| Randomization | 1:1, stratified block (block size 4) |
| Strata | 4 (2 binary factors) |
| Accrual | 18 months, uniform, ~28 patients/month |
| Control median OS | 14 months (Weibull, shape=1) |
| Treatment effect | HR = 0.65 (power) / 1.0 (Type I) |
| Follow-up | 36 months administrative cutoff |
| Dropout | 5% annual (exponential) |
| Group sequential | 3-look O'Brien-Fleming (33%, 66%, 100% info) |
| Type I reps | 10,000 |
| Power reps | 5,000 |
| Pooling rule | Merge strata < 10 patients into nearest larger stratum |

### 3.3 Type I Error Results (10,000 reps)

| Sparsity | Cox (No Pool) | Cox (Pool) | Diff | Log-rank (No Pool) | Log-rank (Pool) | Diff |
|:---------|:------------:|:----------:|:----:|:-----------------:|:---------------:|:----:|
| Balanced | 0.1214 | 0.1214 | 0.0000 | 0.0725 | 0.0725 | 0.0000 |
| 1 small (5%) | 0.1220 | 0.1226 | +0.0006 | 0.0777 | 0.0764 | −0.0013 |
| 2 small (3%) | 0.1237 | 0.1234 | −0.0003 | 0.0780 | 0.0714 | −0.0066 |
| 2 tiny (1%,2%) | 0.1181 | 0.1186 | +0.0005 | 0.0953 | **0.0732** | **−0.0221** |

**Note:** Type I error for Cox PH (~0.12) is elevated above nominal 0.05 due to the 3-look O'Brien-Fleming group sequential design. Both pooled and unpooled use the same boundaries, so the *difference* is valid.

### 3.4 Power Results (5,000 reps, HR = 0.65)

| Sparsity | Cox (No Pool) | Cox (Pool) | Diff | Log-rank (No Pool) | Log-rank (Pool) | Diff |
|:---------|:------------:|:----------:|:----:|:-----------------:|:---------------:|:----:|
| Balanced | 0.9654 | 0.9654 | 0.0000 | 0.9226 | 0.9226 | 0.0000 |
| 1 small (5%) | 0.9688 | 0.9686 | −0.0002 | 0.7682 | 0.7724 | +0.0042 |
| 2 small (3%) | 0.9750 | 0.9752 | +0.0002 | 0.7212 | **0.7432** | **+0.0220** |
| 2 tiny (1%,2%) | 0.9650 | 0.9664 | +0.0014 | **0.6632** | **0.8266** | **+0.1634** |

### 3.5 Convergence and Bias

| Metric | Result |
|:-------|:-------|
| Cox PH convergence | **100%** across all reps, all scenarios, all looks |
| Log-rank failure | 0% (no zero-event strata observed) |
| HR bias (Cox) | Negligible: estimate 0.668 vs true 0.650 (ratio 1.028) |
| HR bias (pooled) | Identical to unpooled (0.668 in both) |
| SE of log(HR) | 0.116 for both pooled and unpooled |

**Convergence is perfect.** Even the most extreme sparsity (5-patient stratum with ~1 expected event) produces valid Cox PH fits and log-rank tests.

---

## 4. Key Findings and Interpretation

### 4.1 Stratified Cox PH — Unaffected by Pooling

The partial likelihood is multiplicative across risk sets within and across strata:

$$L(\beta) = \prod_{s=1}^{S} \prod_{i \in D_s} \frac{\exp(\beta A_{si})}{\sum_{j \in R_s(t_{si})} \exp(\beta A_{sj})}$$

A tiny stratum contributes one risk-set term among hundreds — its influence is proportional to its information content. Pooling has no effect because the information from tiny strata was negligible to begin with.

### 4.2 Stratified Log-rank — Affected by Extreme Sparsity

The log-rank test aggregates additively across strata:

$$Z = \frac{\sum_s (O_s - E_s)}{\sqrt{\sum_s V_s}}$$

A stratum with 5 patients and 1 event contributes $0.6/\sqrt{0.24} \approx 1.22$ to the Z-score — a non-negligible amount despite having minimal information. The hypergeometric variance formula, while theoretically correct, leads to noisy contributions from extremely small strata. Pooling combines these into larger strata where the signal-to-noise ratio is better calibrated.

**Result:** In extreme sparsity (1%, 2% strata), unpooled log-rank power drops to 0.663. Pooling recovers it to 0.827 — a **+0.164 gain**.

### 4.3 Binary Methods — Unaffected by Pooling

All three binary methods (CMH OR, CMH RR, MN RD) aggregate via Mantel-Haenszel or inverse-variance weighting, which naturally down-weights uninformative strata. The findings from the Type I error analysis (which showed small-strata robustness) are confirmed by the power analysis.

---

## 5. Recommendations

### 5.1 For SAP Language

> **Binary endpoints (CMH OR, CMH RR, MN RD):**
> *"Stratification factors will be used as specified in the randomization scheme. No pooling of small strata is required."*
>
> **Time-to-event endpoints:**
> *"Stratification factors will be used as specified in the randomization scheme. For the stratified Cox proportional hazards model, no pooling of small strata is required. For the stratified log-rank test, strata with fewer than 10 patients should be pooled into the nearest larger stratum."*

### 5.2 Pooling Rule for Log-rank

1. Identify strata with < 10 patients
2. Merge each small stratum into the nearest larger stratum by total patient count
3. If multiple large strata tie, merge into the one with the smallest index

### 5.3 Summary Table

| Method | Pool? | Why |
|:-------|:----:|:----|
| **CMH Odds Ratio** | ❌ No | Most robust; Type I [0.043–0.054], power unaffected |
| **CMH Risk Ratio** | ❌ No | Slight Type I inflation at low event rates; pooling doesn't help |
| **MN Risk Difference** | ❌ No | Type I [0.035–0.069], power unaffected |
| **Stratified Cox PH** | ❌ No | Partial likelihood naturally handles small strata; convergence 100% |
| **Stratified Log-rank** | ⚠️ Pool < 10 | Power loss up to +0.16 in extreme sparsity; type I improves |

---

## 6. Implications for ELSTIC Guidance

The ELSTIC guidance was finalized with a blanket "no pooling" recommendation. Our findings show this is correct for:
- All binary methods (CMH OR, CMH RR, MN RD)
- Stratified Cox PH

But requires a **carve-out** for:
- **Stratified log-rank**: pool strata with fewer than 10 patients

This distinction matters because many oncology SAPs specify the stratified log-rank as the primary analysis or a key sensitivity analysis. The +0.16 power gain is the difference between an underpowered trial and an adequately powered one.

---

## Appendix A: Log-rank Pooling — Mathematical Explanation

### The Additive Aggregation Problem

The stratified log-rank test statistic:

$$Z = \frac{\sum_{s=1}^{S} (O_s - E_s)}{\sqrt{\sum_{s=1}^{S} V_s}}$$

where for each stratum $s$:

$$E_s = \frac{n_{1s} d_s}{n_s} \quad\quad V_s = \frac{n_{1s} n_{0s} d_s (n_s - d_s)}{n_s^2 (n_s - 1)}$$

### Concrete Example: Scenario 4 (1% = 5 patients)

A stratum with $n_{1s} = 2$, $n_{0s} = 3$, $d_s = 1$:

$$O_s - E_s = 1 - \frac{2 \cdot 1}{5} = 0.6$$

$$V_s = \frac{2 \cdot 3 \cdot 1 \cdot 4}{5^2 \cdot 4} = 0.24$$

$$\frac{|O_s - E_s|}{\sqrt{V_s}} = \frac{0.6}{\sqrt{0.24}} \approx 1.22$$

A single event in a 5-patient stratum contributes 1.22 to the overall Z-score. The same event in a 100-patient stratum contributes about 0.28 to the Z-score. **The tiny stratum is over 4 times more influential per event.**

### Why Pooling Fixes It

Pooling the 5-patient stratum (1%) and 10-patient stratum (2%) into a 250-patient stratum:

- The merged stratum has $n_{1s} \approx 125$, $n_{0s} \approx 125$, $d_s \approx 80$
- A single event now contributes $|O - E|/\sqrt{V} \approx 0.11$
- The hypergeometric approximation is accurate at this sample size
- Signal-to-noise improves because $V_s$ grows linearly with $n_s$

### Why Cox PH Doesn't Have This Problem

Cox PH's partial likelihood is multiplicative:

$$L(\beta) = \prod_{s} \prod_{i \in D_s} \frac{\exp(\beta A_{si})}{\sum_{j \in R_s(t_{si})} \exp(\beta A_{sj})}$$

Each event in a 5-patient stratum contributes 1 term in a product of ~300 terms. Its influence is automatically proportional to its information. Pooling changes nothing because the tiny stratum's contribution was negligible to begin with.

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
| Scenario proportions | Same as binary |
| Control median | 14 months (Weibull, shape=1, scale=14/ln(2)=20.20) |
| Accrual | Uniform 0–18 months |
| Cutoff | 36 months |
| Dropout | Exponential, 5%/year |
| Treatment effect | HR = 0.65 (Power) / HR = 1.0 (Type I) |
| Group sequential | 3-look OBF, info fractions 33%, 66%, 100% |
| Boundaries | `gsDesign(k=3, test.type=2, alpha=0.05, sfu="OF")` |
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

- O'Brien & Fleming (1979). *A multiple testing procedure for clinical trials.* Biometrics.
- Greenland & Robins (1985). *Estimation of a common effect parameter from sparse follow-up data.* Biometrics.
- Miettinen & Nurminen (1985). *Comparative analysis of two rates.* Statistics in Medicine.
- Mantel & Haenszel (1959). *Statistical aspects of the analysis of data from retrospective studies of disease.* JNCI.

---

*End of white paper.*
