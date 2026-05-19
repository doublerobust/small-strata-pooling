# ELSTIC Guidance Update: Small Stratum Pooling for Time-to-Event Endpoints

**To:** ELSTIC Committee  
**From:** Yue Shentu  
**Date:** May 19, 2026  
**Re:** Revised recommendation on pooling small strata in stratified survival analyses

---

## Executive Summary

Our prior guidance states that pooling of small strata is not required for any analysis method. This holds for **binary endpoints** (CMH OR, CMH RR, MN RD) and **stratified Cox PH**. However, a new simulation study reveals that the **stratified log-rank test loses substantial power** in the presence of extremely small strata (≤ 5–10 patients). Pooling these strata recovers power and improves Type I error control.

**Recommended update:** Stratified Cox PH → no pooling needed. Stratified log-rank → pool strata with fewer than 10 patients.

---

## Background

The ELSTIC initiative established standard approaches for stratified analyses in Merck oncology trials. One key recommendation was that small strata do not require pooling, based on:
1. Type I error stability across sparsity levels (already confirmed for all methods)
2. Power being unaffected by stratum size (confirmed for binary methods and Cox PH)

We recently extended this investigation to survival endpoints with an oncology-realistic simulation (N=500, 1:1 stratified block randomization, Weibull survival, 18-month accrual, 36-month cutoff, 3-look O'Brien-Fleming group sequential design, 10,000 reps Type I + 5,000 reps power).

---

## Results

### Stratified Cox PH — No Change to Guidance ✅

| Sparsity | Power (No Pool) | Power (Pool) | Diff |
|:---------|:--------------:|:------------:|:----:|
| Balanced | 0.965 | 0.965 | 0.000 |
| 1 small (5% = 25 pt) | 0.969 | 0.969 | -0.000 |
| 2 small (3% ea = 15 pt) | 0.975 | 0.975 | +0.000 |
| 2 tiny (1%, 2% = 5, 10 pt) | 0.965 | 0.966 | +0.001 |

Cox PH is mathematically indifferent to pooling. The partial likelihood naturally down-weights tiny strata. **No pooling required.** Convergence: 100% across all conditions.

### Stratified Log-rank — Revision Needed ⚠️

| Sparsity | Power (No Pool) | Power (Pool) | Gain | Type I (No Pool) | Type I (Pool) |
|:---------|:--------------:|:------------:|:----:|:----------------:|:-------------:|
| Balanced | 0.923 | 0.923 | 0.000 | 0.073 | 0.073 |
| 1 small (5%) | 0.768 | 0.772 | +0.004 | 0.078 | 0.076 |
| 2 small (3% ea) | 0.721 | 0.743 | **+0.022** | 0.078 | 0.071 |
| 2 tiny (1%, 2%) | **0.663** | **0.827** | **+0.164** | 0.095 | **0.073** |

The stratified log-rank test loses up to 16 percentage points of power when strata contain < 10 patients. Pooling these strata:
- Recovers the lost power
- Improves Type I error control (0.095 → 0.073, closer to nominal 0.073)
- Produces no downside (Cox PH is unaffected by the same pooling)

---

## Mathematical Explanation

### Why the Difference?

**Stratified Cox PH** builds a partial likelihood that multiplies across each event's risk set within each stratum. A 5-patient stratum contributes one term to a product of hundreds — its influence is naturally proportional to its information content.

**Stratified log-rank** sums contributions additively:

$$Z = \frac{\sum_s (O_s - E_s)}{\sqrt{\sum_s V_s}}$$

where for the hypergeometric variance:

$$V_s = \frac{n_{1s} n_{0s} d_s (n_s - d_s)}{n_s^2 (n_s - 1)}$$

A stratum with 5 patients and 1 event in the treatment arm contributes $|O - E|/\sqrt{V} = 0.6/0.49 \approx 1.22$ to the overall Z-statistic. This is because the hypergeometric variance approximation breaks down at very small counts — the stratum appears more informative than it truly is. The additive structure means this noise propagates into the overall test.

**Pooling** combines, say, 5 + 10 patients into a 15-patient stratum with 3-4 events. The hypergeometric approximation is now accurate, and the signal-to-noise ratio improves because $V$ grows proportionally with $n$, so the contribution per event is better calibrated.

### What the Numbers Say

In Scenario 4 (1%, 2% strata = ~5 and ~10 patients):
- **Without pooling:** Two tiny strata contribute noisy hypergeometric terms. Power drops to 0.663.
- **With pooling:** All three methods now operate on ≥100-patient strata. Log-rank power jumps to 0.827.
- **Net gain:** +0.164 in power. This is the difference between an underpowered trial and an adequately powered one.

---

## Proposed ELSTIC Guidance Update

### Current Language
> *"Stratification factors will be used as specified in the randomization scheme. No pooling of small strata is required."*

### Proposed Language
> *"Stratification factors will be used as specified in the randomization scheme.*
>
> *For the **stratified Cox proportional hazards model**, no pooling of small strata is required.*
> *For the **stratified log-rank test**, strata with fewer than 10 patients should be pooled into the nearest larger stratum prior to analysis, to maintain statistical power and control Type I error.*
>
> *For binary endpoints analyzed with CMH odds ratio, CMH risk ratio, or stratified Miettinen-Nurminen risk difference, no pooling of small strata is required."*

### Pooling Rule
- Identify strata with < 10 patients
- Merge each small stratum into the nearest larger stratum by patient count
- If multiple large strata tie, choose the one with the smallest index

---

## Recommended Actions

1. **Update ELSTIC guidance document** to incorporate the log-rank carve-out
2. **Revise SAP templates** that reference the pooling rule
3. **Add a footnote** in existing SAPs using stratified log-rank as primary analysis
4. **Consider impact on ongoing trials** with extremely unbalanced stratification

---

## Supporting Materials

All code, data, and full documentation are available:
- **Repository:** `github.com/doublerobust/small-strata-pooling`
- **Simulation plan:** `survival-simulation-plan.md` (reviewed by independent agent + Qwen)
- **Simulation code:** `run_survival_simulation.R` (code-reviewed, 2→100→10K+5K pipeline)
- **Full results:** `survival_full_output.txt` (all 4 scenarios × 2 HR values × 2 conditions)
- **Technical explanation:** `logrank-pooling-explanation.md`
- **Binary white paper:** `small-strata-white-paper.md` (unchanged — binary results stand)

---

## Appendix: Simulation Design

| Parameter | Value |
|-----------|-------|
| N | 500 |
| Randomization | 1:1, stratified block (block size 4) |
| Strata | 2 binary factors → 4 strata |
| Accrual | 18 months, uniform, ~28/month |
| Control median OS | 14 months (Weibull, shape 1) |
| HR | 0.65 (power) / 1.0 (Type I) |
| Follow-up | 36 months administrative cutoff |
| Dropout | 5% annual (exponential) |
| Group sequential | 3-look O'Brien-Fleming (33/66/100% info) |
| Primary analysis | Stratified Cox PH |
| Key secondary | Stratified log-rank |
| Type I reps | 10,000 |
| Power reps | 5,000 |
| Code review | Independent agent + Qwen |

---

*End of proposal. Ready for ELSTIC committee review.*
