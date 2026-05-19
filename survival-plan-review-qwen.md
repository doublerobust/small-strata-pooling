# Review: Small Strata Pooling Simulation Plan

**Reviewer:** Qwen  
**Date:** 2026-05-19  
**Verdict:** **Minor Revisions**

---

## 1. Trial Design Parameters

### ✅ Sound
- N = 500, 1:1 randomization with stratified blocks is standard.
- 2 binary stratification factors → 4 strata is the most common oncology scenario.
- Accrual 18 months, uniform, ~28/month is internally consistent (18 × 28 = 504 ≈ 500).
- Administrative cutoff at 36 months from start is a reasonable oncology follow-up window.

### ⚠️ Issues
- **Strata proportions need explicit definition.** The plan references "same 4 scenarios as binary study" but doesn't list the actual proportions. For reproducibility, these must be specified:
  - Balanced: (0.25, 0.25, 0.25, 0.25)
  - 1 small: e.g., (0.50, 0.25, 0.20, 0.05)
  - 2 small: e.g., (0.45, 0.45, 0.03, 0.07)
  - 2 tiny: e.g., (0.49, 0.49, 0.01, 0.01)
- **Target events ~350 needs verification.** With control median = 14 months, HR = 0.65, accrual 18 months, cutoff 36 months, the expected event fraction depends on the Weibull parameterization. The plan should include a quick calculation or reference showing ~350 is achievable. With N = 500, this implies ~70% event rate, which is plausible but should be confirmed via pilot.

---

## 2. O'Brien-Fleming Boundaries with gsDesign

### ✅ Sound
- Using `gsDesign::gsDesign(k=3, test.type=2, alpha=0.05)` for O'Brien-Fleming is correct (`test.type=2` = O'Brien-Fleming spending).
- Information fractions at 33%, 66%, 100% of target events is standard practice (information fraction, not calendar time).
- The reported boundaries (z = 3.004, 2.124, 1.683) match O'Brien-Fleming for k=3.

### ⚠️ Issues
- **Information fraction calculation matters.** The plan says "33%, 66%, 100% of target events" but in practice, information fraction should be computed as (events observed / target events) at each look, not a fixed calendar fraction. This is a subtle but important distinction — if the simulation produces fewer/more events than expected, the looks should still trigger at the correct information fractions.
- **Alpha spending function not specified.** `gsDesign` defaults to O'Brien-Fleming spending, but the plan should explicitly state `test = 2` in the `gsDesign` call for the spending function, not just for the boundary generation.
- **Multiple comparison correction for pooling vs. no-pooling comparison.** The plan measures "power difference (pool − no pool)" but doesn't address how to account for the fact that both analyses use the same data — the comparison itself needs a paired test or confidence interval, not just a raw difference.

---

## 3. Weibull Survival Generation and HR Parameterization

### ⚠️ Critical Issue — Weibull Parameterization

The plan states:
> Control: T ~ Weibull(scale=14, shape=1) → median 14 months

This is **incorrect**. For a Weibull distribution:

- **Median = scale × (ln 2)^(1/shape)**
- For shape = 1 (exponential): median = scale × ln(2) ≈ 0.693 × scale
- So **scale = 14, shape = 1 gives median = 14 × 0.693 = 9.7 months**, NOT 14 months.

To get median = 14 months with Weibull:
- **shape = 1 (exponential):** scale = 14 / ln(2) ≈ **20.20**
- **shape = γ (general):** scale = 14 / (ln 2)^(1/γ)

**The plan should specify a shape parameter γ.** Common choices in oncology:
- γ = 1 (exponential, constant hazard) — simplest, but assumes constant hazard
- γ = 1.5 (monotonically decreasing hazard) — more realistic for many cancers
- γ = 2 (Rayleigh, increasing hazard) — used in some applications

**Recommendation:** Specify γ explicitly and compute scale from the desired median:
```r
scale_control <- median_control / (log(2)^(1/shape))
scale_treatment <- scale_control / HR  # for HR parameterization
```

### ⚠️ Treatment Weibull Parameterization
The plan says:
> Treatment: T ~ Weibull(scale=14/HR, shape=1)

This is **also incorrect** for the same reason. If the control median is 14 months with shape γ, then:
- Control: T ~ Weibull(scale = 14/(ln 2)^(1/γ), shape = γ)
- Treatment: T ~ Weibull(scale = 14/(HR × (ln 2)^(1/γ)), shape = γ)

The HR parameterization via scale ratio is correct in principle (ratio of scales = 1/HR), but the base scale must be computed correctly.

### ✅ Sound
- HR = 0.65 is a realistic oncology effect size.
- Random censoring with C ~ exp(rate = -log(0.95)/12) is correct for 5% annual dropout.

---

## 4. Stratified Cox and Log-Rank Implementation

### ✅ Sound
- `coxph(Surv(time, event) ~ trt + strata(stratum))` is the correct specification for stratified Cox PH.
- `survdiff(Surv(time, event) ~ trt + strata(stratum))` is correct for stratified log-rank.
- Pooling strata < 10 patients into nearest larger stratum is a reasonable pooling strategy.

### ⚠️ Issues
- **Pooling rule is underspecified.** "Nearest larger stratum" needs a definition. Nearest by what metric? Proportion similarity? Random tie-breaking? This matters for reproducibility.
- **Pooling threshold of 10 patients may be too high.** With N = 500 and 4 strata, even the "balanced" scenario has 125 per stratum. But in the "2 tiny" scenario (1%, 2%), the tiny strata have only 5 and 10 patients respectively. The 1% stratum (5 patients) would be pooled, which is fine, but the 2% stratum (10 patients) sits exactly at the threshold. Consider using a proportional threshold (e.g., < 2% of N) rather than absolute.
- **Log-rank test with pooled strata.** The plan says "both with pooling" but the `survdiff` call with `strata()` after pooling would still treat each original stratum level as a separate stratum unless the stratum variable is updated. The pooling step must **reassign** stratum labels, not just merge data rows.

---

## 5. Sample Size, Event Count, and Power Expectations

### ⚠️ Issues
- **10,000 reps for Type I error is fine, but 10,000 reps for power is overkill.** Power estimates stabilize around 1,000–2,000 reps for binary outcomes (power = proportion). The computational cost of 10,000 reps with group sequential looks and multiple analyses (Cox + log-rank, pooled + unpooled) per rep is significant. Consider:
  - 10,000 reps for Type I error (need precision on small α)
  - 2,000–5,000 reps for power (sufficient for ~95% CIs on power estimates)
- **Convergence failure rate needs a precise definition.** "Non-convergent iterations" is vague. Should this check:
  - `coxph` `convergence` code ≠ 0?
  - |gradient| > tolerance?
  - Hessian not positive definite?
- **Zero-count strata handling for log-rank.** The plan mentions this but doesn't specify what happens when a stratum has zero events. `survdiff` will silently drop empty strata, which is fine, but the simulation should verify this behavior and document it.

---

## 6. Potential Issues with Tiny Strata and Cox Convergence

### ⚠️ Critical Issues

1. **Cox model with tiny strata may produce unstable estimates.** With 5 patients in a stratum, the partial likelihood contribution from that stratum is based on very few risk sets. This is the core phenomenon the simulation should investigate — but the plan doesn't explicitly state that this instability is the **hypothesis** being tested.

2. **Stratified log-rank with tiny strata.** If a stratum has 0–2 events, the log-rank statistic for that stratum is essentially noise. The stratified test aggregates across strata, so tiny strata add noise without signal. This is a known issue and the simulation should explicitly test whether pooling improves or worsens power in this regime.

3. **Pool-and-test circularity.** If we pool strata and then test, the pooling itself changes the stratum composition. In the extreme case where all 4 strata are tiny (e.g., N = 100), pooling could merge everything into 1–2 strata, effectively removing stratification. The plan should address this edge case.

4. **No mention of proportional hazards (PH) assumption.** Stratified Cox does not assume a common hazard ratio across strata — that's the whole point of stratification. But the simulation generates data with a **single HR across all strata**, which means PH holds. The plan should explicitly note that this simulation only addresses the **small-strata precision** question, not the PH violation question.

5. **No mention of time-varying effects.** If the true effect varies over time (common in oncology), the Cox model's constant HR assumption may not hold. This is a separate question but worth noting as a limitation.

---

## Summary of Required Revisions

| Priority | Issue | Action |
|----------|-------|--------|
| **Must fix** | Weibull parameterization (scale/median confusion) | Specify shape γ, compute scale from median correctly |
| **Must fix** | Strata proportions not specified | List all 4 scenarios with exact proportions |
| **Must fix** | Information fraction trigger mechanism | Clarify that looks trigger at events observed / target events milestones |
| **Important** | Pooling threshold justification | Justify threshold of 10 or switch to proportional threshold |
| **Important** | "Nearest larger stratum" undefined | Define the pooling rule precisely |
| **Important** | Power vs. Type I error rep counts | Reduce power reps to 2,000–5,000; keep 10,000 for Type I error |
| **Good to add** | PH assumption note | Explicitly state simulation assumes common HR across strata |
| **Good to add** | Edge case analysis | What happens when pooling eliminates all stratification? |
| **Nice to have** | Shape parameter sensitivity | Consider testing multiple γ values (1, 1.5, 2) |

---

## Overall Assessment

This is a **well-structured** simulation plan with clear objectives and a sensible design framework. The primary concern is the **Weibull parameterization error** (items 3 and 6 above), which would produce incorrect survival times and thus biased results if uncorrected. The strata proportions must also be specified for reproducibility.

The statistical approach (stratified Cox + log-rank with O'Brien-Fleming boundaries) is sound for the stated objectives. The main scientific question — whether pooling tiny strata improves power without inflating Type I error — is well-motivated and worth investigating.

**Recommendation: Minor Revisions.** Fix the Weibull parameterization and strata proportions, clarify the pooling rule and information fraction triggers, and the plan is ready for coding.
