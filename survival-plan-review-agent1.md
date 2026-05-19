# Agent Review: Survival Simulation Plan — Small Strata Pooling

**Reviewer:** Agent 1  
**Date:** 2026-05-19  
**Plan file:** `survival-simulation-plan.md`  

---

## Verdict: Major Revisions ⚠️

Three issues need fixing before coding: a material error in the Weibull parameterization (Item 3), an under-specified pooling algorithm (Item 6), and missing details on event monitoring for group sequential looks (Item 2). The other items are sound with minor suggestions.

---

## 1. Trial Design Realism ✅

**Accrual (28 patients/month × 18 months → 500 patients):** Realistic for a moderate-size phase 3 oncology trial. Many registration trials in advanced solid tumors accrue at similar or slightly faster rates.

**Control median survival = 14 months:** Plausible for advanced/metastatic settings (e.g., 2L NSCLC, pancreatic, gastric). No issue.

**HR = 0.65:** Realistic for an active oncology experimental arm. Many contemporary targeted therapies and immunotherapies report HRs in the 0.55–0.75 range.

**Random censoring (5% annual loss to follow-up):** Standard and reasonable.

**Target events (~350 out of 500):** With control median survival 14 months, treatment median ~21.5 months (at HR=0.65), 18-month uniform accrual, and 36-month administrative cutoff, 70% event rate is plausible. However, see Item 3 — the actual median under the stated Weibull parameters gives a *different* event rate that should be rechecked.

**Stratified block randomization (block size 4, 1:1, 2 binary factors → 4 strata):** Standard and realistic for a stratified oncology trial.

---

## 2. Group Sequential Boundaries ⚠️

The stated boundaries from `gsDesign::gsDesign(k=3, test.type=2, alpha=0.05)` are correct for the O'Brien-Fleming design with three **equally spaced** looks at information fractions 1/3, 2/3, 1:

| Look | z-boundary | Nominal α (2-sided) | Correct? |
|------|-----------|---------------------|----------|
| 1    | 3.004     | 0.0027              | ✓        |
| 2    | 2.124     | 0.0337              | ✓        |
| 3    | 1.683     | 0.0923              | ✓        |

**Issue — Event-driven looks not exact fractions:** The plan states looks occur at "33%, 66%, 100% of target events," but this is stochastic. The simulation should specify:

- **Does each rep use the exact same event counts** (e.g., floor(0.33×350) = 115 events at look 1)? Or does the actual event milestone vary per rep?
- **How to handle a rep that never reaches 33% of events within 36 months?** (Unlikely with these parameters, but possible under HR=1.0 with longer survival than expected.)
- **Should looks be calendar-driven** (e.g., at fixed calendar times) **or event-driven** (at specific event counts)? The plan implies event-driven, which is standard — but should be explicit.

**Minor note:** `test.type=2` is two-sided symmetric, testing both superiority and inferiority simultaneously. This is fine, but if the only interest is showing superiority (HR < 1), `test.type=1` with a one-sided α=0.025 is more conventional in oncology.

---

## 3. Weibull Parameterization ❌ **MATERIAL ERROR**

The plan states:

```
Control: T ~ Weibull(scale=14, shape=1) → median 14 months
Treatment: T ~ Weibull(scale=14/HR, shape=1)
```

**This is incorrect for the stated median survival.**

For R's `rweibull(n, shape=k, scale=λ)`, the median is:

```
median = λ × (ln 2)^(1/k)
```

With k=1 (exponential) and λ=14:

```
median = 14 × ln(2) ≈ 14 × 0.693 = **9.7 months**, not 14 months
```

**Impact:** This affects nearly everything downstream:
- Control median will be ~9.7 months instead of 14 months → higher event rate than expected
- Power estimates will be inflated (more events occur before cutoff)
- The "realistic oncology" scenario becomes more aggressive (shorter control survival)
- Treatment median = 14/0.65 × ln(2) ≈ 14.9 months (not ~21.5 months as implied)
- The HR *ratio* is still 0.65 (proportional hazards property holds), but baseline survival is shifted

**Fix:** Either (a) set `scale = 14 / ln(2) ≈ 20.2` with `shape=1`, or (b) if shape ≠ 1 is acceptable, solve `λ × (ln 2)^(1/k) = 14` for any desired k. The simplest fix is option (a).

---

## 4. Output Metrics ✅ (with gaps)

**Well-defined:**
- Power / Type I error by look ✓
- Convergence failure rate ✓
- Power difference (pool − no pool) ✓

**Missing or under-specified:**

| Gap | Why it matters |
|-----|----------------|
| **Cumulative (overall) power** — proportion rejecting at *any* look | Often the primary quantity of interest; "by look" power is useful but the "ever-reject" rate is what regulators care about |
| **Bias in HR estimate** — does pooling bias the estimated treatment effect? | Key for the "pooling" research question; merging strata with different baseline hazards could induce bias |
| **Empirical SE of log(HR)** — does pooling affect precision? | Important for comparing efficiency |
| **HR coverage** — does the 95% CI from the pooled analysis have correct coverage? | Pooling changes variance estimation |
| **Events per stratum** — distribution across reps | Diagnostic for understanding when/why pooling helps or hurts |
| **Proportion of reps with ≥1 empty stratum in one arm** | Critical for assessing when the unpooled Cox model fails |

**Recommendation:** Add cumulative power, HR bias (mean log(HR) − true log(HR)), and empirical SE to the output metrics.

---

## 5. Simulation Pipeline ✅

The 2 → 100 → 10,000 pipeline is standard and sensible:

| Stage | Purpose | Sample size justification |
|-------|---------|--------------------------|
| 2 reps | Syntax check | Minimal, adequate |
| 100 reps | Proof of concept / bug fixes | Enough to surface convergence failures and edge cases |
| 10,000 reps | Final estimation | SE(p̂) ≈ 0.22% for p=0.05, ≈ 0.4% for p=0.80 — fine |

**Computational load:** With 4 scenarios × 2 HR values (0.65, 1.0) × 2 methods (pooled, unpooled) = 16 conditions × 10,000 reps = 160,000 trial simulations. Each trial generates 500 patients × 3 looks of Cox + log-rank. This will take significant compute but is feasible with `furrr` on a multi-core machine (estimate: hours to a day). Consider whether the null (HR=1) needs 10,000 for all 4 scenarios or if fewer reps suffice for some.

**Note:** The plan should specify whether a **common seed scheme** (e.g., `seed = scenario_id * 1e6 + rep_number`) will be used for reproducibility across the pipeline stages. Reproducibility is critical for debugging between the 2-rep and 100-rep stages.

---

## 6. Missing Edge Cases ❌

### 6a. Pooling algorithm under-specified

> "merge strata < 10 patients into nearest larger stratum"

**Unresolved questions:**
- **Distance metric:** "Nearest" by what measure? Proportion? Patient count? Baseline hazard similarity?
- **Recursive merging:** After merging one small stratum, the receiving stratum may become the *next* largest. What if you merge both 1% strata into the same larger one? Does the order of merging matter?
- **All strata small:** In the "2 tiny (1%, 2%)" scenario, you'll have strata of ~5 and ~10 patients (plus two larger ones). After merging the 5-patient stratum, the recipient is still fine. But is there a scenario where after merging, the recipient still has < 10?
- **Pooling across treatment arms:** Does pooling happen on the stratum variable *before* the analysis, i.e., some patients get reassigned to a different stratum label? Or is it purely analytical (like a merged factor level)?

### 6b. Eventless strata in Cox PH

A stratum with 5 patients and block size 4 can easily end up with 2 treatment + 3 control. If neither arm has an event, the Cox model's stratum-specific baseline hazard is undefined. The plan mentions "convergence failure" but doesn't specify how to detect it:

- **Cox convergence detection:** Should check `coxph` return value — `fit$convergence` or `fit$iter` exceeding max, or catching the `singular`/`non-converged` warning via `tryCatch`.
- **What to do with non-convergent reps:** Drop them? Record as failure? Impute? The plan mentions "convergence failure rate" as an output metric, so these reps aren't dropped — but the denominator for power must be clear: is it *all* reps, or only *converged* reps?

### 6c. Log-rank with sparse strata

`survdiff` can produce NA p-values or NaN chi-square statistics when a stratum has events in only one arm or zero total events. The plan mentions "zero-count strata" for log-rank but should also handle the case of **one-event strata** (e.g., 1 event out of 5 patients in a stratum — chi-square with 1 df but highly unstable).

### 6d. Randomization failure with small strata

With block randomization (block size 4) in small strata (e.g., 1% = 5 patients), a single block of 4 fits, leaving 1 patient. That patient gets the next treatment assignment from a partially filled block. This can cause moderate **imbalance** in very small strata (e.g., 4:1 instead of 2:3 in a 5-patient stratum). The simulation should check and possibly report the actual arm sizes per stratum.

### 6e. No null scenario with heterogenous strata effects

The Type I error simulation (HR=1.0) tests the case where treatment has no effect *and* strata have similar baseline hazards (since all Weibull(scale=14)). A more realistic null might include **heterogeneous baseline hazards across strata** (e.g., by adding a stratum-specific frailty or scale multiplier). This would test whether pooling confounded strata inflates Type I error more than not pooling.

### 6f. No mention of multiple-testing correction

The plan evaluates both Cox PH and log-rank on the same data. This is descriptive rather than inferential (comparing a method property), so corrections aren't needed — but it should be stated explicitly that no multiplicity adjustment is applied across methods.

---

## Summary of Required Changes

| # | Severity | Issue | Fix |
|---|----------|-------|-----|
| 1 | **Critical** | Weibull shape=1, scale=14 → median 9.7, not 14 | Change scale to 14/ln(2) ≈ 20.2, or choose a shape that gives median 14 |
| 2 | **High** | Pooling algorithm underspecified | Define distance metric, merging order, and handling of edge cases for "nearest" |
| 3 | **High** | Event-driven looks not fully specified | State exact event targets, boundary handling for reps that don't reach milestones |
| 4 | **Medium** | Cox convergence detection not defined | Specify `tryCatch` logic and convergence flag checking |
| 5 | **Medium** | Missing cumulative (ever-reject) power | Add to output metrics |
| 6 | **Medium** | Missing HR bias / coverage metrics | Add to support the pooling research question |
| 7 | **Low** | No seed scheme for reproducibility | Add a deterministic seed scheme (e.g., `scenario * 1e6 + rep`) |
| 8 | **Low** | test.type=2 vs test.type=1 | Consider whether symmetric two-sided testing is intended or one-sided superiority |

---

## Overall: Major Revisions Required ⚠️

The plan has a clear analytical framework and the overall approach (trial design, group sequential design, simulation pipeline) is sound. However, the Weibull parameterization error is a material problem that will invalidate all downstream results if not fixed. The pooling algorithm needs sharper definition to prevent ambiguous coding. Once items 1–3 are addressed, this plan is ready for implementation.
