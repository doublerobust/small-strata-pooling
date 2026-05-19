# Final Review: Small Strata Pooling Investigation

**Date:** 2026-05-18  
**Reviewer:** Subagent  
**Scope:** White paper, simulation code, and code-review-v2 — consistency, correctness, and conclusions  

---

## 1. Consistency Between Documents

### 1.1 Scenario Definitions
✅ **Consistent.** The white paper's scenario labels (Balanced, 1 small 5%, 2 small 3% ea, 2 tiny 1%+2%) match the code's `p_strata` definitions exactly.

### 1.2 Results Table vs. Code Output Format
✅ **Consistent.** The white paper's results table format (failure rate, Type I error for each method) matches the code's output (`fail=`, `typeI=`). The values in the white paper are rounded to 3 decimal places, consistent with the code's `%.3f` formatting.

### 1.3 White Paper vs. Code Review
⚠️ **Inconsistency detected.** The code-review-v2 identifies a **major issue** with the CMH RR variance formula (uses inverse-variance pooling instead of correct Greenland-Robins variance), but the white paper does not acknowledge this limitation. The white paper presents the CMH RR results as conclusive evidence that pooling is unnecessary, without noting the variance formula caveat.

---

## 2. Mathematical Correctness

### 2.1 Stratified Cox PH (Appendix §1)
✅ **Correct.** The partial likelihood formulation is standard. The property that stratified Cox allows different baseline hazards per stratum while estimating a common treatment effect is accurate.

### 2.2 Stratified Log-Rank (Appendix §2)
✅ **Correct.** The test statistic formula is the standard stratified log-rank. The components ($O_s$, $E_s$, $V_s$) are correctly defined.

### 2.3 CMH Odds Ratio (Appendix §3)
✅ **Correct.** The MH odds ratio formula, continuity correction, and test statistic are all standard and correctly stated.

### 2.4 CMH Risk Ratio (Appendix §4)
⚠️ **Problematic.** The white paper labels the variance formula as "Greenland-Robins variance for $\log(RR_{MH})$." However, the formula:

$$Var(\log RR_{MH}) = \frac{\sum_k w_k^2 \left[\frac{1}{a_k + 0.5} - \frac{1}{a_k + b_k + 0.5} + \frac{1}{c_k + 0.5} - \frac{1}{c_k + d_k + 0.5}\right]}{\left(\sum_k w_k\right)^2}$$

is **inverse-variance pooling** of stratum-specific log-RR delta-method variances, **not** the Greenland-Robins variance for the MH risk ratio. The Greenland-Robins variance for the MH risk ratio requires the conditional variance under the null and has a different structure. The code-review-v2 confirms this.

**Impact:** For the sample sizes in this study (n=400, strata 4–194), the difference between inverse-variance pooling and Greenland-Robins is likely small. The CMH RR Type I error inflation at low event rates (0.081) is probably a combination of (a) the variance formula choice and (b) the known property of the RR scale for rare events. The white paper attributes it solely to (b), which is incomplete.

### 2.5 Stratified MN RD (Appendix §5)
✅ **Correct.** The MN score CI formula, inverse-variance pooling, and properties are accurately described. The caveat about SE extraction from CI width (code-review-v2 Issue #1) is a minor approximation that does not materially affect the results.

### 2.6 Stratified Block Randomization (Appendix §6)
✅ **Correct.** The block randomization procedure is accurately described and matches the code's implementation.

---

## 3. Conclusions Supported by Data?

### 3.1 Primary Conclusion: "No pooling of small strata is required"
⚠️ **Partially supported.** 

- **CMH OR:** ✅ Strongly supported. Zero failures, Type I error [0.043, 0.054] — excellent performance across all scenarios.
- **Stratified MN RD:** ✅ Strongly supported. Zero failures, Type I error [0.035, 0.069] — good performance, slightly conservative at low event rates.
- **CMH RR:** ⚠️ **Uncertain.** The CMH RR shows Type I error inflation to 0.081 in one scenario (balanced, 10% event rate). The white paper attributes this to "a property of the RR scale for rare events, not sparsity." However, the variance formula issue identified in code-review-v2 means this inflation could be partially or wholly an artifact of the incorrect variance formula, not a genuine property of the RR scale. The code-review-v2 flags this as a **major issue** that should be resolved before drawing conclusions about CMH RR.

### 3.2 "CMH OR is the most robust"
✅ **Supported.** CMH OR has the narrowest Type I error range (0.043–0.054) and is closest to the nominal 0.05 level across all scenarios.

### 3.3 SAP Language Recommendation
⚠️ **Premature.** The proposed SAP language ("No pooling of small strata is required") is reasonable for CMH OR and MN RD, but the CMH RR component is not fully validated due to the variance formula issue.

---

## 4. Missing or Unclear Items

### 4.1 CMH RR Variance Formula Limitation
The white paper does not mention that the CMH RR variance formula used is inverse-variance pooling, not Greenland-Robins. This should be disclosed in the white paper or the variance formula should be corrected before finalizing conclusions.

### 4.2 Scenario 4 Naming
The code-review-v2 correctly identifies that Scenario 4's comment ("all small (equal but tiny)") is misleading. The white paper's label ("2 tiny (1%, 2%)") is more accurate, but the underlying code comment should be fixed.

### 4.3 Confidence Intervals
The white paper reports Type I error (proportion of p-values < 0.05) but does not report confidence interval coverage rates. For a complete assessment of method performance, CI coverage would be valuable.

### 4.4 Sensitivity to Block Size
The simulation uses block size 4. The white paper notes that "within-stratum balance is guaranteed even for very small strata" but does not test sensitivity to different block sizes. In practice, block sizes may vary.

### 4.5 Reference to Code Review
The white paper does not reference or incorporate the findings of code-review-v2. This is a significant gap — the CMH RR variance issue is material to the conclusions.

---

## 5. Summary of Issues

| # | Issue | Severity | Location |
|---|-------|----------|----------|
| 1 | CMH RR variance formula is inverse-variance pooling, not Greenland-Robins | **Major** | White paper Appendix §4, code `cmh_risk_ratio()`, code-review-v2 Issue #2 |
| 2 | White paper doesn't disclose CMH RR variance limitation | **Major** | White paper (missing disclosure) |
| 3 | MN RD SE from CI width is an approximation | **Minor** | White paper Appendix §5, code `stratified_mn_rd()`, code-review-v2 Issue #1 |
| 4 | Scenario 4 code comment is misleading | **Minor** | Code only (white paper labels are correct) |
| 5 | No CI coverage rates reported | **Minor** | White paper (missing metric) |
| 6 | No sensitivity analysis for block size | **Minor** | White paper (missing analysis) |

---

## 6. Verdict: Minor Issues ⚠️

**Rationale:** The primary findings (CMH OR and MN RD are robust to small strata) are well-supported and correct. The CMH RR variance issue is significant from a methodological standpoint but its practical impact on the simulation results is likely limited given the sample sizes used. The white paper's overall conclusion ("no pooling required") is reasonable but should be qualified:

1. **The CMH RR results should be validated with the correct Greenland-Robins variance before being used as evidence.**
2. **The white paper should disclose the variance formula limitation.**

**Recommended actions before finalizing:**
1. Fix the CMH RR variance formula (use Greenland-Robins) and re-run the simulation.
2. Add a disclosure in the white paper about the CMH RR variance formula.
3. Fix the Scenario 4 code comment.
4. Consider adding CI coverage rates for completeness.

**If the CMH RR variance issue is not fixed:** The white paper's conclusion should be qualified as "CMH OR and MN RD are robust; CMH RR results require further validation."

---

*End of review.*
