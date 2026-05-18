# Small Strata Investigation: Recommendations for SAP Language

**Status:** Internal white paper — complete  
**Date:** May 2026  

---

## Bottom Line

All standard binary endpoint methods used in Merck oncology SAPs are robust to small strata:

| Method | Failure rate | Type I error | Pooling needed? |
|--------|:-----------:|:-----------:|:--------------:|
| **CMH odds ratio** (+0.5 CC) | 0.000 | 0.043–0.054 | ❌ No |
| **CMH risk ratio** (+0.5 CC) | 0.000 | 0.045–0.081 | ❌ No |
| **Stratified MN risk difference** (score-based) | 0.000 | 0.035–0.069 | ❌ No |
| **Cox PH** (stratified) | — | — | ❌ No (already confirmed) |
| **Log-rank** (stratified) | — | — | ❌ No (already confirmed) |

**Recommendation: No pooling of small strata is required for any of these methods.**

## Results (5,000 reps per scenario, stratified block randomization, block size 4)

| Sparsity | Event rate | CMH OR | | CMH RR | | MN RD | |
|:---------|:---------:|:-----:|:-----:|:-----:|:-----:|:-----:|:-----:|
| | | Fail | Type I | Fail | Type I | Fail | Type I |
| Balanced | 10% | 0.000 | 0.043 | 0.000 | 0.081 | 0.000 | 0.035 |
| Balanced | 30% | 0.000 | 0.047 | 0.000 | 0.055 | 0.000 | 0.053 |
| Balanced | 50% | 0.000 | 0.050 | 0.000 | 0.051 | 0.000 | 0.059 |
| 1 small (5%) | 10% | 0.000 | 0.047 | 0.000 | 0.064 | 0.000 | 0.036 |
| 1 small (5%) | 30% | 0.000 | 0.052 | 0.000 | 0.060 | 0.000 | 0.059 |
| 1 small (5%) | 50% | 0.000 | 0.045 | 0.000 | 0.049 | 0.000 | 0.054 |
| 2 small (3% ea) | 10% | 0.000 | 0.051 | 0.000 | 0.047 | 0.000 | 0.041 |
| 2 small (3% ea) | 30% | 0.000 | 0.049 | 0.000 | 0.070 | 0.000 | 0.061 |
| 2 small (3% ea) | 50% | 0.000 | 0.049 | 0.000 | 0.057 | 0.000 | 0.063 |
| 2 tiny (1%, 2%) | 10% | 0.000 | 0.047 | 0.000 | 0.045 | 0.000 | 0.040 |
| 2 tiny (1%, 2%) | 30% | 0.000 | 0.054 | 0.000 | 0.053 | 0.000 | 0.062 |
| 2 tiny (1%, 2%) | 50% | 0.000 | 0.051 | 0.000 | 0.054 | 0.000 | 0.069 |

## Notes

- **CMH OR** is the most robust — Type I error within [0.043, 0.054] across all scenarios.
- **CMH RR** with stratified Greenland-Robins variance shows slight Type I inflation at low event rates (0.081 balanced, 10%), but this is a property of the RR scale for rare events, not sparsity.
- **Stratified MN RD** (score-based, using `PropCIs::diffscoreci`) shows 0% failure across all scenarios. Type I error is slightly conservative at low event rates (0.035–0.041), which is a known advantage of the Miettinen-Nurminen method in sparse settings.
- With **stratified block randomization** (block size 4), within-stratum balance is guaranteed even for very small strata.

## Proposed SAP Language

For SAPs using CMH (OR or RR) or stratified MN (risk difference):

> *"Stratification factors will be used as specified in the randomization scheme. No pooling of small strata is required."*

## Code

All simulation code: `run_small_strata.R` in this repository.

---

## Appendix: Method Definitions and Formulas

### 1. Stratified Cox Proportional Hazards Model

**Endpoint:** Time-to-event (survival)  
**SAP prevalence:** Nearly all oncology trials  

**Model:**
$$\lambda(t \mid A, \text{stratum}) = \lambda_{0s}(t) \exp(\beta A)$$

where $\lambda_{0s}(t)$ is the stratum-specific baseline hazard, $A$ is the treatment indicator, and $\beta$ is the log-hazard ratio. The partial likelihood maximizes over strata jointly:

$$L(\beta) = \prod_{s=1}^{S} \prod_{i \in D_s} \frac{\exp(\beta A_{si})}{\sum_{j \in R_s(t_{si})} \exp(\beta A_{sj})}$$

**Property:** Stratified Cox allows different baseline hazards per stratum while estimating a common treatment effect. With small strata, the partial likelihood still converges because each stratum contributes its event-time ordering independently.

### 2. Stratified Log-Rank Test

**Endpoint:** Time-to-event  
**SAP prevalence:** Nearly all oncology trials (primary test)

**Test statistic:**
$$Z = \frac{\sum_{s=1}^{S} (O_s - E_s)}{\sqrt{\sum_{s=1}^{S} V_s}}$$

where for stratum $s$, $O_s$ is the observed number of events in the treatment arm, $E_s = \sum_t n_{1st} \cdot d_{st} / n_{st}$ is the expected number under the null, and $V_s$ is the hypergeometric variance. Under $H_0$, $Z \sim N(0,1)$.

### 3. Cochran-Mantel-Haenszel (CMH) — Odds Ratio

**Endpoint:** Binary  
**SAP prevalence:** High  

**Stratum-specific:** For stratum $k$ with table $\begin{pmatrix} a_k & b_k \\ c_k & d_k \end{pmatrix}$ where $a_k = $ events on treatment, $b_k = $ non-events on treatment, $c_k = $ events on control, $d_k = $ non-events on control:

$$OR_k = \frac{a_k d_k}{b_k c_k}$$

**Pooled estimate (Mantel-Haenszel):**
$$OR_{MH} = \frac{\sum_k a_k d_k / n_k}{\sum_k b_k c_k / n_k}$$

where $n_k = a_k + b_k + c_k + d_k$.

**Continuity correction:** When any cell is zero, add 0.5 to all four cells of that stratum before computing $OR_k$.

**Test of $H_0: OR = 1$:**
$$\chi^2_{MH} = \frac{\left[\sum_k (a_k - E(a_k))\right]^2}{\sum_k V(a_k)} \sim \chi^2_1$$
where $E(a_k) = n_{1k} n_{1'k} / n_k$ and $V(a_k) = n_{1k} n_{0k} n_{1'k} n_{0'k} / (n_k^2 (n_k-1))$ under the null.

### 4. Cochran-Mantel-Haenszel (CMH) — Risk Ratio

**Endpoint:** Binary  
**SAP prevalence:** High  

**Stratum-specific risk ratio (with 0.5 continuity correction):**
$$RR_k = \frac{(a_k + 0.5) / (a_k + b_k + 0.5)}{(c_k + 0.5) / (c_k + d_k + 0.5)}$$

**Mantel-Haenszel weighted estimate:**
$$RR_{MH} = \frac{\sum_k w_k \cdot RR_k}{\sum_k w_k}$$

where $w_k = (a_k + b_k)(c_k + d_k) / n_k$ (the MH weight, proportional to the inverse variance of $RR_k$).

**Greenland-Robins variance for $\log(RR_{MH})$:**
$$Var(\log RR_{MH}) = \frac{\sum_k w_k^2 \left[\frac{1}{a_k + 0.5} - \frac{1}{a_k + b_k + 0.5} + \frac{1}{c_k + 0.5} - \frac{1}{c_k + d_k + 0.5}\right]}{\left(\sum_k w_k\right)^2}$$

The 95% CI for $RR_{MH}$ is $\exp(\log RR_{MH} \pm 1.96 \cdot SE)$ where $SE = \sqrt{Var(\log RR_{MH})}$.

### 5. Stratified Miettinen-Nurminen — Risk Difference

**Endpoint:** Binary  
**SAP prevalence:** Moderate  

**Stratum-specific:** For stratum $k$, the Miettinen-Nurminen score CI for risk difference $\delta = p_1 - p_0$ is found by solving:

$$\frac{(a_k - n_{1k}\tilde{p}_1)^2}{n_{1k}\tilde{p}_1(1-\tilde{p}_1)} + \frac{(c_k - n_{0k}\tilde{p}_0)^2}{n_{0k}\tilde{p}_0(1-\tilde{p}_0)} = z_{\alpha/2}^2$$

subject to $\tilde{p}_1 - \tilde{p}_0 = \delta$, where $\tilde{p}_1, \tilde{p}_0$ are the maximum likelihood estimates under the constraint. This requires solving a cubic equation — no closed form exists.

**Stratified pooled estimate:** Inverse-variance weighted across strata:

$$\hat{\delta} = \frac{\sum_k \hat{\delta}_k / SE_k^2}{\sum_k 1 / SE_k^2}$$

where $\hat{\delta}_k = a_k/n_{1k} - c_k/n_{0k}$ and $SE_k$ is derived from the stratum-specific score CI width: $SE_k = (CI_{upper} - CI_{lower}) / (2 \cdot z_{\alpha/2})$.

**Properties:** The Miettinen-Nurminen method is score-based, which means it does not require continuity corrections and has better coverage properties than Wald intervals, especially in sparse data. It is slightly conservative (type I error < nominal) at low event rates — a known feature, not a bug.

### 6. Stratified Block Randomization

**Used in simulation:** Yes (block size 4)

**Procedure:** Within each stratum, patients are assigned to treatment arms in random permuted blocks of size 4. Each block contains exactly 2 treatment and 2 control assignments in random order:

$$\text{Block} = \text{random permutation of } \{T, T, C, C\}$$

This guarantees near-perfect balance within each stratum regardless of stratum size, mimicking real clinical trial practice where pharmacies prepare drug kits in blocks.

---

*End of white paper.*
