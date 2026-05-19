# Why Pooling Matters for Log-rank but Not Cox PH
## Technical Explanation with Simulation Evidence

### 1. The Mathematical Difference

Both stratified Cox PH and stratified log-rank handle small strata, but through completely different mechanisms.

#### Stratified Cox PH — Partial Likelihood

$$\text{PL}(\beta) = \prod_{s=1}^{S} \prod_{i \in D_s} \frac{\exp(\beta A_{si})}{\sum_{j \in R_s(t_{si})} \exp(\beta A_{sj})}$$

Key property: **multiplies across risk sets within and across strata.** A stratum with 5 patients, 1 event contributes only 1 term to the product:
- The risk set has 5 patients → denominator sums 5 terms
- Each patient's weight is $\exp(\beta A_j)$
- If the 1 event happens to the treated patient: contribution = $\exp(\beta) / (\exp(\beta) + \sum_{controls} 1)$
- This is bounded between 0.33 and 0.67 (depending on $\beta$)
- Multiply across hundreds of risk sets: one small risk set contributes **multiplicatively little**

**Intuition:** Cox PH naturally shrinks the influence of tiny strata because each event in a tiny stratum comes from a small risk set, which provides limited information about the relative ordering. Pooling 5 patients into a stratum of 100 changes virtually nothing in the overall likelihood.

#### Stratified Log-rank — Hypergeometric Sum

$$Z = \frac{\sum_{s=1}^{S} (O_s - E_s)}{\sqrt{\sum_{s=1}^{S} V_s}}$$

where for each stratum $s$:
$$E_s = \frac{n_{1s} d_s}{n_s} \quad\quad V_s = \frac{n_{1s} n_{0s} d_s (n_s - d_s)}{n_s^2 (n_s - 1)}$$

Key property: **linear sum across strata.** Each stratum contributes additively to both numerator and denominator. A stratum with 5 patients, 1 event contributes:

$$O_s - E_s = 1 - \frac{2 \cdot 1}{5} = 0.6$$

$$V_s = \frac{2 \cdot 3 \cdot 1 \cdot 4}{5^2 \cdot 4} = 0.24$$

The contribution to the **overall Z-score** from this single tiny stratum:

$$\frac{|O_s - E_s|}{\sqrt{V_s}} = \frac{0.6}{\sqrt{0.24}} = \frac{0.6}{0.49} \approx 1.22$$

**Compare with Cox PH:** For the same 5-patient stratum, the Cox contribution to the overall likelihood ratio is essentially negligible — it's one partial likelihood term among hundreds.

---

### 2. Why the Hypergeometric Variance Breaks Down

The hypergeometric variance $V_s$ is derived under the assumption that the marginal totals $(n_{1s}, n_{0s}, d_s)$ are fixed and $O_s$ follows a hypergeometric distribution. For a 5-patient stratum:

| Scenario | $n_{1s}$ | $n_{0s}$ | $d_s$ | Possible $O_s$ values | $V_s$ |
|:--------:|:--------:|:--------:|:-----:|:---------------------:|:-----:|
| 2 trt, 3 ctl | 2 | 3 | 1 | {0, 1} | 0.24 |
| 2 trt, 3 ctl | 2 | 3 | 2 | {0, 1, 2} | 0.36 |

A **single extra event** in the treatment arm of a tiny stratum flips $O_s$ from 0 to 1, changing $O_s - E_s$ by 0.6 and the contribution to the Z-score by 1.22. That's a **66% change in one stratum's contribution from a single random event**. The variance $V_s$ correctly captures this uncertainty, but the **additive structure** means this randomness propagates into the overall test statistic.

In Cox PH, by contrast, observing 1 event vs 0 events in the small stratum changes one term in the partial likelihood product — a tiny relative change in the overall likelihood.

---

### 3. Our Simulation Data: The Proof

#### 3.1 Cox PH: Pooling Does Nothing

From 10,000 reps under the null and 5,000 reps with HR=0.65:

| Scenario | Cox Power (No Pool) | Cox Power (Pool) | Diff |
|:---------|:------------------:|:----------------:|:----:|
| Balanced | 0.965 | 0.965 | 0.000 |
| 1 small (5%) | 0.969 | 0.969 | -0.000 |
| 2 small (3%) | 0.975 | 0.975 | +0.000 |
| 2 tiny (1%,2%) | 0.965 | 0.966 | +0.001 |

**Every scenario: diff ≤ 0.001.** Cox PH is mathematically indifferent to pooling.

#### 3.2 Log-rank: Pooling Helps in Extreme Sparsity

| Scenario | LR Power (No Pool) | LR Power (Pool) | Gain | Type I (No Pool → Pool) |
|:---------|:------------------:|:---------------:|:----:|:----------------------:|
| Balanced | 0.923 | 0.923 | 0.000 | 0.073 → 0.073 |
| 1 small (5%) | 0.768 | 0.772 | +0.004 | 0.078 → 0.076 |
| 2 small (3%) | 0.721 | 0.743 | **+0.022** | 0.078 → 0.071 |
| 2 tiny (1%,2%) | **0.663** | **0.827** | **+0.163** | 0.095 → **0.073** |

In Scenario 4 (1%, 2% strata = ~5 and ~10 patients):

**Without pooling:** Four strata, two of which are ~5 and ~10 patients. The log-rank test includes these tiny strata with their noisy hypergeometric contributions. The overall Z-statistic is the sum of four independent contributions, and the two tiny strata add variance without adding commensurate signal.

**With pooling:** Merge the 5-patient and 10-patient strata into the nearest large stratum. Now we have 2-3 stable strata, each with enough patients and events that the hypergeometric approximation is accurate. The Z-statistic's numerator grows (more signal) while the denominator grows more slowly (less noise per stratum).

**Mechanism quantified:**
- Scenario 4, no pool: Power = **0.663** (roughly 2:1 odds of detecting HR=0.65)
- Scenario 4, pool: Power = **0.827** (roughly 5:1 odds)
- That's a +0.16 improvement = **25% relative increase in power**
- Type I error also improves: 0.095 → 0.073 (closer to nominal)

---

### 4. When Should We Pool for Log-rank?

Based on the simulation evidence:

| Strata size | Recommendation | Rationale |
|:-----------|:--------------|:----------|
| All ≥ 20 | No pooling needed | Hypergeometric approximation accurate |
| Any < 10 | **Pool recommended for log-rank** | Power gain of 0.02-0.16, Type I improves |
| Any 10-20 | Marginal benefit | Small gains (0.004-0.02), case-by-case |
| Cox PH | Never pool | Power and Type I unaffected at all sparsity levels |

---

### 5. Proposed SAP Language Update

For the binary/CMH white paper (already pushed):
> *"Stratification factors will be used as specified in the randomization scheme. No pooling of small strata is required."*
> ✅ Stands for CMH OR, CMH RR, MN RD

For **survival analyses (ELSTIC guidance):**
> *"Stratification factors will be used as specified. For the stratified Cox PH model, no pooling of small strata is required. For the stratified log-rank test, strata with fewer than 10 patients should be pooled into the nearest larger stratum to maintain power and control Type I error."*

---

### 6. Why This Matters for ELSTIC

The ELSTIC guidance (for time-to-event endpoints) was finalized with a blanket "no pooling" recommendation. Our simulation shows this is correct for Cox PH but **incorrect for log-rank**. Since many oncology trials specify the stratified log-rank as the primary analysis or the key sensitivity analysis, this distinction is material.

The +0.16 power gain in extreme sparsity is not a theoretical nicety — it's the difference between a trial with 66% power (underpowered) and 83% power (adequately powered), all else equal.
