# Internal Investigation: Do Small Strata Require Pooling for Standard SAP Methods?

**Status:** Internal white paper — not a publication  
**Goal:** Evidence-based recommendation for Merck SAP language on stratum pooling  
**Methods:** Only methods currently used in Merck oncology SAPs (no novel methodology)

---

## Background

Stratified randomization creates strata by stratification factors. When some strata have very few patients, SAPs require pre-specified pooling rules. Merck internal investigation found that for Cox PH and log-rank tests, pooling is unnecessary — they converge and produce valid inference regardless of stratum size.

The open question: do the binary endpoint methods used in Merck SAPs have the same robustness?

## Methods to Investigate

| Method | Endpoint | SAP use | Continuity correction |
|--------|----------|:-------:|:---------------------:|
| Cox PH (stratified) | TTE | Nearly all | N/A — already investigated |
| Log-rank (stratified) | TTE | Nearly all | N/A — already investigated |
| CMH (odds ratio) | Binary | High | +0.5 continuity correction |
| CMH (risk ratio) | Binary | High | +0.5 continuity correction |
| MN (risk difference) | Binary | Moderate | None needed (score-based) |

No logistic regression — not used in Merck oncology SAPs for binary endpoints.

## Simulation Design

**Data generation:** Stratified binary outcomes via Bernoulli with stratum-specific response probabilities. Simple DGP — the question is about numerical behavior of standard estimators under sparsity, not model performance.

**Scenarios:**
- 2-4 stratification factors (realistic for Merck oncology)
- N = 200-500 (typical phase 2/3)
- Event rates 10-50%
- Sparsity: some strata with 1-5 patients per arm
- 5,000 replicates per scenario

**Evaluation:**
- **Convergence failure rate** (primary endpoint — drives the recommendation)
- Type I error
- Coverage of 95% CIs
- Bias

## Decision Framework

Estimate failure rates first. Results will determine whether method-specific pooling recommendations are needed.

## Deliverable

Short white paper with:
1. Failure rate table (CMH OR, CMH RR, MN RD) × 3 event rates × 2-3 sparsity patterns
2. Recommendation for SAP language (method-specific if patterns differ)
3. Reproducible code
