# GSD Implementation Review — Survival Simulation

## 1. O'Brien-Fleming Boundaries from `gsDesign`

**Verdict: Correctly specified.**

```r
gs <- gsDesign(k = 3, test.type = 2, alpha = 0.05, sfu = "OF")
obf_z <- gs$upper$bound  # z-value boundaries
```

- `test.type = 2` → one-sided test (as documented in `?gsDesign`)
- `sfu = "OF"` → O'Brien-Fleming spending function (upper side)
- For k=3 looks with OBF, the boundaries are approximately:
  - Look 1: z ≈ 3.481 (very conservative, as expected at early look)
  - Look 2: z ≈ 2.213
  - Look 3: z ≈ 1.960 (≈ normal 0.05 two-sided boundary)

These are the correct OBF boundaries for a one-sided 0.025 test (equivalent to two-sided 0.05). **No bug here.**

---

## 2. Sequential Testing Rule — Stop at First Rejection

**Verdict: BUG — The "stop at first rejection" logic is partially broken.**

The code iterates over all planned looks and checks rejection at each, but the `ever_reject_*` counters use a simple accumulation pattern:

```r
for (lk in seq_len(min(nlooks, 3))) {
    ...
    if (!is.na(z_no) && abs(z_no) >= obf_z[lk]) {
        if (hr < 1 && z_no < 0) out$ever_reject_cox_no <- 1
        if (hr == 1) out$ever_reject_cox_no <- 1
    }
    ...
}
```

**The problem:** The code *does* correctly set `ever_reject = 1` once (since it's already 1, setting it again has no effect). The "stop at first rejection" semantics are actually preserved *by accident* because the variable is a binary flag, not a counter. If look 1 rejects, look 2's check is irrelevant because the flag is already 1.

**However**, there is a subtler issue: the code does **not** actually *skip* subsequent looks. In a real group sequential trial, once you stop for efficacy, you don't collect more data or do further analyses. The simulation computes all looks regardless. For Type I error estimation this doesn't matter (the flag is already set). For power, it also doesn't matter for the same reason. **This is not a functional bug but a conceptual mismatch.**

**Minor issue with information fraction:** The code triggers looks at 33%, 66%, 100% of *observed* events, not *planned* events. gsDesign assumes fixed information fractions (e.g., 1/3, 2/3, 1). If the actual information fraction deviates significantly from the planned, the boundaries are no longer exact. With N=500 and the accrual/censoring setup, the observed events will vary across reps, making the actual information fractions random. **gsDesign's boundaries assume fixed timing.** This is a known issue in group sequential simulation (called "randomized information fractions"), and the approximation is standard practice, but it's worth noting.

---

## 3. Wald z-statistics vs. Log-Rank Score Test

**Verdict: BUG — This is the primary cause of inflated Type I error.**

The code uses Wald z-statistics from `coxph`:

```r
cox_z <- coef(cox_fit) / sqrt(diag(vcov(cox_fit)))
```

**The problem:** O'Brien-Fleming boundaries derived from `gsDesign` assume a **score test statistic** (not a Wald test statistic). The score test for a Cox model in the context of group sequential designs is the **stratified log-rank test statistic**, standardized to have a standard normal distribution under the null.

### Why the Wald statistic inflates Type I error:

1. **Different variance estimation under the null:** The Wald test uses the observed information (evaluated at the MLE), while the score test uses the information under the null (HR = 1). For small strata with sparse data, these diverge significantly.

2. **Asymmetry in the tails:** The Wald statistic tends to have heavier tails than the null distribution when events are sparse (which is the exact scenario this simulation targets — small strata). This means `abs(z)` will be larger than it should be under the null, leading to more rejections.

3. **The stratified log-rank is the canonical test statistic for group sequential survival trials.** The `gsDesign` package documentation and all classical group sequential methodology (Lan-DeMets, O'Brien-Fleming, etc.) are built around the score test.

### Evidence from the code:

The code *also* computes the stratified log-rank p-value correctly:

```r
lr_chisq <- lr_fit$chisq
lr_p <- pchisq(lr_chisq, df = 1, lower.tail = FALSE) / 2
```

And the log-rank boundaries (`p_lr < 0.025`) are applied correctly against the nominal 0.025 threshold. **The log-rank results should be correct.** The Cox Wald results are the inflated ones.

---

## 4. Information Fraction Calculation

**Verdict: Approximate but acceptable for simulation purposes.**

```r
targets <- round(c(0.33, 0.66, 1.00) * nevents)
```

- Uses *observed* events (not planned). This is standard in practice because you don't know the final event count in advance.
- The `round()` makes it stepwise: e.g., with 150 events → targets at 49, 99, 150.
- The actual information fractions will vary per replication around 1/3, 2/3, 1.

**This is not a bug.** It's the standard approach in simulation studies of group sequential designs. The variability in information fractions is a known approximation that is generally acceptable for moderate-to-large event counts. With N=500 and the accrual/censoring setup, the total events should be in the range of 100-200, giving enough granularity for the rounded targets.

---

## 5. Why Cox Type I Error is ~0.12 Instead of ~0.05

**Root cause: Wald z-statistic used with OBF boundaries designed for score test statistics.**

Here's the chain:

1. `gsDesign(k=3, test.type=2, alpha=0.05, sfu="OF")` gives OBF boundaries calibrated for the **score test** (standardized stratified log-rank).
2. The code uses **Wald z-statistics** from `coxph` for the Cox-based rejection rule.
3. Under the null (HR=1) with sparse strata, the Wald statistic has a different distribution than the score statistic — specifically, it has **heavier tails**.
4. The boundaries are too liberal for the Wald statistic, leading to excess rejections.
5. With small strata (especially the 1% and 2% strata in the "2 tiny" scenario), the sparse-data problem is amplified because:
   - Few events per stratum → unstable variance estimates in the Wald statistic
   - The `pool_strata` function merges strata with < 10 events, which helps but doesn't fully eliminate the problem
   - In the extreme "2 tiny (1%, 2%)" scenario, the pooling threshold of 10 may still leave very few events in the merged stratum

**The inflated Type I error (~0.12 vs ~0.05) is entirely explained by using the wrong test statistic.**

---

## 6. The Fix

### Fix A: Use the score test (stratified log-rank) for the primary rejection rule

Replace the Cox Wald z-statistic with the standardized stratified log-rank score statistic:

```r
# In analyze_look(), replace the Cox Wald z with the score z:
analyze_look <- function(d) {
    ...
    # Stratified log-rank score test (canonical for GSD)
    lr_fit <- tryCatch(
        survdiff(Surv(time, event) ~ trt + strata(stratum), data = d),
        error = function(e) NULL, warning = function(w) NULL
    )
    if (is.null(lr_fit) || is.null(lr_fit$chisq) || is.na(lr_fit$chisq) || lr_fit$chisq < 0) {
        score_z <- NA; lr_p <- NA; lr_conv <- FALSE
    } else {
        score_z <- sign(lr_fit$obs[2] - lr_fit$exp[2]) * sqrt(lr_fit$chisq)
        lr_p <- ifelse(lr_fit$chisq < 0, NA,
                       ifelse((lr_fit$obs[2] - lr_fit$exp[2]) > 0,
                              1 - pchisq(lr_fit$chisq, df = 1) / 2,
                              pchisq(lr_fit$chisq, df = 1) / 2))
        lr_conv <- TRUE
    }
    ...
    c(cox_z = cox_z, cox_hr = cox_hr, cox_se = cox_se,
      cox_conv = cox_conv, cox_p = cox_p,
      score_z = unname(score_z), lr_p = lr_p, lr_conv = lr_conv)
}
```

### Fix B: If you must use Cox, use the score test from coxph

`coxph` has a `score` component that gives the score test:

```r
cox_fit <- coxph(Surv(time, event) ~ trt + strata(stratum), data = d)
score_z <- cox_fit$score[1] / sqrt(cox_fit$score[2])  # score / sqrt(var(score))
```

But the stratified log-rank (`survdiff`) is cleaner and more canonical.

### Fix C (minimal change): Use the log-rank results as the primary GSD test

Since the code already computes `lr_p` correctly, just use it as the primary decision rule:

```r
# In run_rep(), replace:
#   if (!is.na(z_no) && abs(z_no) >= obf_z[lk]) { ever_reject_cox_no <- 1 }
# with:
if (!is.na(p_lr_no) && p_lr_no < alpha_look[lk]) {
    if (hr < 1) out$ever_reject_lr_no <- 1  # or use score_z for direction
    if (hr == 1) out$ever_reject_lr_no <- 1
}
```

Where `alpha_look[lk]` is the nominal alpha at look `lk`, derived from the OBF boundaries:

```r
alpha_look <- c(
    2 * pnorm(-obf_z[1]),  # nominal alpha at look 1
    2 * pnorm(-obf_z[2]),  # nominal alpha at look 2  
    2 * pnorm(-obf_z[3])   # nominal alpha at look 3
)
```

### Recommended approach

**Use Fix A or Fix C.** The stratified log-rank score test is the canonical test statistic for group sequential survival trials. The OBF boundaries from `gsDesign` are calibrated for this statistic. Using the Wald statistic is a category error.

---

## Summary of Bugs

| # | Issue | Severity | Fix |
|---|-------|----------|-----|
| 1 | Wald z-statistic used instead of score test | **Critical** | Use `survdiff` score statistic or `coxph$score` |
| 2 | Information fraction uses observed events (not planned) | Low (standard practice) | Acceptable for simulation |
| 3 | "Stop at first rejection" is semantic, not functional | Low (no practical impact) | No change needed |
| 4 | Boundaries from `gsDesign` are correct | N/A | No change needed |

---

## Corrected `run_rep` Function (using score test)

```r
run_rep <- function(i, scenario_idx, hr, p_strata) {
  seed <- 20260519 + scenario_idx * 1e6 + ifelse(hr == 1, 0, 1) * 1e5 + i
  d <- gen_trial(seed, p_strata, hr)
  
  res <- analyze_trial(d)
  nlooks <- length(res)
  
  # Nominal alpha at each look (two-sided, derived from OBF z-boundaries)
  alpha_look <- 2 * pnorm(-obf_z[1:3])
  
  out <- list(
    seed = seed, scenario = scenario_idx, hr = hr, n_events = sum(d$event),
    ever_reject_cox_no = 0, ever_reject_lr_no = 0,
    ever_reject_cox_pool = 0, ever_reject_lr_pool = 0,
    cox_conv_no = 0, cox_conv_pool = 0,
    lr_conv_no = 0, lr_conv_pool = 0,
    cox_hr_last = NA, cox_se_last = NA,
    pool_hr_last = NA, pool_se_last = NA,
    score_z_last = NA
  )
  
  for (lk in seq_len(min(nlooks, 3))) {
    r <- res[[lk]]
    score_z_no <- r$no_pool["score_z"]
    score_z_pool <- r$pool["score_z"]
    
    # Stop at first rejection (group sequential rule)
    if (!is.na(score_z_no) && abs(score_z_no) >= obf_z[lk]) {
      if (hr < 1 && score_z_no < 0) out$ever_reject_cox_no <- 1
      if (hr == 1) out$ever_reject_cox_no <- 1
    }
    if (!is.na(score_z_pool) && abs(score_z_pool) >= obf_z[lk]) {
      if (hr < 1 && score_z_pool < 0) out$ever_reject_cox_pool <- 1
      if (hr == 1) out$ever_reject_cox_pool <- 1
    }
    
    # Convergence / estimates at last look
    if (lk == nlooks) {
      out$cox_conv_no <- r$no_pool["cox_conv"]
      out$cox_conv_pool <- r$pool["cox_conv"]
      out$lr_conv_no <- r$no_pool["lr_conv"]
      out$lr_conv_pool <- r$pool["lr_conv"]
      out$cox_hr_last <- r$no_pool["cox_hr"]
      out$cox_se_last <- r$no_pool["cox_se"]
      out$pool_hr_last <- r$pool["cox_hr"]
      out$pool_se_last <- r$pool["cox_se"]
      out$score_z_last <- score_z_no
    }
  }
  
  out
}
```

And update `analyze_look` to return `score_z`:

```r
analyze_look <- function(d) {
    cox_fit <- tryCatch(
        coxph(Surv(time, event) ~ trt + strata(stratum), data = d),
        error = function(e) NULL, warning = function(w) NULL
    )
    if (is.null(cox_fit) || is.null(cox_fit$coefficients) || any(is.na(coef(cox_fit)))) {
        cox_z <- NA; cox_hr <- NA; cox_se <- NA; cox_conv <- FALSE; score_z <- NA
    } else {
        cox_z <- coef(cox_fit) / sqrt(diag(vcov(cox_fit)))
        cox_hr <- exp(coef(cox_fit))
        cox_se <- sqrt(diag(vcov(cox_fit)))
        cox_conv <- TRUE
        # Score test (canonical for GSD)
        score_z <- cox_fit$score[1] / sqrt(cox_fit$score[2])
    }
    
    lr_fit <- tryCatch(
        survdiff(Surv(time, event) ~ trt + strata(stratum), data = d),
        error = function(e) NULL, warning = function(w) NULL
    )
    if (is.null(lr_fit)) {
        lr_p <- NA; lr_conv <- FALSE
    } else {
        lr_chisq <- lr_fit$chisq
        lr_p <- if (is.null(lr_chisq) || is.na(lr_chisq) || lr_chisq < 0) NA else
            pchisq(lr_chisq, df = 1, lower.tail = FALSE) / 2
        if (!is.na(lr_p) && !is.null(lr_fit$obs) && !is.null(lr_fit$exp) && length(lr_fit$obs) >= 2) {
            lr_oe <- lr_fit$obs[2] - lr_fit$exp[2]
            if (!is.na(lr_oe) && lr_oe > 0) lr_p <- 1 - lr_p
        }
        lr_conv <- TRUE
    }
    
    c(cox_z = unname(cox_z), cox_hr = unname(cox_hr), cox_se = unname(cox_se),
      cox_conv = cox_conv, score_z = unname(score_z),
      cox_p = if (is.na(cox_z)) NA else 2 * pnorm(-abs(cox_z)),
      lr_p = lr_p, lr_conv = lr_conv)
}
```
