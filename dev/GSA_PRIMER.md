# Global Sensitivity Analysis for PBPK Models — Methods, Presentation, Interpretation & Pitfalls

*A methods primer for the `ospgsa` package.*

---

## 1. Why GSA (and not local one-at-a-time) for correlated PBPK parameters

**The question, framed.** Why is a *global* sensitivity analysis necessary instead of the local one-at-a-time (OAT) sensitivity that PK-Sim/`ospsuite` already report, and — given that physiological PBPK parameters are mutually correlated — how should the analysis account for that correlation? This section answers the first half; Section 2 answers the second.

**What local OAT actually computes.** The built-in `ospsuite::SensitivityAnalysis` produces a *normalized local finite-difference* sensitivity of each scalar PK metric to each parameter,

$$
S_{ij} \;=\; \frac{\Delta PK_j / PK_j}{\Delta p_i / p_i} \;\approx\; \frac{\partial PK_j}{\partial p_i}\cdot\frac{p_i}{PK_j},
$$

evaluated by perturbing one parameter at a time around a single reference point (default ±10%, `variationRange = 0.1`, `numberOfSteps = 2`). This is fast and interpretable but has four structural limitations:

1. **It is local.** It characterizes the model only in a small neighborhood of one nominal parameter vector. PBPK responses are nonlinear (saturable Michaelis–Menten clearance, flow- vs capacity-limited regimes, absorption thresholds); a derivative at the nominal point does not describe behavior across the plausible parameter space.
2. **It is one-at-a-time.** OAT holds all other parameters fixed, so it cannot detect or quantify **interactions**. A parameter that matters only through an interaction (large total effect, negligible main effect) is invisible to OAT.
3. **It has no probability model.** OAT uses a fixed ± fraction, not the parameters' distributions, so it cannot apportion *output uncertainty* to inputs — which is precisely the question for uncertainty quantification and for deciding which parameters to estimate vs fix.
4. **It ignores the geometry of the explored region.** OAT samples points along the coordinate axes only, covering a vanishingly small, non-representative slice of a high-dimensional space.

**What GSA adds.** Global methods vary all parameters simultaneously across their full distributions and decompose the resulting output variability. They (a) explore the whole input space, (b) capture nonlinearity and interactions, (c) are defined relative to an explicit input distribution so they answer the uncertainty-apportionment question, and (d) support principled factor fixing (a parameter is fixable only if its **total** effect is negligible). The standard pharmacometric workflow is therefore *Morris screen → variance-based (or moment-independent) quantification on the survivors* (McNally, Cotton & Loizou 2011; Hsieh et al. 2018).

**Why correlation forces the method choice.** Physiology is correlated by construction: organ volumes and blood flows scale with body weight, regional flows sum to cardiac output, enzyme abundances co-vary, and CL/V are linked. The classical variance-based decomposition (Sobol) and the Fourier methods (FAST/eFAST) and Morris all **assume input independence**. Applied to a correlated sample they misattribute variance: first-order indices no longer sum to ≤1 (∑Sᵢ can exceed 1), total indices can fall below the total variance, and a non-influential parameter can look important purely because it is correlated with an influential one. The PBPK case study of Loizou et al. found exactly this — standard Sobol/Morris *overestimated* Vss influence and *underestimated* liver-volume and CYP-clearance effects because correlations were unmodeled. The right approach is therefore not merely "a GSA" but a **correlation-aware** GSA that separates a parameter's own (structural) effect from the effect it inherits through correlation. That separation is the δ two-stage method of Section 2.

---

## 2. The δ (delta) moment-independent index and the two-stage approach for correlated inputs

This is the methodological core of the package and the answer to the correlation half of the question. The method is the **two-stage δ analysis** of De Carlo, Tosca, Melillo & Magni (2023), applied to PBPK by Cuquerella-Gilabert et al. (2026).

### 2.1 The δ index — definition

For output $Y=f(X_1,\dots,X_k)$, let $f_Y(y)$ be the unconditional output density and $f_{Y\mid X_i}(y)$ the density of $Y$ conditional on a fixed value of $X_i$. The **shift** at a fixed $X_i=x_i$ is the area between the two densities,

$$
s(X_i) \;=\; \int \bigl|\, f_Y(y) - f_{Y\mid X_i}(y)\,\bigr|\, dy
$$

(twice the total-variation distance), and Borgonovo's δ is the half of its expectation over the distribution of $X_i$:

$$
\boxed{\;\delta_i \;=\; \tfrac{1}{2}\, \mathbb{E}_{X_i}\!\bigl[\, s(X_i)\,\bigr]
\;=\; \int f_{X_i}(x_i)\!\left(\tfrac{1}{2}\int \bigl| f_Y(y) - f_{Y\mid X_i}(y)\bigr|\,dy\right) dx_i\;}
$$

An equivalent and estimation-friendly **joint-density form** is half the $L^1$ distance between the joint density and the product of marginals:

$$
\delta_i \;=\; \tfrac12 \iint \bigl|\, f_{X_i}(x_i) f_Y(y) - f_{X_i,Y}(x_i,y)\,\bigr|\, dy\,dx_i .
$$

**Properties** (Borgonovo 2007):

- **Range:** $0 \le \delta_i \le 1$ (the ½ normalizes the doubled-area $L^1$ distance, which is at most 2).
- **Moment-independent:** δ responds to changes in *any* feature of the output distribution — shape, skew, tails, multimodality — not only variance. A factor can have Sobol $S_i\approx 0$ yet $\delta_i>0$ if it reshapes the distribution without changing variance. This makes δ well suited to the skewed/heavy-tailed PK outputs (AUC, Cmax) and to eradication/trade-off metrics that variance summarizes poorly.
- **Independence characterization:** $\delta_i = 0 \iff Y$ is independent of $X_i$ **and** $X_i$ is uncorrelated with the other inputs. This exact property is what makes the two-stage decomposition valid.
- **No input-independence assumption:** δ uses only the bivariate law of $(X_i,Y)$, so it is *computable from a correlated sample with no modification* — unlike the Sobol ANOVA terms, which are orthogonal only under independence.
- **Caution — not additive:** the $\delta_i$ do not sum to 1 and do not split cleanly into main/interaction terms. They rank and attribute; they do not partition variance.

### 2.2 The two-stage procedure (correlated vs structural separation)

Because δ on a correlated sample blends a parameter's *structural/causal* effect (through the model $f$) with its *correlation-induced/indirect* effect (inherited from correlated partners), the method computes δ on **two designs** and compares.

- **Stage 1 — independence.** Sample each parameter from its **marginal** with all correlations switched off, and compute $\delta_{1,i}$. Removing correlation isolates the **direct/causal** contribution: a non-zero $\delta_{1,i}$ certifies a genuine structural effect.
- **Stage 2 — full correlated joint.** Sample from the realistic correlated joint distribution and compute $\delta_{2,i}$. These indices reflect importance under the *true* (correlated) output distribution and capture indirect, correlation-transmitted effects.

**Classification rule** (the mechanism that separates correlated from structural contributions):

| Category | Condition | Reading |
|---|---|---|
| Causal / structural | $\delta_{1,i} > 0$ | direct model-driven effect |
| Indirect / correlation-driven only | $\delta_{1,i}=0$ **and** $\delta_{2,i}>0$ | influential *only* via correlation with an influential parameter |
| Both | $\delta_{1,i}>0$ **and** $\delta_{2,i}>0$ | direct effect plus correlation transmission |

The structural part is what survives Stage 1; the correlative part is what is gained (or lost) moving Stage 1 → Stage 2. The canonical example (De Carlo 2023): the potency-related `IC50` has $\delta_1=0$ but $\delta_2>0$ because it is strongly correlated with `k2`, which has a strong direct effect on the threshold concentration $C_T$. The actionable reading: *estimating `k2` accurately is not enough — you must also reduce its correlation with the parameters that ride on it.* This converts the analysis into experimental-design guidance: reduce a parameter's own estimation uncertainty only if it is causal ($\delta_1>0$); if it matters only indirectly, the leverage is in breaking the correlation.

> **An equivalent "full vs independent" framing** (Mara & Tarantola 2012, transposed to δ): define $\delta_i^{\text{full}}$ on the dependent sample (total apparent influence) and $\delta_i^{\text{ind}}$ on a sample where $X_i$ has been orthogonalized/decorrelated from the others via a conditional (Rosenblatt) transform (structural-only). Then $\delta_i^{\text{corr}} = \delta_i^{\text{full}} - \delta_i^{\text{ind}}$. This is an attribution heuristic, not an exact decomposition (δ is non-additive), and the orthogonalization order can affect the split. The De Carlo Stage-1/Stage-2 construction is the practical realization used by the tool.

### 2.3 The given-data estimator (implementation-precise)

δ is estimated from a single input–output sample $\{(X_i^{(n)}, Y^{(n)})\}_{n=1}^N$ — no extra model runs (Plischke, Borgonovo & Smith 2013, "given-data" estimator). The intractable point-conditional density is replaced by conditioning on $X_i$ falling in a **class** (partition cell):

1. **Partition** the range of $X_i$ into $M$ disjoint classes on the **ranks** of $X_i$ (equal-frequency / equiprobable binning), so each class holds ≈ $N/M$ points.
2. **Estimate densities by KDE** (Gaussian kernel, Silverman bandwidth) on a common output grid: the unconditional $\hat f_Y$ from all $Y$, and each class-conditional $\hat f_{Y\mid X_i\in\mathcal M_m}$ from the $Y$ in class $m$.
3. **Aggregate** (Plischke Eqn. 26):

$$
\hat\delta_i \;=\; \frac{1}{2}\sum_{m=1}^{M} \frac{N_m}{N}\, \int \bigl|\, \hat f_Y(y) - \hat f_{Y\mid X_i\in\mathcal{M}_m}(y)\,\bigr|\, dy,
\qquad \sum_m N_m = N.
$$

**Number of classes** grows sub-linearly with $N$ ($M\propto N^{\alpha}$, $\alpha\approx\tfrac12$–$\tfrac23$). SALib's calibration (a faithful operational realization of the paper) is

```
exp = 2.0 / (7.0 + tanh((1500 - N) / 500))   # ≈0.25 small N, ≈0.33 large N
M   = min(ceil(N**exp), 48)                  # capped at 48 classes
```

**Bias correction.** The class-conditioning and finite-$N$ KDE make $\hat\delta_i$ upward-biased ($\mathbb{E}[\hat\delta_i]>0$ even when $\delta_i=0$). Apply the bootstrap bias reduction (Plischke Eqn. 30):

$$
\tilde\delta_i \;=\; 2\,\hat\delta_i \;-\; \bar\delta_{\text{boot}},
$$

i.e. `d = 2.0 * d_hat - d.mean()` over the bootstrap replicates.

**Bootstrap confidence intervals.** δ has no closed-form variance; resample the index set with replacement $B$ times, recompute the (bias-corrected) $\hat\delta_i^{*(b)}$, and report either the **normal half-width** $z_{(1+\gamma)/2}\cdot\mathrm{sd}^*(\hat\delta_i^{*})$ or, preferably for a bounded skewed index, the **percentile** interval $[\text{quantile}_{(1-\gamma)/2},\,\text{quantile}_{(1+\gamma)/2}]$. De Carlo (2023) and Cuquerella-Gilabert (2026) both use **1000 bootstrap resamples**.

**Sampling and distributions used in the source papers.** Inputs modeled as **multivariate lognormal** $\mathrm{LogN}(\theta,\Omega)$ (PK parameters are positive, right-skewed). **N = 100,000** Monte-Carlo draws from the joint, used in both stages; Stage 1 zeroes the off-diagonals of $\Omega$, Stage 2 uses the full $\Omega$. Unknown correlations are filled to a valid positive-definite structure with the R package **`mvLognCorrEst`**; pairs are tagged as correlated ($\rho=\pm0.95, 0.75, 0.6$), uncorrelated, or unknown.

> **R note for the tool.** There is no off-the-shelf R function implementing the exact Plischke class-conditioning + bias correction. `sensitivity::sensiFdiv(fdiv = "TV")` gives Borgonovo's δ up to a constant (total-variation f-divergence, KDE-based, bootstrap CIs via `nboot`/`conf`), but it conditions via continuous KDE divergence rather than equal-frequency classes. For a faithful given-data δ, **port the SALib algorithm** (it is short, pure-numeric); the pseudocode below is the port target.

```
# given-data δ (port target)
calc_delta(Y, X, ygrid, m):                 # m = M+1 rank-bin boundaries
    fY  <- kde(Y, "silverman") on ygrid
    rX  <- rank(X, ties="ordinal")
    dhat <- 0
    for j in 1..(len(m)-1):
        idx <- which(rX > m[j] & rX <= m[j+1])
        if length(idx) < 2: next
        fYc <- kde(Y[idx], "silverman") on ygrid
        dhat <- dhat + (length(idx)/(2*N)) * trapz(abs(fY - fYc), ygrid)
    return dhat

driver(Xi, Y, B, gamma):
    N   <- length(Y); exp <- 2/(7 + tanh((1500-N)/500))
    M   <- min(ceil(N^exp), 48); m <- linspace(0, N, M+1)
    yg  <- linspace(min(Y), max(Y), 100)
    dh  <- calc_delta(Y, Xi, yg, m)
    for b in 1..B: db[b] <- calc_delta(Y[r_b], Xi[r_b], yg, m)   # r_b = resample idx
    dc  <- 2*dh - db                                              # bias reduction
    return mean(dc), qnorm(0.5+gamma/2)*sd(dc)                    # or percentile CI

# Two stages: run driver on the correlated sample (δ2) and on the
# independence sample / orthogonalized Xi (δ1); classify per the table above.
```

The given-data estimator also yields the first-order Sobol $S_1$ for free from the same partition (variance of class-conditional means / total variance), useful as a cross-check.

---

## 3. Companion methods table

All R functions are in the CRAN `sensitivity` package unless noted; `sensobol` is the dedicated Sobol package. "Correlation-safe?" means the index remains interpretable under dependent inputs.

| Method | What it measures | Correlation-safe? | Typical sample size | R implementation | When to use |
|---|---|---|---|---|---|
| **Morris (elementary effects)** | μ\* = overall importance (≈ proxy for total effect); σ = nonlinearity + interactions | **No** (assumes independent box) | $r(k+1)$, $r\approx10$–50 | `morris`, `morrisMultOut` | Cheapest first-pass **screening** to drop clearly inactive factors |
| **Sobol $S_1$ / $S_T$** | $S_1$ = variance share if $X_i$ fixed (prioritization); $S_T$ = main + all interactions (fixing) | **No** (∑Sᵢ breaks under correlation) | $N(k+2)$, $N=2^{10}$–$2^{14}$ (power of 2) | `sensobol::sobol_indices`; `soboljansen`, `sobol2007`, `sobolEff`, `sobolSalt` | Quantitative variance decomposition with **independent** inputs |
| **FAST / eFAST** | Fourier variance at factor frequencies: $S_1$ (FAST), $S_1$+$S_T$ (eFAST) | **No** | $n\cdot k$; eFAST $n\gtrsim65$/factor; RBD-FAST from $N\gtrsim512$ | `fast99` (eFAST), `fast::sensitivity` (FAST) | Efficient variance indices, many factors, independent inputs; eFAST is the PBPK reliability/efficiency sweet spot |
| **SRC / SRRC** | Standardized (rank) regression coefficients; share of variance under (rank-)linear fit | **No** (raw shares break under collinearity) | any MC/LHS, $N\gtrsim10k$ | `src(..., rank=)` | Cheap; **only if model $R^2$ high** and (monotone-)linear |
| **PRCC / PCC** | Partial (rank) correlation of $X_i$ with $Y$, controlling for the others | **Partly — for *linear* correlation** | any MC/LHS | `pcc(..., rank=, semi=)`; `epiR::epi.prcc` | Monotone models with **linearly** correlated inputs |
| **Shapley effects** | Fair allocation of variance (incl. correlation + interaction) among inputs; ≥0, **sum to 1** | **Yes — designed for it** | permutation MC; expensive in $k$ | `shapleyPermEx`, `shapleyPermRand`, `shapleySubsetMc`, `shapleyLinearGaussian` | Proper **variance attribution** under dependence |
| **HSIC (R²-HSIC)** | Kernel dependence; any (incl. non-monotonic) statistical dependence; = 0 ⇔ independent | **Yes** (except `anova` mode → No) | any MC/LHS, a few hundred runs | `sensiHSIC`, `testHSIC` (screening p-values) | General dependence/screening, incl. correlated and non-monotonic |
| **δ (Borgonovo, moment-independent)** | Whole-distribution shift when conditioning on $X_i$; two-stage → structural vs correlative | **Yes** (the tool's primary method) | given-data; several hundred to $10^5$ | `sensiFdiv(fdiv="TV")`; SALib-port given-data estimator | **Correlated** PBPK inputs, skewed/multimodal outputs, structural-vs-correlative split |

**Workflow rule of thumb:** Morris screen → drop inactive factors → on survivors run Sobol/eFAST if inputs are independent, or **Shapley / HSIC / δ** if dependent. For this tool the default on correlated PBPK parameters is the **two-stage δ**; Sobol/eFAST/Morris are offered for independent-input cases and cross-checks.

---

## 4. Result presentation

### 4.1 Tables

Report **one table per output metric** (e.g. one for $\log C_{max}$, one for $\log\mathrm{AUC}$). The tidy estimator output (`sensobol::sobol_indices`, and the analogous δ output) carries: `parameters`, `sensitivity` (index type), `original` (point estimate), `bias`, `std.error`, `low.ci`, `high.ci`.

| Parameter | $S_1$ (95% CI) | $S_{T}$ (95% CI) | $S_T-S_1$ (interaction) | Rank ($S_T$) | Influential? |
|---|---|---|---|---|---|
| CL | 0.41 (0.37–0.45) | 0.49 (0.44–0.55) | 0.08 | 1 | yes |
| V1 | 0.22 (0.19–0.26) | 0.31 (0.27–0.36) | 0.09 | 2 | yes |
| ka | 0.02 (0.00–0.05) | 0.06 (0.03–0.10) | 0.04 | 5 | borderline |
| dummy | 0.00 (−0.01–0.02) | 0.01 (0.00–0.03) | — | — | noise floor |

For the **two-stage δ**, report $\delta_1$ and $\delta_2$ side by side with the classification column:

| Parameter | $\delta_1$ (95% CI) | $\delta_2$ (95% CI) | Effect type |
|---|---|---|---|
| k2 | 0.34 (0.30–0.38) | 0.36 (0.32–0.40) | causal (+ slight indirect) |
| IC50 | 0.00 (0.00–0.02) | 0.21 (0.17–0.25) | indirect only (corr. with k2) |
| Vu10 | 0.05 (0.02–0.08) | 0.04 (0.01–0.07) | causal |

Conventions: always pair every index with a **bootstrap CI**; include a **dummy/phantom-parameter row** as the empirical noise floor (any real parameter whose CI overlaps the dummy's is statistically indistinguishable from zero); **rank by the total/full index** ($S_T$ or $\delta_2$); for time-resolved analyses, summarize each parameter by the **maximum index across time points**. Add a **computational/convergence table** (parameter count, $N$, runtime per method, #influential at each cutoff, max CI half-width).

### 4.2 Key plot types (one-line ggplot recipes)

**(a) Tornado / horizontal bar with CIs** — ranked importance.
```r
ggplot(df, aes(reorder(parameters, original), original)) + geom_col() +
  geom_errorbar(aes(ymin = low.ci, ymax = high.ci), width = .25) +
  geom_hline(yintercept = 0.05, linetype = 2, colour = "red") + coord_flip() + theme_bw()
```

**(b) Grouped $S_1$ vs $S_T$ bars** — the single most informative GSA plot; the gap = interactions.
```r
ggplot(long, aes(reorder(parameters, original), original, fill = sensitivity)) +
  geom_col(position = position_dodge(.6), width = .55) +
  geom_errorbar(aes(ymin = low.ci, ymax = high.ci), position = position_dodge(.6), width = .2) + theme_bw()
```

**(c) Two-stage δ split** — $\delta_1$ (structural) vs $\delta_2$ (full), or stacked structural + correlative per parameter.
```r
ggplot(delta_long, aes(reorder(parameters, d2), value, fill = stage)) +
  geom_col(position = position_dodge(.6)) +
  geom_errorbar(aes(ymin = low.ci, ymax = high.ci), position = position_dodge(.6), width = .2) +
  coord_flip() + theme_bw()
```

**(d) Stacked variance decomposition** — one bar per output, $S_i$ plus a $1-\sum S_i$ interaction remainder.
```r
ggplot(decomp, aes(output, Si, fill = parameters)) + geom_col(position = "stack") + ylim(0,1) + theme_bw()
```

**(e) Morris (μ\*, σ) plane** — screening.
```r
ggplot(morris_df, aes(mu.star, sigma, label = parameters)) + geom_point() +
  geom_abline(slope = 1, linetype = 3) + geom_text(vjust = -.6) + theme_bw()
```

**(f) Time-resolved heatmap (parameters × time)** — how drivers shift over the C–t curve.
```r
ggplot(tr, aes(time, reorder(parameters, value), fill = value)) + geom_tile() +
  scale_fill_viridis_c() + theme_minimal()
```
or `pksensi::heat_check(out, order = "total order")`.

**(g) Stacked-area $S_i(t)$** — absorption vs elimination drivers, concentration overlaid.
```r
ggplot(tr, aes(time, Si, fill = parameters)) + geom_area(position = "fill") + theme_bw()
```

**(h) Convergence vs N** — index with CI ribbon against sample size.
```r
ggplot(conv, aes(N, original, colour = parameters)) + geom_line() +
  geom_ribbon(aes(ymin = low.ci, ymax = high.ci, fill = parameters), alpha = .15, colour = NA) +
  scale_x_log10() + theme_bw()
```

**(i) Scatter / cobweb + output distribution** — raw sanity check before any index.
```r
sensobol::plot_scatter(N = 2000, data = mat, Y = y, params = colnames(mat))
sensobol::plot_uncertainty(Y = y, N = N) + scale_x_log10()
```

### 4.3 Interpretation rules and thresholds

- **"Influential" cutoffs.** Variance-based: $0.05$ is the common cutoff (Zhang et al.; "may not be stringent enough for simple models"), $0.01$ a more inclusive screen. Morris: normalized $\mu^*$ or $\sigma > 0.1$. δ: no universal cutoff — use the dummy floor or a small value (~0.01–0.05) plus CI overlap with the dummy. **Best practice over any fixed cutoff:** use a **dummy parameter** as the empirical significance floor.
- **$S_1$ vs $S_T$.** $S_1\approx S_T$ ⇒ additive (no interactions); $\sum S_i\approx1$ confirms additivity. $S_T\gg S_1$ ⇒ the parameter acts mainly through interactions — do **not** fix it on a small $S_1$. $S_T\approx0$ ⇒ non-influential and **fixable**. Always rank/fix on **$S_T$** (or $\mu^*$, or $\delta_2$/$\delta^{\text{full}}$), never on $S_1$ alone.
- **Morris σ.** Large σ ⇒ effects vary across the space ⇒ strong nonlinearity and/or interactions; use $\mu^*$ (not μ) so opposite-sign effects of non-monotonic factors do not cancel.
- **δ under correlation.** δ measures the average shift in the *entire* output PDF, so it stays meaningful when variance is a poor summary (skew/multimodality). Read the two-stage split: a parameter important *only* through its correlative term is not an independent driver — its apparent importance is borrowed, which determines whether you can perturb/estimate it independently.

---

## 5. Pitfalls checklist (do / don't)

### 5.1 Log-scale vs linear-scale parameters — the #1 PBPK pitfall

PBPK parameters (clearances, $K_m$, permeabilities, partition coefficients, microsomal $CL_{int}$) routinely span 2–4 orders of magnitude and are conventionally **lognormal**. The sampling scale changes the parameter ensemble and therefore the sensitivity ranking.

- **There is no transform-invariant sensitivity index.** Sobol $S_i/S_T$ *and* δ are defined relative to the input marginal — they answer "how much does $Y$ vary *given the assumed spread of $X_i$*." Change the marginal (linear vs log range, or the lognormal $\sigma$) and the indices change, because you changed the variance/uncertainty budget being apportioned. An index is a property of (model + input distribution), not of the model alone.
- **Sample on the right scale.** For a positive parameter with order-of-magnitude uncertainty use **log-uniform** (uniform on $\log X$) or **lognormal**, not linear-uniform. A linear $U(1,1000)$ puts ~90% of mass above 100 and starves the low-clearance regime; log-uniform gives equal weight per decade. Implementation: draw on the unit hypercube, map through the log-scale inverse-CDF (`qlnorm`, or `a*(b/a)^u` for log-uniform), back-transform **before** calling PK-Sim. For lognormal from mean $m$ and CV $c$: `sdlog = sqrt(log(1+c^2)); meanlog = log(m) - sdlog^2/2`.
- **Interpreting log-sampled indices.** When you sample $\log X_i$, $S_i$/$\delta_i$ measure the contribution of *multiplicative* (proportional) variation — usually the physiologically meaningful question for clearances and partition coefficients. An index for a log-sampled clearance is **not comparable** to one for a linearly sampled body weight; report the sampling scale of every input.
- **Transform the OUTPUT too.** AUC and $C_{max}$ are strongly right-skewed and span decades, so $\mathrm{Var}(Y)$ is dominated by a few tail runs and is unstable. Best practice is to run GSA on **$\log\mathrm{AUC}$ / $\log C_{max}$** (Hsieh et al.: "all model outputs were transformed to logarithmic scale"); this stabilizes variance, behaves better in the Sobol decomposition, and matches the lognormal PK likelihood. Note it reframes the question to "drivers of fold-variation."

**Do:** sample order-of-magnitude parameters log-uniform/lognormal and back-transform before simulating; use $\log\mathrm{AUC}/\log C_{max}$ as outputs; report each input's sampling scale.
**Don't:** put a linear-uniform on a decade-spanning parameter; compare a log-sampled index with a linear-sampled one; treat indices as intrinsic to the model.

### 5.2 Correlated inputs break variance-based indices

The Sobol decomposition and $\sum S_i\le1\le\sum S_T$ hold **only under independence**. On a correlated sample, $\sum S_i$ can exceed 1, $\sum S_T$ can fall below the total variance, and a non-influential parameter can look important via its correlated partner. Acute in PBPK because physiology is correlated by construction.
**Do:** encode physiological correlations (sample body weight, then derive organ volumes/flows; enforce ∑ regional flows = cardiac output; use a copula or a virtual-population covariance); use **two-stage δ**, **Shapley**, or **HSIC/PRCC** when inputs are dependent.
**Don't:** report raw Sobol $S_i/S_T$ from a correlated sample as variance fractions; panic at $\sum S_i>1$ — diagnose it as a correlation (or under-convergence) symptom.

### 5.3 Convergence and insufficient N

GSA indices are Monte-Carlo estimates; under-sampling gives noisy, non-reproducible rankings (and the $\sum S_i>1$ / negative-$S_i$ pathologies). Costs: Sobol $N(k+2)$ ($N(2k+2)$ for second order); Morris $r(k+1)$; eFAST $n\cdot k$; δ single-sample (given-data) but needs enough points per class for the conditional KDE.
**Do:** report **bootstrap CIs**; grow $N$ until CIs are tight and rankings stop changing; require **max 95% CI half-width < 0.1** ("convergence index < 0.1"); show a convergence plot (`sensobol::sobol_convergence`); use Sobol quasi-random sequences with $N$ a power of 2.
**Don't:** report point estimates with no CI; assume one big run converged; use second-order Sobol unless you can afford $N(2k+2)$.

### 5.4 Non-monotonicity / nonlinearity

Linear and rank methods silently fail on non-monotone PBPK responses (saturable clearance, U-shaped/threshold effects).
**Do:** report the **SRC $R^2$** (or rank-$R^2$) — SRC is trustworthy only when $R^2$ is high (≳0.7); escalate **SRC → SRRC** (monotone nonlinear) **→ Sobol/δ** (any functional form). δ and Sobol capture non-monotone effects; δ additionally captures distributional change.
**Don't:** quote Pearson/SRC sensitivities without the $R^2$; trust rank methods (SRRC/PRCC) on non-monotone responses.

### 5.5 Output-metric choice

A simulation yields a full C–t curve; the chosen scalar QoI changes the rankings ($C_{max}$ ← absorption/distribution; AUC ← clearance; $T_{max}$ ← rate constants), and even the integration window matters ($AUC_{24}$ vs $AUC_{48}$ vs $AUC_\infty$).
**Do:** pre-specify QoIs and windows ($\log AUC_{0-\tau}$, $\log C_{max}$, $T_{max}$); analyze each separately (reuse the same base sample across QoIs); consider time-resolved indices; for zero-inflated outputs (BLQ, non-absorbing regimes) prefer δ (more robust than variance when a spike at zero dominates $\mathrm{Var}(Y)$).
**Don't:** pick one metric and generalize; ignore zero-inflated outputs.

### 5.6 Failed simulations / unrealistic ranges (NA handling)

The Saltelli pick-freeze design **requires a complete matrix**: a single NaN in an $A/B/A_B^{(i)}$ row corrupts the paired differences and biases multiple indices. Wide/non-physiological bounds push the ODE solver into stiff/crash regions and inflate apparent sensitivity.
**Do:** derive ranges from data/virtual populations and **truncate** (e.g. $\mu_{\log}\pm z\sigma_{\log}$); ensure prior predictions cover the calibration data; **log every failed run and its parameter vector** and report the failure fraction (a high rate is itself a finding); prefer **given-data estimators** (δ, rank/regression) that tolerate dropped rows, or re-draw failed Sobol points.
**Don't:** silently `na.omit` rows out of a Saltelli matrix (breaks its balance and **underestimates** total indices); let NaNs propagate into variance terms.

### 5.7 Factor fixing — and over-fixing

The primary PBPK use of GSA is factor fixing (fix non-influential parameters, estimate the rest). **Fix on the total index $S_T$ (or Morris $\mu^*$, or $\delta_2$), never on $S_i$** — a parameter can have $S_i\approx0$ yet large $S_T$ through interactions. Cutoffs in practice: fix if Sobol total/main < 0.01–0.05, or normalized Morris < 0.1.
**Do:** fix on $S_T/\mu^*/\delta_2$ aggregated as the **max over outputs and time points**, using the **conservative CI bound** (upper bound vs the threshold), and document fixed values; in the two-stage δ, a parameter important only via correlation ($\delta_1=0,\delta_2>0$) is a candidate to **fix** rather than estimate independently.
**Don't:** fix on $S_1$; fix near a threshold when the CI is wide; fix using a single output or a single time point; fix under a too-narrow range or wrong correlation structure (the PBPK study showed correlations caused both over- and under-fixing of enzyme-clearance terms).

*Cross-cutting note on normalization:* variance indices are dimensionless fractions, comparable across parameters but silent about the magnitude of $\mathrm{Var}(Y)$ — report $\mathrm{Var}(\log Y)$ and the mean for context; normalize Morris per output before comparing across QoIs; never equate a **local** elasticity ($\partial\ln Y/\partial\ln X$) with a **global** index.

---

## References

**Primary δ method and PBPK application**
- De Carlo, A., Tosca, E. M., Melillo, N., & Magni, P. (2023). A two-stages global sensitivity analysis by using the δ sensitivity index in presence of correlated inputs: application on a tumor growth inhibition model based on the dynamic energy budget theory. *J Pharmacokinet Pharmacodyn* 50(5):395–409. DOI: https://doi.org/10.1007/s10928-023-09872-w · PMID 37422844 · Open access: https://pmc.ncbi.nlm.nih.gov/articles/PMC10460734/
- Cuquerella-Gilabert, M., De Carlo, A., Sánchez-Herrero, S., Reig-López, J., Merino-Sanjuán, M., Tosca, E. M., Mangas-Sanjuán, V., & Magni, P. (2026). Leveraging Two-Stage δ Global Sensibility Analysis Method to Inform Parameter Estimation in PBPK Models. *Pharmaceutical Statistics* 25(2):e70082. DOI: https://doi.org/10.1002/pst.70082

**δ index theory and estimation**
- Borgonovo, E. (2007). A new uncertainty importance measure. *Reliability Engineering & System Safety* 92(6):771–784. DOI: https://doi.org/10.1016/j.ress.2006.04.015
- Plischke, E., Borgonovo, E., & Smith, C. L. (2013). Global sensitivity measures from given data. *European Journal of Operational Research* 226(3):536–550. DOI: https://doi.org/10.1016/j.ejor.2012.11.047
- Borgonovo, E., & Plischke, E. (2016). A common rationale for global sensitivity measures and their estimation. *Risk Analysis*. DOI: https://doi.org/10.1111/risa.12555
- SALib delta estimator (operational reference): https://salib.readthedocs.io/en/latest/_modules/SALib/analyze/delta.html
- R `sensitivity::sensiFdiv` (TV f-divergence = δ up to a constant): https://search.r-project.org/CRAN/refmans/sensitivity/html/sensiFdiv.html

**Correlated-input decomposition**
- Mara, T. A., & Tarantola, S. (2012). Variance-based sensitivity indices for models with dependent inputs. *Reliability Engineering & System Safety* 107:115–121.
- Mara, T. A., Tarantola, S., & Annoni, P. (2015). Non-parametric methods for global sensitivity analysis of model output with dependent inputs. *Environmental Modelling & Software*. https://par.nsf.gov/servlets/purl/10191596
- Iooss, B., & Prieur, C. (2019). Shapley effects for sensitivity analysis with correlated inputs. *Int. J. Uncertainty Quantification* 9(5). https://arxiv.org/abs/1707.01334

**Companion GSA methods**
- Saltelli, A., et al. (2010). Variance based sensitivity analysis of model output. Design and estimator for the total sensitivity index. *Computer Physics Communications*.
- Puy, A., Lo Piano, S., Saltelli, A., & Levin, S. A. (2022). sensobol: an R package to compute variance-based sensitivity indices. *JSS* 102(5). https://arxiv.org/pdf/2101.10103
- Tarantola, S., Gatelli, D., & Mara, T. A. (2006). Random balance designs for the estimation of first order global sensitivity indices (RBD-FAST). *Reliability Engineering & System Safety*.
- Morris, M. D. (1991). Factorial sampling plans for preliminary computational experiments. Elementary effects: https://en.wikipedia.org/wiki/Elementary_effects_method
- Campolongo, F., Saltelli, A., & Cariboni, J. (2011). From screening to quantitative sensitivity analysis (radial design). *Computer Physics Communications* 182:978.
- Marino, S., et al. (2008). A methodology for performing global uncertainty and sensitivity analysis in systems biology (PRCC). *J Theor Biol*. https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2570191/
- Song, E., Nelson, B. L., & Staum, J. (2016). Shapley effects for global sensitivity analysis: theory and computation.
- Da Veiga, S. (2015). Global sensitivity analysis with dependence measures (HSIC). *J. Stat. Comput. Simul.* 85(7):1283. R: https://search.r-project.org/CRAN/refmans/sensitivity/html/sensiHSIC.html
- CRAN `sensitivity` package manual: https://cran.r-project.org/web/packages/sensitivity/sensitivity.pdf

**Presentation, interpretation, and PBPK GSA practice**
- Hsieh, N.-H., Reisfeld, B., Bois, F. Y., & Chiu, W. A. (2018). Applying a global sensitivity analysis workflow to improve the computational efficiencies in PBPK modeling. *Front. Pharmacol.* https://pmc.ncbi.nlm.nih.gov/articles/PMC6002508/
- Zhang, X.-Y., Trame, M. N., Lesko, L. J., & Schmidt, S. (2015). Sobol sensitivity analysis for systems pharmacology. *CPT:PSP*. https://pmc.ncbi.nlm.nih.gov/articles/PMC5006244/
- McNally, K., Cotton, R., & Loizou, G. D. (2011). A workflow for global sensitivity analysis of PBPK models. *Front. Pharmacol.* https://www.frontiersin.org/articles/10.3389/fphar.2018.00588/full
- Loizou, G. D., et al. Considerations and caveats when applying global sensitivity analysis methods to PBPK models. https://pmc.ncbi.nlm.nih.gov/articles/PMC7367914/
- pksensi APAP-PBPK vignette: https://cran.r-project.org/web/packages/pksensi/vignettes/pbpk_apap.html
- FDA (2018). Physiologically Based Pharmacokinetic Analyses — Format and Content (Guidance for Industry). https://www.fda.gov/files/drugs/published/Physiologically-Based-Pharmacokinetic-Analyses-%E2%80%94-Format-and-Content-Guidance-for-Industry.pdf

**Correlated sampling and existing OSP tooling**
- Iman, R. L., & Conover, W. J. (1982). A distribution-free approach to inducing rank correlation among input variables. *Comm. Stat. – Sim. Comp.* 11(3):311–334.
- R packages: `lhs`, `mc2d::cornode` (Iman–Conover), `copula`, `mvtnorm`, `mvLognCorrEst`.
- Najjar, A., et al. (2024). `ospsuite.globalsensitivity` (Morris/Sobol/eFAST for OSP). *CPT:PSP*. DOI: https://doi.org/10.1002/psp4.13256 · Repo: https://github.com/Open-Systems-Pharmacology/OSPSuite.GlobalSensitivity
- ospsuite (OSP) R API — efficient calculations / SimulationBatch / PK analysis vignettes: https://www.open-systems-pharmacology.org/OSPSuite-R/
