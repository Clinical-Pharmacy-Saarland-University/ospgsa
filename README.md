# ospgsa

> Correlation-aware global sensitivity analysis for Open Systems Pharmacology (PBPK) models.

<!-- badges: start -->

[![R-CMD-check](https://github.com/Clinical-Pharmacy-Saarland-University/ospgsa/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/Clinical-Pharmacy-Saarland-University/ospgsa/actions/workflows/R-CMD-check.yaml)

<!-- badges: end -->

Global sensitivity analysis (GSA) for **Open Systems Pharmacology (OSP)** models.

`ospgsa` performs correlation-aware global sensitivity analysis of OSP models (`.pkml` simulations exported from PK-Sim or MoBi). A _local_, one-at-a-time analysis can be inadequate for PBPK models whose physiological parameters are nonlinear, interacting and mutually correlated.

Its main method is the Borgonovo delta moment-independent index with a two-stage decomposition for correlated inputs (De Carlo et al. 2023, Cuquerella-Gilabert et al. 2026), which separates a parameter's direct (structural) effect from the importance it merely inherits through correlation with other parameters. Sobol, Morris and regression methods are included for cross-checking. The GSA engine itself is model-agnostic. It drives any _evaluator_ that maps an input matrix to an output matrix, with OSP models as the built-in integration.

## Relationship to existing OSP tools

OSP already offers sensitivity analysis, but at a different scope. `ospsuite::SensitivityAnalysis` computes a local, one-at-a-time sensitivity of PK parameters (a normalized finite difference around a single nominal point) and the separate `ospsuite.globalsensitivity` package adds Morris/Sobol/eFAST. `ospgsa` adds a correlation-aware, distribution-based analysis and the delta two-stage method.

|                       | `ospsuite::SensitivityAnalysis` | `ospsuite.globalsensitivity` | `ospgsa`                                               |
| --------------------- | ------------------------------- | ---------------------------- | ------------------------------------------------------ |
| Scope                 | local (one point)               | global                       | global                                                 |
| Parameters varied     | one at a time                   | all jointly                  | all jointly                                            |
| Input distributions   | no (fixed ± %)                  | yes                          | yes (5 marginals)                                      |
| Correlated inputs     | no                              | no                           | yes (delta two-stage, copula or Iman-Conover sampling) |
| Interactions captured | no                              | yes                          | yes                                                    |
| Methods               | local elasticity                | Morris, Sobol, eFAST         | delta (two-stage), Sobol S1/ST, Morris, SRC/PRCC       |
| Outputs               | PK parameters                   | PK parameters                | PK parameters and time-resolved curves                 |

## Correlation-aware two-stage delta

Local one-at-a-time sensitivity describes the model only near one nominal point and cannot capture nonlinearity, interactions or apportion _output_ uncertainty to inputs. Global methods vary all parameters at once across their full distributions and decompose the resulting variability. Because physiological parameters are correlated by construction, variance-based methods such as Sobol and FAST, as well as screening methods such as Morris, violate the common assumption of independent inputs. `ospgsa` defaults to Borgonovo's delta, which is computable from a correlated sample with no modification.

On a correlated sample delta blends a parameter's own _structural_ effect with the importance it _inherits_ through correlation. The two-stage procedure (De Carlo et al. 2023, Cuquerella-Gilabert et al. 2026) computes delta on two designs and compares them.

- **Stage 1, independent design** (correlations switched off). $\delta_1$ isolates the direct, model-driven effect.
- **Stage 2, full correlated design**. $\delta_2$ reflects importance under the true joint distribution, including correlation-transmitted effects.

`delta_classification()` reads the two stages against a `dummy` parameter that fixes the empirical noise floor.

| Class               | Condition                              | Reading                                              |
| ------------------- | -------------------------------------- | ---------------------------------------------------- |
| **causal**          | $\delta_1$ above the floor                | genuine direct, model-driven effect                  |
| **indirect-only**   | $\delta_1$ at the floor, $\delta_2$ above it | influential _only_ through correlation with a driver |
| **both**            | $\delta_1$ and $\delta_2$ above the floor    | direct effect **plus** correlation transmission      |
| **non-influential** | neither above the floor                | fixable                                              |

The [GSA methods primer](https://github.com/Clinical-Pharmacy-Saarland-University/ospgsa/blob/main/dev/gsa-methods.md) gives a short treatment of the methods.

## Installation

```r
# install.packages("remotes")
remotes::install_github("Clinical-Pharmacy-Saarland-University/ospgsa")
```

The OSP evaluator additionally needs **`ospsuite`**, which is not on CRAN. Install it from the [OSP project](https://www.open-systems-pharmacology.org/). The GSA engine and all estimators also work without `ospsuite`, for example on analytic models via `function_evaluator()`.

## Quick start (analytic, no PK-Sim)

This runs with base R only and illustrates the two-stage method. A tiny model where only `k2` drives the output, but the inert `IC50` is strongly correlated with it. A correct, correlation-aware analysis must report `IC50` as influential **only through correlation**, not as a direct driver.

```r
library(ospgsa)

p <- gsa_parameters(
  gsa_parameter("k2",   "uniform", min = 0, max = 1),
  gsa_parameter("IC50", "uniform", min = 0, max = 1) # inert, but correlated with k2
)
R <- gsa_correlation(p, c("k2", "IC50", 0.9))

model <- function_evaluator(function(M) 2 * M[, "k2"], "C_T") # only k2 acts

# gsa() builds both designs (independent and correlated) plus the dummy floor,
# then runs the two-stage delta analysis in one call.
res <- gsa(p, model, method = "delta", n = 3000, correlation = R, boot = 50, seed = 1)

delta_classification(res)       # k2 = causal, IC50 = indirect-only
gsa_table(res, index = "delta") # ranked delta_1 / delta_2 with bootstrap CIs
plot_delta_two_stage(res)       # delta_1 (structural) vs delta_2 (full)
```

## A full two-stage analysis of a PBPK model

A complete pipeline for a real model. It uses an explicit independent and correlated design pair and evaluates the PK-Sim runs with resumable, crash-safe checkpointing, then feeds both designs to the two-stage estimator.

```r
library(ospgsa)
library(ospsuite)

sim <- system.file("extdata", "Aciclovir.pkml", package = "ospsuite")

# 1. Find the parameter paths you want to vary. ospsuite_parameter_paths() lists
#    or searches the model's paths (use fixed = TRUE to avoid escaping the "|").
ospsuite_parameter_paths(sim, "lipophilicity")
ospsuite_parameter_paths(sim, "Fraction unbound", fixed = TRUE, value = TRUE)

# 2. Describe the uncertain parameters as a spec table. ospsuite_parameters()
#    resolves each pattern to a unique path, anchors the marginal at the model's
#    current value and appends a dummy noise-floor parameter.
specs <- data.frame(
  name    = c("Lipophilicity", "FractionUnbound", "Hematocrit"),
  pattern = c(
    "Aciclovir\\|Lipophilicity$",
    "Aciclovir\\|Fraction unbound \\(plasma\\)$",
    "Organism\\|Hematocrit$"
  ),
  dist    = c("normal", "truncnorm", "truncnorm"),
  width   = c(0.3, 0.15, 0.1),
  lower   = c(NA, 1e-4, 0.3),
  upper   = c(NA, 0.999, 0.6)
)
params <- ospsuite_parameters(sim, specs)

# 3. Physiological correlation as a Spearman rank target (a list of triplets).
R <- gsa_correlation(params, list(c("Hematocrit", "FractionUnbound", 0.3)))

# 4. The two designs for the two-stage analysis (distinct seeds).
des <- gsa_two_stage_designs(params, n = 20000, correlation = R, seed = 1)

# 5. A multi-core PK-Sim evaluator that keeps one model loaded.
ev <- ospsuite_evaluator(sim, params, pk_parameters = c("AUC_tEnd", "C_max"), n_cores = parallel::detectCores())

# 6. Crash-safe, resumable evaluation. Re-running the same call resumes from the
#    last finished block. Distinct tags keep the two designs side by side.
ckpt   <- "gsa_out/checkpoints"
Y_corr <- gsa_evaluate(des$full$X, ev, checkpoint_dir = ckpt, block_size = 1000, tag = "corr")$Y
Y_ind  <- gsa_evaluate(des$ind$X, ev, checkpoint_dir = ckpt, block_size = 1000, tag = "indep")$Y

# 7. The two-stage delta estimator on the evaluated designs. AUC and Cmax are
#    skewed, so analyze them on the log scale and set non-positive outputs to NA.
res <- gsa_delta_two_stage(
  X_full = des$full$X, Y_full = gsa_sanitize_positive(Y_corr),
  X_ind  = des$ind$X, Y_ind = gsa_sanitize_positive(Y_ind),
  boot = 500, log_output = TRUE, dummy = "dummy", seed = 1
)

# 8. Read, tabulate and plot.
delta_classification(res)
gsa_table(res, index = "delta")
plot_delta_two_stage(res, output = "Plasma (Peripheral Venous Blood)__AUC_tEnd")
gsa_save_plots(res, dir = "figs")
```

Notes on the checkpointed evaluation.

- **Resume is automatic.** Re-run the identical call (same design, `checkpoint_dir`, `block_size`, `tag`) and finished blocks are read from disk, so only missing blocks are computed. Changing the design (any value, the seed, the marginals or the correlation) changes the content hash and starts a fresh checkpoint set, so stale results are never silently reused.
- **`crash_skip` for models that hard-crash.** A stiff parameter draw can segfault the ODE solver and kill the R _process_, which `tryCatch()` cannot catch. With `crash_skip = TRUE` each block is guarded by a marker file, and a block that crashed the process last time is filled with `NA` rows on the next resume (the estimators drop `NA` rows) instead of crashing forever. Run such jobs in a re-launching loop until they complete. It is off by default, because while on, interrupting a block mid-evaluation also turns that block into `NA`.

For a quick look without checkpointing, `gsa(params, ev, method = "delta", correlation = R, n = 20000)` builds the design pair and runs both stages internally. The explicit pipeline above makes the costly model evaluations crash-safe and resumable.

### Other building blocks

- **Screen first.** `gsa(..., method = "morris")` cheaply drops inactive parameters before a full delta run.
- **Convergence.** `gsa_convergence()` with `plot_convergence()` re-estimates at growing sample sizes to confirm the rankings are stable.
- **Time-resolved drivers.** Build a time-resolved evaluator with `ospsuite_evaluator(..., time_points = ...)` and use `plot_time_heatmap()` to see how importance shifts along the concentration-time curve.
- **Cross-check.** `gsa_sobol()` with `plot_sobol()` gives variance-based S1/ST indices for the independent-input case (the S1 to ST gap reveals interactions).

### Custom output metric

To analyze a metric the built-in PK parameters do not cover (`t_max`, a partial AUC, a slope, a ratio, time above a threshold, ...), pass `ospsuite_evaluator()` a `summary_fn`. It receives each run's simulated curve(s) and the raw results and returns a named numeric vector. Those names become the GSA outputs and flow through the rest of the analysis unchanged.

```r
# t_max and a partial AUC(0-24 h) on the plasma curve.
my_metric <- function(profiles, sr) {
  d <- profiles[["MTX"]] # data.frame(time, value), base units
  early <- d[d$time <= 24 * 60, ]
  c(
    tmax = d$time[which.max(d$value)],
    auc_0_24 = sum(diff(early$time) * (head(early$value, -1) + tail(early$value, -1)) / 2)
  )
}

# Debug it on a single run before launching the full GSA.
ospsuite_test_summary(sim, params, my_metric, outputs = MTX_PATH)

# Then use it like any evaluator.
ev <- ospsuite_evaluator(sim, params, outputs = MTX_PATH, summary_fn = my_metric)
```

`summary_fn` replaces the built-in `pk_parameters` extraction. For PK-Sim's own PK parameters inside it, call `ospsuite::calculatePKAnalyses(sr)`.

## Function reference

**Describe the inputs**

| Function                     | Purpose                                                                                                                                      |
| ---------------------------- | -------------------------------------------------------------------------------------------------------------------------------------------- |
| `gsa_parameter()`            | One uncertain parameter (`uniform`, `loguniform`, `normal`, `lognormal`, `truncnorm`) with optional display unit and log or natural scale    |
| `gsa_parameters()`           | Collect parameters into one object                                                                                                           |
| `gsa_correlation()`          | Build and validate a (rank) correlation matrix from pairwise triplets (or a list of them) or a full matrix                                   |
| `gsa_dummy()`                | Add a phantom parameter as the empirical significance noise floor                                                                            |
| `ospsuite_parameters()`      | One-call setup that builds model-anchored parameters from a spec table (resolve paths, anchor each marginal at the model value, add a dummy) |
| `ospsuite_resolve_paths()`   | Resolve regex patterns to unique PK-Sim/MoBi parameter paths                                                                                 |
| `ospsuite_parameter_paths()` | List or search a model's parameter paths (optionally with their current values) to find the path for a spec                                  |

**Sample and evaluate**

| Function                  | Purpose                                                                                                                                      |
| ------------------------- | -------------------------------------------------------------------------------------------------------------------------------------------- |
| `gsa_sample()`            | Draw a design (Latin hypercube, correlated via Iman-Conover or Gaussian copula)                                                              |
| `ospsuite_evaluator()`    | Run an OSP model many times via `SimulationBatch` (multi-core). Returns PK parameters, time-resolved concentrations, or a custom `summary_fn` metric |
| `ospsuite_test_summary()` | Run one simulation and show what a custom `summary_fn` receives and returns (debug aid)                                                      |
| `function_evaluator()`    | Wrap any R function as an evaluator                                                                                                          |
| `gsa_evaluate()`          | Run a design through an evaluator and collect outputs, optionally checkpointed and resumable for long or crash-prone runs (`checkpoint_dir`) |
| `gsa_two_stage_designs()` | Build the independent and correlated design pair (distinct seeds) for a two-stage delta run                                                  |
| `gsa_sanitize_positive()` | Set non-positive outputs to `NA` before log-scale analysis                                                                                   |

**Compute sensitivity**

| Function                 | Purpose                                                                                                                                                      |
| ------------------------ | ------------------------------------------------------------------------------------------------------------------------------------------------------------ |
| `gsa()`                  | One-call driver that builds the designs, evaluates and computes the requested indices (auto two-stage delta plus dummy floor when a correlation is supplied) |
| `gsa_delta()`            | Borgonovo delta moment-independent index (given-data estimator, bootstrap CIs)                                                                               |
| `gsa_delta_two_stage()`  | delta on an independent versus a correlated design, splitting structural from correlation-driven importance                                                  |
| `delta_classification()` | Label each parameter _causal_, _indirect-only_ or _both_                                                                                                     |
| `gsa_sobol()`            | Sobol first-order (`S1`) and total-order (`ST`) variance indices                                                                                             |
| `gsa_morris()`           | Morris elementary effects (`mu`, `mu*`, `sigma`) for screening                                                                                               |
| `gsa_regression()`       | Standardized or partial (rank) regression (`SRC`, `SRRC`, `PCC`, `PRCC`) with model R²                                                                       |
| `gsa_convergence()`      | Re-estimate at increasing sample sizes to check convergence                                                                                                  |
| `ospgsa_cluster()`       | Build a PSOCK cluster with ospgsa loaded on the workers to parallelize the delta bootstrap (`cl =` in `gsa_delta()`, `gsa_delta_two_stage()` and `gsa()`)    |

**Inspect and present**

| Function                                | Purpose                                                                                 |
| --------------------------------------- | --------------------------------------------------------------------------------------- |
| `print()` / `summary()` / `gsa_table()` | View ranked indices                                                                     |
| `plot()` / `gsa_plot()`                 | Plot method for a result that dispatches to the `plot_*` function for the chosen `type` |
| `plot_indices()`                        | Tornado bar plot of any index with confidence intervals                                 |
| `plot_sobol()`                          | Grouped `S1` versus `ST` bars (the gap indicates interactions)                          |
| `plot_delta_two_stage()`                | $\delta_1$ (structural) versus $\delta_2$ (full) bars                                         |
| `plot_morris()`                         | (mu\*, sigma) screening plane                                                           |
| `plot_time_heatmap()`                   | Index over time (parameters by time)                                                    |
| `plot_convergence()`                    | Index versus sample size                                                                |
| `plot_scatter()`                        | Output versus each input (a raw sanity check)                                           |
| `gsa_save_plots()`                      | Write the standard PNG set (delta two-stage, PRCC/SRC, Sobol, Morris) per output        |

Every estimator returns an `ospgsa_result`, a tidy `data.table` of indices (`estimate`, bootstrap CIs, `significant`, `rank`) plus the evaluated design, the outputs and run metadata.

## Implemented methods

| Method          | What it measures                                           | Correlation-safe  | Typical sample size         |
| --------------- | ---------------------------------------------------------- | ----------------- | --------------------------- |
| delta two-stage | whole-distribution shift, structural vs correlation-driven | **yes** (primary) | hundreds to 10⁴             |
| Sobol `S1`/`ST` | variance shares and interactions                           | no                | `n·(k+2)`, `n = 2¹⁰ to 2¹⁴` |
| Morris `mu*`    | screening importance and nonlinearity                      | no                | `r·(k+1)`                   |
| SRC / PRCC      | (rank) regression strength                                 | partly            | any LHS, `n` above 10³      |

## References

- De Carlo A, Tosca EM, Melillo N, Magni P (2023). _J Pharmacokinet Pharmacodyn_ 50(5):395-409. [doi:10.1007/s10928-023-09872-w](https://doi.org/10.1007/s10928-023-09872-w)
- Cuquerella-Gilabert M et al. (2026). _Pharmaceutical Statistics_ 25(2):e70082. [doi:10.1002/pst.70082](https://doi.org/10.1002/pst.70082)
- Borgonovo E (2007). _Reliab Eng Syst Saf_ 92(6):771-784. [doi:10.1016/j.ress.2006.04.015](https://doi.org/10.1016/j.ress.2006.04.015)
- Plischke E, Borgonovo E, Smith CL (2013). _Eur J Oper Res_ 226(3):536-550. [doi:10.1016/j.ejor.2012.11.047](https://doi.org/10.1016/j.ejor.2012.11.047)

## How to cite

If you use `ospgsa`, please cite the package together with the two-stage delta method papers it implements. Run `citation("ospgsa")` for the formatted entries (and BibTeX via `toBibtex(citation("ospgsa"))`).

```r
citation("ospgsa")
```

- Selzer D, Lehr T, Rüdesheim S, Dette C, Marok F (2026). _ospgsa: Global Sensitivity Analysis for Open Systems Pharmacology Models._ R package version 0.1.0. <https://github.com/Clinical-Pharmacy-Saarland-University/ospgsa>
- De Carlo A, Tosca EM, Melillo N, Magni P (2023). _J Pharmacokinet Pharmacodyn_ 50(5):395-409. [doi:10.1007/s10928-023-09872-w](https://doi.org/10.1007/s10928-023-09872-w)
- Cuquerella-Gilabert M et al. (2026). _Pharmaceutical Statistics_ 25(2):e70082. [doi:10.1002/pst.70082](https://doi.org/10.1002/pst.70082)

## License

GPL-3. Developed at the Chair of Clinical Pharmacy, Saarland University.
