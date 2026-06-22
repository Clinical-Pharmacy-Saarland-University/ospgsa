# End-to-end global sensitivity analysis of a PK-Sim model with ospgsa.
# Uses the Aciclovir example shipped with ospsuite. Run interactively.

library(ospgsa)
library(ospsuite)

sim_path <- system.file("extdata", "Aciclovir.pkml", package = "ospsuite")

## 1. Describe parameter uncertainty ------------------------------------------
# Values are given in display units; ospgsa converts them to base units before
# each run. Order-of-magnitude parameters should use loguniform / lognormal.
p <- gsa_parameters(
  gsa_parameter(
    "Lipophilicity",
    "normal",
    mean = -0.097,
    sd = 0.3,
    path = "Aciclovir|Lipophilicity",
    unit = "Log Units"
  ),
  gsa_parameter(
    "FractionUnbound",
    "truncnorm",
    mean = 0.85,
    sd = 0.1,
    lower = 0.3,
    upper = 0.999,
    path = "Aciclovir|Fraction unbound (plasma)"
  ),
  gsa_parameter(
    "Hematocrit",
    "truncnorm",
    mean = 0.45,
    sd = 0.05,
    lower = 0.3,
    upper = 0.6,
    path = "Organism|Hematocrit"
  )
)

## 2. Optional physiological correlation (Spearman rank target) ---------------
R <- gsa_correlation(p, c("Hematocrit", "FractionUnbound", 0.3))

## 3. Build the PK-Sim evaluator ----------------------------------------------
ev <- ospsuite_evaluator(
  sim_path,
  p,
  pk_parameters = c("AUC_tEnd", "C_max", "t_max"),
  n_cores = max(1L, parallel::detectCores() - 1L)
)

## 4. Two-stage delta GSA on log-transformed PK metrics -----------------------
res <- gsa(
  p,
  ev,
  method = "delta",
  n = 2000,
  correlation = R,
  boot = 500,
  conf = 0.95,
  log_output = TRUE,
  seed = 1
)

print(gsa_table(res, index = "delta"))
print(delta_classification(res)) # causal / indirect-only / both

plot_delta_two_stage(res, output = "Plasma (Peripheral Venous Blood)__AUC_tEnd")

## 5. Cross-check the drivers with Sobol (independent design) ------------------
res_sobol <- gsa_sobol(p, ev, n = 2^11, boot = 200, log_output = TRUE)
print(gsa_table(res_sobol))
plot_sobol(res_sobol, output = "Plasma (Peripheral Venous Blood)__C_max")

## 6. Convergence + standard figure set ---------------------------------------
conv <- gsa_convergence(p, ev, method = "delta", n_seq = c(500, 1000, 2000), boot = 100)
plot_convergence(conv)
gsa_save_plots(res, dir = "aciclovir_figs")

## 7. Time-resolved sensitivity along the concentration-time curve ------------
evT <- ospsuite_evaluator(
  sim_path,
  p,
  pk_parameters = character(0),
  time_points = seq(0, 720, by = 30),
  n_cores = max(1L, parallel::detectCores() - 1L)
)
resT <- gsa(p, evT, method = "delta", n = 2000, correlation = R, boot = 100, seed = 1)
plot_time_heatmap(resT, index = "delta")
