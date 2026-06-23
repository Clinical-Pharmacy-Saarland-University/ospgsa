# Global Sensitivity Analysis for PBPK Models

A short methods primer for the `ospgsa` package.

## Why global instead of local

The built-in `ospsuite::SensitivityAnalysis` reports a normalized local finite-difference sensitivity of each PK metric to each parameter,

$$
S_{ij} = \frac{\Delta PK_j / PK_j}{\Delta p_i / p_i} \approx \frac{\partial PK_j}{\partial p_i}\cdot\frac{p_i}{PK_j},
$$

evaluated by perturbing one parameter at a time around a single reference point. This is fast but local. It describes the model only near one nominal vector, holds all other parameters fixed so it cannot see interactions, and uses a fixed fraction rather than the parameters' distributions so it cannot apportion output uncertainty to inputs.

Global methods vary all parameters at once across their full distributions and decompose the resulting output variability. They explore the whole input space, capture nonlinearity and interactions, and are defined relative to an explicit input distribution. This makes them the right tool when the question is which parameters drive the uncertainty of a PBPK prediction.

Physiology is correlated by construction. Organ volumes and blood flows scale with body weight, regional flows sum to cardiac output, enzyme abundances co-vary and clearance and volume are linked. Sobol, FAST and Morris all assume independent inputs. On a correlated sample they misattribute variance, so first-order indices no longer sum to at most 1 and a non-influential parameter can look important only because it correlates with an influential one. The package therefore uses a correlation-aware method that separates a parameter's own effect from the effect it inherits through correlation.

## The delta moment-independent index

For output $Y = f(X_1, \dots, X_k)$ let $f_Y(y)$ be the unconditional output density and $f_{Y\mid X_i}(y)$ the density of $Y$ conditional on a fixed value of $X_i$. The shift at $X_i = x_i$ is the area between the two densities,

$$
s(X_i) = \int \bigl|\, f_Y(y) - f_{Y\mid X_i}(y)\,\bigr|\, dy,
$$

and Borgonovo's delta is half its expectation over the distribution of $X_i$,

$$
\delta_i = \tfrac{1}{2}\, \mathbb{E}_{X_i}\!\bigl[\, s(X_i)\,\bigr].
$$

Key properties. The index lies in $0 \le \delta_i \le 1$. It is moment-independent, so it responds to any change in the output distribution such as shape, skew, or tails, not only variance. A factor can have Sobol $S_i \approx 0$ yet $\delta_i > 0$ if it reshapes the distribution without changing variance, which suits the skewed PK outputs AUC and Cmax. It uses only the joint law of $(X_i, Y)$, so it is computable from a correlated sample with no modification, unlike the Sobol terms which are orthogonal only under independence. The delta values do not sum to 1 and do not split into clean main and interaction terms, so they rank and attribute rather than partition variance.

## The two-stage approach for correlated inputs

On a correlated sample delta blends a parameter's structural effect through the model $f$ with its correlation-induced effect inherited from correlated partners. The two-stage method of De Carlo et al. (2023) computes delta on two designs and compares.

Stage 1, independence. Each parameter is sampled from its marginal with all correlations switched off, giving $\delta_{1,i}$. Removing correlation isolates the direct causal contribution, so a non-zero $\delta_{1,i}$ certifies a genuine structural effect.

Stage 2, full joint. Parameters are sampled from the realistic correlated joint distribution, giving $\delta_{2,i}$. These indices reflect importance under the true output distribution and capture correlation-transmitted effects.

| Category                | Condition                                 | Reading                                                        |
| ----------------------- | ----------------------------------------- | -------------------------------------------------------------- |
| Causal                  | $\delta_{1,i} > 0$                        | direct model-driven effect                                     |
| Correlation-driven only | $\delta_{1,i} = 0$ and $\delta_{2,i} > 0$ | influential only via correlation with an influential parameter |
| Both                    | $\delta_{1,i} > 0$ and $\delta_{2,i} > 0$ | direct effect plus correlation transmission                    |

The structural part is what survives Stage 1. The correlative part is what is gained or lost moving from Stage 1 to Stage 2. A parameter that matters only through correlation is not an independent driver, which tells you whether you can estimate or perturb it on its own.

## Companion methods

These are offered for independent-input cases and for cross-checks.

Morris elementary effects give a cheap first-pass screen. The mean absolute effect $\mu^*$ proxies overall importance and the spread $\sigma$ flags nonlinearity and interactions.

Sobol $S_1$ and $S_T$ give a variance decomposition under independent inputs. $S_1$ is the variance share from fixing $X_i$ and $S_T$ adds all interactions involving $X_i$.

Regression methods such as SRC, SRRC, PRCC, and PCC are cheap and interpretable when the model is close to linear or monotone, judged by the regression $R^2$.

## References

- De Carlo, A., Tosca, E. M., Melillo, N., Magni, P. (2023). A two-stages global sensitivity analysis by using the delta sensitivity index in presence of correlated inputs. _J Pharmacokinet Pharmacodyn_ 50(5):395-409. https://doi.org/10.1007/s10928-023-09872-w
- Cuquerella-Gilabert, M. et al. (2026). Leveraging Two-Stage delta Global Sensitivity Analysis Method to Inform Parameter Estimation in PBPK Models. _Pharmaceutical Statistics_ 25(2):e70082. https://doi.org/10.1002/pst.70082
- Borgonovo, E. (2007). A new uncertainty importance measure. _Reliability Engineering & System Safety_ 92(6):771-784. https://doi.org/10.1016/j.ress.2006.04.015
- Plischke, E., Borgonovo, E., Smith, C. L. (2013). Global sensitivity measures from given data. _European Journal of Operational Research_ 226(3):536-550. https://doi.org/10.1016/j.ejor.2012.11.047
- Saltelli, A. et al. (2010). Variance based sensitivity analysis of model output. _Computer Physics Communications_ 181(2):259-270.
- Morris, M. D. (1991). Factorial sampling plans for preliminary computational experiments. _Technometrics_ 33(2):161-174.
