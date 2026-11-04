# LFQPCA

**Longitudinal Functional Quantile Principal Component Analysis**

`LFQPCA` analyzes functional observations collected at repeated visits from the same subjects. It separates persistent subject trajectories from temporary visit deviations at a chosen quantile level. The package also fits longitudinal functional quantile principal component regression (LFQPCR), which relates a scalar response to these two components through separate coefficient functions.

The [package manual](LFQPCA_0.1.0.pdf) contains the model definitions, function arguments, and further examples.

## What the model estimates

At each quantile level, a curve is represented by three parts:

| Part | Interpretation |
| --- | --- |
| Population surface | The population quantile curve as a function of visit time. |
| Persistent subject component | A subject's functional intercept and slope in visit time. |
| Transient visit component | The deviation associated with an individual visit, constrained to carry no polynomial trend within a subject. |

The functions are represented by B-splines, with a periodic option for circular domains. Estimation uses check loss and gives every subject the same total weight, including when subjects have different numbers of visits. Missing curve values are allowed on a common functional grid.

For regression, the persistent and transient components enter separately. This allows their associations with a scalar response to differ.

## Installation

Install the package directly from GitHub:

```r
if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}

remotes::install_github("UfukBeyaztas/LFQPCA")

library(LFQPCA)
```

## Preparing the data

Each row corresponds to one subject at one visit. Keep the rows aligned across all inputs.

| Input | Format |
| --- | --- |
| `Y` | Numeric matrix with one curve per row and one functional grid point per column. Use `NA` for missing measurements. |
| `subject` | Subject labels, with one entry per row of `Y`. |
| `visit_time` | Numeric visit times, with one entry per row of `Y`. |
| `grid` | Strictly increasing numeric vector of length `ncol(Y)`, giving the locations along each curve. |
| `V` | For LFQPCR, a numeric scalar response for each row of `Y`. |
| `covariates` | Optional numeric matrix of scalar covariates, with one row per visit. Encode categorical variables before fitting. |

The functional grid and visit time describe different domains. For example, the grid might record time within a daily activity curve, while `visit_time` records years since enrollment. Supply visit times in their original units; the package standardizes them internally.

Every curve needs at least five observed grid points. The population time basis must have full rank within every subject, and estimating subject slopes requires at least two distinct visit times. With visit components and the default `strict_identification = TRUE`, a population time polynomial of degree `p` requires at least `p + 2` suitably distinct visits per subject. Thus, a linear population time effect requires at least three visits, and a quadratic effect requires at least four.

## Fitting LFQPCA

The package includes the simulation designs used in the accompanying study. This example generates 30 subjects with five visits each and fits the median model. The generated curves contain some missing values.

```r
library(LFQPCA)

dat <- lfqpca_simulate("D1", I = 30, J = 5)

fit_F <- lfqpca(
  Y = dat$Y,
  subject = dat$subject,
  visit_time = dat$visit_time,
  grid = dat$grid,
  tau = 0.5,
  K_subject = 2,
  K_visit = 2,
  spline_df = 10,
  mean_time_degree = 2,
  max_iter = 60,
  tol = 1e-5,
  n_starts = 2,
  verbose = FALSE
)

summary(fit_F)

plot(fit_F, type = "components")
plot(fit_F, type = "objective")
plot(fit_F, type = "fitted", curve_indices = 1:4)
```

`K_subject` and `K_visit` set the numbers of persistent and transient components. `spline_df` controls the dimension of the functional spline basis, and `mean_time_degree` controls the population trend in visit time. The choices above match the structure of design D1; they should be reconsidered for other data.

The fitted object contains the population functions, loading functions, and scores:

```r
fit_F$eta_functions       # Population functions on the grid
fit_F$phi0               # Subject intercept loadings
fit_F$phi1               # Subject slope loadings
fit_F$psi                # Visit loadings
head(fit_F$xi)           # Subject scores
head(fit_F$zeta)         # Visit scores

Q_hat <- fitted(fit_F)   # Fitted quantile curves
R_hat <- residuals(fit_F)
```

Check `fit_F$converged`, `fit_F$termination`, and `fit_F$nondegenerate` before interpreting a fit. The criterion is not jointly convex in the functions and scores, so different starts can lead to different solutions. For a substantive analysis, consider more starts and inspect the objective history and rank diagnostics.

### Other quantile levels

`lfqpca_continuation()` fits a target quantile by moving from a reference level through intermediate levels. An existing median fit can be reused:

```r
fit_90 <- lfqpca_continuation(
  Y = dat$Y,
  subject = dat$subject,
  visit_time = dat$visit_time,
  grid = dat$grid,
  tau = 0.9,
  reference_fit = fit_F,
  K_subject = 2,
  K_visit = 2,
  spline_df = 10,
  mean_time_degree = 2,
  max_iter = 60,
  tol = 1e-5,
  n_starts = 2,
  verbose = FALSE
)

fit_90$continuation_tau_path
summary(fit_90)
```

The reference fit must use the same data, grid, and model dimensions. With the default path, the fit moves from 0.5 through 0.7 and 0.8 to 0.9.

## Scalar response quantile regression

LFQPCR uses the fitted persistent and transient components as functional predictors. `tau_F` is the quantile level used to decompose the curves, and `tau_R` is the response quantile of interest. These levels can be chosen separately.

The following example uses median curve components to model the 0.9 quantile of a scalar response. Training and test sets contain different subjects.

```r
sim <- lfqpcr_simulate(
  "Case2", I_train = 30, I_test = 10, J = 5,
  tau_F = 0.5, tau_R = 0.9
)

tr <- sim$train
te <- sim$test

fit_R <- lfqpcr(
  Y = sim$Y[tr, , drop = FALSE],
  V = sim$V[tr],
  subject = sim$subject[tr],
  visit_time = sim$visit_time[tr],
  grid = sim$grid,
  covariates = sim$C[tr, , drop = FALSE],
  tau_F = 0.5,
  tau_R = 0.9,
  K_subject = 2,
  K_visit = 2,
  spline_df = 10,
  mean_time_degree = 2,
  response_time_degree = 2,
  lfqpca_max_iter = 60,
  lfqpca_tol = 1e-5,
  lfqpca_n_starts = 2,
  verbose = FALSE
)

summary(fit_R)
coef(fit_R)

gam <- lfqpcr_coefficient_functions(fit_R)

matplot(
  gam$grid, gam[, c("gamma_B", "gamma_U")],
  type = "l", lty = 1:2, col = c("navy", "firebrick"),
  xlab = "Functional Domain", ylab = "Coefficient"
)
legend(
  "topright", c("Persistent", "Transient"),
  lty = 1:2, col = c("navy", "firebrick"), bty = "n"
)

pred <- predict(
  fit_R,
  newY = sim$Y[te, , drop = FALSE],
  new_subject = sim$subject[te],
  new_visit_time = sim$visit_time[te],
  new_covariates = sim$C[te, , drop = FALSE]
)

head(pred$fitted)
```

`gamma_B` describes the association of the persistent component with the response quantile, and `gamma_U` describes the association of the transient component. Both coefficient functions lie in their respective fitted loading spaces. The first-stage fit is stored in `fit_R$lfqpca_fit`; inspect its convergence and rank diagnostics as well.

Prediction estimates the new subjects' scores from their observed curves while keeping the fitted functions fixed. Supply new curves on the training grid, visit times in the original units, and covariates with the same columns and encoding as in training. When visit components are included, new subjects also need enough distinct visits to satisfy the time-basis constraints. Each predicted visit uses all supplied curves from that subject.

If an LFQPCA model has already been fitted, use `lfqpcr_from_lfqpca()` to fit the response model without repeating the decomposition. This is also the route for using a first-stage fit obtained by continuation.

## Choosing the model and assessing uncertainty

The two cross-validation functions address different questions:

| Function | What is held out | Use |
| --- | --- | --- |
| `lfqpca_cv()` | Observed grid points within curves. All subjects and visits remain in training. | Compare component counts and spline dimensions for curve reconstruction at a fixed quantile level. |
| `lfqpcr_cv()` | Groups of subjects, with their visits kept together. | Evaluate response quantile prediction for a specified model. Compare calls to choose `tau_F`, component counts, or spline dimension. |

`lfqpcr_bootstrap()` resamples whole subjects and refits both stages, keeping each subject's visits together. It returns the fitted bootstrap models from which percentile intervals can be calculated. Inspect the convergence and rank flags of the returned fits before using them for inference. Cross-validation and bootstrap runs require repeated fits and can take substantially longer than the examples above.

See `?lfqpca_cv`, `?lfqpcr_cv`, and `?lfqpcr_bootstrap` for complete examples.

## Simulation designs

`lfqpca_simulate()` generates designs D1 to D7 for the decomposition study. `lfqpcr_simulate()` generates Case1 to Case6 for the regression study. These include quantile-specific scores or loading spaces, unbalanced visits, and departures from the fitted model assumptions.

```r
lfqpca_sim_designs()
```

The simulation functions accept `replication` and `seed` arguments. Their defaults reproduce the first replication; changing `replication` generates another data set. The manual describes the generating models and the settings used in each design.
