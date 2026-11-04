\name{lfqpca}
\alias{lfqpca}
\alias{print.lfqpca_fit}
\alias{summary.lfqpca_fit}
\alias{print.summary.lfqpca_fit}
\alias{fitted.lfqpca_fit}
\alias{residuals.lfqpca_fit}
\title{Longitudinal functional quantile principal component analysis}
\description{
Fits the LFQPCA model at one quantile level to curves recorded at several visits
per subject.
}
\usage{
lfqpca(Y, subject, visit_time, grid = NULL, tau = 0.5, K_subject = 2,
       K_visit = 1, spline_df = 8, spline_degree = 3, periodic = FALSE,
       mean_time_degree = 1, max_iter = 30, tol = 1e-4, n_starts = 2,
       rq_method_functions = c("fnb", "fn", "br"),
       rq_method_scores = c("br", "fnb", "fn"), strict_identification = TRUE,
       objective_increase_tolerance = 1e-7, verbose = TRUE,
       warm_start_fit = NULL, rank_tol = 1e-8)

\method{print}{lfqpca_fit}(x, ...)

\method{summary}{lfqpca_fit}(object, ...)

\method{print}{summary.lfqpca_fit}(x, ...)

\method{fitted}{lfqpca_fit}(object, ...)

\method{residuals}{lfqpca_fit}(object, ...)
}
\arguments{
\item{Y}{numeric matrix with one curve per row and one column per grid point.
Missing values are treated as unobserved; each curve needs at least five
observed values.}
\item{subject}{subject labels, one for each row of \code{Y}.}
\item{visit_time}{visit times, one for each row of \code{Y}, for example years
since the first visit.}
\item{grid}{increasing grid of the functional domain. The default is an equally
spaced grid on [0, 1].}
\item{tau}{quantile level in (0, 1).}
\item{K_subject}{number of subject components.}
\item{K_visit}{number of visit components; 0 gives a model without visit
deviations.}
\item{spline_df}{dimension of the B-spline basis used for all functions.}
\item{spline_degree}{degree of the B-splines; 3 gives cubic splines.}
\item{periodic}{if \code{TRUE}, periodic B-splines are used, for curves on a
circular domain.}
\item{mean_time_degree}{degree \eqn{p}{p} of the polynomial time basis.}
\item{max_iter}{maximum number of iterations for each start.}
\item{tol}{tolerance for the relative change of the objective.}
\item{n_starts}{number of starting values; see Details.}
\item{rq_method_functions, rq_method_scores}{quantile regression algorithms for
the function and the score updates: \code{"fnb"} and \code{"fn"} both use the
Frisch-Newton interior point method and \code{"br"} the simplex method of
Barrodale and Roberts, both from \pkg{quantreg}. If the interior point method
fails, the simplex method is tried.}
\item{strict_identification}{if \code{TRUE} and \code{K_visit > 0}, every
subject needs at least \eqn{p + 2}{p + 2} distinct visit times. With
\code{FALSE}, subjects with fewer visits are kept and get zero visit scores.}
\item{objective_increase_tolerance}{largest increase of the objective accepted
in a block update.}
\item{verbose}{if \code{TRUE}, progress is reported by messages.}
\item{warm_start_fit}{optional \code{lfqpca_fit} for the same data and
dimensions, usually at a nearby quantile level; it is used as the second
start.}
\item{rank_tol}{tolerance for the rank checks of the fitted subject and visit
blocks.}
\item{x, object}{an object of class \code{lfqpca_fit}, or of class
\code{summary.lfqpca_fit} for the print method of the summary.}
\item{...}{not used.}
}
\details{
With \eqn{u_{ij}}{u_ij} the visit time standardized by its mean and standard
deviation and \eqn{H(u) = (1, u, \ldots, u^p)^\top}{H(u) = (1, u, ..., u^p)'},
the fitted quantile curves are
\deqn{Q_{ij}(s) = H(u_{ij})^\top \eta(s) + \xi_i^\top \{\phi_0(s) + u_{ij}\phi_1(s)\} + \zeta_{ij}^\top \psi(s).}{Q_ij(s) = H(u_ij)' eta(s) + xi_i' (phi_0(s) + u_ij phi_1(s)) + zeta_ij' psi(s).}
The visit scores of subject \eqn{i} are \eqn{Z_i = N_i A_i}{Z_i = N_i A_i},
where the columns of \eqn{N_i}{N_i} are an orthonormal basis of the null space
of \eqn{H_i^\top}{H_i'}, so that \eqn{H_i^\top Z_i = 0}{H_i' Z_i = 0}. The
subject scores have mean zero, the pairs
\eqn{(\phi_{0k}, \phi_{1k})}{(phi_0k, phi_1k)} are orthonormal for
\eqn{\int f_0 g_0 + \int f_1 g_1}{int f_0 g_0 + int f_1 g_1}, the
\eqn{\psi_l}{psi_l} are orthonormal, and in each block the scores are
uncorrelated with decreasing variances. Each component is signed so that its
largest loading value is positive.

All functions lie in the span of the B-spline basis, and the estimate minimizes
\deqn{L_\tau = \sum_{i=1}^{I} \sum_{j=1}^{J_i} \sum_{d \in O_{ij}} \frac{\omega_{ijd}}{I J_i}\, \rho_\tau\{Y_{ij}(s_d) - Q_{ij}(s_d)\},}{L_tau = sum_i sum_j sum_(d in O_ij) omega_ijd/(I J_i) rho_tau(Y_ij(s_d) - Q_ij(s_d)),}
where \eqn{\rho_\tau(r) = r\{\tau - 1(r < 0)\}}{rho_tau(r) = r(tau - 1(r < 0))}
is the check function, \eqn{O_{ij}}{O_ij} holds the observed grid points of
curve \eqn{(i, j)}{(i, j)} and the trapezoidal weights
\eqn{\omega_{ijd}}{omega_ijd} sum to one over these points. Every subject has
the same total weight, whatever its number of visits or missing values.

The algorithm cycles through three blocks. Given the scores, the spline
coefficients of all functions solve one weighted linear quantile regression.
Given the functions, the scores \eqn{(\xi_i, A_i)}{(xi_i, A_i)} of each subject
solve a small weighted quantile regression. A rotation then restores the
identifying conditions without changing the fitted curves. An update is accepted
only if it raises the objective by at most \code{objective_increase_tolerance};
the fit has converged when the relative change of the objective is below
\code{tol} in an iteration without rejected updates.

The first start fits the population part alone and takes singular value
decompositions of subject-wise intercept and slope curves of the residuals and
of the remaining visit residuals. With a warm start, the second start is the
fit given in \code{warm_start_fit}; the other starts are random. The fit with the smallest objective
among the starts with full-rank subject and visit blocks is returned, or among
all successful starts if none has full rank.
}
\value{
\code{lfqpca} returns an object of class \code{lfqpca_fit}, a list with
components including
\item{eta_functions}{matrix of the population functions on the grid, one column
per term of \eqn{H(u)}{H(u)}.}
\item{phi0, phi1}{matrices of the subject intercept and slope loadings, one
column per component.}
\item{psi}{matrix of the visit loadings.}
\item{xi}{matrix of subject scores, one row per subject in order of first
appearance in \code{subject}.}
\item{zeta}{matrix of visit scores, one row per curve.}
\item{lambda, nu}{variances of the subject and visit scores; \code{pve_subject}
and \code{pve_visit} give their cumulative proportions.}
\item{fitted}{matrix of fitted quantile curves.}
\item{objective, objective_history}{value of the criterion at the solution and
over the iterations.}
\item{converged, termination, iterations}{convergence flag, reason for stopping
and number of iterations. The reason is one of \code{"converged"},
\code{"max_iter"} and \code{"stalled_rejected_updates"}.}
\item{nondegenerate, rank_diagnostics}{whether the subject and visit blocks
have full rank, and the rank checks behind it.}
\item{constraint_diagnostics}{largest absolute entries of
\eqn{H_i^\top Z_i}{H_i' Z_i} for each subject.}
\item{time_center, time_scale}{mean and standard deviation used to standardize
\code{visit_time}.}
\item{all_start_objectives, start_types, best_start}{objectives and types of the
starts, and the index of the returned start.}
The data, grid, weights, spline and time bases and the settings of the call are
kept as well.

\code{summary} returns the quantile level, convergence information, the weighted
proportion of observed values at or below the fitted quantiles, the score
variances and their cumulative proportions, and the largest constraint errors.
\code{fitted} and \code{residuals} return the fitted curves and \code{Y} minus
the fitted curves.
}
\examples{
\dontrun{
# Experiment I, design D1, first replication
dat <- lfqpca_simulate("D1", I = 60, J = 5, replication = 1)
fit <- lfqpca(dat$Y, dat$subject, dat$visit_time, dat$grid, tau = 0.5,
              K_subject = 2, K_visit = 2, spline_df = 10,
              mean_time_degree = 2, max_iter = 60, tol = 1e-5, n_starts = 5,
              rq_method_functions = "fnb", rq_method_scores = "br",
              verbose = FALSE)
fit
summary(fit)

# Error of the fitted median surface
truth <- lfqpca_sim_truth(dat, tau = 0.5)
sqrt(mean((fitted(fit) - truth$quantile_surface)^2))
}
}
