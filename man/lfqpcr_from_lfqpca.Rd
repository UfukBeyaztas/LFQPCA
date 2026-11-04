\name{lfqpcr_from_lfqpca}
\alias{lfqpcr_from_lfqpca}
\title{Quantile regression on the persistent and transient components}
\description{
Second stage of LFQPCR: regresses a quantile of a scalar response on the
subject and visit components of a fitted LFQPCA model.
}
\usage{
lfqpcr_from_lfqpca(lfqpca_fit, V, covariates = NULL, tau_R = 0.5,
                   response_time_degree = NULL,
                   rq_method = c("fnb", "br", "fn"), space_tol = 1e-10,
                   rank_tol = 1e-7)
}
\arguments{
\item{lfqpca_fit}{an object of class \code{lfqpca_fit}, fitted at level
\eqn{\tau_F}{tau_F}.}
\item{V}{scalar responses, one for each curve of the fit. Visits with a missing
response are left out.}
\item{covariates}{optional numeric matrix or data frame of scalar covariates,
one row for each curve of the fit.}
\item{tau_R}{quantile level of the response.}
\item{response_time_degree}{degree of the time basis \eqn{G(u)}{G(u)}; at least
\code{mean_time_degree} of the fit, which is the default.}
\item{rq_method}{quantile regression algorithm, as in \code{lfqpca}.}
\item{space_tol}{relative tolerance on the singular values when the bases of
the loading spaces are built.}
\item{rank_tol}{tolerance for removing linearly dependent columns of the
design.}
}
\details{
The model for the \eqn{\tau_R}{tau_R}-quantile of the response is
\deqn{Q_{\tau_R}(V_{ij}) = G(u_{ij})^\top \alpha + C_{ij}^\top \beta + \langle \gamma_B, B_{ij} \rangle + \langle \gamma_U, U_{ij} \rangle,}{Q_tauR(V_ij) = G(u_ij)' alpha + C_ij' beta + <gamma_B, B_ij> + <gamma_U, U_ij>,}
where \eqn{B_{ij}}{B_ij} and \eqn{U_{ij}}{U_ij} are the persistent and
transient components of the LFQPCA fit and the time basis is
\eqn{G(u) = (1, u, \ldots, u^{p_G})^\top}{G(u) = (1, u, ..., u^pG)'}. Inner
products \eqn{\langle f, g \rangle = \int f g}{<f, g> = int f g} are computed
with the trapezoidal rule on the grid. Only the parts of \eqn{\gamma_B}{gamma_B} in the
span of the loadings \eqn{\phi_{0k}, \phi_{1k}}{phi_0k, phi_1k} and of
\eqn{\gamma_U}{gamma_U} in the span of the \eqn{\psi_l}{psi_l} enter the model,
so the coefficient functions are taken in these spaces.

Orthonormal bases \eqn{\chi_B}{chi_B} and \eqn{\chi_U}{chi_U} of the two spaces
come from singular value decompositions of the weighted loadings. With the
coordinates
\eqn{x_{B,ij} = \langle \chi_B, B_{ij} \rangle}{x_B,ij = <chi_B, B_ij>} and
\eqn{x_{U,ij} = \langle \chi_U, U_{ij} \rangle}{x_U,ij = <chi_U, U_ij>}, the
estimate minimizes
\deqn{\frac{1}{I_R} \sum_{i} \frac{1}{J_i^R} \sum_{j \in R_i} \rho_{\tau_R}\left(V_{ij} - G(u_{ij})^\top \alpha - C_{ij}^\top \beta - x_{B,ij}^\top a_B - x_{U,ij}^\top a_U\right),}{(1/I_R) sum_i (1/J_i^R) sum_(j in R_i) rho_tauR(V_ij - G(u_ij)' alpha - C_ij' beta - x_B,ij' a_B - x_U,ij' a_U),}
where \eqn{R_i}{R_i} holds the \eqn{J_i^R}{J_i^R} visits of subject \eqn{i} with
an observed response and \eqn{I_R}{I_R} is the number of such subjects. Then
\eqn{\gamma_B = a_B^\top \chi_B}{gamma_B = a_B' chi_B} and
\eqn{\gamma_U = a_U^\top \chi_U}{gamma_U = a_U' chi_U}. Linearly dependent
columns of the design are removed before the fit, and their coefficients are
set to zero.

In terms of the scores,
\eqn{\langle \gamma_B, B_{ij} \rangle = \sum_k \xi_{ik}(\theta_{0k} + u_{ij}\theta_{1k})}{<gamma_B, B_ij> = sum_k xi_ik (theta_0k + u_ij theta_1k)}
and
\eqn{\langle \gamma_U, U_{ij} \rangle = \sum_l \zeta_{ijl}\delta_l}{<gamma_U, U_ij> = sum_l zeta_ijl delta_l},
with \eqn{\theta_{0k} = \langle \gamma_B, \phi_{0k} \rangle}{theta_0k = <gamma_B, phi_0k>},
\eqn{\theta_{1k} = \langle \gamma_B, \phi_{1k} \rangle}{theta_1k = <gamma_B, phi_1k>}
and \eqn{\delta_l = \langle \gamma_U, \psi_l \rangle}{delta_l = <gamma_U, psi_l>}.
}
\value{
An object of class \code{lfqpcr_fit}, a list with components including
\item{tau_F, tau_R}{quantile levels of the curves and of the response.}
\item{coefficients}{named vector of the regression coefficients: time terms,
covariates and coordinates.}
\item{a_alpha, beta, a_B, a_U}{the blocks of \code{coefficients}.}
\item{gamma_B, gamma_U}{coefficient functions on the grid.}
\item{theta0, theta1, delta}{coefficients of the score form.}
\item{fitted, residuals}{fitted \eqn{\tau_R}{tau_R}-quantiles and residuals of
\code{V}.}
\item{objective}{value of the criterion.}
\item{B_hat, U_hat, x_B, x_U}{fitted components and their coordinates.}
\item{chi_B, chi_U, r_B, r_U}{bases of the loading spaces and their
dimensions.}
\item{design_rank, dropped_columns}{rank of the design and the removed
columns.}
\item{I_R, J_R, response_weights}{number of subjects with a response, their
numbers of visits and the weights of the criterion.}
\item{lfqpca_fit, lfqpca_converged, lfqpca_nondegenerate}{the first-stage fit
and its convergence and full-rank flags. The regression is only as reliable as
this fit.}
}
\examples{
\dontrun{
# Experiment II, Case 2, with tau_F = tau_R = 0.9
sim <- lfqpcr_simulate("Case2", I_train = 60, I_test = 20, J = 5,
                       tau_F = 0.9, tau_R = 0.9, replication = 1)
tr <- sim$train

# First stage at tau_F = 0.9, reached from the median
fit_F <- lfqpca_continuation(sim$Y[tr, ], sim$subject[tr], sim$visit_time[tr],
                             sim$grid, tau = 0.9, K_subject = 2, K_visit = 2,
                             spline_df = 10, mean_time_degree = 2,
                             max_iter = 60, tol = 1e-5,
                             rq_method_functions = "fnb",
                             rq_method_scores = "br")

# Second stage at tau_R = 0.9
fit_R <- lfqpcr_from_lfqpca(fit_F, V = sim$V[tr], covariates = sim$C[tr, ],
                            tau_R = 0.9, response_time_degree = 2,
                            rq_method = "fnb")
fit_R
coef(fit_R)
matplot(sim$grid, cbind(fit_R$gamma_B, sim$truth$gamma_B), type = "l",
        lty = 1:2, col = 1, xlab = "s", ylab = "gamma_B")
}
}
