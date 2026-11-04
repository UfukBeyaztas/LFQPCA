\name{lfqpcr_coefficient_functions}
\alias{lfqpcr_coefficient_functions}
\alias{lfqpcr_alpha}
\title{Coefficient functions and time effect of an LFQPCR fit}
\description{
\code{lfqpcr_coefficient_functions} returns the coefficient functions
\eqn{\gamma_B}{gamma_B} and \eqn{\gamma_U}{gamma_U} on the grid, and
\code{lfqpcr_alpha} evaluates the time effect
\eqn{G(u)^\top \alpha}{G(u)' alpha} at given visit times.
}
\usage{
lfqpcr_coefficient_functions(object)

lfqpcr_alpha(object, visit_time)
}
\arguments{
\item{object}{an object of class \code{lfqpcr_fit}.}
\item{visit_time}{visit times on the scale of the original \code{visit_time}.}
}
\details{
If the persistent component of a visit changes by a function \eqn{h}{h}, the
fitted \eqn{\tau_R}{tau_R}-quantile of the response changes by
\eqn{\int \gamma_B h}{int gamma_B h}; \eqn{\gamma_U}{gamma_U} plays the same
role for the transient component. Both functions lie in the spans of the fitted
loadings.
}
\value{
\code{lfqpcr_coefficient_functions} returns a data frame with columns
\code{grid}, \code{gamma_B} and \code{gamma_U}. \code{lfqpcr_alpha} returns a
numeric vector with one value per visit time.
}
\examples{
\dontrun{
sim <- lfqpcr_simulate("Case2", I_train = 60, I_test = 20, J = 5,
                       tau_F = 0.5, tau_R = 0.9, replication = 1)
tr <- sim$train
te <- sim$test
fit <- lfqpcr(sim$Y[tr, ], sim$V[tr], sim$subject[tr], sim$visit_time[tr],
              sim$grid, covariates = sim$C[tr, ], tau_F = 0.5, tau_R = 0.9,
              K_subject = 2, K_visit = 2, spline_df = 10,
              mean_time_degree = 2, response_time_degree = 2,
              lfqpca_max_iter = 60, lfqpca_tol = 1e-5, lfqpca_n_starts = 5,
              lfqpca_rq_method_functions = "fnb",
              lfqpca_rq_method_scores = "br", verbose = FALSE)
gam <- lfqpcr_coefficient_functions(fit)
matplot(gam$grid, cbind(gam$gamma_B, sim$truth$gamma_B), type = "l",
        lty = 1:2, col = 1, xlab = "s", ylab = "gamma_B")
lfqpcr_alpha(fit, visit_time = c(-1, 0, 1))
}
}
