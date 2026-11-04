\name{lfqpca_continuation}
\alias{lfqpca_continuation}
\alias{lfqpca_continuation_path}
\title{LFQPCA along a path of quantile levels}
\description{
Fits LFQPCA at a target quantile level by moving from a reference level in
steps, each fit starting from the previous one.
}
\usage{
lfqpca_continuation(Y, subject, visit_time, grid = NULL, tau,
                    reference_tau = 0.5, fractions = c(0, 0.5, 0.75, 1),
                    n_starts = 5, path_starts = 2, reference_fit = NULL,
                    verbose = FALSE, ...)

lfqpca_continuation_path(tau, reference_tau = 0.5,
                         fractions = c(0, 0.5, 0.75, 1))
}
\arguments{
\item{Y, subject, visit_time, grid}{the data, as in \code{lfqpca}.}
\item{tau}{target quantile level.}
\item{reference_tau}{level at which the path starts.}
\item{fractions}{fractions of the distance from \code{reference_tau} to
\code{tau} at which the model is fitted; 0 and 1 are always included.}
\item{n_starts}{number of starts of the reference fit and of the fit at
\code{tau}.}
\item{path_starts}{number of starts of the intermediate fits; at least 2, since
one of them is the warm start.}
\item{reference_fit}{optional fit at \code{reference_tau}. It replaces the
reference fit, which saves time when several levels share one reference.}
\item{verbose}{passed to \code{lfqpca}.}
\item{...}{further arguments of \code{lfqpca}, for example \code{K_subject},
\code{K_visit}, \code{spline_df}, \code{mean_time_degree}, \code{max_iter} or
\code{tol}.}
}
\details{
The criterion is not convex in the functions and scores jointly, so the block
algorithm finds a local minimum that depends on its start. The continuation fits
the model at \eqn{\tau_0 + f(\tau - \tau_0)}{tau_0 + f (tau - tau_0)} for each
fraction \eqn{f}{f}, where \eqn{\tau_0}{tau_0} is \code{reference_tau}, and
passes each fit to the next one as \code{warm_start_fit}. The structured and the
random starts of \code{lfqpca} stay available at every step. With the default
fractions, a fit at 0.9 goes through 0.5, 0.7 and 0.8.

A failed fit at the reference level or at an intermediate level is recorded and
skipped, and the next fit starts from the last successful one. Only a failure
at \code{tau} is an error.
}
\value{
The fit at \code{tau} is returned, with two further components: the messages
of failed fits along the path, \code{continuation_errors}, and the levels of the
path, \code{continuation_tau_path}. The path function returns these levels
only.
}
\examples{
\dontrun{
# Levels visited on the way from 0.5 to 0.9
lfqpca_continuation_path(tau = 0.9)

# Experiment I, design D3, at tau = 0.9 as in the simulations
dat <- lfqpca_simulate("D3", I = 60, J = 5, replication = 1)
fit90 <- lfqpca_continuation(dat$Y, dat$subject, dat$visit_time, dat$grid,
                             tau = 0.9, K_subject = 2, K_visit = 2,
                             spline_df = 10, mean_time_degree = 2,
                             max_iter = 60, tol = 1e-5,
                             rq_method_functions = "fnb",
                             rq_method_scores = "br")
fit90$continuation_tau_path
fit90$continuation_errors
truth90 <- lfqpca_sim_truth(dat, tau = 0.9)
sqrt(mean((fitted(fit90) - truth90$quantile_surface)^2))

# Lower tail, starting from an existing median fit
fit50 <- lfqpca(dat$Y, dat$subject, dat$visit_time, dat$grid, tau = 0.5,
                K_subject = 2, K_visit = 2, spline_df = 10,
                mean_time_degree = 2, max_iter = 60, tol = 1e-5, n_starts = 5,
                rq_method_functions = "fnb", rq_method_scores = "br",
                verbose = FALSE)
fit10 <- lfqpca_continuation(dat$Y, dat$subject, dat$visit_time, dat$grid,
                             tau = 0.1, reference_fit = fit50,
                             K_subject = 2, K_visit = 2, spline_df = 10,
                             mean_time_degree = 2, max_iter = 60, tol = 1e-5,
                             rq_method_functions = "fnb",
                             rq_method_scores = "br")
fit10$continuation_tau_path
}
}
