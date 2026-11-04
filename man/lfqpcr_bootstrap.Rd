\name{lfqpcr_bootstrap}
\alias{lfqpcr_bootstrap}
\title{Subject bootstrap for LFQPCR}
\description{
Refits LFQPCR to samples of subjects drawn with replacement.
}
\usage{
lfqpcr_bootstrap(Y, V, subject, visit_time, grid = NULL, covariates = NULL,
                 B = 200, seed = 1, ..., verbose = TRUE)
}
\arguments{
\item{Y, V, subject, visit_time, grid, covariates}{the data, as in
\code{lfqpcr}.}
\item{B}{number of bootstrap samples.}
\item{seed}{seed for the resampling. The random number state of the session is
restored afterwards.}
\item{...}{further arguments of \code{lfqpcr}, such as \code{tau_F},
\code{tau_R} and the dimensions.}
\item{verbose}{if \code{TRUE}, progress is reported by messages.}
}
\details{
Each bootstrap sample draws as many subjects as the data contain, with
replacement, and keeps all their visits; a subject drawn twice enters as two
subjects. Fits that fail are kept as error objects. Percentile intervals for
the coefficient functions or other quantities follow from the successful fits.
A fit counts as successful when \code{lfqpcr} returns; the components
\code{lfqpca_converged} and \code{lfqpca_nondegenerate} of each fit show whether
its first stage converged with full-rank blocks.
}
\value{
A list with components \code{fits} (list of \code{lfqpcr_fit} objects and error
objects), \code{successful} (logical vector), \code{n_successful},
\code{n_failed} and \code{seed}.
}
\examples{
\dontrun{
sim <- lfqpcr_simulate("Case2", I_train = 60, I_test = 20, J = 5,
                       tau_F = 0.5, tau_R = 0.9, replication = 1)
tr <- sim$train
boot <- lfqpcr_bootstrap(sim$Y[tr, ], sim$V[tr], sim$subject[tr],
                         sim$visit_time[tr], sim$grid,
                         covariates = sim$C[tr, ], B = 50, seed = 1,
                         tau_F = 0.5, tau_R = 0.9, K_subject = 2,
                         K_visit = 2, spline_df = 10, mean_time_degree = 2,
                         response_time_degree = 2, lfqpca_max_iter = 60,
                         lfqpca_tol = 1e-5, lfqpca_n_starts = 5)
boot$n_successful

# Pointwise 95\% percentile band for gamma_B
gamma_B <- sapply(boot$fits[boot$successful], function(f) f$gamma_B)
band <- apply(gamma_B, 1, quantile, probs = c(0.025, 0.975))
matplot(sim$grid, cbind(t(band), sim$truth$gamma_B), type = "l",
        lty = c(2, 2, 1), col = 1, xlab = "s", ylab = "gamma_B")
}
}
