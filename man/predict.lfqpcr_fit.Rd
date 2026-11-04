\name{predict.lfqpcr_fit}
\alias{predict.lfqpcr_fit}
\title{Response quantiles of new subjects}
\description{
Predicts the \eqn{\tau_R}{tau_R}-quantile of the response of new subjects from
their curves, times and covariates.
}
\usage{
\method{predict}{lfqpcr_fit}(object, newY, new_subject, new_visit_time,
  new_covariates = NULL, trajectory_only = FALSE, ...)
}
\arguments{
\item{object}{an object of class \code{lfqpcr_fit}.}
\item{newY, new_subject, new_visit_time}{curves, subject labels and visit times
of the new subjects, as in the predict method for \code{lfqpca_fit} objects.}
\item{new_covariates}{covariates of the rows of \code{newY}, with the same
columns as in the fit.}
\item{trajectory_only}{if \code{TRUE}, the coordinates of the transient
component are set to zero.}
\item{...}{not used.}
}
\details{
The LFQPCA scores of the new subjects are estimated from their curves with the
functions of the first stage held fixed. The persistent and transient components
are rebuilt from these scores, and their coordinates enter the fitted
regression. With \code{trajectory_only = TRUE} only the persistent component is
used, so the prediction depends on the curves through the subject trajectory
alone.
}
\value{
A list with components \code{fitted} (predicted quantiles), \code{x_B},
\code{x_U}, \code{B_hat}, \code{U_hat} (coordinates and components),
\code{xi}, \code{zeta}, \code{time_std}, \code{subject_id} and
\code{trajectory_only}.
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

# Prediction for the 20 test subjects and empirical coverage
pred <- predict(fit, newY = sim$Y[te, ], new_subject = sim$subject[te],
                new_visit_time = sim$visit_time[te],
                new_covariates = sim$C[te, ])
mean(sim$V[te] <= pred$fitted)

# Persistent component only
pred_B <- predict(fit, newY = sim$Y[te, ], new_subject = sim$subject[te],
                  new_visit_time = sim$visit_time[te],
                  new_covariates = sim$C[te, ], trajectory_only = TRUE)
mean(sim$V[te] <= pred_B$fitted)
}
}
