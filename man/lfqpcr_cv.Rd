\name{lfqpcr_cv}
\alias{lfqpcr_cv}
\title{Subject-level cross-validation for LFQPCR}
\description{
Estimates the prediction loss of LFQPCR by leaving out groups of subjects.
}
\usage{
lfqpcr_cv(Y, V, subject, visit_time, grid = NULL, covariates = NULL,
          tau_F = 0.5, tau_R = 0.5, K_subject = 2, K_visit = 1, spline_df = 8,
          spline_degree = 3, periodic = FALSE, mean_time_degree = 1,
          response_time_degree = NULL, n_folds = 5, seed = 1,
          lfqpca_max_iter = 30, lfqpca_tol = 1e-4, lfqpca_n_starts = 2,
          lfqpca_rq_method_functions = c("fnb", "fn", "br"),
          lfqpca_rq_method_scores = c("br", "fnb", "fn"),
          strict_identification = TRUE, objective_increase_tolerance = 1e-7,
          lfqpcr_rq_method = c("fnb", "br", "fn"), space_tol = 1e-10,
          rank_tol = 1e-7, verbose = TRUE)
}
\arguments{
\item{Y, V, subject, visit_time, grid, covariates}{the data, as in
\code{lfqpcr}.}
\item{tau_F, tau_R, K_subject, K_visit, spline_df, spline_degree, periodic}{
settings of \code{lfqpcr}.}
\item{mean_time_degree, response_time_degree, strict_identification}{settings
of \code{lfqpcr}.}
\item{lfqpca_max_iter, lfqpca_tol, lfqpca_n_starts}{settings of the first
stage.}
\item{lfqpca_rq_method_functions, lfqpca_rq_method_scores, lfqpcr_rq_method}{
quantile regression algorithms of the two stages.}
\item{objective_increase_tolerance, space_tol, rank_tol}{tolerances, as in
\code{lfqpcr}.}
\item{n_folds}{number of folds, between 2 and the number of subjects.}
\item{seed}{seed for the allocation of subjects to folds. The random number
state of the session is restored afterwards.}
\item{verbose}{if \code{TRUE}, progress is reported by messages.}
}
\details{
The subjects are allocated at random to \code{n_folds} folds. For each fold,
\code{lfqpcr} is fitted to the other subjects, and the responses of the
held-out subjects are predicted from their curves with the predict method.
The criterion \code{cv_objective} is the subject-balanced check loss at level
\code{tau_R} over all held-out predictions. Comparing it over \code{tau_F}, the
numbers of components or \code{spline_df} gives a choice of these settings.
}
\value{
A list with components \code{tau_F}, \code{tau_R}, \code{n_folds},
\code{fold_id} (fold of each row), \code{prediction} (held-out predictions),
\code{eligible} (rows used in the criterion), \code{response_weights} and
\code{cv_objective}.
}
\examples{
\dontrun{
# Choice of tau_F by five-fold cross-validation in Case 3
sim <- lfqpcr_simulate("Case3", I_train = 60, I_test = 20, J = 5,
                       tau_F = 0.9, tau_R = 0.9, replication = 1)
tr <- sim$train
tau_F_grid <- c(0.1, 0.5, 0.9)
cv_loss <- sapply(tau_F_grid, function(tf) {
  lfqpcr_cv(sim$Y[tr, ], sim$V[tr], sim$subject[tr], sim$visit_time[tr],
            sim$grid, covariates = sim$C[tr, ], tau_F = tf, tau_R = 0.9,
            K_subject = 2, K_visit = 2, spline_df = 10, mean_time_degree = 2,
            response_time_degree = 2, n_folds = 5, seed = 1,
            lfqpca_max_iter = 60, lfqpca_tol = 1e-5, lfqpca_n_starts = 5,
            verbose = FALSE)$cv_objective
})
tau_F_grid[which.min(cv_loss)]
}
}
