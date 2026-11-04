\name{lfqpca_cv}
\alias{lfqpca_cv}
\title{Grid-point cross-validation for LFQPCA}
\description{
Compares numbers of components and spline dimensions by holding out grid points
within the curves.
}
\usage{
lfqpca_cv(Y, subject, visit_time, grid = NULL,
          tuning_grid = expand.grid(K_subject = c(1, 2), K_visit = c(1, 2),
                                    spline_df = c(6, 8),
                                    KEEP.OUT.ATTRS = FALSE),
          tau = 0.5, n_folds = 3, max_iter = 20, n_starts = 1, verbose = TRUE,
          ...)
}
\arguments{
\item{Y, subject, visit_time, grid}{the data, as in \code{lfqpca}.}
\item{tuning_grid}{data frame with columns \code{K_subject}, \code{K_visit} and
\code{spline_df}, one row per candidate model.}
\item{tau}{quantile level.}
\item{n_folds}{number of folds.}
\item{max_iter, n_starts}{passed to \code{lfqpca}.}
\item{verbose}{if \code{TRUE}, progress is reported by messages.}
\item{...}{further arguments of \code{lfqpca}, for example
\code{mean_time_degree}.}
}
\details{
The observed points of every curve are split at random into \code{n_folds}
folds of nearly equal size. For each candidate and each fold, the model is
fitted with the points of the fold removed, and the held-out loss is the check
loss at these points with the same subject-balanced weights as the criterion of
\code{lfqpca}. Since only grid points are held out, every subject and visit
takes part in every fit. Each curve needs at least \code{n_folds + 5} observed
points. The folds are drawn with the current random number generator, so
\code{set.seed} makes the result reproducible.
}
\value{
A list with components
\item{fold_results}{data frame with the held-out loss of every candidate and
fold.}
\item{summary}{mean and standard deviation of the held-out loss of each
candidate, ordered by the mean.}
\item{best}{first row of \code{summary}.}
\item{fold_id}{matrix with the fold of every observed grid point.}
}
\examples{
\dontrun{
dat <- lfqpca_simulate("D1", I = 60, J = 5, replication = 1)
tuning <- expand.grid(K_subject = 1:3, K_visit = 1:2, spline_df = c(8, 10),
                      KEEP.OUT.ATTRS = FALSE)
set.seed(1)
cv <- lfqpca_cv(dat$Y, dat$subject, dat$visit_time, dat$grid,
                tuning_grid = tuning, tau = 0.5, n_folds = 3, max_iter = 20,
                verbose = FALSE, mean_time_degree = 2)
cv$summary
cv$best
}
}
