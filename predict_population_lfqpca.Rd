\name{predict_population_lfqpca}
\alias{predict_population_lfqpca}
\title{Population quantile curves at given visit times}
\description{
Evaluates the population part \eqn{H(u)^\top \eta(s)}{H(u)' eta(s)} of a fitted
LFQPCA model at given visit times.
}
\usage{
predict_population_lfqpca(object, new_visit_time)
}
\arguments{
\item{object}{an object of class \code{lfqpca_fit}.}
\item{new_visit_time}{visit times on the scale of the original
\code{visit_time}.}
}
\value{
A matrix with one row per visit time and one column per grid point.
}
\examples{
\dontrun{
dat <- lfqpca_simulate("D1", I = 60, J = 5, replication = 1)
fit <- lfqpca(dat$Y, dat$subject, dat$visit_time, dat$grid, tau = 0.5,
              K_subject = 2, K_visit = 2, spline_df = 10,
              mean_time_degree = 2, max_iter = 60, tol = 1e-5, n_starts = 5,
              rq_method_functions = "fnb", rq_method_scores = "br",
              verbose = FALSE)

# Population medians at visit times -1, 0 and 1
pop <- predict_population_lfqpca(fit, new_visit_time = c(-1, 0, 1))
matplot(dat$grid, t(pop), type = "l", lty = 1:3, col = 1,
        xlab = "s", ylab = "population median")
}
}
