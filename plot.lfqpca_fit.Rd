\name{plot.lfqpca_fit}
\alias{plot.lfqpca_fit}
\title{Plots of a fitted LFQPCA model}
\description{
Draws the loadings, the objective over the iterations, observed and fitted
curves, or population quantile curves of a fitted LFQPCA model.
}
\usage{
\method{plot}{lfqpca_fit}(x,
  type = c("components", "objective", "fitted", "mean_surface"),
  curve_indices = NULL, mean_times = c(-1, 0, 1), ...)
}
\arguments{
\item{x}{an object of class \code{lfqpca_fit}.}
\item{type}{what to draw; see Details.}
\item{curve_indices}{rows of \code{Y} drawn for \code{type = "fitted"}; by
default the first six.}
\item{mean_times}{standardized visit times at which the population curves are
drawn for \code{type = "mean_surface"}; 0 is the mean visit time.}
\item{...}{not used.}
}
\details{
\describe{
\item{\code{"components"}}{the subject intercept, subject slope and visit
loadings;}
\item{\code{"objective"}}{the objective over the iterations;}
\item{\code{"fitted"}}{observed and fitted curves of the rows in
\code{curve_indices};}
\item{\code{"mean_surface"}}{population quantile curves at the standardized
times \code{mean_times}.}
}
}
\value{
\code{x}, invisibly.
}
\examples{
\dontrun{
dat <- lfqpca_simulate("D1", I = 60, J = 5, replication = 1)
fit <- lfqpca(dat$Y, dat$subject, dat$visit_time, dat$grid, tau = 0.5,
              K_subject = 2, K_visit = 2, spline_df = 10,
              mean_time_degree = 2, max_iter = 60, tol = 1e-5, n_starts = 5,
              rq_method_functions = "fnb", rq_method_scores = "br",
              verbose = FALSE)
plot(fit, type = "components")
plot(fit, type = "objective")
plot(fit, type = "fitted", curve_indices = 1:4)
plot(fit, type = "mean_surface")
}
}
