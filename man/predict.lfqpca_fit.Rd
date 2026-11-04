\name{predict.lfqpca_fit}
\alias{predict.lfqpca_fit}
\title{Scores and quantile curves of new subjects}
\description{
Estimates the scores of new subjects from their curves, with the functions of a
fitted LFQPCA model held fixed, and returns the fitted quantile curves.
}
\usage{
\method{predict}{lfqpca_fit}(object, newY, new_subject, new_visit_time, ...)
}
\arguments{
\item{object}{an object of class \code{lfqpca_fit}.}
\item{newY}{matrix of curves of the new subjects on the grid of the fit, one per
row; missing values are allowed.}
\item{new_subject}{subject labels of the rows of \code{newY}.}
\item{new_visit_time}{visit times of the rows of \code{newY}, on the scale of
the original \code{visit_time}.}
\item{...}{not used.}
}
\details{
The visit times are standardized with the mean and standard deviation of the
training times. For each new subject, the scores \eqn{(\xi, A)}{(xi, A)}
minimize the weighted check loss of its curves with
\eqn{\eta, \phi_0, \phi_1, \psi}{eta, phi_0, phi_1, psi} fixed, and the visit
scores are \eqn{N_0 A}{N_0 A}, with \eqn{N_0}{N_0} a basis of the null space of
the transposed time design of the subject. Each reconstructed curve therefore
uses all curves of its subject. When the model has visit components, every new
subject needs at least \eqn{p + 2}{p + 2} distinct visit times.
}
\value{
A list with components \code{fitted} (matrix of quantile curves), \code{xi}
(subject scores, one row per new subject), \code{zeta} (visit scores),
\code{subject_id}, \code{subject_levels}, \code{time_std} (standardized visit
times), \code{objective} (criterion of the new curves) and
\code{constraint_diagnostics}.
}
\examples{
\dontrun{
dat <- lfqpca_simulate("D1", I = 60, J = 5, replication = 1)
fit <- lfqpca(dat$Y, dat$subject, dat$visit_time, dat$grid, tau = 0.5,
              K_subject = 2, K_visit = 2, spline_df = 10,
              mean_time_degree = 2, max_iter = 60, tol = 1e-5, n_starts = 5,
              rq_method_functions = "fnb", rq_method_scores = "br",
              verbose = FALSE)

# Twenty new subjects from the same design
new <- lfqpca_simulate("D1", I = 20, J = 5, replication = 2)
pred <- predict(fit, newY = new$Y, new_subject = new$subject,
                new_visit_time = new$visit_time)
dim(pred$fitted)
head(pred$xi)
}
}
