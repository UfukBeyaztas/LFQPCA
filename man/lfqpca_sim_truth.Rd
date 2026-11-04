\name{lfqpca_sim_truth}
\alias{lfqpca_sim_truth}
\title{True components of an Experiment I data set}
\description{
Returns the true population surface, components and quantile curves of a
simulated data set at a given quantile level.
}
\usage{
lfqpca_sim_truth(data, tau, settings = lfqpca_sim_settings())
}
\arguments{
\item{data}{output of \code{lfqpca_simulate} or \code{lfqpca_sim_data}.}
\item{tau}{quantile level.}
\item{settings}{design constants, from \code{lfqpca_sim_settings}.}
}
\details{
The \eqn{\tau}{tau}-quantile is the location plus
\eqn{q = \Phi^{-1}(\tau) - \Phi^{-1}(0.5)}{q = qnorm(tau) - qnorm(0.5)} times
the scale, split into population, subject and visit parts. In D3 the scores
move with \eqn{q}{q}, in D4 the loadings do, and in the other designs only the
population part changes. In D7, \code{B} contains the quadratic subject term
and \code{B_intercept_slope} does not. The true functions need not satisfy the
normalization of \code{lfqpca}, so a fit is best compared with the surfaces
\code{population}, \code{B}, \code{U} and \code{quantile_surface}, or through
the spans of the loadings.
}
\value{
A list with the level \code{tau}, the shift \code{q_shift}, the functions
\code{eta}, \code{phi0}, \code{phi1} and \code{psi}, the scores \code{xi} and
\code{zeta}, and the surfaces \code{population}, \code{B_intercept_slope},
\code{B}, \code{U} and \code{quantile_surface}, which have one row per curve
and one column per grid point.
}
\examples{
\dontrun{
dat <- lfqpca_simulate("D3", I = 60, J = 5, replication = 1)
truth <- lfqpca_sim_truth(dat, tau = 0.9)
matplot(dat$grid, t(truth$quantile_surface[dat$subject == 1, ]), type = "l",
        lty = 1, xlab = "s", ylab = "0.9-quantile, subject 1")

# Share of the subject part that is due to the quadratic term in D7
dat7 <- lfqpca_simulate("D7", I = 60, J = 5, replication = 1)
truth7 <- lfqpca_sim_truth(dat7, tau = 0.5)
sum((truth7$B - truth7$B_intercept_slope)^2) / sum(truth7$B^2)
}
}
