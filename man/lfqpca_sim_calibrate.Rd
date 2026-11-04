\name{lfqpca_sim_calibrate}
\alias{lfqpca_sim_calibrate}
\title{Calibration of the D4 scale intercept}
\description{
Pilot computation of the intercept of the scale surface in design D4.
}
\usage{
lfqpca_sim_calibrate(n_rep = 500, I = 120, J = c(5, 10), start = 0.75,
                     margin = 0.25, rounding = 0.1, seed = NULL,
                     settings = lfqpca_sim_settings())
}
\arguments{
\item{n_rep}{number of probe draws for each number of visits.}
\item{I}{number of subjects in the probe draws.}
\item{J}{numbers of visits; one pilot is run for each.}
\item{start}{intercept used in the probe draws.}
\item{margin}{amount added to the largest required intercept.}
\item{rounding}{the result is rounded up to a multiple of \code{rounding}.}
\item{seed}{base seed of the probe draws; by default
\code{seed_experiment1 + 1e6}.}
\item{settings}{design constants, from \code{lfqpca_sim_settings}.}
}
\details{
For each probe draw, the intercept needed to keep the D4 scale surface above
\code{min_scale} is \code{start} plus the shortfall of the surface. The value
for each \code{J} is the largest of these intercepts plus \code{margin}, rounded
up, and the result is the largest value over \code{J}. With the defaults the
largest required intercepts are 2.677 for five visits and 2.394 for ten
visits, which gives the value 3 used in the paper. The reserve-seed rule of
\code{lfqpca_simulate} covers the rare draws that still fall below the bound.
}
\value{
A list with the calibrated \code{intercept} and a data frame \code{by_J} with
the largest required intercept and the rounded value for each \code{J}.
}
\examples{
\dontrun{
cal <- lfqpca_sim_calibrate(n_rep = 500, I = 120, J = c(5, 10))
cal$by_J
cal$intercept
}
}
