\name{lfqpca_simulate}
\alias{lfqpca_simulate}
\alias{lfqpca_sim_data}
\title{Simulation designs of Experiment I}
\description{
Generates data from the seven designs used to study LFQPCA in the paper.
}
\usage{
lfqpca_simulate(design = paste0("D", 1:7), I = 60, J = 5, replication = 1,
                seed = NULL, missing_rate = 0.05,
                settings = lfqpca_sim_settings())

lfqpca_sim_data(design = paste0("D", 1:7), seed, I = 60, J = 5,
                missing_rate = 0.05, settings = lfqpca_sim_settings(),
                scale_only = FALSE)
}
\arguments{
\item{design}{one of \code{"D1"} to \code{"D7"}; see Details.}
\item{I}{number of subjects.}
\item{J}{number of visits per subject; ignored in design D6.}
\item{replication}{replication number, which sets the seed.}
\item{seed}{seed of the data set. In \code{lfqpca_simulate} it overrides
\code{replication}.}
\item{missing_rate}{proportion of curve values removed at random.}
\item{settings}{design constants, from \code{lfqpca_sim_settings}.}
\item{scale_only}{if \code{TRUE}, only the scale surface is generated and its
minimum is returned. This is used by the reserve-seed rule and by
\code{lfqpca_sim_calibrate}.}
}
\details{
Every subject has \code{J} visits (in D6, a number drawn from 4 to 9), with
times drawn uniformly on [-1.2, 1.2] and sorted. The curves are observed at 61
equally spaced points of [0, 1] and follow the location-scale model
\deqn{Y_{ij}(s) = L_{ij}(s) + S_{ij}(s)\,\varepsilon_{ij}(s),}{Y_ij(s) = L_ij(s) + S_ij(s) eps_ij(s),}
with independent standard normal errors, so that the
\eqn{\tau}{tau}-quantile is \eqn{L_{ij} + \Phi^{-1}(\tau) S_{ij}}{L_ij + qnorm(tau) S_ij}.
The location \eqn{L}{L} has the LFQPCA form with a quadratic time basis, two
subject components and two visit components, orthonormal trigonometric
loadings, centered normal subject scores with standard deviations 1.25 and 0.70,
and visit scores \eqn{N_i A_i}{N_i A_i} with standard deviations 0.75 and 0.40.
The population functions are
\eqn{1.2 + 0.8\sin(2\pi s) + 0.2\cos(4\pi s)}{1.2 + 0.8 sin(2 pi s) + 0.2 cos(4 pi s)}
for the intercept,
\eqn{0.35\cos(2\pi s) + 0.2(s - 0.5)}{0.35 cos(2 pi s) + 0.2 (s - 0.5)} for
\eqn{u}{u} and \eqn{-0.18\sin(\pi s)}{-0.18 sin(pi s)} for \eqn{u^2}{u^2}.
\describe{
\item{D1}{time effects, with subject and visit components common to all
quantile levels; \eqn{S(s) = 0.75 + 0.15\cos(2\pi s)}{S(s) = 0.75 + 0.15 cos(2 pi s)}.}
\item{D2}{as D1 without time effects: the population functions of
\eqn{u}{u} and \eqn{u^2}{u^2} and the slope loadings are zero.}
\item{D3}{quantile-specific scores. The scale has its own first subject and
visit scores on the location loadings, so the scores at level
\eqn{\tau}{tau} are \eqn{\xi_i + \Phi^{-1}(\tau)\xi_i^S}{xi_i + qnorm(tau) xi_i^S},
and similarly for the visit scores. The first location visit score is
doubled, and the scale intercept keeps \eqn{S}{S} above \code{min_scale} when
\eqn{|u| \le 2.2}{|u| <= 2.2}.}
\item{D4}{quantile-specific loadings. The scale uses the location scores with
new loading directions orthogonal to the location loadings, and the intercept
3 from \code{lfqpca_sim_calibrate}.}
\item{D5}{visit scores drawn independently, which violates
\eqn{H_i^\top Z_i = 0}{H_i' Z_i = 0}.}
\item{D6}{unbalanced visits.}
\item{D7}{subject curvature outside the model: the term
\eqn{\xi_i^\top \phi_2(s) u_{ij}^2}{xi_i' phi_2(s) u_ij^2} with
\eqn{\phi_2 = 1.13\,\phi_1}{phi_2 = 1.13 phi_1} is added.}
}
\code{lfqpca_simulate} uses the seed \code{seed_experiment1 + replication}
for every design. In D3 and D4, a draw
whose scale surface falls below \code{min_scale} is replaced by the draw with
seed \code{seed + k * reserve_seed_offset} for the smallest k that gives an
admissible surface, with k at most \code{reserve_seeds}; \code{reserve_draws}
records k. \code{lfqpca_sim_data} generates one data set
from a given seed, without this rule.
}
\value{
A list with components \code{design}, \code{seed}, \code{Y} (curves with
missing values), \code{Y_complete}, \code{subject}, \code{visit_time},
\code{time_std}, \code{grid}, \code{weights} (trapezoidal weights),
\code{visits} (visits per subject) and \code{truth} (generating functions,
scores and surfaces of the location and scale parts). \code{lfqpca_simulate}
adds \code{reserve_draws}. With \code{scale_only = TRUE}, a short list with the
scale intercept and the minimum of the scale surface is returned.
}
\examples{
\dontrun{
# Design D3, replication 1
dat <- lfqpca_simulate("D3", I = 60, J = 5, replication = 1)
dim(dat$Y)
mean(is.na(dat$Y))

# One replication of each design
sims <- lapply(paste0("D", 1:7), lfqpca_simulate, I = 60, J = 5,
               replication = 1)
names(sims) <- paste0("D", 1:7)
sapply(sims, function(x) nrow(x$Y))

# Unbalanced visits, and the same data from lfqpca_sim_data
dat6 <- lfqpca_simulate("D6", I = 120, replication = 3)
table(dat6$visits)
seed3 <- lfqpca_sim_settings()$seed_experiment1 + 3
dat6b <- lfqpca_sim_data("D6", seed = seed3, I = 120)
identical(dat6$Y, dat6b$Y)
}
}
