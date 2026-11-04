\name{lfqpca_sim_settings}
\alias{lfqpca_sim_settings}
\alias{lfqpca_sim_designs}
\title{Constants and list of the simulation designs}
\description{
\code{lfqpca_sim_settings} returns the constants of the simulation designs and
\code{lfqpca_sim_designs} lists the designs of both experiments.
}
\usage{
lfqpca_sim_settings()

lfqpca_sim_designs()
}
\details{
The main constants are the number of grid points \code{D} (61), the numbers of
subject and visit components (2 each), the degree of the time basis (2), the
spline dimension used in the fits (10), the range of the visit times
(\code{time_range}), the numbers of visits in the unbalanced designs
(\code{unbalanced_visits}, 4 to 9), the scale intercepts of the designs with
a common scale and of D4 (\code{base_scale_intercept},
\code{loading_scale_intercept}) and the lower
bound \code{min_scale} of the scale surfaces.

In D3 and Case 3, the scale scores are chosen so that the first subject scores
at the reference and the calibration levels (\code{reference_tau},
\code{calibration_tau}) have the correlation \code{score_correlation}; the first visit score is multiplied by
\code{visit_signal}, and the scale intercept keeps the scale above \code{min_scale} for
\eqn{|u| \le}{|u| <=} \code{max_abs_time} and at most \code{bound_visits}
visits. \code{loading_norms} are the norms of the D4 scale loadings and
\code{curvature} is the factor of the D7 quadratic term. In Case 3 the response
depends on the components at the level \code{response_tau_star}. The seeds of the two experiments are \code{seed_experiment1} and
\code{seed_experiment2}; \code{iid_seed_offset}, \code{reserve_seed_offset} and
\code{reserve_seeds} set the seeds of the independent visit scores and of the
reserve draws.
}
\value{
\code{lfqpca_sim_settings} returns a named list. \code{lfqpca_sim_designs}
returns a data frame with columns \code{experiment}, \code{code} and
\code{description}.
}
\examples{
\dontrun{
lfqpca_sim_designs()
settings <- lfqpca_sim_settings()
str(settings[c("D", "K_subject", "K_visit", "spline_df", "min_scale")])

# Curves on a coarser grid
settings$D <- 31
dat <- lfqpca_simulate("D1", I = 30, settings = settings)
dim(dat$Y)
}
}
