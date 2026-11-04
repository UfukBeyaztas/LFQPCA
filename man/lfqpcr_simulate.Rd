\name{lfqpcr_simulate}
\alias{lfqpcr_simulate}
\alias{lfqpcr_sim_data}
\title{Simulation designs of Experiment II}
\description{
Generates curves, scalar responses and covariates from the six cases used to
study LFQPCR in the paper.
}
\usage{
lfqpcr_simulate(case = paste0("Case", 1:6), I_train = 60, I_test = 20, J = 5,
                tau_F = 0.5, tau_R = 0.9, replication = 1, seed = NULL,
                settings = lfqpca_sim_settings())

lfqpcr_sim_data(case = paste0("Case", 1:6), seed, I_train = 60, I_test = 20,
                J = 5, tau_F = 0.5, tau_R = 0.9,
                settings = lfqpca_sim_settings())
}
\arguments{
\item{case}{one of \code{"Case1"} to \code{"Case6"}; see Details.}
\item{I_train, I_test}{numbers of training and test subjects.}
\item{J}{number of visits per subject; ignored in Case 6.}
\item{tau_F}{quantile level of the curves at which the true components are
returned.}
\item{tau_R}{quantile level of the response.}
\item{replication}{replication number, which sets the seed.}
\item{seed}{seed of the data set. In \code{lfqpcr_simulate} it overrides
\code{replication}.}
\item{settings}{design constants, from \code{lfqpca_sim_settings}.}
}
\details{
Visit times and grid are as in Experiment I, the curves are complete, and time
is standardized with the visits of the training subjects. The location part of
the curves has the form used in Experiment I, with subject scores centered over
the training subjects. Except in Case 3, the curves are
\deqn{Y_{ij}(s) = L_{ij}(s) + \sigma_{ij}(s)\{E_{ij}(s) - 1\}, \qquad \sigma_{ij}(s) = (0.55 + 0.2u_{ij}^2)(0.25 + 0.45s),}{Y_ij(s) = L_ij(s) + sigma_ij(s) (E_ij(s) - 1),   sigma_ij(s) = (0.55 + 0.2 u_ij^2)(0.25 + 0.45 s),}
with standard exponential errors, so the persistent and transient components
are the same at every level \eqn{\tau_F}{tau_F}. The response is
\deqn{V_{ij} = \alpha(u_{ij}) + C_{ij}^\top \beta + \langle \gamma_B, B_{ij} \rangle + \langle \gamma_U, U_{ij} \rangle + \sigma_R(u_{ij})\{E_{ij} - q_E(\tau_R)\},}{V_ij = alpha(u_ij) + C_ij' beta + <gamma_B, B_ij> + <gamma_U, U_ij> + sigma_R(u_ij) (E_ij - qexp(tau_R)),}
with \eqn{\alpha(u) = 0.5 + 0.4u - 0.15u^2}{alpha(u) = 0.5 + 0.4 u - 0.15 u^2},
a subject-level and a visit-level standard normal covariate with
\eqn{\beta = (0.5, -0.35)}{beta = (0.5, -0.35)},
\eqn{\sigma_R(u) = 0.55 + 0.15|u|}{sigma_R(u) = 0.55 + 0.15 |u|}, a standard
exponential error \eqn{E_{ij}}{E_ij} and its quantile function
\eqn{q_E}{q_E}. The \eqn{\tau_R}{tau_R}-quantile of \eqn{V_{ij}}{V_ij} is thus
the linear predictor.
\describe{
\item{Case1}{one coefficient function acting on \eqn{B + U}{B + U}; its
projections on the two loading spaces give \eqn{\gamma_B}{gamma_B} and
\eqn{\gamma_U}{gamma_U}.}
\item{Case2}{separate coefficient functions for \eqn{B}{B} and \eqn{U}{U}.}
\item{Case3}{normal errors with quantile-specific scores as in design D3; the
components change with \eqn{\tau_F}{tau_F}, and the response depends on the
components at level 0.9.}
\item{Case4}{no time effect in the curves, the error scales or
\eqn{\alpha}{alpha}.}
\item{Case5}{visit scores drawn independently, which violates
\eqn{H_i^\top Z_i = 0}{H_i' Z_i = 0}.}
\item{Case6}{unbalanced visits, 4 to 9 per subject.}
}
Cases 2 to 6 share the coefficient functions. \code{lfqpcr_simulate} uses the
seed \code{seed_experiment2 + replication - 1}; the curves and responses do not
depend on \code{tau_F}.
}
\value{
A list with components \code{case}, \code{seed}, \code{Y}, \code{V}, \code{C}
(covariates), \code{subject}, \code{visit_time}, \code{grid}, \code{visits},
\code{train} and \code{test} (logical indicators of the rows), \code{tau_F},
\code{tau_R} and \code{truth}. The list \code{truth} holds the standardization
of time, the true functions and scores at \code{tau_F}, the true
\eqn{\tau_F}{tau_F}-quantile curves (\code{functional_quantile}), the
components that enter the response (\code{B_response}, \code{U_response}), the
true \code{alpha}, \code{beta}, \code{gamma_B} and \code{gamma_U}, and the true
response quantiles.
}
\examples{
\dontrun{
# Case 2, replication 1
sim <- lfqpcr_simulate("Case2", I_train = 60, I_test = 20, J = 5,
                       tau_F = 0.5, tau_R = 0.9, replication = 1)
table(train = sim$train)
str(sim$truth[c("beta", "gamma_B", "gamma_U")])

# One replication of each case at tau_F = 0.9
sims <- lapply(paste0("Case", 1:6), lfqpcr_simulate, tau_F = 0.9,
               tau_R = 0.9, replication = 1)
names(sims) <- paste0("Case", 1:6)
sapply(sims, function(x) length(x$V))
}
}
