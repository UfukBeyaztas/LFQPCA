\name{lfqpcr}
\alias{lfqpcr}
\alias{print.lfqpcr_fit}
\alias{summary.lfqpcr_fit}
\alias{print.summary.lfqpcr_fit}
\alias{coef.lfqpcr_fit}
\alias{fitted.lfqpcr_fit}
\alias{residuals.lfqpcr_fit}
\title{Longitudinal functional quantile principal component regression}
\description{
Fits LFQPCA to the curves at level \eqn{\tau_F}{tau_F} and then the quantile
regression of a scalar response at level \eqn{\tau_R}{tau_R} on the persistent
and transient components.
}
\usage{
lfqpcr(Y, V, subject, visit_time, grid = NULL, covariates = NULL, tau_F = 0.5,
       tau_R = 0.5, K_subject = 2, K_visit = 1, spline_df = 8,
       spline_degree = 3, periodic = FALSE, mean_time_degree = 1,
       response_time_degree = NULL, lfqpca_max_iter = 30, lfqpca_tol = 1e-4,
       lfqpca_n_starts = 2, lfqpca_rq_method_functions = c("fnb", "fn", "br"),
       lfqpca_rq_method_scores = c("br", "fnb", "fn"),
       strict_identification = TRUE, objective_increase_tolerance = 1e-7,
       lfqpcr_rq_method = c("fnb", "br", "fn"), space_tol = 1e-10,
       rank_tol = 1e-7, verbose = TRUE)

\method{print}{lfqpcr_fit}(x, ...)

\method{summary}{lfqpcr_fit}(object, ...)

\method{print}{summary.lfqpcr_fit}(x, ...)

\method{coef}{lfqpcr_fit}(object, ...)

\method{fitted}{lfqpcr_fit}(object, ...)

\method{residuals}{lfqpcr_fit}(object, ...)
}
\arguments{
\item{Y, subject, visit_time, grid}{the curves and their design, as in
\code{lfqpca}.}
\item{V}{scalar responses, one for each row of \code{Y}; missing values are
allowed.}
\item{covariates}{optional numeric matrix or data frame of scalar covariates,
one row for each row of \code{Y}.}
\item{tau_F}{quantile level of the curves.}
\item{tau_R}{quantile level of the response.}
\item{K_subject, K_visit, spline_df, spline_degree, periodic, mean_time_degree}{
dimensions and bases of the first stage, as in \code{lfqpca}.}
\item{response_time_degree, space_tol, rank_tol}{as in
\code{lfqpcr_from_lfqpca}.}
\item{lfqpca_max_iter, lfqpca_tol, lfqpca_n_starts}{\code{max_iter}, \code{tol}
and \code{n_starts} of the first stage.}
\item{lfqpca_rq_method_functions, lfqpca_rq_method_scores}{quantile regression
algorithms of the first stage.}
\item{strict_identification, objective_increase_tolerance}{as in
\code{lfqpca}.}
\item{lfqpcr_rq_method}{quantile regression algorithm of the second stage.}
\item{verbose}{if \code{TRUE}, progress of the first stage is reported.}
\item{x, object}{an object of class \code{lfqpcr_fit}, or of class
\code{summary.lfqpcr_fit} for the print method of the summary.}
\item{...}{not used.}
}
\details{
The first stage calls \code{lfqpca} directly, without quantile continuation.
For levels \eqn{\tau_F}{tau_F} far from the median it can be better to fit the
first stage with \code{lfqpca_continuation} and pass the result to
\code{lfqpcr_from_lfqpca}, which is how the simulations of the paper were run.
The model and the estimator of the second stage are described with
\code{lfqpcr_from_lfqpca}.
}
\value{
An object of class \code{lfqpcr_fit}, as returned by \code{lfqpcr_from_lfqpca}.
\code{coef}, \code{fitted} and \code{residuals} return the regression
coefficients, the fitted \eqn{\tau_R}{tau_R}-quantiles and the residuals of
\code{V}. \code{summary} returns the quantile levels, the number of subjects
with a response, the dimensions of the loading spaces, the rank of the design,
the number of observations, the criterion, the subject-balanced proportion of
responses at or below the fitted quantiles and the coefficients.
}
\examples{
\dontrun{
# Experiment II, Case 2, with tau_F = 0.5 and tau_R = 0.9
sim <- lfqpcr_simulate("Case2", I_train = 60, I_test = 20, J = 5,
                       tau_F = 0.5, tau_R = 0.9, replication = 1)
tr <- sim$train
te <- sim$test
fit <- lfqpcr(sim$Y[tr, ], sim$V[tr], sim$subject[tr], sim$visit_time[tr],
              sim$grid, covariates = sim$C[tr, ], tau_F = 0.5, tau_R = 0.9,
              K_subject = 2, K_visit = 2, spline_df = 10,
              mean_time_degree = 2, response_time_degree = 2,
              lfqpca_max_iter = 60, lfqpca_tol = 1e-5, lfqpca_n_starts = 5,
              lfqpca_rq_method_functions = "fnb",
              lfqpca_rq_method_scores = "br", verbose = FALSE)
fit
summary(fit)
head(fitted(fit))
}
}
