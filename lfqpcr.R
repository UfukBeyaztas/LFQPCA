lfqpcr <- function(
    Y,
    V,
    subject,
    visit_time,
    grid = NULL,
    covariates = NULL,
    tau_F = 0.5,
    tau_R = 0.5,
    K_subject = 2,
    K_visit = 1,
    spline_df = 8,
    spline_degree = 3,
    periodic = FALSE,
    mean_time_degree = 1,
    response_time_degree = NULL,
    lfqpca_max_iter = 30,
    lfqpca_tol = 1e-4,
    lfqpca_n_starts = 2,
    lfqpca_rq_method_functions = c("fnb", "fn", "br"),
    lfqpca_rq_method_scores = c("br", "fnb", "fn"),
    strict_identification = TRUE,
    objective_increase_tolerance = 1e-7,
    lfqpcr_rq_method = c("fnb", "br", "fn"),
    space_tol = 1e-10,
    rank_tol = 1e-7,
    verbose = TRUE) {
  
  lfqpca_rq_method_functions <- match.arg(lfqpca_rq_method_functions)
  lfqpca_rq_method_scores <- match.arg(lfqpca_rq_method_scores)
  lfqpcr_rq_method <- match.arg(lfqpcr_rq_method)
  
  fit_F <- lfqpca(
    Y = Y,
    subject = subject,
    visit_time = visit_time,
    grid = grid,
    tau = tau_F,
    K_subject = K_subject,
    K_visit = K_visit,
    spline_df = spline_df,
    spline_degree = spline_degree,
    periodic = periodic,
    mean_time_degree = mean_time_degree,
    max_iter = lfqpca_max_iter,
    tol = lfqpca_tol,
    n_starts = lfqpca_n_starts,
    rq_method_functions = lfqpca_rq_method_functions,
    rq_method_scores = lfqpca_rq_method_scores,
    strict_identification = strict_identification,
    objective_increase_tolerance = objective_increase_tolerance,
    verbose = verbose
  )
  
  lfqpcr_from_lfqpca(
    lfqpca_fit = fit_F,
    V = V,
    covariates = covariates,
    tau_R = tau_R,
    response_time_degree = response_time_degree,
    rq_method = lfqpcr_rq_method,
    space_tol = space_tol,
    rank_tol = rank_tol
  )
}