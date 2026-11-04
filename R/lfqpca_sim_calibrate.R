lfqpca_sim_calibrate <- function(n_rep = 500, I = 120, J = c(5, 10),
                                 start = 0.75, margin = 0.25,
                                 rounding = 0.1, seed = NULL,
                                 settings = lfqpca_sim_settings()) {
  if (is.null(seed)) seed <- settings$seed_experiment1 + 1e6
  s <- settings
  s$loading_scale_intercept <- start
  by_J <- do.call(rbind, lapply(J, function(J_k) {
    required <- vapply(seq_len(n_rep), function(b) {
      probe <- lfqpca_sim_data("D4", seed + 10000 * J_k + b - 1, I = I,
                               J = J_k, missing_rate = 0, settings = s,
                               scale_only = TRUE)
      probe$required_intercept
    }, numeric(1))
    frozen <- ceiling((max(required) + margin) / rounding - 1e-12) * rounding
    data.frame(J = J_k, max_required = max(required),
               intercept = max(start, frozen))
  }))
  list(intercept = max(by_J$intercept), by_J = by_J)
}
