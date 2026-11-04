lfqpca_simulate <- function(design = paste0("D", 1:7), I = 60, J = 5,
                            replication = 1, seed = NULL, missing_rate = 0.05,
                            settings = lfqpca_sim_settings()) {
  design <- match.arg(design)
  if (is.null(seed)) seed <- settings$seed_experiment1 + replication
  if (!(design %in% c("D3", "D4"))) {
    dat <- lfqpca_sim_data(design, seed, I = I, J = J,
                           missing_rate = missing_rate, settings = settings)
    dat$reserve_draws <- 0
    return(dat)
  }
  
  for (k in 0:settings$reserve_seeds) {
    seed_k <- seed + k * settings$reserve_seed_offset
    probe <- lfqpca_sim_data(design, seed_k, I = I, J = J,
                             settings = settings, scale_only = TRUE)
    if (is.finite(probe$scale_min) &&
        probe$scale_min >= settings$min_scale - 1e-10) {
      dat <- lfqpca_sim_data(design, seed_k, I = I, J = J,
                             missing_rate = missing_rate, settings = settings)
      dat$reserve_draws <- k
      return(dat)
    }
  }
  stop("No admissible ", design, " draw among ", settings$reserve_seeds,
       " reserve seeds.", call. = FALSE)
}
