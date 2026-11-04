lfqpca_sim_with_seed <- function(seed, expr) {
  restore <- lfqpcr_rng_state_saver()
  on.exit(restore(), add = TRUE)
  set.seed(seed)
  expr
}