lfqpca_sim_antithetic <- function(n, scale) {
  n_pair <- n %/% 2
  out <- numeric(0)
  if (n_pair > 0) {
    v <- sample(c(-1, 1), n_pair, replace = TRUE)
    out <- c(v, -v)
  }
  if (n %% 2 == 1) out <- c(out, 0)
  sample(out, length(out), replace = FALSE) * scale
}