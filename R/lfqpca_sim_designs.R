lfqpca_sim_designs <- function() {
  data.frame(
    experiment = rep(c("I", "II"), c(7, 6)),
    code = c(paste0("D", 1:7), paste0("Case", 1:6)),
    description = c(
      "Time effects; quantile-invariant subject and visit parts",
      "No time effect in the population and subject parts",
      "Quantile-specific scores",
      "Quantile-specific loading spaces",
      "Visit scores violating H_i' Z_i = 0",
      "Unbalanced visits, J_i uniform on 4, ..., 9",
      "Subject curvature outside the model",
      "One functional coefficient acting on B + U",
      "Separate coefficients for B and U",
      "Quantile-specific scores, response at tau* = 0.9",
      "No time effect",
      "Visit scores violating H_i' Z_i = 0",
      "Unbalanced visits, J_i uniform on 4, ..., 9"
    ),
    stringsAsFactors = FALSE
  )
}
