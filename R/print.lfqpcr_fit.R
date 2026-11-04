print.lfqpcr_fit <- function(x, ...) {
  cat("Longitudinal Functional Quantile Principal Component Regression\n")
  cat("Predictor quantile tau_F:", x$tau_F, "\n")
  cat("Response quantile tau_R:", x$tau_R, "\n")
  cat("Contributing subjects I_R:", x$I_R, "\n")
  cat("Subject-space dimension r_B:", x$r_B, "\n")
  cat("Visit-space dimension r_U:", x$r_U, "\n")
  cat("Stage-2 design rank:", x$design_rank, "\n")
  cat("LFQPCR objective:", format(x$objective, digits = 8), "\n")
  
  if (length(x$dropped_columns)) {
    cat(
      "Exact dependent columns removed:",
      paste(colnames(x$design)[x$dropped_columns], collapse = ", "),
      "\n"
    )
  }
  
  invisible(x)
}
