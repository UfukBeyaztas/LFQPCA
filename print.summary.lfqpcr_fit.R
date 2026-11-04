print.summary.lfqpcr_fit <- function(x, ...) {
  cat("LFQPCR summary\n")
  cat("tau_F:", x$tau_F, "  tau_R:", x$tau_R, "\n")
  cat("Contributing subjects:", x$I_R, "\n")
  cat("Contributing visits:", x$n_stage2, "\n")
  cat("r_B:", x$r_B, "  r_U:", x$r_U, "\n")
  cat("Stage-2 design rank:", x$design_rank, "\n")
  cat("LFQPCR objective:", format(x$objective, digits = 8), "\n")
  cat(
    "Subject-balanced empirical coverage:",
    format(x$weighted_coverage, digits = 6),
    "\n"
  )
  cat(
    "B-space orthonormality error:",
    format(x$B_orthonormality_error, scientific = TRUE),
    "\n"
  )
  cat(
    "U-space orthonormality error:",
    format(x$U_orthonormality_error, scientific = TRUE),
    "\n"
  )
  
  if (length(x$dropped_columns)) {
    cat(
      "Exact dependent columns removed:",
      paste(x$dropped_columns, collapse = ", "),
      "\n"
    )
  }
  
  invisible(x)
}