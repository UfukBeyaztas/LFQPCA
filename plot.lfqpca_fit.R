plot.lfqpca_fit <- function(
    x,
    type = c("components", "objective", "fitted", "mean_surface"),
    curve_indices = NULL,
    mean_times = c(-1, 0, 1),
    ...) {
  
  type <- match.arg(type)
  
  if (type == "objective") {
    plot(
      seq_along(x$objective_history) - 1,
      x$objective_history,
      type = "b",
      pch = 19,
      xlab = "Iteration",
      ylab = "LFQPCA objective",
      main = "LFQPCA objective history"
    )
    return(invisible(x))
  }
  
  if (type == "components") {
    old <- par(no.readonly = TRUE)
    on.exit(par(old), add = TRUE)
    columns <- max(x$K_subject, x$K_visit, 1)
    par(mfrow = c(3, columns), mar = c(3.4, 3.6, 2.5, 0.8))
    
    for (k in seq_len(columns)) {
      if (k <= x$K_subject) {
        plot(
          x$grid, x$phi0[, k], type = "l", lwd = 2,
          xlab = "Functional domain", ylab = "Loading",
          main = paste0("Subject intercept ", k)
        )
      } else plot.new()
    }
    for (k in seq_len(columns)) {
      if (k <= x$K_subject) {
        plot(
          x$grid, x$phi1[, k], type = "l", lwd = 2,
          xlab = "Functional domain", ylab = "Loading",
          main = paste0("Subject slope ", k)
        )
      } else plot.new()
    }
    for (k in seq_len(columns)) {
      if (k <= x$K_visit) {
        plot(
          x$grid, x$psi[, k], type = "l", lwd = 2,
          xlab = "Functional domain", ylab = "Loading",
          main = paste0("Visit component ", k)
        )
      } else plot.new()
    }
    return(invisible(x))
  }
  
  if (type == "fitted") {
    if (is.null(curve_indices)) curve_indices <- seq_len(min(6, nrow(x$Y)))
    old <- par(no.readonly = TRUE)
    on.exit(par(old), add = TRUE)
    par(mfrow = c(ceiling(length(curve_indices) / 2), 2),
        mar = c(3.4, 3.6, 2.5, 0.8))
    
    for (r in curve_indices) {
      yrange <- range(x$Y[r, ], x$fitted[r, ], na.rm = TRUE)
      plot(
        x$grid, x$Y[r, ],
        pch = 16, cex = 0.45, ylim = yrange,
        xlab = "Functional domain", ylab = "Y",
        main = paste0(
          "Curve ", r, ", subject ", x$subject_levels[x$subject_id[r]],
          ", u=", round(x$time_std[r], 2)
        )
      )
      lines(x$grid, x$fitted[r, ], lwd = 2)
    }
    return(invisible(x))
  }
  
  u <- as.numeric(mean_times)
  Hnew <- lfqpca_time_basis(u, degree = x$mean_time_degree)
  eta <- lfqpca_mean_surface(Hnew, x$eta_functions)
  matplot(
    x$grid,
    t(eta),
    type = "l",
    lty = seq_along(u),
    lwd = 2,
    xlab = "Functional domain",
    ylab = "Population quantile",
    main = paste0("Population quantile surface, tau = ", x$tau)
  )
  legend(
    "topright",
    legend = paste0("u = ", u),
    lty = seq_along(u),
    lwd = 2,
    bty = "n"
  )
  invisible(x)
}