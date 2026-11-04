lfqpca_function_basis <- function(grid, df = 8, degree = 3,
                                  periodic = FALSE) {
  grid <- as.numeric(grid)
  df <- as.integer(df)
  degree <- as.integer(degree)
  
  if (df < degree + 1) {
    stop("'spline_df' must be at least spline_degree + 1.", call. = FALSE)
  }
  if (df > length(grid)) {
    stop("'spline_df' cannot exceed the number of functional grid points.",
         call. = FALSE)
  }
  
  if (periodic) {
    if (!requireNamespace("pbs", quietly = TRUE)) {
      stop("Package 'pbs' is needed for periodic = TRUE.", call. = FALSE)
    }
    B <- pbs(
      grid,
      df = df,
      degree = degree,
      intercept = TRUE,
      periodic = TRUE,
      Boundary.knots = range(grid)
    )
  } else {
    B <- bs(
      grid,
      df = df,
      degree = degree,
      intercept = TRUE,
      Boundary.knots = range(grid)
    )
  }
  
  unname(as.matrix(B))
}