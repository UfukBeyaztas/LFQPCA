lfqpca_validate_data <- function(Y, subject, visit_time, grid = NULL,
                                 K_visit = 1,
                                 strict_identification = TRUE,
                                 min_curve_points = 5) {
  Y <- as.matrix(Y)
  storage.mode(Y) <- "double"
  n <- nrow(Y)
  D <- ncol(Y)
  
  if (n < 2 || D < min_curve_points) {
    stop(
      "Y must contain at least two curves and at least ", min_curve_points,
      " functional grid points.",
      call. = FALSE
    )
  }
  if (length(subject) != n || length(visit_time) != n) {
    stop("Lengths of 'subject' and 'visit_time' must equal nrow(Y).",
         call. = FALSE)
  }
  if (anyNA(subject)) {
    stop("'subject' cannot contain missing values.", call. = FALSE)
  }
  if (any(!is.finite(visit_time))) {
    stop("'visit_time' must contain only finite values.", call. = FALSE)
  }
  if (all(is.na(Y))) {
    stop("Y contains no observed values.", call. = FALSE)
  }
  if (any(rowSums(is.finite(Y)) < min_curve_points)) {
    stop(
      "Every curve must contain at least ", min_curve_points,
      " observed grid points.",
      call. = FALSE
    )
  }
  
  if (is.null(grid)) grid <- seq(0, 1, length.out = D)
  grid <- as.numeric(grid)
  if (length(grid) != D) {
    stop("length(grid) must equal ncol(Y).", call. = FALSE)
  }
  if (any(!is.finite(grid)) || is.unsorted(grid, strictly = TRUE)) {
    stop("'grid' must be a strictly increasing finite vector.", call. = FALSE)
  }
  
  subject_factor <- factor(subject, levels = unique(subject))
  subject_id <- as.integer(subject_factor)
  subject_levels <- levels(subject_factor)
  I <- length(subject_levels)
  counts <- tabulate(subject_id, nbins = I)
  
  if (I < 2) stop("At least two subjects are required.", call. = FALSE)
  
  list(
    Y = Y,
    subject_id = subject_id,
    subject_levels = subject_levels,
    visit_time = as.numeric(visit_time),
    grid = grid,
    n = n,
    D = D,
    I = I,
    counts = counts
  )
}