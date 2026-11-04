lfqpca_fill_curve_na <- function(x, grid) {
  x <- as.numeric(x)
  good <- is.finite(x)
  if (!sum(good)) return(rep(0, length(x)))
  if (sum(good) == 1) return(rep(x[good], length(x)))
  approx(grid[good], x[good], xout = grid, rule = 2)$y
}