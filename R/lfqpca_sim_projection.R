lfqpca_sim_projection <- function(f, basis, w, tol = 1e-10) {
  f <- as.numeric(f)
  basis <- as.matrix(basis)
  sqrt_w <- sqrt(w)
  qr_b <- qr(sweep(basis, 1, sqrt_w, "*"), tol = tol, LAPACK = FALSE)
  coef <- qr.coef(qr_b, f * sqrt_w)
  coef[!is.finite(coef)] <- 0
  as.numeric(basis %*% coef)
}
