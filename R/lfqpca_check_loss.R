lfqpca_check_loss <- function(residual, tau, mean = FALSE) {
  residual <- residual[is.finite(residual)]
  if (!length(residual)) return(NA_real_)
  out <- sum(lfqpca_check_loss_values(residual, tau))
  if (mean) out / length(residual) else out
}
