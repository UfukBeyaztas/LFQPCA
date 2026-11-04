lfqpca_check_loss_values <- function(residual, tau) {
  residual * (tau - (residual < 0))
}