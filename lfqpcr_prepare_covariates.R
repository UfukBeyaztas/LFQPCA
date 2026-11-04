lfqpcr_prepare_covariates <- function(covariates, n) {
  if (is.null(covariates)) {
    return(matrix(0, nrow = n, ncol = 0))
  }
  
  C <- as.matrix(covariates)
  storage.mode(C) <- "double"
  
  if (nrow(C) != n) {
    stop(
      "nrow(covariates) must equal the number of subject-visit curves.",
      call. = FALSE
    )
  }
  
  if (is.null(colnames(C))) {
    colnames(C) <- paste0("C", seq_len(ncol(C)))
  }
  
  C
}