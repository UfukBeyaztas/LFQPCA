lfqpca_rotate_block <- function(scores, loadings, weights,
                                center_scores = TRUE,
                                eps = 1e-10) {
  scores <- as.matrix(scores)
  loadings <- as.matrix(loadings)
  K <- ncol(loadings)
  
  if (!K) {
    return(list(
      scores = scores,
      loadings = loadings,
      score_mean = numeric(0),
      variances = numeric(0)
    ))
  }
  if (ncol(scores) != K) {
    stop("scores and loadings must have the same number of columns.",
         call. = FALSE)
  }
  
  score_mean <- if (center_scores) colMeans(scores) else rep(0, K)
  S0 <- if (center_scores) sweep(scores, 2, score_mean, "-") else scores
  
  G <- lfqpca_weighted_gram(loadings, weights)
  G_half <- lfqpca_sym_power(G, 0.5, eps = eps)
  G_inv_half <- lfqpca_sym_power(G, -0.5, eps = eps)
  
  L1 <- loadings %*% G_inv_half
  S1 <- S0 %*% G_half
  
  if (nrow(S1) > 1) {
    score_cov <- cov(S1)
  } else {
    score_cov <- crossprod(S1)
  }
  score_cov <- (score_cov + t(score_cov)) / 2
  es <- eigen(score_cov, symmetric = TRUE)
  ord <- order(es$values, decreasing = TRUE)
  V <- es$vectors[, ord, drop = FALSE]
  
  L2 <- L1 %*% V
  S2 <- S1 %*% V
  variances <- pmax(es$values[ord], 0)
  
  peak_index <- apply(abs(L2), 2, which.max)
  signs <- sign(L2[cbind(peak_index, seq_len(K))])
  signs[signs == 0] <- 1
  L2 <- sweep(L2, 2, signs, "*")
  S2 <- sweep(S2, 2, signs, "*")
  
  list(
    scores = S2,
    loadings = L2,
    score_mean = score_mean,
    variances = variances
  )
}