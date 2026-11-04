lfqpca_block_rank_info <- function(scores, loadings, weights,
                                   rank_tol = 1e-8,
                                   orthonormality_tol = 1e-6) {
  scores <- as.matrix(scores)
  loadings <- as.matrix(loadings)
  K <- ncol(loadings)
  
  if (!K) {
    return(list(
      K = 0,
      loading_rank = 0,
      score_rank = 0,
      orthonormality_error = 0,
      gram_eigenvalues = numeric(0),
      score_variances = numeric(0),
      full_rank = TRUE
    ))
  }
  
  G <- lfqpca_weighted_gram(loadings, weights)
  G <- (G + t(G)) / 2
  ev_G <- if (all(is.finite(G))) {
    eigen(G, symmetric = TRUE, only.values = TRUE)$values
  } else {
    rep(NA_real_, K)
  }
  
  S_c <- sweep(scores, 2, colMeans(scores), "-")
  C <- if (nrow(S_c) > 1) crossprod(S_c) / (nrow(S_c) - 1) else crossprod(S_c)
  C <- (C + t(C)) / 2
  ev_S <- if (all(is.finite(C))) {
    eigen(C, symmetric = TRUE, only.values = TRUE)$values
  } else {
    rep(NA_real_, K)
  }
  
  finite <- all(is.finite(ev_G)) && all(is.finite(ev_S))
  loading_rank <- if (finite) sum(ev_G > rank_tol * max(ev_G, 0)) else 0
  score_rank <- if (finite) sum(ev_S > rank_tol * max(ev_S, 0)) else 0
  orthonormality_error <- if (all(is.finite(G))) max(abs(G - diag(K))) else Inf
  
  list(
    K = as.integer(K),
    loading_rank = as.integer(loading_rank),
    score_rank = as.integer(score_rank),
    orthonormality_error = orthonormality_error,
    gram_eigenvalues = ev_G,
    score_variances = ev_S,
    full_rank = finite && loading_rank == K && score_rank == K &&
      orthonormality_error <= orthonormality_tol
  )
}
