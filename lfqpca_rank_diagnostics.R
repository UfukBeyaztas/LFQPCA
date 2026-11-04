lfqpca_rank_diagnostics <- function(xi, phi0, phi1, zeta, psi, weights,
                                    rank_tol = 1e-8,
                                    orthonormality_tol = 1e-6) {
  subject <- lfqpca_block_rank_info(
    scores = xi,
    loadings = rbind(phi0, phi1),
    weights = rep(weights, 2),
    rank_tol = rank_tol,
    orthonormality_tol = orthonormality_tol
  )
  visit <- lfqpca_block_rank_info(
    scores = zeta,
    loadings = psi,
    weights = weights,
    rank_tol = rank_tol,
    orthonormality_tol = orthonormality_tol
  )
  
  list(
    subject = subject,
    visit = visit,
    rank_tol = rank_tol,
    orthonormality_tol = orthonormality_tol,
    nondegenerate = isTRUE(subject$full_rank) && isTRUE(visit$full_rank)
  )
}