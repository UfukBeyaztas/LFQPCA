lfqpcr_split_coefficients <- function(coefficients, block_sizes) {
  coefficients <- as.numeric(coefficients)
  block_sizes <- as.integer(block_sizes)
  
  if (length(block_sizes) != 4) {
    stop(
      "block_sizes must contain exactly four blocks: alpha, beta, B, and U.",
      call. = FALSE
    )
  }
  
  names(block_sizes) <- c("alpha", "beta", "B", "U")
  
  if (anyNA(block_sizes) || any(block_sizes < 0)) {
    stop(
      "block_sizes must contain non-negative finite integer sizes.",
      call. = FALSE
    )
  }
  
  if (sum(block_sizes) != length(coefficients)) {
    stop(
      "The coefficient vector length does not match the sum of block_sizes.",
      call. = FALSE
    )
  }
  
  starts <- c(1, head(cumsum(block_sizes), -1) + 1)
  names(starts) <- names(block_sizes)
  
  ends <- cumsum(block_sizes)
  names(ends) <- names(block_sizes)
  
  take <- function(block) {
    m <- unname(block_sizes[[block]])
    
    if (m == 0) {
      return(numeric(0))
    }
    
    first <- unname(starts[[block]])
    last <- unname(ends[[block]])
    
    coefficients[seq.int(first, last)]
  }
  
  list(
    a_alpha = take("alpha"),
    beta = take("beta"),
    a_B = take("B"),
    a_U = take("U")
  )
}
