lfqpca_rq_solve <- function(x, y, tau, which_method) {
  messages <- character(0)
  
  fit <- tryCatch(
    withCallingHandlers(
      switch(
        which_method,
        # quantreg has no rq.fit.fn(); rq(method = "fn") uses rq.fit.fnb().
        fnb = rq.fit.fnb(x = x, y = y, tau = tau),
        fn = rq.fit.fnb(x = x, y = y, tau = tau),
        br = rq.fit.br(x = x, y = y, tau = tau),
        stop("Unknown quantile-regression method: ", which_method,
             call. = FALSE)
      ),
      warning = function(w) {
        messages <<- c(messages, conditionMessage(w))
        invokeRestart("muffleWarning")
      }
    ),
    error = function(e) e
  )
  
  if (inherits(fit, "error")) {
    return(list(
      ok = FALSE,
      coefficients = NULL,
      message = conditionMessage(fit),
      warnings = messages
    ))
  }
  
  beta <- as.numeric(fit$coefficients)
  failure_warning <- grepl(
    "Error info|singular|Premature end|conditioning",
    messages,
    ignore.case = TRUE
  )
  problems <- c(
    if (length(beta) != ncol(x)) "wrong number of returned coefficients",
    if (length(beta) &&
        !all(is.finite(beta))) "non-finite returned coefficients",
    messages[failure_warning]
  )
  
  list(
    ok = !length(problems),
    coefficients = beta,
    message = paste(problems, collapse = "; "),
    warnings = messages
  )
}