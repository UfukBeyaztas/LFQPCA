lfqpcr_subject_weights <- function(subject_id, eligible) {
  subject_id <- as.integer(subject_id)
  eligible <- as.logical(eligible)
  
  if (length(subject_id) != length(eligible)) {
    stop("subject_id and eligible must have the same length.", call. = FALSE)
  }
  
  if (anyNA(subject_id) || any(subject_id < 1)) {
    stop("subject_id must contain positive integer labels.", call. = FALSE)
  }
  
  active_subjects <- sort(unique(subject_id[eligible]))
  
  if (!length(active_subjects)) {
    stop("No subject has an eligible scalar-response visit.", call. = FALSE)
  }
  
  I_R <- length(active_subjects)
  weights <- numeric(length(subject_id))
  J_R <- setNames(integer(I_R), as.character(active_subjects))
  
  for (ii in seq_along(active_subjects)) {
    i <- active_subjects[ii]
    rows_i <- which(eligible & subject_id == i)
    J_i_R <- length(rows_i)
    J_R[ii] <- J_i_R
    weights[rows_i] <- (1 / I_R) * (1 / J_i_R)
  }
  
  list(
    weights = weights,
    active_subjects = active_subjects,
    I_R = I_R,
    J_R = J_R
  )
}
