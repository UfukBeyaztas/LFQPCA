lfqpcr_sim_data <- function(case = paste0("Case", 1:6), seed, I_train = 60,
                            I_test = 20, J = 5,
                            tau_F = 0.5, tau_R = 0.9,
                            settings = lfqpca_sim_settings()) {
  case <- match.arg(case)
  s <- settings
  D <- s$D
  K_subject <- s$K_subject
  K_visit <- s$K_visit
  set.seed(seed)
  
  I_total <- I_train + I_test
  train_subjects <- seq_len(I_train)
  test_subjects <- I_train + seq_len(I_test)
  visits <- if (case == "Case6") {
    sample(s$unbalanced_visits, I_total, replace = TRUE)
  } else {
    rep.int(J, I_total)
  }
  subject <- rep(seq_len(I_total), visits)
  n <- length(subject)
  train <- subject %in% train_subjects
  test <- subject %in% test_subjects
  grid <- seq(0, 1, length.out = D)
  w_norm <- lfqpca_trapezoid_weights(grid)
  w <- lfqpca_trapezoid_weights(grid, normalize = FALSE)
  visit_time <- numeric(n)
  pos <- 1
  for (i in seq_len(I_total)) {
    J_i <- visits[i]
    visit_time[pos:(pos + J_i - 1)] <-
      sort(runif(J_i, min = s$time_range[1], max = s$time_range[2]))
    pos <- pos + J_i
  }
  time_center <- mean(visit_time[train])
  time_scale <- sd(visit_time[train])
  time_std <- (visit_time - time_center) / time_scale
  H <- lfqpca_time_basis(time_std, degree = s$mean_time_degree)
  
  no_time <- case == "Case4"
  eta0 <- 1.2 + 0.8 * sin(2 * pi * grid) + 0.20 * cos(4 * pi * grid)
  eta1 <- if (no_time) {
    rep(0, D)
  } else {
    0.35 * cos(2 * pi * grid) + 0.20 * (grid - 0.5)
  }
  eta2 <- if (no_time) rep(0, D) else -0.18 * sin(pi * grid)
  eta <- cbind(eta0, eta1, eta2)
  raw_pair <- cbind(
    c(sin(2 * pi * grid), 0.55 * cos(2 * pi * grid)),
    c(cos(pi * grid) + 0.20 * sin(3 * pi * grid), -0.45 * cos(3 * pi * grid))
  )[, seq_len(K_subject), drop = FALSE]
  pair <- lfqpca_weighted_orthonormalize(raw_pair, rep(w_norm, 2))
  phi0 <- pair[seq_len(D), , drop = FALSE]
  phi1 <- pair[D + seq_len(D), , drop = FALSE]
  raw_psi <- cbind(phi0[, 1],
                   0.7 * sin(4 * pi * grid) + 0.3 * phi1[, min(2, K_subject)])
  psi <- lfqpca_weighted_orthonormalize(
    raw_psi[, seq_len(K_visit), drop = FALSE], w_norm
  )
  if (no_time) phi1[, ] <- 0
  
  xi_sd <- c(1.25, 0.70)[seq_len(K_subject)]
  xi <- matrix(rnorm(I_total * K_subject), nrow = I_total, ncol = K_subject)
  xi <- sweep(xi, 2, xi_sd, "*")
  xi <- sweep(xi, 2, colMeans(xi[train_subjects, , drop = FALSE]), "-")
  
  zeta_sd <- c(0.75, 0.40)[seq_len(K_visit)]
  zeta <- matrix(0, nrow = n, ncol = K_visit)
  null_bases <- vector("list", I_total)
  for (i in seq_len(I_total)) {
    rows <- which(subject == i)
    N_i <- lfqpca_constraint_null_basis(H[rows, , drop = FALSE])
    null_bases[[i]] <- N_i
    if (ncol(N_i) > 0) {
      A_i <- matrix(rnorm(ncol(N_i) * K_visit), nrow = ncol(N_i),
                    ncol = K_visit)
      A_i <- sweep(A_i, 2, zeta_sd * sqrt(length(rows) / ncol(N_i)), "*")
      zeta[rows, ] <- N_i %*% A_i
    }
  }
  if (case == "Case5") {
    zeta <- lfqpca_sim_with_seed(seed + s$iid_seed_offset, {
      z <- matrix(rnorm(n * K_visit), nrow = n, ncol = K_visit)
      sweep(z, 2, zeta_sd, "*")
    })
  }
  heterogeneous <- case == "Case3"
  if (heterogeneous) zeta[, 1] <- s$visit_signal * zeta[, 1]
  
  location_population <- lfqpca_mean_surface(H, eta)
  location_subject <- lfqpca_subject_contribution(subject, time_std, xi,
                                                  phi0, phi1)
  location_visit <- lfqpca_visit_contribution(zeta, psi)
  location_surface <- location_population + location_subject + location_visit
  error_exp <- matrix(rexp(n * D, rate = 1), nrow = n, ncol = D)
  
  xi_tau <- xi
  zeta_tau <- zeta
  B_response <- NULL
  U_response <- NULL
  if (!heterogeneous) {
    B <- location_subject
    U <- location_visit
    sigma_time <- if (no_time) rep(0.55, n) else 0.55 + 0.20 * time_std^2
    sigma_surface <- outer(sigma_time, 0.25 + 0.45 * grid, "*")
    k_tau <- qexp(tau_F, rate = 1) - 1
    population <- location_population + k_tau * sigma_surface
    functional_quantile <- location_surface + k_tau * sigma_surface
    Y <- location_surface + sigma_surface * (error_exp - 1)
  } else {
    error_normal <- matrix(rnorm(n * D, mean = 0, sd = 1), nrow = n, ncol = D)
    if (max(abs(time_std)) > s$max_abs_time + 1e-12) {
      stop("Case 3 requires max |u| <= ", s$max_abs_time, ".", call. = FALSE)
    }
    multiplier <- lfqpca_sim_score_multiplier(s)
    q_bound <- s$bound_visits - (s$mean_time_degree + 1)
    scale_intercept <- lfqpca_sim_scale_intercept(
      phi0[, 1], phi1[, 1], psi[, 1],
      subject_coef = multiplier * xi_sd[1],
      visit_coef = s$visit_signal * multiplier * zeta_sd[1],
      J_design = max(visits), settings = s
    )
    scale_eta <- cbind(rep(scale_intercept, D), rep(0, D), rep(0, D))
    xi_scale <- matrix(0, nrow = I_total, ncol = K_subject)
    xi_scale[train_subjects, 1] <-
      lfqpca_sim_antithetic(length(train_subjects), xi_sd[1])
    xi_scale[test_subjects, 1] <-
      lfqpca_sim_antithetic(length(test_subjects), xi_sd[1])
    zeta_scale <- matrix(0, nrow = n, ncol = K_visit)
    for (i in seq_len(I_total)) {
      rows <- which(subject == i)
      N_i <- null_bases[[i]]
      if (ncol(N_i) > 0) {
        if (ncol(N_i) > q_bound) {
          stop("A subject has more than ", s$bound_visits, " visits.",
               call. = FALSE)
        }
        A_i <- matrix(0, nrow = ncol(N_i), ncol = K_visit)
        g <- rnorm(ncol(N_i))
        A_i[, 1] <- s$visit_signal * g / sqrt(sum(g^2)) * sqrt(length(rows)) *
          zeta_sd[1]
        zeta_scale[rows, ] <- N_i %*% A_i
      }
    }
    xi_scale <- multiplier * xi_scale
    zeta_scale <- multiplier * zeta_scale
    scale_population <- lfqpca_mean_surface(H, scale_eta)
    scale_surface <- scale_population +
      lfqpca_subject_contribution(subject, time_std, xi_scale, phi0, phi1) +
      lfqpca_visit_contribution(zeta_scale, psi)
    if (!is.finite(min(scale_surface)) ||
        min(scale_surface) < s$min_scale - 1e-10) {
      stop("The scale surface falls below ", s$min_scale, ".", call. = FALSE)
    }
    q_ref <- qnorm(s$reference_tau, mean = 0, sd = 1)
    q_shift <- qnorm(tau_F, mean = 0, sd = 1) - q_ref
    eta <- eta + q_shift * scale_eta
    xi_tau <- xi + q_shift * xi_scale
    zeta_tau <- zeta + q_shift * zeta_scale
    population <- location_population + q_shift * scale_population
    B <- lfqpca_subject_contribution(subject, time_std, xi_tau, phi0, phi1)
    U <- lfqpca_visit_contribution(zeta_tau, psi)
    functional_quantile <- population + B + U
    Y <- location_surface + scale_surface * (error_normal - q_ref)
  }
  
  c_subject <- rnorm(I_total)
  C <- cbind(c_subject = c_subject[subject], c_visit = rnorm(n))
  beta <- c(c_subject = 0.50, c_visit = -0.35)
  alpha <- if (no_time) {
    rep(0.50, n)
  } else {
    0.50 + 0.40 * time_std - 0.15 * time_std^2
  }
  gamma_single <- rep(NA_real_, D)
  if (case == "Case1") {
    gamma_single <- lfqpca_sim_l2_normalize(
      0.90 * phi0[, 1] - 0.55 * phi1[, 1] + 0.40 * phi0[, min(2, K_subject)] +
        0.35 * psi[, min(2, K_visit)],
      w, 0.80
    )
    gamma_B <- lfqpca_sim_projection(gamma_single, cbind(phi0, phi1), w)
    gamma_U <- lfqpca_sim_projection(gamma_single, psi, w)
    B_effect <- rowSums(sweep(B, 2, gamma_single * w, "*"))
    U_effect <- rowSums(sweep(U, 2, gamma_single * w, "*"))
  } else {
    if (heterogeneous && (K_subject < 2 || K_visit < 2)) {
      stop("Case 3 requires K_subject >= 2 and K_visit >= 2.", call. = FALSE)
    }
    gamma_B <- lfqpca_sim_l2_normalize(
      0.90 * phi0[, 1] - 0.60 * phi1[, 1] + 0.50 * phi0[, min(2, K_subject)] +
        0.30 * phi1[, min(2, K_subject)],
      w, 0.80
    )
    gamma_U <- lfqpca_sim_l2_normalize(
      0.70 * psi[, 1] - 0.50 * psi[, min(2, K_visit)], w, 0.60
    )
    if (heterogeneous) {
      # The response depends on the components at the level tau* = 0.9.
      q_star <- qnorm(s$response_tau_star) - qnorm(s$reference_tau)
      B_response <- lfqpca_subject_contribution(
        subject, time_std, xi + q_star * xi_scale, phi0, phi1
      )
      U_response <- lfqpca_visit_contribution(zeta + q_star * zeta_scale, psi)
      B_effect <- rowSums(sweep(B_response, 2, gamma_B * w, "*"))
      U_effect <- rowSums(sweep(U_response, 2, gamma_U * w, "*"))
    } else {
      B_effect <- rowSums(sweep(B, 2, gamma_B * w, "*"))
      U_effect <- rowSums(sweep(U, 2, gamma_U * w, "*"))
    }
  }
  response_quantile <- alpha + as.numeric(C %*% beta) + B_effect + U_effect
  sigma_R <- if (no_time) rep(0.55, n) else 0.55 + 0.15 * abs(time_std)
  V <- response_quantile + sigma_R * (rexp(n, rate = 1) - qexp(tau_R, rate = 1))
  
  list(
    case = case,
    seed = seed,
    Y = Y,
    V = V,
    C = C,
    subject = subject,
    visit_time = visit_time,
    grid = grid,
    visits = visits,
    train = train,
    test = test,
    tau_F = tau_F,
    tau_R = tau_R,
    truth = list(
      time_center = time_center,
      time_scale = time_scale,
      time_std = time_std,
      eta = eta,
      phi0 = phi0,
      phi1 = phi1,
      psi = psi,
      xi = xi_tau,
      zeta = zeta_tau,
      population = population,
      B = B,
      U = U,
      functional_quantile = functional_quantile,
      B_response = if (is.null(B_response)) B else B_response,
      U_response = if (is.null(U_response)) U else U_response,
      alpha = alpha,
      beta = beta,
      gamma_single = gamma_single,
      gamma_B = gamma_B,
      gamma_U = gamma_U,
      B_effect = B_effect,
      U_effect = U_effect,
      response_quantile = response_quantile
    )
  )
}
