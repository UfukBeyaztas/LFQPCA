lfqpca_sim_data <- function(design = paste0("D", 1:7), seed, I = 60, J = 5,
                            missing_rate = 0.05,
                            settings = lfqpca_sim_settings(),
                            scale_only = FALSE) {
  design <- match.arg(design)
  s <- settings
  D <- s$D
  K_subject <- s$K_subject
  K_visit <- s$K_visit
  set.seed(seed)
  
  visits <- if (design == "D6") {
    sample(s$unbalanced_visits, I, replace = TRUE)
  } else {
    rep.int(J, I)
  }
  subject <- rep(seq_len(I), visits)
  n <- length(subject)
  grid <- seq(0, 1, length.out = D)
  w <- lfqpca_trapezoid_weights(grid, normalize = TRUE)
  visit_time <- numeric(n)
  pos <- 1
  for (i in seq_len(I)) {
    J_i <- visits[i]
    visit_time[pos:(pos + J_i - 1)] <-
      sort(runif(J_i, min = s$time_range[1], max = s$time_range[2]))
    pos <- pos + J_i
  }
  time_std <- (visit_time - mean(visit_time)) / sd(visit_time)
  H <- lfqpca_time_basis(time_std, degree = s$mean_time_degree)
  
  time_on <- design != "D2"
  eta0 <- 1.2 + 0.8 * sin(2 * pi * grid) + 0.20 * cos(4 * pi * grid)
  eta1 <- if (time_on) {
    0.35 * cos(2 * pi * grid) + 0.20 * (grid - 0.5)
  } else {
    rep(0, D)
  }
  eta2 <- if (time_on) -0.18 * sin(pi * grid) else rep(0, D)
  eta <- cbind(eta0, eta1, eta2)
  
  raw_phi0 <- cbind(sin(2 * pi * grid),
                    cos(pi * grid) + 0.20 * sin(3 * pi * grid))
  raw_phi1 <- cbind(0.55 * cos(2 * pi * grid), -0.45 * cos(3 * pi * grid))
  pair <- lfqpca_weighted_orthonormalize(
    rbind(raw_phi0[, seq_len(K_subject), drop = FALSE],
          raw_phi1[, seq_len(K_subject), drop = FALSE]),
    rep(w, 2)
  )
  phi0 <- pair[seq_len(D), , drop = FALSE]
  phi1 <- if (time_on) {
    pair[D + seq_len(D), , drop = FALSE]
  } else {
    matrix(0, D, K_subject)
  }
  raw_psi <- cbind(phi0[, 1],
                   0.7 * sin(4 * pi * grid) + 0.3 * cos(3 * pi * grid))
  psi <- lfqpca_weighted_orthonormalize(
    raw_psi[, seq_len(K_visit), drop = FALSE], w
  )
  
  xi_sd <- c(1.25, 0.70)[seq_len(K_subject)]
  xi <- matrix(rnorm(I * K_subject), I, K_subject)
  xi <- sweep(xi, 2, xi_sd, "*")
  xi <- scale(xi, center = TRUE, scale = FALSE)
  
  zeta_sd <- c(0.75, 0.40)[seq_len(K_visit)]
  zeta <- matrix(0, n, K_visit)
  null_bases <- vector("list", I)
  for (i in seq_len(I)) {
    rows <- which(subject == i)
    N_i <- lfqpca_constraint_null_basis(H[rows, , drop = FALSE])
    null_bases[[i]] <- N_i
    if (ncol(N_i)) {
      A_i <- matrix(rnorm(ncol(N_i) * K_visit), nrow = ncol(N_i),
                    ncol = K_visit)
      A_i <- sweep(A_i, 2, zeta_sd * sqrt(length(rows) / ncol(N_i)), "*")
      zeta[rows, ] <- N_i %*% A_i
    }
  }
  if (design == "D5") {
    zeta <- lfqpca_sim_with_seed(seed + s$iid_seed_offset, {
      z <- matrix(rnorm(n * K_visit), nrow = n, ncol = K_visit)
      sweep(z, 2, zeta_sd, "*")
    })
  }
  if (design == "D3") zeta[, 1] <- s$visit_signal * zeta[, 1]
  
  location_population <- lfqpca_mean_surface(H, eta)
  location_subject_intercept_slope <-
    lfqpca_subject_contribution(subject, time_std, xi, phi0, phi1)
  phi2 <- matrix(0, D, K_subject)
  location_subject_quadratic <- matrix(0, n, D)
  if (design == "D7") {
    phi2 <- s$curvature * phi1
    quadratic_scores <- sweep(xi[subject, , drop = FALSE], 1, time_std^2, "*")
    location_subject_quadratic <- quadratic_scores %*% t(phi2)
  }
  location_subject <-
    location_subject_intercept_slope + location_subject_quadratic
  location_visit <- lfqpca_visit_contribution(zeta, psi)
  location_surface <- location_population + location_subject + location_visit
  
  q_bound <- s$bound_visits - (s$mean_time_degree + 1)
  multiplier <- lfqpca_sim_score_multiplier(s)
  xi_scale <- matrix(0, nrow = I, ncol = K_subject)
  xi_scale[, 1] <- lfqpca_sim_antithetic(I, xi_sd[1])
  zeta_scale <- matrix(0, nrow = n, ncol = K_visit)
  for (i in seq_len(I)) {
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
  
  scale_phi0 <- matrix(0, D, K_subject)
  scale_phi1 <- matrix(0, D, K_subject)
  scale_psi <- matrix(0, D, K_visit)
  scale_subject <- matrix(0, n, D)
  scale_visit <- matrix(0, n, D)
  if (design == "D3") {
    if (max(abs(time_std)) > s$max_abs_time + 1e-12) {
      if (scale_only) return(list(design = design, scale_min = -Inf))
      stop("D3 requires max |u| <= ", s$max_abs_time, ".", call. = FALSE)
    }
    scale_intercept <- lfqpca_sim_scale_intercept(
      phi0[, 1], phi1[, 1], psi[, 1],
      subject_coef = multiplier * xi_sd[1],
      visit_coef = s$visit_signal * multiplier * zeta_sd[1],
      J_design = max(visits), settings = s
    )
    scale_eta <- cbind(rep(scale_intercept, D), rep(0, D), rep(0, D))
    scale_phi0 <- phi0
    scale_phi1 <- phi1
    scale_psi <- psi
    scale_subject <- lfqpca_subject_contribution(subject, time_std, xi_scale,
                                                 scale_phi0, scale_phi1)
    scale_visit <- lfqpca_visit_contribution(zeta_scale, scale_psi)
  } else {
    scale_intercept <- if (design == "D4") {
      s$loading_scale_intercept
    } else {
      s$base_scale_intercept
    }
    scale_eta <- cbind(scale_intercept + 0.15 * cos(2 * pi * grid),
                       rep(0, D), rep(0, D))
  }
  if (design == "D4") {
    raw_scale_phi0 <- cbind(
      0.70 * cos(3 * pi * grid) + 0.20 * sin(5 * pi * grid),
      sin(4 * pi * grid) - 0.30 * cos(2 * pi * grid)
    )[, seq_len(K_subject), drop = FALSE]
    raw_scale_phi1 <- cbind(
      sin(3 * pi * grid) + 0.15 * cos(5 * pi * grid),
      -cos(4 * pi * grid) + 0.20 * sin(2 * pi * grid)
    )[, seq_len(K_subject), drop = FALSE]
    raw_scale_psi <- cbind(
      cos(3 * pi * grid) + 0.25 * sin(2 * pi * grid),
      sin(5 * pi * grid) - 0.20 * cos(3 * pi * grid)
    )[, seq_len(K_visit), drop = FALSE]
    spline_basis <- lfqpca_function_basis(grid = grid, df = s$spline_df,
                                          degree = s$spline_degree,
                                          periodic = s$periodic)
    scale_B <- lfqpca_sim_new_directions(
      raw = lfqpca_sim_project(cbind(raw_scale_phi0, raw_scale_phi1),
                               spline_basis, w),
      base_space = lfqpca_sim_project(cbind(phi0, phi1), spline_basis, w),
      norms = c(s$loading_norms$phi0[seq_len(K_subject)],
                s$loading_norms$phi1[seq_len(K_subject)]),
      w = w
    )
    scale_phi0 <- scale_B[, seq_len(K_subject), drop = FALSE]
    scale_phi1 <- scale_B[, K_subject + seq_len(K_subject), drop = FALSE]
    scale_psi <- lfqpca_sim_new_directions(
      raw = lfqpca_sim_project(raw_scale_psi, spline_basis, w),
      base_space = lfqpca_sim_project(psi, spline_basis, w),
      norms = s$loading_norms$psi[seq_len(K_visit)],
      w = w
    )
    scale_subject <- lfqpca_subject_contribution(subject, time_std, xi,
                                                 scale_phi0, scale_phi1)
    scale_visit <- lfqpca_visit_contribution(zeta, scale_psi)
  }
  scale_population <- lfqpca_mean_surface(H, scale_eta)
  scale_surface <- scale_population + scale_subject + scale_visit
  scale_min <- min(scale_surface)
  if (scale_only) {
    return(list(design = design, scale_intercept = scale_intercept,
                scale_min = scale_min,
                required_intercept = scale_intercept +
                  max(0, s$min_scale - scale_min)))
  }
  if (!is.finite(scale_min) || scale_min < s$min_scale - 1e-10) {
    stop("The scale surface falls below ", s$min_scale, ".", call. = FALSE)
  }
  
  Y_complete <- location_surface +
    scale_surface * matrix(rnorm(n * D, mean = 0, sd = 1), nrow = n, ncol = D)
  Y <- Y_complete
  if (missing_rate > 0) {
    Y[sample(seq_along(Y), floor(missing_rate * length(Y)))] <- NA
  }
  
  list(
    design = design,
    seed = seed,
    Y = Y,
    Y_complete = Y_complete,
    subject = subject,
    visit_time = visit_time,
    time_std = time_std,
    grid = grid,
    weights = w,
    visits = visits,
    truth = list(
      heterogeneity = switch(design, D3 = "score", D4 = "loading", "none"),
      eta = eta,
      phi0 = phi0,
      phi1 = phi1,
      phi2 = phi2,
      psi = psi,
      xi = xi,
      zeta = zeta,
      xi_scale = xi_scale,
      zeta_scale = zeta_scale,
      location_population = location_population,
      location_subject_intercept_slope = location_subject_intercept_slope,
      location_subject_quadratic = location_subject_quadratic,
      location_visit = location_visit,
      scale_intercept = scale_intercept,
      scale_eta = scale_eta,
      scale_phi0 = scale_phi0,
      scale_phi1 = scale_phi1,
      scale_psi = scale_psi,
      scale_population = scale_population,
      scale_subject = scale_subject,
      scale_visit = scale_visit
    )
  )
}