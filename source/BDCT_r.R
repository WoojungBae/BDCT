# ==============================================================================
# Define functions -------------------------------------------------------
# ==============================================================================

# ------------------------------------------------------------------------------
# Sign convention and design-stage priors
# ------------------------------------------------------------------------------
#   h(t | z, x) = h0(t) exp(-z * eta_z(t) - x * eta_x).
#
# Therefore:
#   eta_z(t) > 0 means treatment is beneficial.
#   eta_z(t) = 0 is the exact-null reference value.
#   eta_z(t) < 0 means treatment is harmful.
#   eta_x > 0 means larger x lowers the hazard.
#   The interval-specific treatment hazard ratio is exp(-eta_z(t)).
#
# Current decision criterion:
#   delta_star = 0.00 month.
#
# Common favorable temporal profile shapes:
#   constant   : (1.00, 1.00, 1.00, 1.00, 1.00)
#   increasing : (0.50, 0.75, 1.00, 1.25, 1.50)
#   waning     : (1.50, 1.25, 1.00, 0.75, 0.50)
#
# All three shapes have interval mean equal to 1, so kappa_s is the
# interval-average eta value of the profile center kappa_s * c_s.
#
# Design-stage simulation priors:
#   prior_id = 1: exact-null point-mass design prior
#                 omega_D = 0, H_D = 0
#                 eta_j = 0, j = 1, ..., 5
#
#   prior_id = 2: continuous no-benefit design prior
#                 omega_D = 0, H_D = 0
#                 eta_j ~ TN(0, 0.30^2; -0.25, 0),
#                 independently across intervals
#
#   prior_id = 3: continuous positive-benefit design prior
#                 omega_D = 1, H_D = 1
#                 eta_j ~ TN(kappa_s * c_{s,j}, 0.10^2; 0, Inf),
#                 independently across intervals
#
#   prior_id = 4: more dispersed continuous positive-benefit design prior
#                 omega_D = 1, H_D = 1
#                 eta_j ~ TN(kappa_s * c_{s,j}, 0.30^2; 0, Inf),
#                 independently across intervals
#
# For Priors 3 and 4, kappa_s is calibrated separately for each temporal
# profile so that
#
#   E_eta[Delta_RMST(24)] = 2.00 months.
#
# Benefit-profile input:
#   "constant", "increasing", or "waning".
# Legacy aliases "delay" and "wane" are accepted and normalized.
#
# Priors 1 and 2 do not depend on the temporal benefit profile and should each
# be simulated only once per design. Priors 3 and 4 are simulated separately
# for the constant, increasing, and waning profiles.
#
# Intermediate mixture weights are NOT simulated separately. Predictive
# operating characteristics for 0 < omega_D < 1 are obtained by
# post-processing component-specific operating characteristics.
#
# Final-look rule:
#   final_cutoff = max(tau, K2 event time).
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# Basic utilities
# ------------------------------------------------------------------------------

collapse_num <- function(x, digits = 4) {
  paste(round(as.numeric(x), digits), collapse = ", ")
}

safe_max <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0) return(NA_real_)
  max(x)
}

safe_min <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0) return(NA_real_)
  min(x)
}

safe_mean <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) == 0) return(NA_real_)
  mean(x)
}

safe_rmse <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) == 0) return(NA_real_)
  sqrt(mean(x))
}

draw_truncnorm <- function(n = 1L, location = 0, sd = 1, lower = -Inf, upper = Inf) {
  if (!is.finite(sd) || sd <= 0) stop("`sd` must be positive.")
  if (lower >= upper) stop("`lower` must be smaller than `upper`.")
  
  p_lower <- pnorm(lower, mean = location, sd = sd)
  p_upper <- pnorm(upper, mean = location, sd = sd)
  
  if (!is.finite(p_lower) || !is.finite(p_upper) || p_upper <= p_lower) stop("Invalid truncated-normal probability interval.")
  
  u <- runif(n, min = p_lower, max = p_upper)
  qnorm(u, mean = location, sd = sd)
}

# ------------------------------------------------------------------------------
# Favorable temporal profile shapes and calibrated scaling constants
# ------------------------------------------------------------------------------

get_eta_z_profile_shapes <- function() {
  list(
    constant = c(1.00, 1.00, 1.00, 1.00, 1.00),
    increasing = c(0.50, 0.75, 1.00, 1.25, 1.50),
    waning = c(1.50, 1.25, 1.00, 0.75, 0.50)
  )
}

get_eta_z_kappas <- function(prior_id) {
  prior_id <- as.integer(prior_id)
  
  if (
    length(prior_id) != 1L ||
    is.na(prior_id) ||
    !(prior_id %in% c(3L, 4L))
  ) {
    stop(
      "`prior_id` must be 3 or 4 when requesting positive-benefit profile kappa values."
    )
  }
  
  if (prior_id == 3L) {
    # Prior 3:
    # eta_j ~ TN(kappa_s * c_{s,j}, 0.10^2; 0, Inf), independently,
    # with E_eta[Delta_RMST(24)] = 2.00.
    return(c(
      constant = 0.2881833829,
      increasing = 0.3985155474,
      waning = 0.2265051174
    ))
  }
  
  # Prior 4:
  # eta_j ~ TN(kappa_s * c_{s,j}, 0.30^2; 0, Inf), independently,
  # with E_eta[Delta_RMST(24)] = 2.00.
  c(
    constant = 0.1469563494,
    increasing = 0.1990671295,
    waning = 0.1149792825
  )
}

get_eta_z_centers <- function(prior_id = 3L) {
  shapes <- get_eta_z_profile_shapes()
  kappas <- get_eta_z_kappas(prior_id)
  
  list(
    constant = kappas[["constant"]] * shapes$constant,
    increasing = kappas[["increasing"]] * shapes$increasing,
    waning = kappas[["waning"]] * shapes$waning
  )
}

normalize_benefit_profile <- function(profile) {
  if (length(profile) != 1L || is.na(profile)) {
    stop("`profile` must be a single non-missing value.")
  }
  
  key <- tolower(gsub("[_ ]", "-", trimws(profile)))
  
  map <- c(
    "constant" = "constant",
    "delay" = "increasing",
    "increasing" = "increasing",
    "increasing-effect" = "increasing",
    "wane" = "waning",
    "waning" = "waning"
  )
  
  if (!(key %in% names(map))) {
    stop("`profile` must be one of constant, increasing, or waning.")
  }
  
  unname(map[[key]])
}

get_eta_z_calibration_report <- function(delta_star = 0.00) {
  shapes <- get_eta_z_profile_shapes()
  profiles <- c("constant", "increasing", "waning")
  
  shape_text <- vapply(
    profiles,
    function(p) collapse_num(shapes[[p]], digits = 2),
    character(1)
  )
  
  centers_p3 <- get_eta_z_centers(prior_id = 3L)
  centers_p4 <- get_eta_z_centers(prior_id = 4L)
  kappas_p3 <- get_eta_z_kappas(prior_id = 3L)
  kappas_p4 <- get_eta_z_kappas(prior_id = 4L)
  
  location_rmst_p3 <- c(
    constant = 2.0213294007,
    increasing = 2.0046265329,
    waning = 2.0102601115
  )
  
  location_rmst_p4 <- c(
    constant = 1.0460907852,
    increasing = 1.0171928555,
    waning = 1.0425245934
  )
  
  make_rows <- function(prior_id, sigma_star, centers, kappas, location_delta) {
    do.call(
      rbind,
      lapply(seq_along(profiles), function(ii) {
        p <- profiles[ii]
        eta_location <- centers[[p]]
        
        data.frame(
          prior_id = prior_id,
          profile = p,
          base_shape = shape_text[ii],
          sigma_star = sigma_star,
          kappa_value = unname(kappas[[p]]),
          eta_z_center = collapse_num(eta_location),
          HR_center = collapse_num(exp(-eta_location)),
          calibration_target = "E_eta[Delta_RMST(24)] = 2",
          achieved_target_delta_rmst = 2.00,
          center_delta_rmst = unname(location_delta[[p]]),
          delta_star = delta_star,
          H_D_center = as.integer(location_delta[[p]] > delta_star),
          stringsAsFactors = FALSE
        )
      })
    )
  }
  
  rbind(
    make_rows(
      prior_id = 3L,
      sigma_star = 0.10,
      centers = centers_p3,
      kappas = kappas_p3,
      location_delta = location_rmst_p3
    ),
    make_rows(
      prior_id = 4L,
      sigma_star = 0.30,
      centers = centers_p4,
      kappas = kappas_p4,
      location_delta = location_rmst_p4
    )
  )
}

# ------------------------------------------------------------------------------
# True RMST and true median calculations for a given eta_z profile
# ------------------------------------------------------------------------------

true_delta_rmst_from_eta <- function(eta_z, breaks, lambda0, eta_x, tau, x_std) {
  J <- length(lambda0)
  
  rmst_arm <- function(z_val) {
    total <- rep(0, length(x_std))
    cumhaz <- rep(0, length(x_std))
    
    for (j in seq_len(J)) {
      left <- breaks[j]
      right <- breaks[j + 1]
      
      if (tau <= left) break
      
      width_j <- min(tau, right) - left
      
      if (width_j > 0) {
        haz_j <- lambda0[j] * exp(-z_val * eta_z[j] - x_std * eta_x)
        surv_left <- exp(-cumhaz)
        
        total <- total +
          surv_left * (1 - exp(-haz_j * width_j)) / haz_j
        
        cumhaz <- cumhaz + haz_j * width_j
      }
      
      if (tau <= right) break
    }
    
    if (tau > breaks[J + 1]) {
      haz_last <- lambda0[J] * exp(-z_val * eta_z[J] - x_std * eta_x)
      surv_left <- exp(-cumhaz)
      width_last <- tau - breaks[J + 1]
      
      total <- total +
        surv_left * (1 - exp(-haz_last * width_last)) / haz_last
    }
    
    mean(total)
  }
  
  rmst_arm(1) - rmst_arm(0)
}

calc_true_marginal_surv_at_t <- function(t_val, z_val, eta_z, breaks, lambda0, eta_x, x_std) {
  J <- length(lambda0)
  cumhaz <- rep(0, length(x_std))
  
  for (j in seq_len(J)) {
    left <- breaks[j]
    right <- breaks[j + 1]
    
    if (t_val <= left) break
    
    width_j <- min(t_val, right) - left
    
    if (width_j > 0) {
      haz_j <- lambda0[j] * exp(-z_val * eta_z[j] - x_std * eta_x)
      cumhaz <- cumhaz + haz_j * width_j
    }
    
    if (t_val <= right) break
  }
  
  if (t_val > breaks[J + 1]) {
    haz_last <- lambda0[J] * exp(-z_val * eta_z[J] - x_std * eta_x)
    width_last <- t_val - breaks[J + 1]
    cumhaz <- cumhaz + haz_last * width_last
  }
  
  mean(exp(-cumhaz))
}

true_delta_med_from_eta <- function(eta_z, breaks, lambda0, eta_x, x_std, q = 0.5) {
  target <- 1 - q
  max_t <- 500
  
  find_med <- function(z_val) {
    if (calc_true_marginal_surv_at_t(max_t, z_val, eta_z, breaks, lambda0, eta_x, x_std) > target) return(NA_real_)
    
    obj_fun <- function(t_val) {
      calc_true_marginal_surv_at_t(t_val, z_val, eta_z, breaks, lambda0, eta_x, x_std) - target
    }
    
    uniroot(obj_fun, lower = 0, upper = max_t, extendInt = "no")$root
  }
  
  med1 <- find_med(1)
  med0 <- find_med(0)
  
  if (is.na(med1) || is.na(med0)) return(NA_real_)
  
  med1 - med0
}

# ------------------------------------------------------------------------------
# Design-stage priors
# ------------------------------------------------------------------------------

get_design_prior_spec <- function(prior_id) {
  prior_id <- as.integer(prior_id)
  
  if (
    length(prior_id) != 1L ||
    is.na(prior_id) ||
    !(prior_id %in% 1:4)
  ) {
    stop("`prior_id` must be one of 1, 2, 3, or 4.")
  }
  
  specs <- data.frame(
    prior_id = 1:4,
    prior_label = c(
      "exact_null_point_mass",
      "no_benefit_continuous",
      "benefit_positive_truncnorm_sd010",
      "benefit_positive_truncnorm_sd030"
    ),
    prior_component = c(
      "no_benefit",
      "no_benefit",
      "benefit",
      "benefit"
    ),
    prior_family = c(
      "point_mass",
      "continuous",
      "continuous",
      "continuous"
    ),
    omega_D = c(0, 0, 1, 1),
    H_D = c(0L, 0L, 1L, 1L),
    stringsAsFactors = FALSE
  )
  
  specs[specs$prior_id == prior_id, , drop = FALSE]
}

get_design_omega <- function(prior_id) {
  get_design_prior_spec(prior_id)$omega_D[[1]]
}

get_design_prior_label <- function(prior_id) {
  get_design_prior_spec(prior_id)$prior_label[[1]]
}

# ------------------------------------------------------------------------------
# Draw eta_z from the final design-stage priors
# ------------------------------------------------------------------------------

draw_eta_z_from_prior <- function(
    prior_id,
    profile = NULL
) {
  prior_id <- as.integer(prior_id)
  spec <- get_design_prior_spec(prior_id)
  
  if (prior_id == 1L) {
    # --------------------------------------------------------------------------
    # Prior 1: exact-null point mass.
    # eta_j = 0 for all five intervals.
    # --------------------------------------------------------------------------
    profile_out <- "none"
    eta_draw <- rep(0, 5L)
    
  } else if (prior_id == 2L) {
    # --------------------------------------------------------------------------
    # Prior 2: continuous no-benefit prior.
    # eta_j ~ TN(0, 0.30^2; -0.25, 0), independently across intervals.
    # --------------------------------------------------------------------------
    profile_out <- "none"
    
    eta_draw <- vapply(
      seq_len(5L),
      function(j) {
        draw_truncnorm(
          n = 1L,
          location = 0.00,
          sd = 0.30,
          lower = -0.25,
          upper = 0.00
        )
      },
      numeric(1)
    )
    
  } else {
    # --------------------------------------------------------------------------
    # Priors 3 and 4: continuous positive-benefit priors.
    #
    # Prior 3:
    #   eta_j ~ TN(kappa_s * c_{s,j}, 0.10^2; 0, Inf)
    #
    # Prior 4:
    #   eta_j ~ TN(kappa_s * c_{s,j}, 0.30^2; 0, Inf)
    #
    # Draws are independent across intervals. The profile-specific kappas are
    # calibrated so that E_eta[Delta_RMST(24)] = 2.00.
    # --------------------------------------------------------------------------
    if (is.null(profile)) {
      stop("A benefit profile is required when `prior_id` is 3 or 4.")
    }
    
    profile_out <- normalize_benefit_profile(profile)
    centers <- get_eta_z_centers(prior_id = prior_id)
    eta_location <- centers[[profile_out]]
    sigma_star <- if (prior_id == 3L) 0.10 else 0.30
    
    eta_draw <- vapply(
      seq_along(eta_location),
      function(j) {
        draw_truncnorm(
          n = 1L,
          location = eta_location[j],
          sd = sigma_star,
          lower = 0,
          upper = Inf
        )
      },
      numeric(1)
    )
  }
  
  list(
    eta_z = eta_draw,
    prior_id = prior_id,
    scenario_id = prior_id,
    prior_label = spec$prior_label[[1]],
    design_prior_type = spec$prior_label[[1]],
    prior_component = spec$prior_component[[1]],
    prior_family = spec$prior_family[[1]],
    profile = profile_out,
    omega_D = spec$omega_D[[1]],
    H_D = spec$H_D[[1]],
    component_id = spec$H_D[[1]]
  )
}

# ------------------------------------------------------------------------------
# Full event-driven data generation function
# ------------------------------------------------------------------------------

generate_data <- function(
    prior_id = 1,
    profile = "constant",
    N = NULL,
    K1 = NULL,
    K2 = 360,
    target_censor = 0.40,
    breaks = c(0, 4.8, 9.6, 14.4, 19.2, 24),
    lambda0 = c(0.0734, 0.0661, 0.0587, 0.0514, 0.0440),
    eta_x = 0.2,
    pi_z = 0.5,
    tau = 24,
    delta_star = 0.00,
    dropout = FALSE,
    dropout_rate = 0,
    truth_n = 50000,
    return_truth = TRUE,
    verbose = TRUE
) {
  prior_id <- as.integer(prior_id)
  stopifnot(prior_id %in% 1:4)
  stopifnot(length(breaks) - 1 == length(lambda0))
  stopifnot(target_censor > 0, target_censor < 1)
  stopifnot(K2 > 0)
  
  if (is.null(N)) {
    N <- ceiling(K2 / (1 - target_censor))
  }
  
  if (is.null(K1)) {
    K1 <- floor(0.6 * K2)
  }
  
  stopifnot(K1 > 0, K1 < K2, K2 <= N)
  
  J <- length(lambda0)
  
  eta_draw_obj <- draw_eta_z_from_prior(
    prior_id = prior_id,
    profile = profile
  )
  
  eta_z <- eta_draw_obj$eta_z
  prior_label <- eta_draw_obj$prior_label
  design_prior_type <- eta_draw_obj$design_prior_type
  prior_component <- eta_draw_obj$prior_component
  prior_family <- eta_draw_obj$prior_family
  component_id <- eta_draw_obj$component_id
  scenario_id <- eta_draw_obj$scenario_id
  omega_D <- eta_draw_obj$omega_D
  H_D <- eta_draw_obj$H_D
  profile <- eta_draw_obj$profile
  
  stopifnot(length(eta_z) == J)
  
  if (verbose) {
    message(
      "Generating prior ", prior_id,
      " [", design_prior_type, "]",
      ", omega_D = ", omega_D,
      ", H_D = ", H_D,
      ", profile = ", profile,
      ", eta_z = (", paste(round(eta_z, 4), collapse = ", "), ")."
    )
  }
  
  sim_one_pwexp <- function(xi, zi) {
    u <- runif(1)
    target <- -log(u)
    cumhaz <- 0
    
    for (j in seq_len(J)) {
      haz_j <- lambda0[j] * exp(-zi * eta_z[j] - xi * eta_x)
      width_j <- breaks[j + 1] - breaks[j]
      inc_j <- haz_j * width_j
      
      if (cumhaz + inc_j >= target) {
        return(breaks[j] + (target - cumhaz) / haz_j)
      }
      
      cumhaz <- cumhaz + inc_j
    }
    
    haz_last <- lambda0[J] * exp(-zi * eta_z[J] - xi * eta_x)
    breaks[J + 1] + (target - cumhaz) / haz_last
  }
  
  x <- rnorm(N, mean = 0, sd = 1)
  z <- rbinom(N, size = 1, prob = pi_z)
  
  T_event <- vapply(
    seq_len(N),
    function(i) sim_one_pwexp(xi = x[i], zi = z[i]),
    numeric(1)
  )
  
  if (dropout) {
    if (dropout_rate <= 0) {
      stop("If dropout = TRUE, `dropout_rate` must be positive.")
    }
    C_loss <- rexp(N, rate = dropout_rate)
  } else {
    C_loss <- rep(Inf, N)
  }
  
  observable_event <- T_event <= C_loss
  
  if (sum(observable_event) < K2) {
    stop(
      "Fewer than K2 observable events were generated among N = ", N,
      ". Increase N, reduce K2, or reduce dropout."
    )
  }
  
  event_times_observed <- sort(T_event[observable_event])
  
  y_K1 <- event_times_observed[K1]
  y_K2_event <- event_times_observed[K2]
  y_final <- max(tau, y_K2_event)
  
  make_look_data <- function(y_cut, K_target, look_name) {
    time <- pmin(T_event, C_loss, y_cut)
    status <- as.integer(T_event <= C_loss & T_event <= y_cut)
    
    data.frame(
      id = seq_len(N),
      x = x,
      z = z,
      time = time,
      status = status,
      T_event = T_event,
      C_loss = C_loss,
      C_admin = y_cut,
      look = look_name,
      look_cutoff = y_cut,
      K_target = K_target,
      dropout_censor = as.integer(C_loss < T_event & C_loss <= y_cut),
      admin_censor = as.integer(T_event > y_cut & C_loss > y_cut),
      event_after_tau = as.integer(status == 1 & T_event > tau),
      stringsAsFactors = FALSE
    )
  }
  
  O1 <- make_look_data(
    y_cut = y_K1,
    K_target = K1,
    look_name = "K1"
  )
  
  O2 <- make_look_data(
    y_cut = y_final,
    K_target = K2,
    look_name = "K2"
  )
  
  full_data <- data.frame(
    id = seq_len(N),
    x = x,
    z = z,
    T_event = T_event,
    C_loss = C_loss,
    observable_event = observable_event,
    stringsAsFactors = FALSE
  )
  
  summarize_look <- function(O) {
    make_row <- function(dd, arm_label) {
      data.frame(
        arm = arm_label,
        N = nrow(dd),
        events = sum(dd$status == 1),
        censored = sum(dd$status == 0),
        censor_rate = mean(dd$status == 0),
        dropout_censor = sum(dd$dropout_censor == 1),
        dropout_censor_rate = mean(dd$dropout_censor == 1),
        admin_censor = sum(dd$admin_censor == 1),
        admin_censor_rate = mean(dd$admin_censor == 1),
        events_after_tau = sum(dd$event_after_tau == 1),
        event_after_tau_rate = mean(dd$event_after_tau == 1),
        look_cutoff = unique(dd$look_cutoff),
        K_target = unique(dd$K_target),
        stringsAsFactors = FALSE
      )
    }
    
    rbind(
      make_row(O, "overall"),
      make_row(O[O$z == 0, ], "Z=0"),
      make_row(O[O$z == 1, ], "Z=1")
    )
  }
  
  truth <- NULL
  
  if (return_truth) {
    x_truth <- rnorm(truth_n, mean = 0, sd = 1)
    
    delta_rmst_tau <- true_delta_rmst_from_eta(
      eta_z = eta_z,
      breaks = breaks,
      lambda0 = lambda0,
      eta_x = eta_x,
      tau = tau,
      x_std = x_truth
    )
    
    delta_med_rep <- true_delta_med_from_eta(
      eta_z = eta_z,
      breaks = breaks,
      lambda0 = lambda0,
      eta_x = eta_x,
      x_std = x_truth
    )
    
    H_D_rmst <- as.integer(delta_rmst_tau > delta_star)
    
    if (!identical(H_D_rmst, as.integer(component_id))) {
      warning(
        "RMST truth classification did not match the design-stage component. ",
        "This may indicate Monte Carlo error in the truth calculation."
      )
    }
    
    truth <- list(
      delta_rmst_tau = delta_rmst_tau,
      delta_med_rep = delta_med_rep,
      H_D = as.integer(component_id),
      H_D_rmst_check = H_D_rmst
    )
  }
  
  out <- list(
    prior_id = prior_id,
    scenario_id = scenario_id,
    prior_label = prior_label,
    design_prior_type = design_prior_type,
    prior_family = prior_family,
    omega_D = omega_D,
    profile = profile,
    prior_component = prior_component,
    component_id = component_id,
    H_D = H_D,
    N = N,
    K1 = K1,
    K2 = K2,
    target_censor = target_censor,
    implied_final_censor_rate = 1 - K2 / N,
    y_K1 = y_K1,
    y_K2_event = y_K2_event,
    y_final = y_final,
    tau = tau,
    delta_star = delta_star,
    x = x,
    z = z,
    X = cbind(z = z, x = x),
    T_event = T_event,
    C_loss = C_loss,
    full_data = full_data,
    data = O2,
    O1 = O1,
    O2 = O2,
    OL = O2,
    time = O2$time,
    status = O2$status,
    breaks = breaks,
    lambda0 = lambda0,
    eta_x = eta_x,
    eta_z = eta_z,
    eta_z_summary = data.frame(
      interval = seq_along(eta_z),
      eta_z = eta_z,
      HR = exp(-eta_z),
      stringsAsFactors = FALSE
    ),
    max_eta_z = max(eta_z),
    min_eta_z = min(eta_z),
    min_HR = min(exp(-eta_z)),
    max_HR = max(exp(-eta_z)),
    pi_z = pi_z,
    dropout = dropout,
    dropout_rate = dropout_rate,
    censoring_O1 = summarize_look(O1),
    censoring_O2 = summarize_look(O2),
    censoring_OL = summarize_look(O2),
    event_after_tau_at_final = sum(O2$event_after_tau == 1),
    event_after_tau_rate_at_final = mean(O2$event_after_tau == 1)
  )
  
  if (return_truth) {
    out$truth <- truth
  }
  
  out
}

# ------------------------------------------------------------------------------
# Read saved truth values
# ------------------------------------------------------------------------------

read_section_table <- function(file, section_title) {
  lines <- readLines(file, warn = FALSE)
  start <- grep(section_title, lines, fixed = TRUE)
  
  if (length(start) == 0) stop("Section not found: ", section_title, " in ", file)
  
  start <- start[1] + 1
  end_candidates <- which(seq_along(lines) > start & trimws(lines) == "")
  end <- end_candidates[end_candidates > start][1]
  
  if (is.na(end)) end <- length(lines) + 1
  
  tab_lines <- lines[start:(end - 1)]
  tab_lines <- tab_lines[nchar(trimws(tab_lines)) > 0]
  
  read.table(text = paste(tab_lines, collapse = "\n"), header = TRUE, sep = "\t", stringsAsFactors = FALSE)
}

read_truth_lookup <- function(truth_dir = ".", prior_ids = 1:4, require_all = TRUE) {
  out <- list()
  
  for (pid in prior_ids) {
    file <- file.path(truth_dir, paste0("Truth_Check_Scenario", pid, ".txt"))
    
    if (!file.exists(file)) {
      msg <- paste0("Truth file not found: ", file)
      if (require_all) stop(msg) else {
        warning(msg)
        next
      }
    }
    
    tab <- read_section_table(file = file, section_title = "1. Truth summaries by profile")
    
    tab <- tab[tab$estimand == "RMST_contrast", ]
    
    if ("scenario_id" %in% names(tab)) {
      names(tab)[names(tab) == "scenario_id"] <- "prior_id"
    }
    
    tab$prior_id <- as.integer(tab$prior_id)
    
    out[[as.character(pid)]] <- tab[, c("prior_id", "profile", "mean")]
  }
  
  if (length(out) == 0) stop("No truth lookup tables were read from truth_dir: ", truth_dir)
  
  truth_lookup <- do.call(rbind, out)
  names(truth_lookup)[names(truth_lookup) == "mean"] <- "true_delta_rmst_lookup"
  row.names(truth_lookup) <- NULL
  
  truth_lookup
}

read_truth_lookup_from_true_dir <- function(true_dir, prior_ids = 1:4, require_all = TRUE) {
  old_wd <- getwd()
  on.exit(setwd(old_wd), add = TRUE)
  
  setwd(true_dir)
  
  read_truth_lookup(truth_dir = ".", prior_ids = prior_ids, require_all = require_all)
}

# ------------------------------------------------------------------------------
# Look-specific analysis datasets
# ------------------------------------------------------------------------------

make_analysis_datasets <- function(dat) {
  O1 <- dat$O1
  O2 <- dat$O2
  
  attr(O1, "X") <- as.matrix(dat$X)
  attr(O2, "X") <- as.matrix(dat$X)
  
  attr(O1, "n_enrolled") <- nrow(O1)
  attr(O2, "n_enrolled") <- nrow(O2)
  
  attr(O1, "events") <- sum(O1$status == 1)
  attr(O2, "events") <- sum(O2$status == 1)
  
  attr(O1, "censor_rate") <- mean(O1$status == 0)
  attr(O2, "censor_rate") <- mean(O2$status == 0)
  
  attr(O1, "calendar_cutoff") <- unique(O1$look_cutoff)
  attr(O2, "calendar_cutoff") <- unique(O2$look_cutoff)
  
  list(O1 = O1, O2 = O2)
}

# ------------------------------------------------------------------------------
# Gamma-process PH model on prespecified grid
# ------------------------------------------------------------------------------

make_gp_grid_data <- function(look_data, X, B0_fun, a0 = 2, grid_width = 1, tau = 24) {
  X <- as.matrix(X)
  time <- as.numeric(look_data$time)
  status <- as.integer(look_data$status)
  
  N <- length(time)
  P <- ncol(X)
  
  grid_end <- ceiling(max(tau, max(time, na.rm = TRUE)) / grid_width) * grid_width
  t_grid <- seq(0, grid_end, by = grid_width)
  
  J <- length(t_grid) - 1
  widths <- diff(t_grid)
  
  dN <- integer(J)
  
  for (j in seq_len(J)) {
    dN[j] <- sum(status == 1 & time > t_grid[j] & time <= t_grid[j + 1])
  }
  
  dB0 <- diff(B0_fun(t_grid))
  dB0 <- pmax(dB0, 1e-8)
  
  list(
    N = N,
    P = P,
    J = J,
    X = X,
    y = time,
    D = status,
    t_grid = t_grid,
    widths = widths,
    dN = as.integer(dN),
    mu_eta = rep(0, P),
    sigma_eta = rep(2.5, P),
    a0 = a0,
    dB0 = dB0
  )
}
gp_grid_stan_code <- '
data {
  int<lower=1> N;
  int<lower=1> P;
  int<lower=1> J;
  matrix[N, P] X;
  vector<lower=0>[N] y;
  array[N] int<lower=0, upper=1> D;
  vector<lower=0>[J + 1] t_grid;
  vector<lower=0>[J] widths;
  array[J] int<lower=0> dN;
  vector[P] mu_eta;
  vector<lower=0>[P] sigma_eta;
  real<lower=0> a0;
  vector<lower=0>[J] dB0;
}

parameters {
  vector[P] eta_coef;
  vector<lower=1e-12>[J] dLambda0;
}

model {
  vector[N] xb = -(X * eta_coef);
  vector[N] exb = exp(xb);
  vector[J] risk_exposure;

  eta_coef ~ normal(mu_eta, sigma_eta);
  dLambda0 ~ gamma(a0 * dB0, a0);

  for (j in 1:J) {
    real left = t_grid[j];
    real width = widths[j];
    risk_exposure[j] = 0;

    for (i in 1:N) {
      real dt = fmin(fmax(y[i] - left, 0), width);
      risk_exposure[j] += dt * exb[i] / width;
    }
  }

  target += dot_product(to_vector(dN), log(dLambda0 ./ widths));
  target += dot_product(to_vector(D), xb);
  target += -dot_product(dLambda0, risk_exposure);
}
'
compile_gp_model <- function() {
  if (!requireNamespace("cmdstanr", quietly = TRUE)) stop("Package 'cmdstanr' is required.")
  
  gp_file <- cmdstanr::write_stan_file(gp_grid_stan_code)
  cmdstanr::cmdstan_model(gp_file)
}

# ------------------------------------------------------------------------------
# Posterior RMST and PH coefficient calculation
# ------------------------------------------------------------------------------

extract_gp_draws <- function(fit) {
  eta_draws <- fit$draws(variables = "eta_coef", format = "matrix")
  dL_draws <- fit$draws(variables = "dLambda0", format = "matrix")
  
  eta_cols <- grep("^eta_coef\\[", colnames(eta_draws), value = TRUE)
  dL_cols <- grep("^dLambda0\\[", colnames(dL_draws), value = TRUE)
  
  list(
    eta = as.matrix(eta_draws[, eta_cols, drop = FALSE]),
    dLambda0 = as.matrix(dL_draws[, dL_cols, drop = FALSE])
  )
}

summarize_eta_PH_trt_gp <- function(fit, X, suffix = "1") {
  draws <- extract_gp_draws(fit)
  eta_mat <- draws$eta
  
  X <- as.matrix(X)
  
  if (is.null(colnames(X))) stop("X must have column names so that the treatment column can be identified.")
  
  if (!("z" %in% colnames(X))) stop("Treatment column `z` was not found in X.")
  
  z_col <- which(colnames(X) == "z")[1]
  
  if (z_col > ncol(eta_mat)) stop("Treatment column index is larger than the number of posterior eta columns.")
  
  eta_trt_draw <- as.numeric(eta_mat[, z_col])
  
  out <- data.frame(
    post_eta_PH_trt_mean = mean(eta_trt_draw, na.rm = TRUE),
    post_eta_PH_trt_q025 = unname(quantile(eta_trt_draw, 0.025, na.rm = TRUE)),
    post_eta_PH_trt_q500 = unname(quantile(eta_trt_draw, 0.500, na.rm = TRUE)),
    post_eta_PH_trt_q975 = unname(quantile(eta_trt_draw, 0.975, na.rm = TRUE)),
    stringsAsFactors = FALSE
  )
  
  names(out) <- paste0(names(out), "_", suffix)
  
  out
}

make_X_standardized <- function(X) {
  X <- as.matrix(X)
  
  if (!("z" %in% colnames(X))) {
    colnames(X)[1] <- "z"
  }
  
  X1 <- X
  X0 <- X
  
  X1[, "z"] <- 1
  X0[, "z"] <- 0
  
  storage.mode(X1) <- "double"
  storage.mode(X0) <- "double"
  
  list(X1 = X1, X0 = X0)
}

rmst_gp_grid_marginal <- function(tau, Xfixed, eta, dLambda0, t_grid) {
  Xfixed <- as.matrix(Xfixed)
  eta <- as.numeric(eta)
  dLambda0 <- as.numeric(dLambda0)
  
  mu <- exp(-drop(Xfixed %*% eta))
  
  J <- length(dLambda0)
  lefts <- t_grid[-length(t_grid)]
  rights <- t_grid[-1]
  widths <- rights - lefts
  
  lens <- pmax(0, pmin(tau, rights) - lefts)
  
  lambda_j <- dLambda0 / widths
  H_prev <- c(0, cumsum(dLambda0))[seq_len(J)]
  
  S_left <- exp(-outer(mu, H_prev))
  denom <- outer(mu, lambda_j)
  len_mat <- matrix(lens, nrow = nrow(Xfixed), ncol = J, byrow = TRUE)
  
  factor <- ifelse(denom > 0, (1 - exp(-denom * len_mat)) / denom, len_mat)
  
  mean(rowSums(S_left * factor))
}

summarize_draws <- function(x) {
  c(
    mean = mean(x, na.rm = TRUE),
    q025 = unname(quantile(x, 0.025, na.rm = TRUE)),
    q500 = unname(quantile(x, 0.500, na.rm = TRUE)),
    q975 = unname(quantile(x, 0.975, na.rm = TRUE)),
    width = unname(quantile(x, 0.975, na.rm = TRUE) - quantile(x, 0.025, na.rm = TRUE))
  )
}

posterior_rmst_summary_gp <- function(fit, gp_data, X, tau = 24, delta_star = 0.00) {
  draws <- extract_gp_draws(fit)
  
  eta_mat <- draws$eta
  dL_mat <- draws$dLambda0
  
  Xstd <- make_X_standardized(X)
  
  rmst_diff <- numeric(nrow(eta_mat))
  
  for (m in seq_len(nrow(eta_mat))) {
    rmst1 <- rmst_gp_grid_marginal(tau = tau, Xfixed = Xstd$X1, eta = eta_mat[m,], dLambda0 = dL_mat[m,], t_grid = gp_data$t_grid)
    
    rmst0 <- rmst_gp_grid_marginal(tau = tau, Xfixed = Xstd$X0, eta = eta_mat[m,], dLambda0 = dL_mat[m,], t_grid = gp_data$t_grid)
    
    rmst_diff[m] <- rmst1 - rmst0
  }
  
  ss <- summarize_draws(rmst_diff)
  
  list(
    draws = rmst_diff,
    Pi = unname(mean(rmst_diff > delta_star, na.rm = TRUE)),
    post_mean = unname(ss["mean"]),
    post_q025 = unname(ss["q025"]),
    post_q500 = unname(ss["q500"]),
    post_q975 = unname(ss["q975"]),
    post_width = unname(ss["width"])
  )
}

extract_mcmc_diagnostics <- function(fit, variables = c("eta_coef", "dLambda0")) {
  ss <- fit$summary(variables = variables)
  sampler_diag <- fit$sampler_diagnostics(format = "df")
  
  data.frame(
    max_Rhat = safe_max(ss$rhat),
    min_ESS = safe_min(ss$ess_bulk),
    n_divergent = if ("divergent__" %in% names(sampler_diag)) {
      sum(sampler_diag$divergent__, na.rm = TRUE)
    } else {
      NA_integer_
    }
  )
}

# ------------------------------------------------------------------------------
# Fit one look
# ------------------------------------------------------------------------------

fit_gp_look <- function(
    mod_gp,
    look_data,
    X,
    B0_fun,
    a0_gp = 2,
    grid_width_gp = 1,
    tau = 24,
    cmdstan_seed = 1,
    chains = 4,
    parallel_chains = 1,
    iter_warmup = 5000,
    iter_sampling = 2500,
    refresh = 200,
    thin = 1,
    adapt_delta = 0.95,
    max_treedepth = 12
) {
  gp_data <- make_gp_grid_data(look_data = look_data, X = X, B0_fun = B0_fun, a0 = a0_gp, grid_width = grid_width_gp, tau = tau)
  
  init_fun_gp <- function() {
    list(eta_coef = rep(0, gp_data$P), dLambda0 = pmax(gp_data$dB0, 1e-4))
  }
  
  fit <- mod_gp$sample(
    data = gp_data,
    seed = as.integer(cmdstan_seed),
    chains = chains,
    parallel_chains = parallel_chains,
    iter_warmup = iter_warmup,
    iter_sampling = iter_sampling,
    refresh = refresh,
    init = init_fun_gp,
    thin = thin,
    adapt_delta = adapt_delta,
    max_treedepth = max_treedepth
  )
  
  list(fit = fit, gp_data = gp_data)
}

# ------------------------------------------------------------------------------
# MCMC function for BDCT
# ------------------------------------------------------------------------------

BDCT_MCMC <- function(
    dat,
    rep_id = NA_integer_,
    mod_gp,
    tau = NULL,
    delta_star = NULL,
    B0_fun,
    B0_median_pfs = NA_real_,
    a0_gp = 2,
    grid_width_gp = 1,
    chains = 4,
    parallel_chains = 1,
    iter_warmup = 5000,
    iter_sampling = 2500,
    refresh = 200,
    thin = 1,
    adapt_delta = 0.95,
    max_treedepth = 12,
    cmdstan_seed = 1,
    truth_n = 50000
) {
  cmdstan_seed <- as.integer(cmdstan_seed)
  
  if (is.null(tau)) {
    tau <- dat$tau
  }
  
  if (is.null(delta_star)) {
    delta_star <- dat$delta_star
  }
  
  if (is.null(dat$breaks)) stop("dat$breaks is required.")
  if (is.null(dat$lambda0)) stop("dat$lambda0 is required.")
  if (is.null(dat$eta_x)) stop("dat$eta_x is required.")
  if (is.null(dat$eta_z)) stop("dat$eta_z is required.")
  if (is.null(dat$prior_id)) stop("dat$prior_id is required.")
  if (is.null(dat$prior_label)) stop("dat$prior_label is required.")
  if (is.null(dat$profile)) stop("dat$profile is required.")
  if (is.null(dat$N)) stop("dat$N is required.")
  if (is.null(dat$K1)) stop("dat$K1 is required.")
  if (is.null(dat$K2)) stop("dat$K2 is required.")
  
  true_delta_rep <- if (!is.null(dat$truth) && !is.null(dat$truth$delta_rmst_tau)) {
    dat$truth$delta_rmst_tau
  } else {
    true_delta_rmst_from_eta(
      eta_z = dat$eta_z,
      breaks = dat$breaks,
      lambda0 = dat$lambda0,
      eta_x = dat$eta_x,
      tau = tau,
      x_std = rnorm(truth_n, mean = 0, sd = 1)
    )
  }
  
  true_delta_med_rep <- if (!is.null(dat$truth) && !is.null(dat$truth$delta_med_rep)) {
    dat$truth$delta_med_rep
  } else {
    true_delta_med_from_eta(
      eta_z = dat$eta_z,
      breaks = dat$breaks,
      lambda0 = dat$lambda0,
      eta_x = dat$eta_x,
      x_std = rnorm(truth_n, mean = 0, sd = 1)
    )
  }
  
  looks <- make_analysis_datasets(dat)
  
  O1 <- looks$O1
  O2 <- looks$O2
  
  X1 <- attr(O1, "X")
  X2 <- attr(O2, "X")
  
  fit1_obj <- fit_gp_look(
    mod_gp = mod_gp,
    look_data = O1,
    X = X1,
    B0_fun = B0_fun,
    a0_gp = a0_gp,
    grid_width_gp = grid_width_gp,
    tau = tau,
    cmdstan_seed = cmdstan_seed,
    chains = chains,
    parallel_chains = parallel_chains,
    iter_warmup = iter_warmup,
    iter_sampling = iter_sampling,
    refresh = refresh,
    thin = thin,
    adapt_delta = adapt_delta,
    max_treedepth = max_treedepth
  )
  
  fit2_obj <- fit_gp_look(
    mod_gp = mod_gp,
    look_data = O2,
    X = X2,
    B0_fun = B0_fun,
    a0_gp = a0_gp,
    grid_width_gp = grid_width_gp,
    tau = tau,
    cmdstan_seed = cmdstan_seed,
    chains = chains,
    parallel_chains = parallel_chains,
    iter_warmup = iter_warmup,
    iter_sampling = iter_sampling,
    refresh = refresh,
    thin = thin,
    adapt_delta = adapt_delta,
    max_treedepth = max_treedepth
  )
  
  constants <- list(
    rep_id = rep_id,
    prior_id = dat$prior_id,
    scenario_id = if (!is.null(dat$scenario_id)) dat$scenario_id else dat$prior_id,
    prior_label = dat$prior_label,
    design_prior_type = if (!is.null(dat$design_prior_type)) dat$design_prior_type else dat$prior_label,
    prior_family = if (!is.null(dat$prior_family)) dat$prior_family else NA_character_,
    omega_D = if (!is.null(dat$omega_D)) dat$omega_D else get_design_omega(dat$prior_id),
    profile = dat$profile,
    component = dat$prior_component,
    component_id = dat$component_id,
    H_D = if (!is.null(dat$truth) && !is.null(dat$truth$H_D)) {
      as.integer(dat$truth$H_D)
    } else if (!is.null(dat$H_D)) {
      as.integer(dat$H_D)
    } else {
      as.integer(dat$component_id)
    },
    true_delta_rmst_rep = true_delta_rep,
    true_delta_med_rep = true_delta_med_rep,
    N_planned = dat$N,
    K1 = dat$K1,
    K2 = dat$K2,
    n_interim = attr(O1, "n_enrolled"),
    n_final = attr(O2, "n_enrolled"),
    events_interim = attr(O1, "events"),
    events_final = attr(O2, "events"),
    censor_rate_interim = attr(O1, "censor_rate"),
    censor_rate_final = attr(O2, "censor_rate"),
    calendar_cutoff_interim = attr(O1, "calendar_cutoff"),
    calendar_cutoff_final = attr(O2, "calendar_cutoff"),
    tau = tau,
    delta_star = delta_star,
    a0_gp = a0_gp,
    grid_width_gp = grid_width_gp,
    B0_median_pfs = B0_median_pfs,
    chains = chains,
    parallel_chains = parallel_chains,
    iter_warmup = iter_warmup,
    iter_sampling = iter_sampling,
    thin = thin,
    adapt_delta = adapt_delta,
    max_treedepth = max_treedepth,
    n_saved_draws = chains * floor(iter_sampling / thin),
    cmdstan_seed = cmdstan_seed
  )
  
  dat_info <- list(
    scenario_id = if (!is.null(dat$scenario_id)) dat$scenario_id else dat$prior_id,
    design_prior_type = if (!is.null(dat$design_prior_type)) dat$design_prior_type else dat$prior_label,
    prior_family = if (!is.null(dat$prior_family)) dat$prior_family else NA_character_,
    omega_D = if (!is.null(dat$omega_D)) dat$omega_D else get_design_omega(dat$prior_id),
    breaks = dat$breaks,
    lambda0 = dat$lambda0,
    eta_x = dat$eta_x,
    eta_z = dat$eta_z,
    eta_z_summary = dat$eta_z_summary,
    truth = dat$truth
  )
  
  analysis_data <- list(O1 = O1, O2 = O2, X1 = X1, X2 = X2)
  
  MCMCposteriors <- list(fit1 = fit1_obj$fit, fit2 = fit2_obj$fit, gp_data1 = fit1_obj$gp_data, gp_data2 = fit2_obj$gp_data)
  
  MCMCresult <- list(constants = constants, dat_info = dat_info, analysis_data = analysis_data, MCMCposteriors = MCMCposteriors)
  
  class(MCMCresult) <- c("BDCT_MCMC", class(MCMCresult))
  
  MCMCresult
}

# ------------------------------------------------------------------------------
# POST function for BDCT (Integrated with Median Survival Contrast)
# ------------------------------------------------------------------------------

BDCT_POST <- function(object, truth_lookup = NULL, tau = NULL, delta_star = NULL) {
  if (is.null(object$constants)) stop("object$constants is required.")
  if (is.null(object$analysis_data)) stop("object$analysis_data is required.")
  if (is.null(object$MCMCposteriors)) stop("object$MCMCposteriors is required.")
  
  constants <- object$constants
  
  if (is.null(tau)) {
    tau <- constants$tau
  }
  
  if (is.null(delta_star)) {
    delta_star <- constants$delta_star
  }
  
  O1 <- object$analysis_data$O1
  O2 <- object$analysis_data$O2
  X1 <- object$analysis_data$X1
  X2 <- object$analysis_data$X2
  
  fit1 <- object$MCMCposteriors$fit1
  fit2 <- object$MCMCposteriors$fit2
  gp_data1 <- object$MCMCposteriors$gp_data1
  gp_data2 <- object$MCMCposteriors$gp_data2
  
  post1 <- posterior_rmst_summary_gp(fit = fit1, gp_data = gp_data1, X = X1, tau = tau, delta_star = delta_star)
  
  post2 <- posterior_rmst_summary_gp(fit = fit2, gp_data = gp_data2, X = X2, tau = tau, delta_star = delta_star)
  
  post_med1 <- posterior_median_surv_summary_gp(fit = fit1, gp_data = gp_data1, X = X1, q = 0.5, delta_star = delta_star)
  
  post_med2 <- posterior_median_surv_summary_gp(fit = fit2, gp_data = gp_data2, X = X2, q = 0.5, delta_star = delta_star)
  
  eta_PH_trt_1 <- summarize_eta_PH_trt_gp(fit = fit1, X = X1, suffix = "1")
  
  eta_PH_trt_2 <- summarize_eta_PH_trt_gp(fit = fit2, X = X2, suffix = "2")
  
  diag1 <- extract_mcmc_diagnostics(fit1)
  diag2 <- extract_mcmc_diagnostics(fit2)
  
  true_delta_lookup <- NA_real_
  
  if (!is.null(truth_lookup)) {
    truth_row <- truth_lookup[truth_lookup$prior_id == constants$prior_id & truth_lookup$profile == constants$profile,]
    
    true_delta_lookup <- if (nrow(truth_row) == 1) {
      truth_row$true_delta_rmst_lookup
    } else {
      NA_real_
    }
  }
  
  true_delta_rep <- constants$true_delta_rmst_rep
  true_delta_med_rep <- constants$true_delta_med_rep
  H_D_current <- as.integer(constants$H_D)
  
  POSTresult <- data.frame(
    rep = constants$rep_id,
    prior_id = constants$prior_id,
    scenario_id = constants$scenario_id,
    prior_label = constants$prior_label,
    design_prior_type = constants$design_prior_type,
    prior_family = constants$prior_family,
    omega_D = constants$omega_D,
    profile = constants$profile,
    component = constants$component,
    component_id = constants$component_id,
    H_D = H_D_current,
    true_delta_rmst_lookup = true_delta_lookup,
    true_delta_rmst_rep = true_delta_rep,
    true_delta_med_rep = true_delta_med_rep,
    N_planned = constants$N_planned,
    K1 = constants$K1,
    K2 = constants$K2,
    n_interim = constants$n_interim,
    n_final = constants$n_final,
    events_interim = constants$events_interim,
    events_final = constants$events_final,
    censor_rate_interim = constants$censor_rate_interim,
    censor_rate_final = constants$censor_rate_final,
    calendar_cutoff_interim = constants$calendar_cutoff_interim,
    calendar_cutoff_final = constants$calendar_cutoff_final,
    tau = tau,
    delta_star = delta_star,
    a0_gp = constants$a0_gp,
    grid_width_gp = constants$grid_width_gp,
    B0_median_pfs = constants$B0_median_pfs,
    chains = constants$chains,
    parallel_chains = constants$parallel_chains,
    iter_warmup = constants$iter_warmup,
    iter_sampling = constants$iter_sampling,
    thin = constants$thin,
    adapt_delta = constants$adapt_delta,
    max_treedepth = constants$max_treedepth,
    n_saved_draws = constants$n_saved_draws,
    cmdstan_seed = constants$cmdstan_seed,
    
    Pi1 = post1$Pi,
    Pi2_cf = post2$Pi,
    post_mean_1 = post1$post_mean,
    post_q025_1 = post1$post_q025,
    post_q500_1 = post1$post_q500,
    post_q975_1 = post1$post_q975,
    post_width_1 = post1$post_width,
    post_mean_2 = post2$post_mean,
    post_q025_2 = post2$post_q025,
    post_q500_2 = post2$post_q500,
    post_q975_2 = post2$post_q975,
    post_width_2 = post2$post_width,
    
    Pi_beneficial_med_1 = post_med1$Pi_beneficial,
    post_mean_med_1     = post_med1$post_mean,
    post_q025_med_1     = post_med1$post_q025,
    post_q500_med_1     = post_med1$post_q500,
    post_q975_med_1     = post_med1$post_q975,
    post_width_med_1    = post_med1$post_width,
    na_rate_med_1       = post_med1$na_rate,
    
    Pi_beneficial_med_2 = post_med2$Pi_beneficial,
    post_mean_med_2     = post_med2$post_mean,
    post_q025_med_2     = post_med2$post_q025,
    post_q500_med_2     = post_med2$post_q500,
    post_q975_med_2     = post_med2$post_q975,
    post_width_med_2    = post_med2$post_width,
    na_rate_med_2       = post_med2$na_rate,
    
    cover_lookup_1 = if (is.na(true_delta_lookup)) {
      NA_integer_
    } else {
      as.integer(post1$post_q025 <= true_delta_lookup & true_delta_lookup <= post1$post_q975)
    },
    bias_lookup_1 = if (is.na(true_delta_lookup)) {
      NA_real_
    } else {
      post1$post_mean - true_delta_lookup
    },
    sqerr_lookup_1 = if (is.na(true_delta_lookup)) {
      NA_real_
    } else {
      (post1$post_mean - true_delta_lookup)^2
    },
    
    cover_lookup_2 = if (is.na(true_delta_lookup)) {
      NA_integer_
    } else {
      as.integer(post2$post_q025 <= true_delta_lookup & true_delta_lookup <= post2$post_q975)
    },
    bias_lookup_2 = if (is.na(true_delta_lookup)) {
      NA_real_
    } else {
      post2$post_mean - true_delta_lookup
    },
    sqerr_lookup_2 = if (is.na(true_delta_lookup)) {
      NA_real_
    } else {
      (post2$post_mean - true_delta_lookup)^2
    },
    
    cover_rep_1 = as.integer(post1$post_q025 <= true_delta_rep & true_delta_rep <= post1$post_q975),
    bias_rep_1 = post1$post_mean - true_delta_rep,
    sqerr_rep_1 = (post1$post_mean - true_delta_rep)^2,
    
    cover_rep_2 = as.integer(post2$post_q025 <= true_delta_rep & true_delta_rep <= post2$post_q975),
    bias_rep_2 = post2$post_mean - true_delta_rep,
    sqerr_rep_2 = (post2$post_mean - true_delta_rep)^2,
    
    cover_med_rep_1 = if (is.finite(post_med1$post_q025) && is.finite(post_med1$post_q975) && is.finite(true_delta_med_rep)) {
      as.integer(post_med1$post_q025 <= true_delta_med_rep && true_delta_med_rep <= post_med1$post_q975)
    } else {
      NA_integer_
    },
    bias_med_rep_1 = if (is.finite(post_med1$post_mean) && is.finite(true_delta_med_rep)) {
      post_med1$post_mean - true_delta_med_rep
    } else {
      NA_real_
    },
    sqerr_med_rep_1 = if (is.finite(post_med1$post_mean) && is.finite(true_delta_med_rep)) {
      (post_med1$post_mean - true_delta_med_rep)^2
    } else {
      NA_real_
    },
    
    cover_med_rep_2 = if (is.finite(post_med2$post_q025) && is.finite(post_med2$post_q975) && is.finite(true_delta_med_rep)) {
      as.integer(post_med2$post_q025 <= true_delta_med_rep && true_delta_med_rep <= post_med2$post_q975)
    } else {
      NA_integer_
    },
    bias_med_rep_2 = if (is.finite(post_med2$post_mean) && is.finite(true_delta_med_rep)) {
      post_med2$post_mean - true_delta_med_rep
    } else {
      NA_real_
    },
    sqerr_med_rep_2 = if (is.finite(post_med2$post_mean) && is.finite(true_delta_med_rep)) {
      (post_med2$post_mean - true_delta_med_rep)^2
    } else {
      NA_real_
    },
    max_Rhat_1 = diag1$max_Rhat,
    min_ESS_1 = diag1$min_ESS,
    n_divergent_1 = diag1$n_divergent,
    max_Rhat_2 = diag2$max_Rhat,
    min_ESS_2 = diag2$min_ESS,
    n_divergent_2 = diag2$n_divergent,
    stringsAsFactors = FALSE
  )
  
  POSTresult <- cbind(POSTresult, eta_PH_trt_1, eta_PH_trt_2)
  
  POSTresult
}

# ------------------------------------------------------------------------------
# GP-PH aliases used by the simulation driver
# ------------------------------------------------------------------------------

BDCT_MCMC_ph <- BDCT_MCMC
BDCT_POST_ph <- BDCT_POST

# ------------------------------------------------------------------------------
# Optional wrapper for backward compatibility
# ------------------------------------------------------------------------------

BDCT <- function(
    dat,
    rep_id = NA_integer_,
    mod_gp,
    truth_lookup = NULL,
    tau = NULL,
    delta_star = NULL,
    B0_fun,
    B0_median_pfs = NA_real_,
    a0_gp = 2,
    grid_width_gp = 1,
    chains = 4,
    parallel_chains = 1,
    iter_warmup = 5000,
    iter_sampling = 2500,
    refresh = 200,
    thin = 1,
    adapt_delta = 0.95,
    max_treedepth = 12,
    cmdstan_seed = 1,
    truth_n = 50000
) {
  mcmc_object <- BDCT_MCMC(
    dat = dat,
    rep_id = rep_id,
    mod_gp = mod_gp,
    tau = tau,
    delta_star = delta_star,
    B0_fun = B0_fun,
    B0_median_pfs = B0_median_pfs,
    a0_gp = a0_gp,
    grid_width_gp = grid_width_gp,
    chains = chains,
    parallel_chains = parallel_chains,
    iter_warmup = iter_warmup,
    iter_sampling = iter_sampling,
    refresh = refresh,
    thin = thin,
    adapt_delta = adapt_delta,
    max_treedepth = max_treedepth,
    cmdstan_seed = cmdstan_seed,
    truth_n = truth_n
  )
  
  post_object <- BDCT_POST(object = mcmc_object, truth_lookup = truth_lookup, tau = tau, delta_star = delta_star)
  
  post_object
}

# ------------------------------------------------------------------------------
# Threshold-grid post-processing
# ------------------------------------------------------------------------------

apply_threshold_grid <- function(res_prob, threshold_grid) {
  out <- do.call(
    rbind,
    lapply(seq_len(nrow(threshold_grid)), function(ii) {
      nu1_E <- threshold_grid$nu1_E[ii]
      nu2_E <- threshold_grid$nu2_E[ii]
      
      tmp <- res_prob
      
      tmp$nu1_E <- nu1_E
      tmp$nu2_E <- nu2_E
      
      # 1. RMST-based decisions
      tmp$E1_success <- as.integer(tmp$Pi1 >= nu1_E)
      tmp$Final_success <- as.integer(tmp$Pi2_cf >= nu2_E)
      tmp$Terminal_final_success <- as.integer(tmp$E1_success == 0 & tmp$Final_success == 1)
      tmp$S_R <- as.integer(tmp$E1_success == 1 | tmp$Terminal_final_success == 1)
      
      # 2. Median-based decisions (delta_star = 0 is implicit in Pi_beneficial_med)
      tmp$E1_success_med <- as.integer(tmp$Pi_beneficial_med_1 >= nu1_E)
      tmp$Final_success_med <- as.integer(tmp$Pi_beneficial_med_2 >= nu2_E)
      tmp$Terminal_final_success_med <- as.integer(tmp$E1_success_med == 0 & tmp$Final_success_med == 1)
      tmp$S_R_med <- as.integer(tmp$E1_success_med == 1 | tmp$Terminal_final_success_med == 1)
      
      tmp
    })
  )
  
  row.names(out) <- NULL
  out
}

# ------------------------------------------------------------------------------
# Operating-characteristic summaries
# ------------------------------------------------------------------------------

mcse_mean <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) <= 1) return(NA_real_)
  sd(as.numeric(x)) / sqrt(length(x))
}

summarize_oc_one_group <- function(dd) {
  dd <- dd[!is.na(dd$Pi1) & !is.na(dd$Pi2_cf), , drop = FALSE]
  
  if (nrow(dd) == 0) return(data.frame(M = 0))
  
  M <- nrow(dd)
  H_values <- unique(dd$H_D[!is.na(dd$H_D)])
  
  if (length(H_values) != 1L) {
    stop(
      "Each simulated design prior must contain a single H_D value. ",
      "Intermediate mixtures are now obtained by post-processing rather than direct simulation."
    )
  }
  
  H_group <- as.integer(H_values[[1]])
  design_prior_type <- if ("design_prior_type" %in% names(dd)) {
    unique(dd$design_prior_type)
  } else if ("prior_label" %in% names(dd)) {
    unique(dd$prior_label)
  } else {
    NA_character_
  }
  
  # --------------------------------------------------------------------------
  # RMST decision operating characteristics
  # --------------------------------------------------------------------------
  
  PES <- mean(dd$E1_success == 1, na.rm = TRUE)
  PFS <- mean(dd$Terminal_final_success == 1, na.rm = TRUE)
  success_prob <- mean(dd$S_R == 1, na.rm = TRUE)
  
  K1 <- unique(dd$K1)
  K2 <- unique(dd$K2)
  
  ENE <- if (length(K1) == 1L && length(K2) == 1L) {
    K1 * PES + K2 * (1 - PES)
  } else {
    NA_real_
  }
  
  BP <- if (H_group == 1L) success_prob else NA_real_
  FPR <- if (H_group == 0L) success_prob else NA_real_
  
  BP_mcse <- if (H_group == 1L) mcse_mean(dd$S_R == 1) else NA_real_
  FPR_mcse <- if (H_group == 0L) mcse_mean(dd$S_R == 1) else NA_real_
  
  pointwise_type1 <- if (
    H_group == 0L &&
    length(design_prior_type) == 1L &&
    identical(design_prior_type, "exact_null_point_mass")
  ) {
    FPR
  } else {
    NA_real_
  }
  
  # --------------------------------------------------------------------------
  # Median decision operating characteristics
  # --------------------------------------------------------------------------
  
  PES_med <- mean(dd$E1_success_med == 1, na.rm = TRUE)
  PFS_med <- mean(dd$Terminal_final_success_med == 1, na.rm = TRUE)
  success_prob_med <- mean(dd$S_R_med == 1, na.rm = TRUE)
  
  ENE_med <- if (length(K1) == 1L && length(K2) == 1L) {
    K1 * PES_med + K2 * (1 - PES_med)
  } else {
    NA_real_
  }
  
  BP_med <- if (H_group == 1L) success_prob_med else NA_real_
  FPR_med <- if (H_group == 0L) success_prob_med else NA_real_
  
  BP_med_mcse <- if (H_group == 1L) mcse_mean(dd$S_R_med == 1) else NA_real_
  FPR_med_mcse <- if (H_group == 0L) mcse_mean(dd$S_R_med == 1) else NA_real_
  
  pointwise_type1_med <- if (
    H_group == 0L &&
    length(design_prior_type) == 1L &&
    identical(design_prior_type, "exact_null_point_mass")
  ) {
    FPR_med
  } else {
    NA_real_
  }
  
  # --------------------------------------------------------------------------
  # Estimation operating characteristics
  # --------------------------------------------------------------------------
  
  interim_bias <- safe_mean(dd$bias_rep_1)
  final_bias_cf <- safe_mean(dd$bias_rep_2)
  final_rmse_cf <- safe_rmse(dd$sqerr_rep_2)
  final_coverage_cf <- safe_mean(dd$cover_rep_2)
  
  interim_bias_med <- safe_mean(dd$bias_med_rep_1)
  final_bias_cf_med <- safe_mean(dd$bias_med_rep_2)
  final_rmse_cf_med <- safe_rmse(dd$sqerr_med_rep_2)
  final_coverage_cf_med <- safe_mean(dd$cover_med_rep_2)
  
  data.frame(
    M = M,
    H_D = H_group,
    
    # RMST decision OCs
    PES = PES,
    PFS = PFS,
    ENE = ENE,
    BP = BP,
    BP_mcse = BP_mcse,
    FPR = FPR,
    FPR_mcse = FPR_mcse,
    Pointwise_TypeI = pointwise_type1,
    
    # RMST estimation OCs
    Interim_bias = interim_bias,
    Final_bias_cf = final_bias_cf,
    Final_RMSE_cf = final_rmse_cf,
    Final_coverage_cf = final_coverage_cf,
    
    # Median decision OCs
    PES_med = PES_med,
    PFS_med = PFS_med,
    ENE_med = ENE_med,
    BP_med = BP_med,
    BP_med_mcse = BP_med_mcse,
    FPR_med = FPR_med,
    FPR_med_mcse = FPR_med_mcse,
    Pointwise_TypeI_med = pointwise_type1_med,
    
    # Median estimation OCs
    Interim_bias_med = interim_bias_med,
    Final_bias_cf_med = final_bias_cf_med,
    Final_RMSE_cf_med = final_rmse_cf_med,
    Final_coverage_cf_med = final_coverage_cf_med,
    
    # Posterior-width and numerical-availability summaries
    mean_post_width_1 = mean(dd$post_width_1, na.rm = TRUE),
    mean_post_width_2 = mean(dd$post_width_2, na.rm = TRUE),
    mean_post_width_med_1 = mean(dd$post_width_med_1, na.rm = TRUE),
    mean_post_width_med_2 = mean(dd$post_width_med_2, na.rm = TRUE),
    mean_na_rate_med_1 = mean(dd$na_rate_med_1, na.rm = TRUE),
    mean_na_rate_med_2 = mean(dd$na_rate_med_2, na.rm = TRUE),
    
    stringsAsFactors = FALSE
  )
}

summarize_oc_by_threshold <- function(res_grid) {
  split_list <- split(
    res_grid,
    list(
      res_grid$prior_id,
      res_grid$profile,
      res_grid$nu1_E,
      res_grid$nu2_E
    ),
    drop = TRUE
  )
  
  out <- do.call(
    rbind,
    lapply(split_list, function(dd) {
      ss <- summarize_oc_one_group(dd)
      
      data.frame(
        prior_id = unique(dd$prior_id),
        scenario_id = if ("scenario_id" %in% names(dd)) {
          unique(dd$scenario_id)
        } else {
          unique(dd$prior_id)
        },
        prior_label = if ("prior_label" %in% names(dd)) {
          unique(dd$prior_label)
        } else {
          get_design_prior_label(unique(dd$prior_id))
        },
        design_prior_type = if ("design_prior_type" %in% names(dd)) {
          unique(dd$design_prior_type)
        } else {
          get_design_prior_label(unique(dd$prior_id))
        },
        omega_D = if ("omega_D" %in% names(dd)) {
          unique(dd$omega_D)
        } else {
          get_design_omega(unique(dd$prior_id))
        },
        profile = unique(dd$profile),
        nu1_E = unique(dd$nu1_E),
        nu2_E = unique(dd$nu2_E),
        ss,
        stringsAsFactors = FALSE
      )
    })
  )
  
  row.names(out) <- NULL
  out
}

# ------------------------------------------------------------------------------
# DPPS, PPV, and NPV from the benefit and continuous no-benefit components
# ------------------------------------------------------------------------------

predictive_ocs_from_components <- function(
    BP,
    FPR,
    omega_D = c(0.8, 0.6, 0.4, 0.2)
) {
  if (length(BP) != 1L || !is.finite(BP) || BP < 0 || BP > 1) {
    stop("`BP` must be a single probability in [0, 1].")
  }
  
  if (length(FPR) != 1L || !is.finite(FPR) || FPR < 0 || FPR > 1) {
    stop("`FPR` must be a single probability in [0, 1].")
  }
  
  if (any(!is.finite(omega_D)) || any(omega_D < 0) || any(omega_D > 1)) {
    stop("All `omega_D` values must lie in [0, 1].")
  }
  
  DPPS <- omega_D * BP + (1 - omega_D) * FPR
  
  PPV_den <- DPPS
  NPV_den <- (1 - omega_D) * (1 - FPR) + omega_D * (1 - BP)
  
  PPV <- ifelse(PPV_den > 0, omega_D * BP / PPV_den, NA_real_)
  NPV <- ifelse(
    NPV_den > 0,
    (1 - omega_D) * (1 - FPR) / NPV_den,
    NA_real_
  )
  
  data.frame(
    omega_D = omega_D,
    BP = BP,
    FPR = FPR,
    DPPS = DPPS,
    PPV = PPV,
    NPV = NPV,
    stringsAsFactors = FALSE
  )
}

# ------------------------------------------------------------------------------
# Exact Continuous-Time Marginal Survival & Quantile Extraction
# ------------------------------------------------------------------------------

calc_marginal_surv_at_t <- function(t_val, Xfixed, eta, dLambda0, t_grid) {
  Xfixed <- as.matrix(Xfixed)
  
  lefts <- t_grid[-length(t_grid)]
  rights <- t_grid[-1]
  widths <- rights - lefts
  lambda_j <- dLambda0 / widths
  
  # Calculate exposure time up to t_val for each interval
  dt <- pmax(0, pmin(t_val, rights) - lefts)
  
  # Cumulative baseline hazard up to t_val
  H0_t <- sum(lambda_j * dt)
  
  # Hazard ratio for each subject
  mu <- exp(-drop(Xfixed %*% as.numeric(eta)))
  
  # Return marginal survival probability
  mean(exp(-mu * H0_t))
}

extract_survival_quantile_uniroot <- function(Xfixed, eta, dLambda0, t_grid, q = 0.5) {
  target <- 1 - q
  max_t <- max(t_grid)
  
  # Check survival probability at the maximum observed/grid time
  S_max <- calc_marginal_surv_at_t(max_t, Xfixed, eta, dLambda0, t_grid)
  
  # Return NA if the target survival probability is not reached within the horizon
  if (S_max > target) return(NA_real_)
  
  # Objective function for root finding
  obj_fun <- function(t_val) {
    calc_marginal_surv_at_t(t_val, Xfixed, eta, dLambda0, t_grid) - target
  }
  
  # Find the exact continuous time using uniroot
  res <- uniroot(obj_fun, lower = 0, upper = max_t, extendInt = "no")
  
  return(res$root)
}

# ------------------------------------------------------------------------------
# Posterior Median Survival Time Contrast Calculation (Exact Root-Finding)
# ------------------------------------------------------------------------------

posterior_median_surv_summary_gp <- function(fit, gp_data, X, q = 0.5, delta_star = 0.00) {
  draws <- extract_gp_draws(fit)
  
  eta_mat <- draws$eta
  dL_mat <- draws$dLambda0
  
  Xstd <- make_X_standardized(X)
  t_grid <- gp_data$t_grid
  
  M <- nrow(eta_mat)
  diff_draws <- rep(NA_real_, M)
  
  for (m in seq_len(M)) {
    med1 <- extract_survival_quantile_uniroot(Xfixed = Xstd$X1, eta = eta_mat[m,], dLambda0 = dL_mat[m,], t_grid = t_grid, q = q)
    
    med0 <- extract_survival_quantile_uniroot(Xfixed = Xstd$X0, eta = eta_mat[m,], dLambda0 = dL_mat[m,], t_grid = t_grid, q = q)
    
    if (is.finite(med1) && is.finite(med0)) {
      diff_draws[m] <- med1 - med0
    }
  }
  
  valid_draws <- diff_draws[is.finite(diff_draws)]
  
  if (length(valid_draws) == 0L) {
    return(list(
      draws = diff_draws,
      Pi_beneficial = NA_real_,
      post_mean = NA_real_,
      post_q025 = NA_real_,
      post_q500 = NA_real_,
      post_q975 = NA_real_,
      post_width = NA_real_,
      na_rate = 1
    ))
  }
  
  ss <- summarize_draws(valid_draws)
  
  list(
    draws = diff_draws,
    Pi_beneficial = unname(mean(valid_draws > delta_star)),
    post_mean = unname(ss["mean"]),
    post_q025 = unname(ss["q025"]),
    post_q500 = unname(ss["q500"]),
    post_q975 = unname(ss["q975"]),
    post_width = unname(ss["width"]),
    na_rate = mean(!is.finite(diff_draws))
  )
}

# ==============================================================================
# INDEPENDENT GAMMA-PROCESS EXTENSION
#
# Purpose:
#   Fits separate Gamma-process baseline hazard models in the treatment and
#   control arms. This avoids imposing a common proportional-hazards treatment
#   coefficient and allows the treatment effect to vary over time.
#
#   Marginal RMST and median contrasts are standardized over the same empirical
#   covariate distribution in both treatment arms.
# ==============================================================================

# ------------------------------------------------------------------------------
# Stan code for a single-arm Gamma-process model
# ------------------------------------------------------------------------------
gp_single_arm_stan_code <- '
data {
  int<lower=1> N;
  int<lower=1> P;
  int<lower=1> J;
  matrix[N, P] X;
  vector<lower=0>[N] y;
  array[N] int<lower=0, upper=1> D;
  vector<lower=0>[J + 1] t_grid;
  vector<lower=0>[J] widths;
  array[J] int<lower=0> dN;
  vector[P] mu_eta;
  vector<lower=0>[P] sigma_eta;
  real<lower=0> a0;
  vector<lower=0>[J] dB0;
}

parameters {
  vector[P] eta_coef;
  vector<lower=1e-12>[J] dLambda0;
}

model {
  vector[N] xb;
  vector[N] exb;
  vector[J] risk_exposure;

  xb = -(X * eta_coef);
  exb = exp(xb);

  eta_coef ~ normal(mu_eta, sigma_eta);
  dLambda0 ~ gamma(a0 * dB0, a0);

  for (j in 1:J) {
    real left = t_grid[j];
    real width = widths[j];
    risk_exposure[j] = 0;

    for (i in 1:N) {
      real dt = fmin(fmax(y[i] - left, 0), width);
      risk_exposure[j] += dt * exb[i] / width;
    }
  }

  target += dot_product(to_vector(dN), log(dLambda0 ./ widths));
  target += dot_product(to_vector(D), xb);
  target += -dot_product(dLambda0, risk_exposure);
}
'
compile_independent_gp_model <- function() {
  if (!requireNamespace("cmdstanr", quietly = TRUE)) stop("Package 'cmdstanr' is required.")
  
  gp_file <- cmdstanr::write_stan_file(gp_single_arm_stan_code)
  cmdstanr::cmdstan_model(gp_file)
}

# ------------------------------------------------------------------------------
# Fit separate Gamma-process models to the treatment and control arms
# ------------------------------------------------------------------------------

fit_gp_independent_look <- function(
    mod_single_gp,
    look_data,
    X,
    B0_fun,
    a0_gp = 2,
    grid_width_gp = 1,
    tau = 24,
    cmdstan_seed = 1,
    chains = 4,
    parallel_chains = 1,
    iter_warmup = 5000,
    iter_sampling = 2500,
    refresh = 200,
    thin = 1,
    adapt_delta = 0.99,
    max_treedepth = 14
) {
  if (!("z" %in% names(look_data))) stop("look_data must contain the treatment indicator `z`.")
  
  idx_trt <- which(look_data$z == 1)
  idx_ctrl <- which(look_data$z == 0)
  
  if (length(idx_trt) == 0 || length(idx_ctrl) == 0) stop("Both treatment arms must contain at least one individual.")
  
  look_trt <- look_data[idx_trt, , drop = FALSE]
  look_ctrl <- look_data[idx_ctrl, , drop = FALSE]
  
  X_covars <- as.matrix(X)
  
  if (is.null(colnames(X_covars))) stop("X must have column names so that treatment `z` can be removed.")
  
  if ("z" %in% colnames(X_covars)) {
    X_covars <- X_covars[, colnames(X_covars) != "z", drop = FALSE]
  }
  
  if (ncol(X_covars) == 0) {
    stop(
      "The current single-arm Stan model requires at least one baseline covariate. ",
      "No covariates remained after removing treatment `z`."
    )
  }
  
  storage.mode(X_covars) <- "double"
  
  X_trt <- X_covars[idx_trt, , drop = FALSE]
  X_ctrl <- X_covars[idx_ctrl, , drop = FALSE]
  
  # Common empirical covariate distribution used to standardize both arms.
  X_common <- X_covars
  
  gp_data_trt <- make_gp_grid_data(look_data = look_trt, X = X_trt, B0_fun = B0_fun, a0 = a0_gp, grid_width = grid_width_gp, tau = tau)
  
  gp_data_ctrl <- make_gp_grid_data(look_data = look_ctrl, X = X_ctrl, B0_fun = B0_fun, a0 = a0_gp, grid_width = grid_width_gp, tau = tau)
  
  init_fun_trt <- function() {
    list(eta_coef = rep(0, gp_data_trt$P), dLambda0 = pmax(gp_data_trt$dB0, 1e-4))
  }
  
  init_fun_ctrl <- function() {
    list(eta_coef = rep(0, gp_data_ctrl$P), dLambda0 = pmax(gp_data_ctrl$dB0, 1e-4))
  }
  
  fit_trt <- mod_single_gp$sample(
    data = gp_data_trt,
    seed = as.integer(cmdstan_seed),
    chains = chains,
    parallel_chains = parallel_chains,
    iter_warmup = iter_warmup,
    iter_sampling = iter_sampling,
    refresh = refresh,
    init = init_fun_trt,
    thin = thin,
    adapt_delta = adapt_delta,
    max_treedepth = max_treedepth
  )
  
  fit_ctrl <- mod_single_gp$sample(
    data = gp_data_ctrl,
    seed = as.integer(cmdstan_seed),
    chains = chains,
    parallel_chains = parallel_chains,
    iter_warmup = iter_warmup,
    iter_sampling = iter_sampling,
    refresh = refresh,
    init = init_fun_ctrl,
    thin = thin,
    adapt_delta = adapt_delta,
    max_treedepth = max_treedepth
  )
  
  list(
    fit_trt = fit_trt,
    gp_data_trt = gp_data_trt,
    fit_ctrl = fit_ctrl,
    gp_data_ctrl = gp_data_ctrl,
    X_trt = X_trt,
    X_ctrl = X_ctrl,
    X_common = X_common
  )
}

# ------------------------------------------------------------------------------
# Posterior RMST contrast for the Independent GP model
# ------------------------------------------------------------------------------

posterior_rmst_summary_independent_gp <- function(indep_fit_obj, tau = 24, delta_star = 0.00) {
  if (is.null(indep_fit_obj$X_common)) stop("indep_fit_obj$X_common is required for common-X standardization.")
  
  draws_trt <- extract_gp_draws(indep_fit_obj$fit_trt)
  draws_ctrl <- extract_gp_draws(indep_fit_obj$fit_ctrl)
  
  dL_trt <- draws_trt$dLambda0
  dL_ctrl <- draws_ctrl$dLambda0
  eta_trt <- draws_trt$eta
  eta_ctrl <- draws_ctrl$eta
  
  X_common <- as.matrix(indep_fit_obj$X_common)
  storage.mode(X_common) <- "double"
  
  M <- min(nrow(dL_trt), nrow(dL_ctrl), nrow(eta_trt), nrow(eta_ctrl))
  
  if (!is.finite(M) || M < 1) stop("No posterior draws were available for the Independent GP contrast.")
  
  rmst_diff <- numeric(M)
  
  for (m in seq_len(M)) {
    rmst1 <- rmst_gp_grid_marginal(
      tau = tau,
      Xfixed = X_common,
      eta = eta_trt[m, ],
      dLambda0 = dL_trt[m, ],
      t_grid = indep_fit_obj$gp_data_trt$t_grid
    )
    
    rmst0 <- rmst_gp_grid_marginal(
      tau = tau,
      Xfixed = X_common,
      eta = eta_ctrl[m, ],
      dLambda0 = dL_ctrl[m, ],
      t_grid = indep_fit_obj$gp_data_ctrl$t_grid
    )
    
    rmst_diff[m] <- rmst1 - rmst0
  }
  
  ss <- summarize_draws(rmst_diff)
  
  list(
    draws = rmst_diff,
    Pi = unname(mean(rmst_diff > delta_star, na.rm = TRUE)),
    post_mean = unname(ss["mean"]),
    post_q025 = unname(ss["q025"]),
    post_q500 = unname(ss["q500"]),
    post_q975 = unname(ss["q975"]),
    post_width = unname(ss["width"])
  )
}

# ------------------------------------------------------------------------------
# Posterior median-survival contrast for the Independent GP model
# ------------------------------------------------------------------------------

posterior_median_surv_summary_independent_gp <- function(indep_fit_obj, q = 0.5, delta_star = 0.00) {
  if (is.null(indep_fit_obj$X_common)) stop("indep_fit_obj$X_common is required for common-X standardization.")
  
  draws_trt <- extract_gp_draws(indep_fit_obj$fit_trt)
  draws_ctrl <- extract_gp_draws(indep_fit_obj$fit_ctrl)
  
  dL_trt <- draws_trt$dLambda0
  dL_ctrl <- draws_ctrl$dLambda0
  eta_trt <- draws_trt$eta
  eta_ctrl <- draws_ctrl$eta
  
  X_common <- as.matrix(indep_fit_obj$X_common)
  storage.mode(X_common) <- "double"
  
  M <- min(nrow(dL_trt), nrow(dL_ctrl), nrow(eta_trt), nrow(eta_ctrl))
  
  if (!is.finite(M) || M < 1) stop("No posterior draws were available for the Independent GP contrast.")
  
  diff_draws <- rep(NA_real_, M)
  
  for (m in seq_len(M)) {
    med1 <- extract_survival_quantile_uniroot(Xfixed = X_common, eta = eta_trt[m,], dLambda0 = dL_trt[m,], t_grid = indep_fit_obj$gp_data_trt$t_grid, q = q)
    
    med0 <- extract_survival_quantile_uniroot(
      Xfixed = X_common,
      eta = eta_ctrl[m, ],
      dLambda0 = dL_ctrl[m, ],
      t_grid = indep_fit_obj$gp_data_ctrl$t_grid,
      q = q
    )
    
    if (is.finite(med1) && is.finite(med0)) {
      diff_draws[m] <- med1 - med0
    }
  }
  
  valid_draws <- diff_draws[is.finite(diff_draws)]
  
  if (length(valid_draws) == 0) {
    return(list(
      draws = diff_draws,
      Pi_beneficial = NA_real_,
      post_mean = NA_real_,
      post_q025 = NA_real_,
      post_q500 = NA_real_,
      post_q975 = NA_real_,
      post_width = NA_real_,
      na_rate = 1
    ))
  }
  
  ss <- summarize_draws(valid_draws)
  
  list(
    draws = diff_draws,
    Pi_beneficial = unname(mean(valid_draws > delta_star)),
    post_mean = unname(ss["mean"]),
    post_q025 = unname(ss["q025"]),
    post_q500 = unname(ss["q500"]),
    post_q975 = unname(ss["q975"]),
    post_width = unname(ss["width"]),
    na_rate = mean(!is.finite(diff_draws))
  )
}

# ==============================================================================
# INDEPENDENT GP WRAPPERS FOR MCMC AND POST-PROCESSING
# ==============================================================================

BDCT_MCMC_indep <- function(
    dat,
    rep_id = NA_integer_,
    mod_single_gp,
    tau = NULL,
    delta_star = NULL,
    B0_fun,
    B0_median_pfs = NA_real_,
    a0_gp = 2,
    grid_width_gp = 1,
    chains = 4,
    parallel_chains = 1,
    iter_warmup = 5000,
    iter_sampling = 2500,
    refresh = 200,
    thin = 1,
    adapt_delta = 0.99,
    max_treedepth = 14,
    cmdstan_seed = 1,
    truth_n = 50000
) {
  cmdstan_seed <- as.integer(cmdstan_seed)
  
  if (is.null(tau)) {
    tau <- dat$tau
  }
  
  if (is.null(delta_star)) {
    delta_star <- dat$delta_star
  }
  
  if (is.null(dat$breaks)) stop("dat$breaks is required.")
  if (is.null(dat$lambda0)) stop("dat$lambda0 is required.")
  if (is.null(dat$eta_x)) stop("dat$eta_x is required.")
  if (is.null(dat$eta_z)) stop("dat$eta_z is required.")
  if (is.null(dat$prior_id)) stop("dat$prior_id is required.")
  if (is.null(dat$prior_label)) stop("dat$prior_label is required.")
  if (is.null(dat$profile)) stop("dat$profile is required.")
  if (is.null(dat$N)) stop("dat$N is required.")
  if (is.null(dat$K1)) stop("dat$K1 is required.")
  if (is.null(dat$K2)) stop("dat$K2 is required.")
  
  true_delta_rep <- if (!is.null(dat$truth) && !is.null(dat$truth$delta_rmst_tau)) {
    dat$truth$delta_rmst_tau
  } else {
    true_delta_rmst_from_eta(
      eta_z = dat$eta_z,
      breaks = dat$breaks,
      lambda0 = dat$lambda0,
      eta_x = dat$eta_x,
      tau = tau,
      x_std = rnorm(truth_n, mean = 0, sd = 1)
    )
  }
  
  true_delta_med_rep <- if (!is.null(dat$truth) && !is.null(dat$truth$delta_med_rep)) {
    dat$truth$delta_med_rep
  } else {
    true_delta_med_from_eta(
      eta_z = dat$eta_z,
      breaks = dat$breaks,
      lambda0 = dat$lambda0,
      eta_x = dat$eta_x,
      x_std = rnorm(truth_n, mean = 0, sd = 1)
    )
  }
  
  looks <- make_analysis_datasets(dat)
  
  O1 <- looks$O1
  O2 <- looks$O2
  X1 <- attr(O1, "X")
  X2 <- attr(O2, "X")
  
  fit1_indep <- fit_gp_independent_look(
    mod_single_gp = mod_single_gp,
    look_data = O1,
    X = X1,
    B0_fun = B0_fun,
    a0_gp = a0_gp,
    grid_width_gp = grid_width_gp,
    tau = tau,
    cmdstan_seed = cmdstan_seed,
    chains = chains,
    parallel_chains = parallel_chains,
    iter_warmup = iter_warmup,
    iter_sampling = iter_sampling,
    refresh = refresh,
    thin = thin,
    adapt_delta = adapt_delta,
    max_treedepth = max_treedepth
  )
  
  fit2_indep <- fit_gp_independent_look(
    mod_single_gp = mod_single_gp,
    look_data = O2,
    X = X2,
    B0_fun = B0_fun,
    a0_gp = a0_gp,
    grid_width_gp = grid_width_gp,
    tau = tau,
    cmdstan_seed = cmdstan_seed,
    chains = chains,
    parallel_chains = parallel_chains,
    iter_warmup = iter_warmup,
    iter_sampling = iter_sampling,
    refresh = refresh,
    thin = thin,
    adapt_delta = adapt_delta,
    max_treedepth = max_treedepth
  )
  
  constants <- list(
    rep_id = rep_id,
    prior_id = dat$prior_id,
    scenario_id = if (!is.null(dat$scenario_id)) dat$scenario_id else dat$prior_id,
    prior_label = dat$prior_label,
    design_prior_type = if (!is.null(dat$design_prior_type)) dat$design_prior_type else dat$prior_label,
    prior_family = if (!is.null(dat$prior_family)) dat$prior_family else NA_character_,
    omega_D = if (!is.null(dat$omega_D)) dat$omega_D else get_design_omega(dat$prior_id),
    profile = dat$profile,
    component = dat$prior_component,
    component_id = dat$component_id,
    H_D = if (!is.null(dat$truth) && !is.null(dat$truth$H_D)) {
      as.integer(dat$truth$H_D)
    } else if (!is.null(dat$H_D)) {
      as.integer(dat$H_D)
    } else {
      as.integer(dat$component_id)
    },
    true_delta_rmst_rep = true_delta_rep,
    true_delta_med_rep = true_delta_med_rep,
    N_planned = dat$N,
    K1 = dat$K1,
    K2 = dat$K2,
    n_interim = attr(O1, "n_enrolled"),
    n_final = attr(O2, "n_enrolled"),
    events_interim = attr(O1, "events"),
    events_final = attr(O2, "events"),
    censor_rate_interim = attr(O1, "censor_rate"),
    censor_rate_final = attr(O2, "censor_rate"),
    calendar_cutoff_interim = attr(O1, "calendar_cutoff"),
    calendar_cutoff_final = attr(O2, "calendar_cutoff"),
    tau = tau,
    delta_star = delta_star,
    a0_gp = a0_gp,
    grid_width_gp = grid_width_gp,
    B0_median_pfs = B0_median_pfs,
    chains = chains,
    parallel_chains = parallel_chains,
    iter_warmup = iter_warmup,
    iter_sampling = iter_sampling,
    thin = thin,
    adapt_delta = adapt_delta,
    max_treedepth = max_treedepth,
    n_saved_draws = chains * floor(iter_sampling / thin),
    cmdstan_seed = cmdstan_seed
  )
  
  dat_info <- list(
    scenario_id = if (!is.null(dat$scenario_id)) dat$scenario_id else dat$prior_id,
    design_prior_type = if (!is.null(dat$design_prior_type)) dat$design_prior_type else dat$prior_label,
    prior_family = if (!is.null(dat$prior_family)) dat$prior_family else NA_character_,
    omega_D = if (!is.null(dat$omega_D)) dat$omega_D else get_design_omega(dat$prior_id),
    breaks = dat$breaks,
    lambda0 = dat$lambda0,
    eta_x = dat$eta_x,
    eta_z = dat$eta_z,
    eta_z_summary = dat$eta_z_summary,
    truth = dat$truth
  )
  
  analysis_data <- list(O1 = O1, O2 = O2, X1 = X1, X2 = X2)
  
  MCMCposteriors <- list(fit1_indep = fit1_indep, fit2_indep = fit2_indep)
  
  MCMCresult <- list(constants = constants, dat_info = dat_info, analysis_data = analysis_data, MCMCposteriors = MCMCposteriors)
  
  class(MCMCresult) <- c("BDCT_MCMC_indep", class(MCMCresult))
  
  MCMCresult
}

BDCT_POST_indep <- function(object, truth_lookup = NULL, tau = NULL, delta_star = NULL) {
  if (is.null(object$constants)) stop("object$constants is required.")
  if (is.null(object$analysis_data)) stop("object$analysis_data is required.")
  if (is.null(object$MCMCposteriors)) stop("object$MCMCposteriors is required.")
  
  constants <- object$constants
  
  if (is.null(tau)) {
    tau <- constants$tau
  }
  
  if (is.null(delta_star)) {
    delta_star <- constants$delta_star
  }
  
  fit1_indep <- object$MCMCposteriors$fit1_indep
  fit2_indep <- object$MCMCposteriors$fit2_indep
  
  post1 <- posterior_rmst_summary_independent_gp(indep_fit_obj = fit1_indep, tau = tau, delta_star = delta_star)
  
  post2 <- posterior_rmst_summary_independent_gp(indep_fit_obj = fit2_indep, tau = tau, delta_star = delta_star)
  
  post_med1 <- posterior_median_surv_summary_independent_gp(indep_fit_obj = fit1_indep, q = 0.5, delta_star = delta_star)
  
  post_med2 <- posterior_median_surv_summary_independent_gp(indep_fit_obj = fit2_indep, q = 0.5, delta_star = delta_star)
  
  # No common treatment coefficient exists in the Independent GP model.
  # These columns are retained as NA to match the GP-PH output structure.
  eta_PH_trt_1 <- data.frame(
    post_eta_PH_trt_mean_1 = NA_real_,
    post_eta_PH_trt_q025_1 = NA_real_,
    post_eta_PH_trt_q500_1 = NA_real_,
    post_eta_PH_trt_q975_1 = NA_real_
  )
  
  eta_PH_trt_2 <- data.frame(
    post_eta_PH_trt_mean_2 = NA_real_,
    post_eta_PH_trt_q025_2 = NA_real_,
    post_eta_PH_trt_q500_2 = NA_real_,
    post_eta_PH_trt_q975_2 = NA_real_
  )
  
  # Calculate each set of MCMC diagnostics once.
  diag1_trt <- extract_mcmc_diagnostics(fit1_indep$fit_trt)
  diag1_ctrl <- extract_mcmc_diagnostics(fit1_indep$fit_ctrl)
  diag2_trt <- extract_mcmc_diagnostics(fit2_indep$fit_trt)
  diag2_ctrl <- extract_mcmc_diagnostics(fit2_indep$fit_ctrl)
  
  true_delta_lookup <- NA_real_
  
  if (!is.null(truth_lookup)) {
    truth_row <- truth_lookup[truth_lookup$prior_id == constants$prior_id & truth_lookup$profile == constants$profile,, drop = FALSE]
    
    true_delta_lookup <- if (nrow(truth_row) == 1) {
      truth_row$true_delta_rmst_lookup
    } else {
      NA_real_
    }
  }
  
  true_delta_rep <- constants$true_delta_rmst_rep
  true_delta_med_rep <- constants$true_delta_med_rep
  H_D_current <- as.integer(constants$H_D)
  
  POSTresult <- data.frame(
    rep = constants$rep_id,
    prior_id = constants$prior_id,
    scenario_id = constants$scenario_id,
    prior_label = constants$prior_label,
    design_prior_type = constants$design_prior_type,
    prior_family = constants$prior_family,
    omega_D = constants$omega_D,
    profile = constants$profile,
    component = constants$component,
    component_id = constants$component_id,
    H_D = H_D_current,
    true_delta_rmst_lookup = true_delta_lookup,
    true_delta_rmst_rep = true_delta_rep,
    true_delta_med_rep = true_delta_med_rep,
    N_planned = constants$N_planned,
    K1 = constants$K1,
    K2 = constants$K2,
    n_interim = constants$n_interim,
    n_final = constants$n_final,
    events_interim = constants$events_interim,
    events_final = constants$events_final,
    censor_rate_interim = constants$censor_rate_interim,
    censor_rate_final = constants$censor_rate_final,
    calendar_cutoff_interim = constants$calendar_cutoff_interim,
    calendar_cutoff_final = constants$calendar_cutoff_final,
    tau = tau,
    delta_star = delta_star,
    a0_gp = constants$a0_gp,
    grid_width_gp = constants$grid_width_gp,
    B0_median_pfs = constants$B0_median_pfs,
    chains = constants$chains,
    parallel_chains = constants$parallel_chains,
    iter_warmup = constants$iter_warmup,
    iter_sampling = constants$iter_sampling,
    thin = constants$thin,
    adapt_delta = constants$adapt_delta,
    max_treedepth = constants$max_treedepth,
    n_saved_draws = constants$n_saved_draws,
    cmdstan_seed = constants$cmdstan_seed,
    
    Pi1 = post1$Pi,
    Pi2_cf = post2$Pi,
    post_mean_1 = post1$post_mean,
    post_q025_1 = post1$post_q025,
    post_q500_1 = post1$post_q500,
    post_q975_1 = post1$post_q975,
    post_width_1 = post1$post_width,
    post_mean_2 = post2$post_mean,
    post_q025_2 = post2$post_q025,
    post_q500_2 = post2$post_q500,
    post_q975_2 = post2$post_q975,
    post_width_2 = post2$post_width,
    
    Pi_beneficial_med_1 = post_med1$Pi_beneficial,
    post_mean_med_1 = post_med1$post_mean,
    post_q025_med_1 = post_med1$post_q025,
    post_q500_med_1 = post_med1$post_q500,
    post_q975_med_1 = post_med1$post_q975,
    post_width_med_1 = post_med1$post_width,
    na_rate_med_1 = post_med1$na_rate,
    
    Pi_beneficial_med_2 = post_med2$Pi_beneficial,
    post_mean_med_2 = post_med2$post_mean,
    post_q025_med_2 = post_med2$post_q025,
    post_q500_med_2 = post_med2$post_q500,
    post_q975_med_2 = post_med2$post_q975,
    post_width_med_2 = post_med2$post_width,
    na_rate_med_2 = post_med2$na_rate,
    
    cover_lookup_1 = if (is.na(true_delta_lookup)) {
      NA_integer_
    } else {
      as.integer(post1$post_q025 <= true_delta_lookup & true_delta_lookup <= post1$post_q975)
    },
    bias_lookup_1 = if (is.na(true_delta_lookup)) {
      NA_real_
    } else {
      post1$post_mean - true_delta_lookup
    },
    sqerr_lookup_1 = if (is.na(true_delta_lookup)) {
      NA_real_
    } else {
      (post1$post_mean - true_delta_lookup)^2
    },
    
    cover_lookup_2 = if (is.na(true_delta_lookup)) {
      NA_integer_
    } else {
      as.integer(post2$post_q025 <= true_delta_lookup & true_delta_lookup <= post2$post_q975)
    },
    bias_lookup_2 = if (is.na(true_delta_lookup)) {
      NA_real_
    } else {
      post2$post_mean - true_delta_lookup
    },
    sqerr_lookup_2 = if (is.na(true_delta_lookup)) {
      NA_real_
    } else {
      (post2$post_mean - true_delta_lookup)^2
    },
    
    cover_rep_1 = as.integer(post1$post_q025 <= true_delta_rep & true_delta_rep <= post1$post_q975),
    bias_rep_1 = post1$post_mean - true_delta_rep,
    sqerr_rep_1 = (post1$post_mean - true_delta_rep)^2,
    
    cover_rep_2 = as.integer(post2$post_q025 <= true_delta_rep & true_delta_rep <= post2$post_q975),
    bias_rep_2 = post2$post_mean - true_delta_rep,
    sqerr_rep_2 = (post2$post_mean - true_delta_rep)^2,
    
    cover_med_rep_1 = if (is.finite(post_med1$post_q025) && is.finite(post_med1$post_q975) && is.finite(true_delta_med_rep)) {
      as.integer(post_med1$post_q025 <= true_delta_med_rep && true_delta_med_rep <= post_med1$post_q975)
    } else {
      NA_integer_
    },
    bias_med_rep_1 = if (is.finite(post_med1$post_mean) && is.finite(true_delta_med_rep)) {
      post_med1$post_mean - true_delta_med_rep
    } else {
      NA_real_
    },
    sqerr_med_rep_1 = if (is.finite(post_med1$post_mean) && is.finite(true_delta_med_rep)) {
      (post_med1$post_mean - true_delta_med_rep)^2
    } else {
      NA_real_
    },
    
    cover_med_rep_2 = if (is.finite(post_med2$post_q025) && is.finite(post_med2$post_q975) && is.finite(true_delta_med_rep)) {
      as.integer(post_med2$post_q025 <= true_delta_med_rep && true_delta_med_rep <= post_med2$post_q975)
    } else {
      NA_integer_
    },
    bias_med_rep_2 = if (is.finite(post_med2$post_mean) && is.finite(true_delta_med_rep)) {
      post_med2$post_mean - true_delta_med_rep
    } else {
      NA_real_
    },
    sqerr_med_rep_2 = if (is.finite(post_med2$post_mean) && is.finite(true_delta_med_rep)) {
      (post_med2$post_mean - true_delta_med_rep)^2
    } else {
      NA_real_
    },
    
    max_Rhat_1 = safe_max(c(diag1_trt$max_Rhat, diag1_ctrl$max_Rhat)),
    min_ESS_1 = safe_min(c(diag1_trt$min_ESS, diag1_ctrl$min_ESS)),
    n_divergent_1 = sum(diag1_trt$n_divergent, diag1_ctrl$n_divergent, na.rm = TRUE),
    max_Rhat_2 = safe_max(c(diag2_trt$max_Rhat, diag2_ctrl$max_Rhat)),
    min_ESS_2 = safe_min(c(diag2_trt$min_ESS, diag2_ctrl$min_ESS)),
    n_divergent_2 = sum(diag2_trt$n_divergent, diag2_ctrl$n_divergent, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
  
  POSTresult <- cbind(POSTresult, eta_PH_trt_1, eta_PH_trt_2)
  
  POSTresult
}