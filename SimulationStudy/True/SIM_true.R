# ==============================================================================
# Calibrate positive-benefit design priors for BDCT
#
# FINAL DESIGN-PRIOR NUMBERING
#
# Prior 1:
#   Exact-null point mass
#
#       eta_{j,D}^Z = 0,  j = 1, ..., 5.
#
#   No kappa calibration is required.
#
#
# Prior 2:
#   Continuous no-benefit design prior
#
#       eta_{j,D}^Z ~ TN(0, 0.30^2; -0.25, 0),
#
#   independently across j = 1, ..., 5.
#
#   No favorable-profile kappa calibration is required.
#
#
# Prior 3:
#   Continuous positive-benefit design prior
#
#       eta_{j,D,+,s}^Z ~ TN(
#           kappa_s * c_{s,j},
#           0.10^2;
#           0,
#           Inf
#       ),
#
#   independently across j = 1, ..., 5.
#
#
# Prior 4:
#   More dispersed continuous positive-benefit design prior
#
#       eta_{j,D,+,s}^Z ~ TN(
#           kappa_s * c_{s,j},
#           0.30^2;
#           0,
#           Inf
#       ),
#
#   independently across j = 1, ..., 5.
#
#
# For Priors 3 and 4, kappa_s is calibrated separately for each temporal
# profile s so that
#
#       E_eta[
#           Delta_RMST(
#               24;
#               eta_D^Z
#           )
#       ] = 2.00.
#
# The expectation is taken jointly over
#
#       X ~ N(0,1)
#
# and the five mutually independent positive-truncated-normal treatment
# effects.
#
#
# Common normalized temporal profile shapes:
#
#   constant:
#       (1.00, 1.00, 1.00, 1.00, 1.00)
#
#   increasing:
#       (0.50, 0.75, 1.00, 1.25, 1.50)
#
#   waning:
#       (1.50, 1.25, 1.00, 0.75, 0.50)
#
# All three shapes have interval mean equal to 1. Therefore
#
#       mean_j(kappa_s * c_{s,j}) = kappa_s.
#
# ==============================================================================


# ------------------------------------------------------------------------------
# Source current BDCT functions
# ------------------------------------------------------------------------------

source_dir <- "/Users/woojung/Documents/Rproject/BDCT/source"
source(file.path(source_dir, "BDCT_r.R"))


# ------------------------------------------------------------------------------
# Common data-generating settings
# ------------------------------------------------------------------------------

breaks <- c(0, 4.8, 9.6, 14.4, 19.2, 24)
lambda0 <- c(0.0734, 0.0661, 0.0587, 0.0514, 0.0440)

eta_x <- 0.2
tau <- 24
target_delta_rmst <- 2.00


# ------------------------------------------------------------------------------
# Positive-benefit design priors to calibrate
# ------------------------------------------------------------------------------

positive_benefit_priors <- data.frame(
  prior_id = c(3L, 4L),
  prior_label = c(
    "benefit_positive_truncnorm_sd010",
    "benefit_positive_truncnorm_sd030"
  ),
  sigma_star = c(0.10, 0.30),
  stringsAsFactors = FALSE
)


# ------------------------------------------------------------------------------
# Favorable temporal profile shapes
# ------------------------------------------------------------------------------

profile_shapes <- list(
  constant = c(1.00, 1.00, 1.00, 1.00, 1.00),
  increasing = c(0.50, 0.75, 1.00, 1.25, 1.50),
  waning = c(1.50, 1.25, 1.00, 0.75, 0.50)
)

profile_means <- vapply(profile_shapes, mean, numeric(1))
stopifnot(all(abs(profile_means - 1) < 1e-12))


# ==============================================================================
# FIXED JOINT INTEGRATION SAMPLE
#
# A fixed stratified integration sample is reused for every kappa, profile,
# and prior. This keeps the calibration objective deterministic and makes the
# Prior-3 versus Prior-4 comparison depend on the specified sigma rather than
# on changing Monte Carlo draws.
# ==============================================================================

n_joint_calibration <- 200000L
seed_calibration <- 870001L


# ------------------------------------------------------------------------------
# Deterministic standard-normal integration points for X
# ------------------------------------------------------------------------------

x_std_calibration <- qnorm(
  (seq_len(n_joint_calibration) - 0.5) / n_joint_calibration
)


# ------------------------------------------------------------------------------
# Fixed stratified Uniform(0,1) points for the five eta intervals
#
# Each column uses the same stratified marginal grid in an independent random
# permutation. Thus each marginal distribution is represented evenly without
# forcing the five interval-specific eta values to move together.
# ------------------------------------------------------------------------------

set.seed(seed_calibration)

u_base <- (seq_len(n_joint_calibration) - 0.5) / n_joint_calibration

u_eta <- matrix(
  NA_real_,
  nrow = n_joint_calibration,
  ncol = 5L
)

for (jj in seq_len(5L)) {
  u_eta[, jj] <- sample(
    u_base,
    size = n_joint_calibration,
    replace = FALSE
  )
}


# ==============================================================================
# POSITIVE-TRUNCATED-NORMAL UTILITIES
# ==============================================================================


# ------------------------------------------------------------------------------
# Inverse CDF for TN(mu, sigma^2; 0, Inf)
# ------------------------------------------------------------------------------

qtnorm_positive_from_u <- function(u, mu, sigma) {
  if (length(mu) != 1L || !is.finite(mu)) {
    stop("mu must be one finite scalar.")
  }
  
  if (length(sigma) != 1L || !is.finite(sigma) || sigma <= 0) {
    stop("sigma must be one positive finite scalar.")
  }
  
  p_lower <- pnorm(0, mean = mu, sd = sigma)
  p <- p_lower + u * (1 - p_lower)
  
  p <- pmin(
    pmax(p, .Machine$double.eps),
    1 - .Machine$double.eps
  )
  
  qnorm(p, mean = mu, sd = sigma)
}


# ------------------------------------------------------------------------------
# Exact moments of TN(mu, sigma^2; 0, Inf)
#
# These are used for reporting/checking the marginal eta distributions.
# ------------------------------------------------------------------------------

positive_tnorm_moments <- function(mu, sigma) {
  alpha <- -mu / sigma
  tail_prob <- pnorm(alpha, lower.tail = FALSE)
  mills <- dnorm(alpha) / tail_prob
  
  mean_value <- mu + sigma * mills
  variance_value <- sigma^2 * (1 + alpha * mills - mills^2)
  
  c(
    mean = mean_value,
    variance = variance_value,
    sd = sqrt(variance_value)
  )
}


# ==============================================================================
# ROW-WISE RMST UNDER THE PIECEWISE-EXPONENTIAL DGP
# ==============================================================================


# ------------------------------------------------------------------------------
# eta_z is either
#
#   - an n x 5 matrix for the treatment arm, or
#   - NULL for the control arm, corresponding to eta_z = 0.
#
# Each row corresponds to one joint draw of X and eta_z.
# ------------------------------------------------------------------------------

rmst_piecewise_rows <- function(x, eta_z = NULL, breaks, lambda0, eta_x, tau) {
  x <- as.numeric(x)
  
  n <- length(x)
  J <- length(lambda0)
  
  if (length(breaks) != J + 1L) {
    stop("length(breaks) must equal length(lambda0) + 1.")
  }
  
  if (is.null(eta_z)) {
    eta_z <- matrix(0, nrow = n, ncol = J)
  } else {
    eta_z <- as.matrix(eta_z)
    
    if (nrow(eta_z) != n || ncol(eta_z) != J) {
      stop("eta_z must have n rows and length(lambda0) columns.")
    }
  }
  
  interval_width <- pmax(
    0,
    pmin(tau, breaks[-1]) - breaks[-length(breaks)]
  )
  
  rmst <- numeric(n)
  survival_start <- rep(1, n)
  
  for (jj in seq_len(J)) {
    width_j <- interval_width[jj]
    
    if (width_j <= 0) {
      next
    }
    
    hazard_j <- lambda0[jj] * exp(
      -eta_z[, jj] - eta_x * x
    )
    
    interval_rmst <- survival_start *
      (-expm1(-hazard_j * width_j)) / hazard_j
    
    rmst <- rmst + interval_rmst
    survival_start <- survival_start * exp(-hazard_j * width_j)
  }
  
  rmst
}


# ------------------------------------------------------------------------------
# Control-arm RMST
#
# This does not depend on kappa, profile, or positive-benefit prior.
# ------------------------------------------------------------------------------

rmst_control <- rmst_piecewise_rows(
  x = x_std_calibration,
  eta_z = NULL,
  breaks = breaks,
  lambda0 = lambda0,
  eta_x = eta_x,
  tau = tau
)

mean_rmst_control <- mean(rmst_control)


# ==============================================================================
# IMPLEMENTATION CHECK
#
# Validate the row-wise RMST calculation against the current BDCT truth
# function using a fixed eta vector.
# ==============================================================================

eta_validation <- rep(0.2850356, 5L)

eta_validation_matrix <- matrix(
  eta_validation,
  nrow = n_joint_calibration,
  ncol = 5L,
  byrow = TRUE
)

rmst_treatment_validation <- mean(
  rmst_piecewise_rows(
    x = x_std_calibration,
    eta_z = eta_validation_matrix,
    breaks = breaks,
    lambda0 = lambda0,
    eta_x = eta_x,
    tau = tau
  )
)

delta_validation_rowwise <- rmst_treatment_validation - mean_rmst_control

delta_validation_existing <- true_delta_rmst_from_eta(
  eta_z = eta_validation,
  breaks = breaks,
  lambda0 = lambda0,
  eta_x = eta_x,
  tau = tau,
  x_std = x_std_calibration
)

cat(
  "\n",
  "============================================================\n",
  "RMST IMPLEMENTATION CHECK\n",
  "============================================================\n",
  "Existing truth function = ", sprintf("%.10f", delta_validation_existing), "\n",
  "Row-wise calculation    = ", sprintf("%.10f", delta_validation_rowwise), "\n",
  "Absolute difference     = ",
  sprintf("%.12f", abs(delta_validation_existing - delta_validation_rowwise)),
  "\n",
  sep = ""
)

if (abs(delta_validation_existing - delta_validation_rowwise) > 1e-6) {
  stop(
    "Row-wise RMST implementation does not agree with ",
    "true_delta_rmst_from_eta()."
  )
}


# ==============================================================================
# GENERIC POSITIVE-BENEFIT CALIBRATION FUNCTIONS
# ==============================================================================


# ------------------------------------------------------------------------------
# Construct eta draws for a given kappa, temporal shape, and sigma
# ------------------------------------------------------------------------------

get_eta_draws_positive_tn <- function(kappa, shape, sigma) {
  if (length(shape) != 5L) {
    stop("shape must contain five interval-specific values.")
  }
  
  eta_location <- kappa * shape
  
  eta_draws <- matrix(
    NA_real_,
    nrow = n_joint_calibration,
    ncol = 5L
  )
  
  for (jj in seq_len(5L)) {
    eta_draws[, jj] <- qtnorm_positive_from_u(
      u = u_eta[, jj],
      mu = eta_location[jj],
      sigma = sigma
    )
  }
  
  colnames(eta_draws) <- paste0("eta_", 1:5)
  eta_draws
}


# ------------------------------------------------------------------------------
# Expected RMST contrast under a positive-truncated-normal benefit prior
# ------------------------------------------------------------------------------

mean_delta_rmst_positive_tn <- function(kappa, shape, sigma) {
  eta_draws <- get_eta_draws_positive_tn(
    kappa = kappa,
    shape = shape,
    sigma = sigma
  )
  
  rmst_treatment <- rmst_piecewise_rows(
    x = x_std_calibration,
    eta_z = eta_draws,
    breaks = breaks,
    lambda0 = lambda0,
    eta_x = eta_x,
    tau = tau
  )
  
  mean(rmst_treatment - rmst_control)
}


# ------------------------------------------------------------------------------
# Calibrate one prior/profile combination
#
# Important for the SD = 0.30 prior:
# even kappa = 0 does not imply eta = 0 after positive truncation.
# Therefore the code explicitly evaluates the objective at kappa = 0 and
# confirms that the target can be bracketed before calling uniroot().
# ------------------------------------------------------------------------------

calibrate_one_positive_tn_profile <- function(
    prior_id,
    prior_label,
    sigma,
    profile,
    shape,
    target = 2.00,
    root_upper = 2.00
) {
  objective <- function(kappa) {
    mean_delta_rmst_positive_tn(
      kappa = kappa,
      shape = shape,
      sigma = sigma
    ) - target
  }
  
  objective_lower <- objective(0)
  objective_upper <- objective(root_upper)
  
  if (objective_lower > 0) {
    stop(
      "For prior ", prior_id, ", profile ", profile,
      ", E_eta[Delta_RMST(24)] already exceeds the target at kappa = 0.\n",
      "objective(0) = ", objective_lower,
      "\nA nonnegative kappa cannot attain the requested target."
    )
  }
  
  if (objective_upper < 0) {
    stop(
      "Could not bracket the kappa root for prior ", prior_id,
      ", profile ", profile,
      ".\nobjective(0) = ", objective_lower,
      "\nobjective(", root_upper, ") = ", objective_upper
    )
  }
  
  kappa_hat <- uniroot(
    f = objective,
    interval = c(0, root_upper),
    tol = 1e-9
  )$root
  
  eta_location <- kappa_hat * shape
  
  eta_draws <- get_eta_draws_positive_tn(
    kappa = kappa_hat,
    shape = shape,
    sigma = sigma
  )
  
  eta_draw_mean <- colMeans(eta_draws)
  
  marginal_moments <- t(
    vapply(
      eta_location,
      function(mu_j) {
        positive_tnorm_moments(
          mu = mu_j,
          sigma = sigma
        )
      },
      numeric(3)
    )
  )
  
  achieved_mean_rmst <- mean_delta_rmst_positive_tn(
    kappa = kappa_hat,
    shape = shape,
    sigma = sigma
  )
  
  location_rmst <- true_delta_rmst_from_eta(
    eta_z = eta_location,
    breaks = breaks,
    lambda0 = lambda0,
    eta_x = eta_x,
    tau = tau,
    x_std = x_std_calibration
  )
  
  location_median <- true_delta_med_from_eta(
    eta_z = eta_location,
    breaks = breaks,
    lambda0 = lambda0,
    eta_x = eta_x,
    x_std = x_std_calibration
  )
  
  data.frame(
    prior_id = prior_id,
    prior_label = prior_label,
    profile = profile,
    
    sigma_star = sigma,
    sigma2_star = sigma^2,
    
    kappa = kappa_hat,
    
    objective_at_kappa_0 = objective_lower,
    
    location_eta_1 = eta_location[1],
    location_eta_2 = eta_location[2],
    location_eta_3 = eta_location[3],
    location_eta_4 = eta_location[4],
    location_eta_5 = eta_location[5],
    
    mean_location_eta = mean(eta_location),
    
    exact_mean_eta_1 = marginal_moments[1, "mean"],
    exact_mean_eta_2 = marginal_moments[2, "mean"],
    exact_mean_eta_3 = marginal_moments[3, "mean"],
    exact_mean_eta_4 = marginal_moments[4, "mean"],
    exact_mean_eta_5 = marginal_moments[5, "mean"],
    
    exact_var_eta_1 = marginal_moments[1, "variance"],
    exact_var_eta_2 = marginal_moments[2, "variance"],
    exact_var_eta_3 = marginal_moments[3, "variance"],
    exact_var_eta_4 = marginal_moments[4, "variance"],
    exact_var_eta_5 = marginal_moments[5, "variance"],
    
    exact_sd_eta_1 = marginal_moments[1, "sd"],
    exact_sd_eta_2 = marginal_moments[2, "sd"],
    exact_sd_eta_3 = marginal_moments[3, "sd"],
    exact_sd_eta_4 = marginal_moments[4, "sd"],
    exact_sd_eta_5 = marginal_moments[5, "sd"],
    
    mc_mean_eta_1 = eta_draw_mean[1],
    mc_mean_eta_2 = eta_draw_mean[2],
    mc_mean_eta_3 = eta_draw_mean[3],
    mc_mean_eta_4 = eta_draw_mean[4],
    mc_mean_eta_5 = eta_draw_mean[5],
    
    mean_exact_eta = mean(marginal_moments[, "mean"]),
    
    Mean_Delta_RMST_24 = achieved_mean_rmst,
    Location_Delta_RMST_24 = location_rmst,
    Location_Delta_Median = location_median,
    
    stringsAsFactors = FALSE
  )
}


# ==============================================================================
# RUN PRIORS 3 AND 4
# ==============================================================================

calibration_results <- do.call(
  rbind,
  lapply(
    seq_len(nrow(positive_benefit_priors)),
    function(ii) {
      prior_id <- positive_benefit_priors$prior_id[ii]
      prior_label <- positive_benefit_priors$prior_label[ii]
      sigma <- positive_benefit_priors$sigma_star[ii]
      
      do.call(
        rbind,
        lapply(
          names(profile_shapes),
          function(profile) {
            message(
              "Calibrating prior ", prior_id,
              " (SD = ", sprintf("%.2f", sigma),
              "), profile: ", profile
            )
            
            calibrate_one_positive_tn_profile(
              prior_id = prior_id,
              prior_label = prior_label,
              sigma = sigma,
              profile = profile,
              shape = profile_shapes[[profile]],
              target = target_delta_rmst,
              root_upper = 2.00
            )
          }
        )
      )
    }
  )
)

row.names(calibration_results) <- NULL

calibration_results_p3 <- calibration_results[
  calibration_results$prior_id == 3L,
  ,
  drop = FALSE
]

calibration_results_p4 <- calibration_results[
  calibration_results$prior_id == 4L,
  ,
  drop = FALSE
]


# ==============================================================================
# COMPLETE RESULTS
# ==============================================================================

cat(
  "\n",
  "============================================================\n",
  "COMPLETE CALIBRATION RESULTS\n",
  "============================================================\n",
  "Prior 3: eta_j ~ TN(kappa_s * c_sj, 0.10^2; 0, Inf)\n",
  "Prior 4: eta_j ~ TN(kappa_s * c_sj, 0.30^2; 0, Inf)\n",
  "Target for both: E_eta[Delta_RMST(24)] = 2.00\n",
  "============================================================\n",
  sep = ""
)

print(
  calibration_results,
  digits = 10,
  row.names = FALSE
)


# ==============================================================================
# CALIBRATED KAPPA VALUES
# ==============================================================================

cat(
  "\n",
  "============================================================\n",
  "CALIBRATED KAPPA VALUES\n",
  "============================================================\n",
  sep = ""
)

for (ii in seq_len(nrow(calibration_results))) {
  cat(
    "prior ", calibration_results$prior_id[ii],
    " | SD = ", sprintf("%.2f", calibration_results$sigma_star[ii]),
    " | ", calibration_results$profile[ii],
    " | kappa = ", sprintf("%.10f", calibration_results$kappa[ii]),
    "\n",
    sep = ""
  )
}


# ==============================================================================
# LOCATION ETA VECTORS
# ==============================================================================

cat(
  "\n",
  "============================================================\n",
  "LOCATION ETA VECTORS: kappa_s * c_s\n",
  "============================================================\n",
  sep = ""
)

for (ii in seq_len(nrow(calibration_results))) {
  eta_vec <- as.numeric(
    calibration_results[
      ii,
      paste0("location_eta_", 1:5)
    ]
  )
  
  cat(
    "prior ", calibration_results$prior_id[ii],
    " | SD = ", sprintf("%.2f", calibration_results$sigma_star[ii]),
    " | ", calibration_results$profile[ii],
    ": (",
    paste(sprintf("%.10f", eta_vec), collapse = ", "),
    ")\n",
    sep = ""
  )
}


# ==============================================================================
# EXACT MARGINAL MEAN ETA VECTORS AFTER POSITIVE TRUNCATION
# ==============================================================================

cat(
  "\n",
  "============================================================\n",
  "EXACT MARGINAL MEAN ETA VECTORS AFTER TRUNCATION\n",
  "============================================================\n",
  sep = ""
)

for (ii in seq_len(nrow(calibration_results))) {
  eta_mean_vec <- as.numeric(
    calibration_results[
      ii,
      paste0("exact_mean_eta_", 1:5)
    ]
  )
  
  cat(
    "prior ", calibration_results$prior_id[ii],
    " | SD = ", sprintf("%.2f", calibration_results$sigma_star[ii]),
    " | ", calibration_results$profile[ii],
    ": (",
    paste(sprintf("%.10f", eta_mean_vec), collapse = ", "),
    ")\n",
    sep = ""
  )
}


# ==============================================================================
# CALIBRATION CHECK
# ==============================================================================

cat(
  "\n",
  "============================================================\n",
  "CALIBRATION CHECK\n",
  "============================================================\n",
  sep = ""
)

for (ii in seq_len(nrow(calibration_results))) {
  cat(
    "prior ", calibration_results$prior_id[ii],
    " | ", calibration_results$profile[ii], "\n",
    "  sigma_star                      = ",
    sprintf("%.10f", calibration_results$sigma_star[ii]), "\n",
    "  kappa                           = ",
    sprintf("%.10f", calibration_results$kappa[ii]), "\n",
    "  E_eta[Delta_RMST(24)] at k=0    = ",
    sprintf(
      "%.10f",
      calibration_results$objective_at_kappa_0[ii] + target_delta_rmst
    ), "\n",
    "  mean location eta               = ",
    sprintf("%.10f", calibration_results$mean_location_eta[ii]), "\n",
    "  mean exact eta after truncation = ",
    sprintf("%.10f", calibration_results$mean_exact_eta[ii]), "\n",
    "  E_eta[Delta_RMST(24)]           = ",
    sprintf("%.10f", calibration_results$Mean_Delta_RMST_24[ii]), "\n",
    "  Delta_RMST at location vector   = ",
    sprintf("%.10f", calibration_results$Location_Delta_RMST_24[ii]), "\n",
    "  Delta_Median at location vector = ",
    sprintf("%.10f", calibration_results$Location_Delta_Median[ii]), "\n\n",
    sep = ""
  )
}


# ==============================================================================
# PRIOR 3 VS PRIOR 4 COMPARISON
# ==============================================================================

prior_3_vs_4 <- merge(
  calibration_results_p3[
    ,
    c(
      "profile",
      "sigma_star",
      "kappa",
      "mean_location_eta",
      "mean_exact_eta",
      "Mean_Delta_RMST_24",
      "Location_Delta_RMST_24"
    )
  ],
  calibration_results_p4[
    ,
    c(
      "profile",
      "sigma_star",
      "kappa",
      "mean_location_eta",
      "mean_exact_eta",
      "Mean_Delta_RMST_24",
      "Location_Delta_RMST_24"
    )
  ],
  by = "profile",
  suffixes = c("_p3", "_p4"),
  sort = FALSE
)

prior_3_vs_4$kappa_difference <- prior_3_vs_4$kappa_p4 - prior_3_vs_4$kappa_p3
prior_3_vs_4$kappa_ratio <- prior_3_vs_4$kappa_p4 / prior_3_vs_4$kappa_p3
prior_3_vs_4$mean_exact_eta_difference <-
  prior_3_vs_4$mean_exact_eta_p4 - prior_3_vs_4$mean_exact_eta_p3

cat(
  "\n",
  "============================================================\n",
  "PRIOR 3 VS PRIOR 4 COMPARISON\n",
  "============================================================\n",
  sep = ""
)

print(
  prior_3_vs_4,
  digits = 10,
  row.names = FALSE
)


# ==============================================================================
# FINAL COMPACT SUMMARY
# ==============================================================================

final_kappa_summary <- calibration_results[
  ,
  c(
    "prior_id",
    "prior_label",
    "profile",
    "sigma_star",
    "kappa",
    "mean_location_eta",
    "mean_exact_eta",
    "Mean_Delta_RMST_24",
    "Location_Delta_RMST_24",
    "Location_Delta_Median"
  )
]

cat(
  "\n",
  "============================================================\n",
  "FINAL KAPPA SUMMARY: PRIORS 3 AND 4\n",
  "============================================================\n",
  sep = ""
)

print(
  final_kappa_summary,
  digits = 10,
  row.names = FALSE
)


# ==============================================================================
# VALUES TO COPY INTO BDCT_r.R
# ==============================================================================

cat(
  "\n",
  "============================================================\n",
  "KAPPA VALUES TO COPY INTO BDCT_r.R\n",
  "============================================================\n",
  sep = ""
)

for (prior_id in c(3L, 4L)) {
  tmp <- calibration_results[
    calibration_results$prior_id == prior_id,
    ,
    drop = FALSE
  ]
  
  tmp <- tmp[
    match(
      c("constant", "increasing", "waning"),
      tmp$profile
    ),
    ,
    drop = FALSE
  ]
  
  cat(
    "\n",
    "prior ", prior_id,
    " (SD = ", sprintf("%.2f", unique(tmp$sigma_star)), ")\n",
    "kappa_constant   <- ", sprintf("%.10f", tmp$kappa[tmp$profile == "constant"]), "\n",
    "kappa_increasing <- ", sprintf("%.10f", tmp$kappa[tmp$profile == "increasing"]), "\n",
    "kappa_waning     <- ", sprintf("%.10f", tmp$kappa[tmp$profile == "waning"]), "\n",
    sep = ""
  )
}