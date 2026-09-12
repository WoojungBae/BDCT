# ==============================================================================
# BDCT FULL SMOKE TEST -- FINAL FOUR DESIGN-PRIOR VERSION
#
# Bayesian Design for Clinical Trial (BDCT)
#
# Final design-stage priors:
#
#   Prior 1:
#     eta_j = 0, j = 1, ..., 5.
#
#   Prior 2:
#     eta_j ~ TN(0, 0.30^2; -0.25, 0), independently across intervals.
#
#   Prior 3:
#     eta_j ~ TN(kappa_s * c_{s,j}, 0.10^2; 0, Inf),
#     independently across intervals.
#
#   Prior 4:
#     eta_j ~ TN(kappa_s * c_{s,j}, 0.30^2; 0, Inf),
#     independently across intervals.
#
# Priors 3 and 4 are calibrated separately by profile so that
#
#     E_eta[Delta_RMST(24)] = 2.00 months.
#
# Checks:
#   1. BDCT naming and required functions.
#   2. Final prior_id = 1, ..., 4 mapping.
#   3. Temporal profile shapes and calibrated kappa values.
#   4. Calibration-report contents.
#   5. Direct distribution-level RMST calibration for Priors 3 and 4.
#   6. Design-prior draw distributions for all four priors.
#   7. generate_data() for all eight DGP evaluation settings.
#   8. H_D agrees with generated RMST truth.
#   9. Event-driven interim/final datasets are valid.
#  10. GP-PH and Independent GP Stan models compile.
#  11. GP-PH and Independent GP short MCMC/post-processing run for one
#      representative trial from each prior.
#  12. Current post_eta_PH_* output names are present and old BGSD/gamma/beta
#      treatment-coefficient names are absent.
#  13. Threshold-grid post-processing works on a synthetic example.
#  14. DPPS/PPV/NPV identities are correct.
#  15. Schoenfeld event calculations reproduce the current anchors.
#
# This is a smoke test, not the production simulation.
# ============================================================================== 


Sys.setenv(
  cmdstanr_no_ver_check = "TRUE",
  MAKEFLAGS = "-j1"
)

suppressPackageStartupMessages({
  library(cmdstanr)
})


# ==============================================================================
# SOURCE
# ==============================================================================

source_dir <- "/Users/woojung/Documents/Rproject/BDCT/source"

if (!dir.exists(source_dir)) {
  stop("source_dir was not found: ", source_dir)
}

source_file <- file.path(source_dir, "BDCT_r.R")

if (!file.exists(source_file)) {
  stop("BDCT_r.R was not found in source_dir: ", source_dir)
}

# Read source text before sourcing so naming checks are based on the file itself,
# not on objects that may remain from an older R session.
source_text <- readLines(source_file, warn = FALSE)
source(source_file)

cat("\n")
cat("============================================================\n")
cat("BDCT FULL SMOKE TEST -- FINAL FOUR DESIGN PRIORS\n")
cat("============================================================\n")


# ==============================================================================
# HELPERS
# ==============================================================================

check_true <- function(condition, message) {
  if (!isTRUE(condition)) {
    stop("\nCHECK FAILED:\n", message, call. = FALSE)
  }
}

check_equal <- function(x, y, message, tolerance = 1e-8) {
  ok <- isTRUE(
    all.equal(
      x,
      y,
      tolerance = tolerance,
      check.attributes = FALSE
    )
  )
  
  if (!ok) {
    stop(
      "\nCHECK FAILED:\n",
      message,
      "\nObserved: ", paste(x, collapse = ", "),
      "\nExpected: ", paste(y, collapse = ", "),
      call. = FALSE
    )
  }
}

check_between <- function(x, lower, upper, message) {
  check_true(
    all(is.finite(x)) && all(x >= lower) && all(x <= upper),
    message
  )
}

benefit_profiles <- c("constant", "increasing", "waning")


# ==============================================================================
# 1. BDCT NAMING, REQUIRED FUNCTIONS, AND FUNCTION INTERFACES
# ==============================================================================

cat("\n")
cat("------------------------------------------------------------\n")
cat("[1] CHECKING BDCT NAMING AND FUNCTION INTERFACES\n")
cat("------------------------------------------------------------\n")

# The current source file itself should not contain the old project name.
check_true(
  !any(grepl("BGSD", source_text, fixed = TRUE)),
  "BDCT_r.R still contains the old project name `BGSD`."
)

# Old uniform-prior machinery should be absent from the final four-prior source.
old_uniform_tokens <- c(
  "benefit_uniform",
  "no_benefit_uniform",
  "benefit_amplitude_U",
  "uniform_benefit",
  "uniform_no_benefit"
)

for (token in old_uniform_tokens) {
  check_true(
    !any(grepl(token, source_text, fixed = TRUE)),
    paste0("BDCT_r.R still contains obsolete uniform-prior token: ", token)
  )
}

required_functions <- c(
  "get_eta_z_profile_shapes",
  "get_eta_z_kappas",
  "get_eta_z_centers",
  "get_eta_z_calibration_report",
  "get_design_prior_spec",
  "get_design_omega",
  "get_design_prior_label",
  "draw_eta_z_from_prior",
  "generate_data",
  "true_delta_rmst_from_eta",
  "true_delta_med_from_eta",
  "compile_gp_model",
  "compile_independent_gp_model",
  "BDCT_MCMC",
  "BDCT_POST",
  "BDCT_MCMC_ph",
  "BDCT_POST_ph",
  "BDCT_MCMC_indep",
  "BDCT_POST_indep",
  "apply_threshold_grid",
  "summarize_oc_one_group",
  "summarize_oc_by_threshold",
  "predictive_ocs_from_components"
)

missing_functions <- required_functions[
  !vapply(
    required_functions,
    exists,
    logical(1),
    mode = "function",
    inherits = TRUE
  )
]

check_true(
  length(missing_functions) == 0L,
  paste0(
    "Missing required functions: ",
    paste(missing_functions, collapse = ", ")
  )
)

# Old project-level function names should not be defined by BDCT_r.R.
old_function_patterns <- c(
  "BGSD_MCMC",
  "BGSD_POST",
  "BGSD_MCMC_ph",
  "BGSD_POST_ph",
  "BGSD_MCMC_indep",
  "BGSD_POST_indep"
)

for (token in old_function_patterns) {
  check_true(
    !any(grepl(token, source_text, fixed = TRUE)),
    paste0("BDCT_r.R still contains obsolete function name: ", token)
  )
}

# Old gamma-based DGP function names should also be absent.
old_eta_function_tokens <- c(
  "get_gamma_z_centers",
  "draw_gamma_z_from_prior",
  "true_delta_rmst_from_gamma",
  "true_delta_med_from_gamma"
)

for (token in old_eta_function_tokens) {
  check_true(
    !any(grepl(token, source_text, fixed = TRUE)),
    paste0("BDCT_r.R still contains obsolete gamma-based name: ", token)
  )
}

# Final draw function should expose only prior_id and profile.
draw_formals <- names(formals(draw_eta_z_from_prior))
check_true(
  identical(draw_formals, c("prior_id", "profile")),
  paste0(
    "draw_eta_z_from_prior() has unexpected arguments: ",
    paste(draw_formals, collapse = ", ")
  )
)

# generate_data() should not expose prior-distribution tuning arguments anymore.
generate_formals <- names(formals(generate_data))

required_generate_args <- c(
  "prior_id",
  "profile",
  "N",
  "K1",
  "K2",
  "target_censor",
  "breaks",
  "lambda0",
  "eta_x",
  "pi_z",
  "tau",
  "delta_star",
  "dropout",
  "dropout_rate",
  "truth_n",
  "return_truth",
  "verbose"
)

check_true(
  all(required_generate_args %in% generate_formals),
  paste0(
    "generate_data() is missing arguments: ",
    paste(setdiff(required_generate_args, generate_formals), collapse = ", ")
  )
)

obsolete_generate_args <- c(
  "benefit_sd",
  "benefit_width",
  "benefit_min",
  "benefit_max",
  "positive_benefit_sd",
  "no_benefit_location",
  "no_benefit_sd",
  "no_benefit_min",
  "no_benefit_max",
  "uniform_benefit_u_min",
  "uniform_benefit_u_max",
  "uniform_no_benefit_min",
  "uniform_no_benefit_max"
)

check_true(
  !any(obsolete_generate_args %in% generate_formals),
  paste0(
    "generate_data() still exposes obsolete prior-tuning arguments: ",
    paste(intersect(obsolete_generate_args, generate_formals), collapse = ", ")
  )
)

cat("[PASS] BDCT naming and final function interfaces are correct.\n")


# ==============================================================================
# 2. FINAL FOUR-PRIOR MAPPING
# ==============================================================================

cat("\n")
cat("------------------------------------------------------------\n")
cat("[2] CHECKING FINAL PRIOR MAPPING\n")
cat("------------------------------------------------------------\n")

expected_prior_specs <- data.frame(
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

for (ii in seq_len(nrow(expected_prior_specs))) {
  pid <- expected_prior_specs$prior_id[ii]
  observed <- get_design_prior_spec(pid)
  
  check_true(
    nrow(observed) == 1L,
    paste0("Unexpected number of rows for prior_id = ", pid)
  )
  
  check_equal(
    observed$prior_id,
    pid,
    paste0("prior_id is incorrect for prior_id = ", pid)
  )
  
  check_true(
    identical(
      as.character(observed$prior_label),
      expected_prior_specs$prior_label[ii]
    ),
    paste0("prior_label is incorrect for prior_id = ", pid)
  )
  
  check_true(
    identical(
      as.character(observed$prior_component),
      expected_prior_specs$prior_component[ii]
    ),
    paste0("prior_component is incorrect for prior_id = ", pid)
  )
  
  check_true(
    identical(
      as.character(observed$prior_family),
      expected_prior_specs$prior_family[ii]
    ),
    paste0("prior_family is incorrect for prior_id = ", pid)
  )
  
  check_equal(
    observed$omega_D,
    expected_prior_specs$omega_D[ii],
    paste0("omega_D is incorrect for prior_id = ", pid)
  )
  
  check_equal(
    observed$H_D,
    expected_prior_specs$H_D[ii],
    paste0("H_D is incorrect for prior_id = ", pid)
  )
  
  check_equal(
    get_design_omega(pid),
    expected_prior_specs$omega_D[ii],
    paste0("get_design_omega() is incorrect for prior_id = ", pid)
  )
  
  check_true(
    identical(
      get_design_prior_label(pid),
      expected_prior_specs$prior_label[ii]
    ),
    paste0("get_design_prior_label() is incorrect for prior_id = ", pid)
  )
}

invalid_prior_ids <- c(-1L, 0L, 5L, 7L)

for (pid in invalid_prior_ids) {
  err <- try(get_design_prior_spec(pid), silent = TRUE)
  check_true(
    inherits(err, "try-error"),
    paste0("prior_id = ", pid, " should be invalid.")
  )
}

cat("[PASS] prior_id = 1, ..., 4 mapping is correct.\n")


# ==============================================================================
# 3. TEMPORAL PROFILE SHAPES, KAPPAS, AND LOCATION VECTORS
# ==============================================================================

cat("\n")
cat("------------------------------------------------------------\n")
cat("[3] CHECKING TEMPORAL PROFILES AND KAPPAS\n")
cat("------------------------------------------------------------\n")

expected_shapes <- list(
  constant = c(1.00, 1.00, 1.00, 1.00, 1.00),
  increasing = c(0.50, 0.75, 1.00, 1.25, 1.50),
  waning = c(1.50, 1.25, 1.00, 0.75, 0.50)
)

observed_shapes <- get_eta_z_profile_shapes()

check_true(
  identical(names(observed_shapes), names(expected_shapes)),
  "Unexpected benefit-profile shape names."
)

for (pp in benefit_profiles) {
  check_equal(
    observed_shapes[[pp]],
    expected_shapes[[pp]],
    paste0("Temporal profile shape is incorrect for ", pp),
    tolerance = 1e-12
  )
  
  check_equal(
    mean(observed_shapes[[pp]]),
    1,
    paste0("Temporal profile shape does not have mean 1 for ", pp),
    tolerance = 1e-12
  )
}

expected_kappas_p3 <- c(
  constant = 0.2881833829,
  increasing = 0.3985155474,
  waning = 0.2265051174
)

expected_kappas_p4 <- c(
  constant = 0.1469563494,
  increasing = 0.1990671295,
  waning = 0.1149792825
)

observed_kappas_p3 <- get_eta_z_kappas(prior_id = 3L)
observed_kappas_p4 <- get_eta_z_kappas(prior_id = 4L)

check_equal(
  observed_kappas_p3,
  expected_kappas_p3,
  "Prior-3 kappa values are incorrect.",
  tolerance = 1e-10
)

check_equal(
  observed_kappas_p4,
  expected_kappas_p4,
  "Prior-4 kappa values are incorrect.",
  tolerance = 1e-10
)

# Priors 1 and 2 do not have favorable-profile kappas.
for (pid in c(1L, 2L)) {
  err <- try(get_eta_z_kappas(pid), silent = TRUE)
  check_true(
    inherits(err, "try-error"),
    paste0("get_eta_z_kappas() should reject prior_id = ", pid)
  )
}

expected_centers_p3 <- list(
  constant = expected_kappas_p3[["constant"]] * expected_shapes$constant,
  increasing = expected_kappas_p3[["increasing"]] * expected_shapes$increasing,
  waning = expected_kappas_p3[["waning"]] * expected_shapes$waning
)

expected_centers_p4 <- list(
  constant = expected_kappas_p4[["constant"]] * expected_shapes$constant,
  increasing = expected_kappas_p4[["increasing"]] * expected_shapes$increasing,
  waning = expected_kappas_p4[["waning"]] * expected_shapes$waning
)

observed_centers_p3 <- get_eta_z_centers(prior_id = 3L)
observed_centers_p4 <- get_eta_z_centers(prior_id = 4L)

for (pp in benefit_profiles) {
  check_equal(
    observed_centers_p3[[pp]],
    expected_centers_p3[[pp]],
    paste0("Prior-3 location vector is incorrect for ", pp),
    tolerance = 1e-10
  )
  
  check_equal(
    observed_centers_p4[[pp]],
    expected_centers_p4[[pp]],
    paste0("Prior-4 location vector is incorrect for ", pp),
    tolerance = 1e-10
  )
  
  check_equal(
    mean(observed_centers_p3[[pp]]),
    expected_kappas_p3[[pp]],
    paste0("Prior-3 mean location eta does not equal kappa for ", pp),
    tolerance = 1e-10
  )
  
  check_equal(
    mean(observed_centers_p4[[pp]]),
    expected_kappas_p4[[pp]],
    paste0("Prior-4 mean location eta does not equal kappa for ", pp),
    tolerance = 1e-10
  )
}

cat("[PASS] Profile shapes, kappas, and location vectors are correct.\n")


# ==============================================================================
# 4. CALIBRATION REPORT CONTENTS
# ==============================================================================

cat("\n")
cat("------------------------------------------------------------\n")
cat("[4] CHECKING CALIBRATION REPORT\n")
cat("------------------------------------------------------------\n")

calibration_report <- get_eta_z_calibration_report(delta_star = 0)

check_true(
  nrow(calibration_report) == 6L,
  "Calibration report should contain exactly six rows: 2 priors x 3 profiles."
)

check_true(
  setequal(unique(calibration_report$prior_id), c(3L, 4L)),
  "Calibration report should contain only Priors 3 and 4."
)

check_true(
  setequal(unique(calibration_report$profile), benefit_profiles),
  "Calibration report has incorrect profile labels."
)

check_true(
  all(calibration_report$achieved_target_delta_rmst == 2),
  "Calibration report target RMST values are not all 2.00."
)

check_equal(
  unique(calibration_report$sigma_star[calibration_report$prior_id == 3L]),
  0.10,
  "Prior-3 sigma_star in calibration report is incorrect."
)

check_equal(
  unique(calibration_report$sigma_star[calibration_report$prior_id == 4L]),
  0.30,
  "Prior-4 sigma_star in calibration report is incorrect."
)

expected_location_rmst_p3 <- c(
  constant = 2.0213294007,
  increasing = 2.0046265329,
  waning = 2.0102601115
)

expected_location_rmst_p4 <- c(
  constant = 1.0460907852,
  increasing = 1.0171928555,
  waning = 1.0425245934
)

for (pp in benefit_profiles) {
  obs3 <- calibration_report$center_delta_rmst[
    calibration_report$prior_id == 3L & calibration_report$profile == pp
  ]
  
  obs4 <- calibration_report$center_delta_rmst[
    calibration_report$prior_id == 4L & calibration_report$profile == pp
  ]
  
  check_equal(
    obs3,
    expected_location_rmst_p3[[pp]],
    paste0("Prior-3 location RMST is incorrect for ", pp),
    tolerance = 1e-9
  )
  
  check_equal(
    obs4,
    expected_location_rmst_p4[[pp]],
    paste0("Prior-4 location RMST is incorrect for ", pp),
    tolerance = 1e-9
  )
}

cat("[PASS] Calibration report is consistent with final Prior 3/4 definitions.\n")


# ==============================================================================
# 5. DIRECT DISTRIBUTION-LEVEL RMST CALIBRATION FOR PRIORS 3 AND 4
# ==============================================================================

cat("\n")
cat("------------------------------------------------------------\n")
cat("[5] DIRECTLY VERIFYING E_eta[Delta_RMST(24)] = 2\n")
cat("------------------------------------------------------------\n")

breaks_cal <- c(0, 4.8, 9.6, 14.4, 19.2, 24)
lambda0_cal <- c(0.0734, 0.0661, 0.0587, 0.0514, 0.0440)
eta_x_cal <- 0.20
tau_cal <- 24
target_delta_rmst_cal <- 2.00

# Use the same deterministic integration construction used for kappa calibration.
n_joint_calibration_check <- 200000L
seed_calibration_check <- 870001L

x_std_cal <- qnorm(
  (seq_len(n_joint_calibration_check) - 0.5) /
    n_joint_calibration_check
)

u_base_check <-
  (seq_len(n_joint_calibration_check) - 0.5) /
  n_joint_calibration_check

set.seed(seed_calibration_check)

u_eta_check <- matrix(
  NA_real_,
  nrow = n_joint_calibration_check,
  ncol = 5L
)

for (jj in seq_len(5L)) {
  u_eta_check[, jj] <- sample(
    u_base_check,
    size = n_joint_calibration_check,
    replace = FALSE
  )
}

qtnorm_positive_from_u_check <- function(u, mu, sigma) {
  p_lower <- pnorm(0, mean = mu, sd = sigma)
  p <- p_lower + u * (1 - p_lower)
  
  p <- pmin(
    pmax(p, .Machine$double.eps),
    1 - .Machine$double.eps
  )
  
  qnorm(p, mean = mu, sd = sigma)
}

rmst_piecewise_rows_check <- function(x, eta_z = NULL, breaks, lambda0, eta_x, tau) {
  x <- as.numeric(x)
  n <- length(x)
  J <- length(lambda0)
  
  check_true(
    length(breaks) == J + 1L,
    "Calibration check has incompatible breaks/lambda0."
  )
  
  if (is.null(eta_z)) {
    eta_z <- matrix(0, nrow = n, ncol = J)
  } else {
    eta_z <- as.matrix(eta_z)
    
    check_true(
      nrow(eta_z) == n && ncol(eta_z) == J,
      "Calibration eta matrix has incorrect dimensions."
    )
  }
  
  interval_width <- pmax(
    0,
    pmin(tau, breaks[-1]) - breaks[-length(breaks)]
  )
  
  rmst <- numeric(n)
  survival_start <- rep(1, n)
  
  for (jj in seq_len(J)) {
    width_j <- interval_width[jj]
    
    if (width_j <= 0) next
    
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

rmst_control_check <- rmst_piecewise_rows_check(
  x = x_std_cal,
  eta_z = NULL,
  breaks = breaks_cal,
  lambda0 = lambda0_cal,
  eta_x = eta_x_cal,
  tau = tau_cal
)

# First validate the row-wise implementation against the existing truth function.
eta_validation <- rep(0.2850356, 5L)
eta_validation_matrix <- matrix(
  eta_validation,
  nrow = n_joint_calibration_check,
  ncol = 5L,
  byrow = TRUE
)

delta_validation_rows <- mean(
  rmst_piecewise_rows_check(
    x = x_std_cal,
    eta_z = eta_validation_matrix,
    breaks = breaks_cal,
    lambda0 = lambda0_cal,
    eta_x = eta_x_cal,
    tau = tau_cal
  ) - rmst_control_check
)

delta_validation_truth <- true_delta_rmst_from_eta(
  eta_z = eta_validation,
  breaks = breaks_cal,
  lambda0 = lambda0_cal,
  eta_x = eta_x_cal,
  tau = tau_cal,
  x_std = x_std_cal
)

check_equal(
  delta_validation_rows,
  delta_validation_truth,
  "Row-wise RMST implementation does not match true_delta_rmst_from_eta().",
  tolerance = 1e-8
)

calibration_check_rows <- list()
row_id <- 1L

for (pid in c(3L, 4L)) {
  sigma_star <- if (pid == 3L) 0.10 else 0.30
  centers_current <- get_eta_z_centers(pid)
  kappas_current <- get_eta_z_kappas(pid)
  
  for (profile_run in benefit_profiles) {
    eta_location <- centers_current[[profile_run]]
    
    eta_draws <- matrix(
      NA_real_,
      nrow = n_joint_calibration_check,
      ncol = 5L
    )
    
    for (jj in seq_len(5L)) {
      eta_draws[, jj] <- qtnorm_positive_from_u_check(
        u = u_eta_check[, jj],
        mu = eta_location[jj],
        sigma = sigma_star
      )
    }
    
    mean_delta <- mean(
      rmst_piecewise_rows_check(
        x = x_std_cal,
        eta_z = eta_draws,
        breaks = breaks_cal,
        lambda0 = lambda0_cal,
        eta_x = eta_x_cal,
        tau = tau_cal
      ) - rmst_control_check
    )
    
    location_delta <- true_delta_rmst_from_eta(
      eta_z = eta_location,
      breaks = breaks_cal,
      lambda0 = lambda0_cal,
      eta_x = eta_x_cal,
      tau = tau_cal,
      x_std = x_std_cal
    )
    
    check_equal(
      mean_delta,
      target_delta_rmst_cal,
      paste0(
        "Prior-", pid,
        " distribution-level RMST calibration failed for ",
        profile_run
      ),
      tolerance = 1e-6
    )
    
    calibration_check_rows[[row_id]] <- data.frame(
      prior_id = pid,
      profile = profile_run,
      sigma_star = sigma_star,
      kappa = kappas_current[[profile_run]],
      achieved_mean_delta_RMST = mean_delta,
      delta_RMST_at_location = location_delta,
      target_delta_RMST = target_delta_rmst_cal,
      stringsAsFactors = FALSE
    )
    
    row_id <- row_id + 1L
  }
}

calibration_check <- do.call(rbind, calibration_check_rows)

cat("\nDistribution-level RMST calibration check:\n")
print(calibration_check, row.names = FALSE, digits = 10)

cat("\n[PASS] Prior-3 and Prior-4 distribution-level RMST calibrations were verified directly.\n")


# ==============================================================================
# 6. DESIGN-PRIOR DRAW CHECKS
# ==============================================================================

cat("\n")
cat("------------------------------------------------------------\n")
cat("[6] CHECKING DESIGN-PRIOR DRAWS\n")
cat("------------------------------------------------------------\n")

B_prior_check <- 2000L
prior_check_rows <- list()
row_id <- 1L

# ------------------------------------------------------------------------------
# Prior 1: exact-null point mass
# ------------------------------------------------------------------------------

set.seed(810001)

draws_p1 <- replicate(
  100L,
  draw_eta_z_from_prior(prior_id = 1L, profile = NULL),
  simplify = FALSE
)

eta_p1 <- do.call(rbind, lapply(draws_p1, function(x) x$eta_z))

check_true(all(eta_p1 == 0), "Prior 1 did not generate eta_z = 0.")
check_true(
  all(vapply(draws_p1, function(x) x$profile, character(1)) == "none"),
  "Prior 1 should store profile = 'none'."
)
check_true(
  all(vapply(draws_p1, function(x) x$component_id, integer(1)) == 0L),
  "Prior 1 should have component_id = 0."
)
check_true(
  all(vapply(draws_p1, function(x) x$omega_D, numeric(1)) == 0),
  "Prior 1 should have omega_D = 0."
)
check_true(
  all(vapply(draws_p1, function(x) x$prior_family, character(1)) == "point_mass"),
  "Prior 1 should be labeled point_mass."
)

prior_check_rows[[row_id]] <- data.frame(
  prior_id = 1L,
  design_prior = "exact_null_point_mass",
  profile = "none",
  mean_eta_1 = 0,
  sd_eta_1 = 0,
  max_abs_corr = NA_real_,
  stringsAsFactors = FALSE
)
row_id <- row_id + 1L

# ------------------------------------------------------------------------------
# Prior 2: TN(0, 0.30^2; -0.25, 0), independently across intervals
# ------------------------------------------------------------------------------

set.seed(820001)

draws_p2 <- replicate(
  B_prior_check,
  draw_eta_z_from_prior(prior_id = 2L, profile = NULL),
  simplify = FALSE
)

eta_p2 <- do.call(rbind, lapply(draws_p2, function(x) x$eta_z))

check_between(
  eta_p2,
  -0.25,
  0,
  "Prior 2 generated eta outside [-0.25, 0]."
)

check_true(
  all(apply(eta_p2, 2, sd) > 0),
  "Prior 2 contains a degenerate interval coefficient."
)

check_true(
  all(vapply(draws_p2, function(x) x$profile, character(1)) == "none"),
  "Prior 2 should store profile = 'none'."
)

check_true(
  all(vapply(draws_p2, function(x) x$component_id, integer(1)) == 0L),
  "Prior 2 should have component_id = 0."
)

cor_p2 <- cor(eta_p2)
max_abs_corr_p2 <- max(abs(cor_p2[upper.tri(cor_p2)]))

check_true(
  max_abs_corr_p2 < 0.15,
  "Prior 2 shows unexpectedly strong interval dependence."
)

# Empirical mean should agree with the intended truncated-normal marginal.
alpha_p2 <- (-0.25 - 0) / 0.30
beta_p2 <- (0 - 0) / 0.30
mean_p2_exact <- 0 + 0.30 *
  (dnorm(alpha_p2) - dnorm(beta_p2)) /
  (pnorm(beta_p2) - pnorm(alpha_p2))

check_true(
  max(abs(colMeans(eta_p2) - mean_p2_exact)) < 0.015,
  "Prior-2 empirical means are inconsistent with TN(0,0.30^2;-0.25,0)."
)

prior_check_rows[[row_id]] <- data.frame(
  prior_id = 2L,
  design_prior = "no_benefit_continuous",
  profile = "none",
  mean_eta_1 = mean(eta_p2[, 1]),
  sd_eta_1 = sd(eta_p2[, 1]),
  max_abs_corr = max_abs_corr_p2,
  stringsAsFactors = FALSE
)
row_id <- row_id + 1L

# ------------------------------------------------------------------------------
# Priors 3 and 4: positive truncated normal, independently across intervals
# ------------------------------------------------------------------------------

positive_prior_info <- data.frame(
  prior_id = c(3L, 4L),
  sigma_star = c(0.10, 0.30),
  design_prior = c(
    "benefit_positive_truncnorm_sd010",
    "benefit_positive_truncnorm_sd030"
  ),
  stringsAsFactors = FALSE
)

for (kk in seq_len(nrow(positive_prior_info))) {
  pid <- positive_prior_info$prior_id[kk]
  sigma_star <- positive_prior_info$sigma_star[kk]
  design_prior <- positive_prior_info$design_prior[kk]
  centers_current <- get_eta_z_centers(pid)
  
  for (profile_run in benefit_profiles) {
    eta_location <- centers_current[[profile_run]]
    
    set.seed(
      830000 +
        pid * 10000 +
        match(profile_run, benefit_profiles) * 1000
    )
    
    draws <- replicate(
      B_prior_check,
      draw_eta_z_from_prior(
        prior_id = pid,
        profile = profile_run
      ),
      simplify = FALSE
    )
    
    eta_mat <- do.call(rbind, lapply(draws, function(x) x$eta_z))
    
    check_true(
      all(is.finite(eta_mat)),
      paste0("Prior ", pid, " generated non-finite eta for ", profile_run)
    )
    
    check_true(
      all(eta_mat >= 0),
      paste0("Prior ", pid, " generated negative eta for ", profile_run)
    )
    
    check_true(
      all(apply(eta_mat, 2, sd) > 0),
      paste0("Prior ", pid, " contains a degenerate interval for ", profile_run)
    )
    
    check_true(
      all(vapply(draws, function(x) x$component_id, integer(1)) == 1L),
      paste0("Prior ", pid, " should have component_id = 1 for ", profile_run)
    )
    
    check_true(
      all(vapply(draws, function(x) x$omega_D, numeric(1)) == 1),
      paste0("Prior ", pid, " should have omega_D = 1 for ", profile_run)
    )
    
    check_true(
      all(vapply(draws, function(x) x$profile, character(1)) == profile_run),
      paste0("Prior ", pid, " stored incorrect profile for ", profile_run)
    )
    
    cor_current <- cor(eta_mat)
    max_abs_corr <- max(abs(cor_current[upper.tri(cor_current)]))
    
    check_true(
      max_abs_corr < 0.15,
      paste0(
        "Prior ", pid,
        " shows unexpectedly strong interval dependence for ",
        profile_run
      )
    )
    
    analytic_mean <- vapply(
      eta_location,
      function(mu_j) {
        alpha_j <- -mu_j / sigma_star
        mu_j + sigma_star * dnorm(alpha_j) / pnorm(alpha_j, lower.tail = FALSE)
      },
      numeric(1)
    )
    
    mean_tolerance <- if (pid == 3L) 0.015 else 0.030
    
    check_true(
      max(abs(colMeans(eta_mat) - analytic_mean)) < mean_tolerance,
      paste0(
        "Prior-", pid,
        " empirical eta means are inconsistent with the intended positive ",
        "truncated-normal marginals for ",
        profile_run
      )
    )
    
    prior_check_rows[[row_id]] <- data.frame(
      prior_id = pid,
      design_prior = design_prior,
      profile = profile_run,
      mean_eta_1 = mean(eta_mat[, 1]),
      sd_eta_1 = sd(eta_mat[, 1]),
      max_abs_corr = max_abs_corr,
      stringsAsFactors = FALSE
    )
    
    row_id <- row_id + 1L
  }
}

prior_check <- do.call(rbind, prior_check_rows)

cat("\nDesign-prior draw check:\n")
print(prior_check, row.names = FALSE)

cat("\n[PASS] All four final design-prior constructions passed.\n")


# ==============================================================================
# 7. generate_data() FOR ALL EIGHT DGP EVALUATION SETTINGS
# ==============================================================================

cat("\n")
cat("------------------------------------------------------------\n")
cat("[7] CHECKING generate_data() FOR ALL EIGHT SETTINGS\n")
cat("------------------------------------------------------------\n")

simulation_settings <- data.frame(
  setting_id = 1:8,
  prior_id = c(
    1L,
    2L,
    3L, 3L, 3L,
    4L, 4L, 4L
  ),
  profile = c(
    NA_character_,
    NA_character_,
    "constant", "increasing", "waning",
    "constant", "increasing", "waning"
  ),
  stringsAsFactors = FALSE
)

tau <- 24
delta_star <- 0
K2 <- 30L
K1 <- floor(0.6 * K2)
target_censor <- 0.40
N <- ceiling(K2 / (1 - target_censor))

breaks <- c(0, 4.8, 9.6, 14.4, 19.2, 24)
lambda0 <- c(0.0734, 0.0661, 0.0587, 0.0514, 0.0440)
eta_x <- 0.20

data_check_rows <- list()

for (ii in seq_len(nrow(simulation_settings))) {
  prior_id_run <- simulation_settings$prior_id[ii]
  profile_value <- simulation_settings$profile[ii]
  
  profile_run <- if (is.na(profile_value)) {
    NULL
  } else {
    as.character(profile_value)
  }
  
  set.seed(910000 + ii * 100)
  
  dat_one <- generate_data(
    prior_id = prior_id_run,
    profile = profile_run,
    N = N,
    K1 = K1,
    K2 = K2,
    target_censor = target_censor,
    breaks = breaks,
    lambda0 = lambda0,
    eta_x = eta_x,
    tau = tau,
    delta_star = delta_star,
    truth_n = 10000L,
    return_truth = TRUE,
    verbose = FALSE
  )
  
  expected_spec <- get_design_prior_spec(prior_id_run)
  expected_profile <- if (prior_id_run %in% c(3L, 4L)) profile_run else "none"
  
  check_equal(
    dat_one$prior_id,
    prior_id_run,
    paste0("dat_one$prior_id is incorrect for setting ", ii)
  )
  
  check_equal(
    dat_one$scenario_id,
    prior_id_run,
    paste0("dat_one$scenario_id is incorrect for setting ", ii)
  )
  
  check_true(
    identical(
      as.character(dat_one$design_prior_type),
      as.character(expected_spec$prior_label[[1]])
    ),
    paste0("design_prior_type is incorrect for setting ", ii)
  )
  
  check_equal(
    dat_one$omega_D,
    expected_spec$omega_D[[1]],
    paste0("omega_D is incorrect for setting ", ii)
  )
  
  check_equal(
    dat_one$H_D,
    expected_spec$H_D[[1]],
    paste0("H_D is incorrect for setting ", ii)
  )
  
  check_true(
    identical(as.character(dat_one$profile), as.character(expected_profile)),
    paste0("profile is incorrect for setting ", ii)
  )
  
  check_true(
    nrow(dat_one$O1) == N && nrow(dat_one$O2) == N,
    paste0("O1/O2 do not contain all randomized participants for setting ", ii)
  )
  
  check_true(
    sum(dat_one$O1$status) == K1,
    paste0("Interim event count does not equal K1 for setting ", ii)
  )
  
  check_true(
    sum(dat_one$O2$status) >= K2,
    paste0("Final event count is below K2 for setting ", ii)
  )
  
  check_true(
    dat_one$y_final >= tau,
    paste0("Final administrative cutoff is below tau for setting ", ii)
  )
  
  check_true(
    length(dat_one$eta_z) == 5L && all(is.finite(dat_one$eta_z)),
    paste0("eta_z is invalid for setting ", ii)
  )
  
  if (prior_id_run == 1L) {
    check_true(all(dat_one$eta_z == 0), "Prior 1 did not generate eta_z = 0.")
    check_equal(
      dat_one$truth$delta_rmst_tau,
      0,
      "Prior-1 exact-null RMST truth is not zero.",
      tolerance = 1e-12
    )
  }
  
  if (prior_id_run == 2L) {
    check_between(
      dat_one$eta_z,
      -0.25,
      0,
      "Prior-2 eta_z lies outside [-0.25, 0]."
    )
    
    check_true(
      dat_one$truth$delta_rmst_tau <= delta_star,
      "Prior 2 generated a positive RMST truth."
    )
  }
  
  if (prior_id_run %in% c(3L, 4L)) {
    check_true(
      all(dat_one$eta_z >= 0),
      paste0("Prior ", prior_id_run, " generated negative eta_z.")
    )
    
    check_true(
      dat_one$truth$delta_rmst_tau > delta_star,
      paste0("Prior ", prior_id_run, " generated a non-positive RMST truth.")
    )
  }
  
  check_equal(
    dat_one$truth$H_D,
    dat_one$component_id,
    paste0("truth$H_D does not equal component_id for setting ", ii)
  )
  
  check_equal(
    dat_one$truth$H_D_rmst_check,
    dat_one$component_id,
    paste0("RMST truth classification does not equal component_id for setting ", ii)
  )
  
  data_check_rows[[ii]] <- data.frame(
    setting = ii,
    prior_id = prior_id_run,
    design_prior_type = dat_one$design_prior_type,
    profile = dat_one$profile,
    omega_D = dat_one$omega_D,
    H_D = dat_one$H_D,
    eta_1 = dat_one$eta_z[1],
    eta_2 = dat_one$eta_z[2],
    eta_3 = dat_one$eta_z[3],
    eta_4 = dat_one$eta_z[4],
    eta_5 = dat_one$eta_z[5],
    delta_RMST = dat_one$truth$delta_rmst_tau,
    delta_median = dat_one$truth$delta_med_rep,
    events_K1 = sum(dat_one$O1$status),
    events_final = sum(dat_one$O2$status),
    final_cutoff = dat_one$y_final,
    stringsAsFactors = FALSE
  )
}

data_check <- do.call(rbind, data_check_rows)

cat("\nData-generation check:\n")
print(data_check, row.names = FALSE)

cat("\n[PASS] generate_data() passed for all eight DGP evaluation settings.\n")


# ==============================================================================
# 8. COMPILE BOTH ANALYSIS MODELS
# ==============================================================================

cat("\n")
cat("------------------------------------------------------------\n")
cat("[8] COMPILING BOTH ANALYSIS MODELS\n")
cat("------------------------------------------------------------\n")

mod_gp <- compile_gp_model()

check_true(
  file.exists(mod_gp$exe_file()),
  "GP-PH model executable was not created."
)

cat("[PASS] GP-PH model compiled.\n")

mod_single_gp <- compile_independent_gp_model()

check_true(
  file.exists(mod_single_gp$exe_file()),
  "Independent GP model executable was not created."
)

cat("[PASS] Independent GP model compiled.\n")


# ==============================================================================
# 9. SHORT GP-PH + INDEPENDENT GP MCMC SMOKE TEST FOR ALL FOUR PRIORS
# ==============================================================================

cat("\n")
cat("------------------------------------------------------------\n")
cat("[9] SHORT GP-PH + INDEPENDENT GP MCMC SMOKE TEST\n")
cat("------------------------------------------------------------\n")

# Representative positive-benefit setting: constant profile for Priors 3 and 4.
# Increasing/waning have already been tested through generate_data() above.
mcmc_settings <- data.frame(
  prior_id = 1:4,
  profile = c(
    NA_character_,
    NA_character_,
    "constant",
    "constant"
  ),
  stringsAsFactors = FALSE
)

a0_gp <- 2
grid_width_gp <- 1
B0_median_pfs <- 10
B0_rate <- log(2) / B0_median_pfs
B0_fun <- function(t) B0_rate * pmax(as.numeric(t), 0)

# Execution-only MCMC settings. These are intentionally short and are NOT
# suitable for convergence assessment or operating-characteristic estimation.
chains_run <- 2L
parallel_chains_run <- 1L
iter_warmup_run <- 200L
iter_sampling_run <- 100L
thin_run <- 1L
refresh_run <- 100L

adapt_delta_gp <- 0.95
max_treedepth_gp <- 12L

adapt_delta_indep <- 0.99
max_treedepth_indep <- 14L

required_common_output <- c(
  "prior_id",
  "scenario_id",
  "prior_label",
  "design_prior_type",
  "prior_family",
  "omega_D",
  "profile",
  "component",
  "component_id",
  "H_D",
  "true_delta_rmst_rep",
  "true_delta_med_rep",
  "Pi1",
  "Pi2_cf",
  "post_mean_1",
  "post_q025_1",
  "post_q975_1",
  "post_mean_2",
  "post_q025_2",
  "post_q975_2",
  "Pi_beneficial_med_1",
  "Pi_beneficial_med_2",
  "max_Rhat_1",
  "min_ESS_1",
  "n_divergent_1",
  "max_Rhat_2",
  "min_ESS_2",
  "n_divergent_2"
)

required_eta_output <- c(
  "post_eta_PH_trt_mean_1",
  "post_eta_PH_trt_q025_1",
  "post_eta_PH_trt_q500_1",
  "post_eta_PH_trt_q975_1",
  "post_eta_PH_trt_mean_2",
  "post_eta_PH_trt_q025_2",
  "post_eta_PH_trt_q500_2",
  "post_eta_PH_trt_q975_2"
)

mcmc_rows <- list()

for (ii in seq_len(nrow(mcmc_settings))) {
  prior_id_run <- mcmc_settings$prior_id[ii]
  profile_value <- mcmc_settings$profile[ii]
  
  profile_run <- if (is.na(profile_value)) {
    NULL
  } else {
    as.character(profile_value)
  }
  
  cat("\n")
  cat("============================================================\n")
  cat(
    "MCMC smoke: prior_id = ", prior_id_run,
    ", profile = ", if (is.null(profile_run)) "none" else profile_run,
    "\n",
    sep = ""
  )
  cat("============================================================\n")
  
  set.seed(101000 + prior_id_run)
  
  dat_one <- generate_data(
    prior_id = prior_id_run,
    profile = profile_run,
    N = N,
    K1 = K1,
    K2 = K2,
    target_censor = target_censor,
    breaks = breaks,
    lambda0 = lambda0,
    eta_x = eta_x,
    tau = tau,
    delta_star = delta_star,
    truth_n = 10000L,
    return_truth = TRUE,
    verbose = TRUE
  )
  
  # --------------------------------------------------------------------------
  # GP-PH
  # --------------------------------------------------------------------------
  
  mcmc_gp <- BDCT_MCMC_ph(
    dat = dat_one,
    rep_id = 1L,
    mod_gp = mod_gp,
    tau = tau,
    delta_star = delta_star,
    B0_fun = B0_fun,
    B0_median_pfs = B0_median_pfs,
    a0_gp = a0_gp,
    grid_width_gp = grid_width_gp,
    chains = chains_run,
    parallel_chains = parallel_chains_run,
    iter_warmup = iter_warmup_run,
    iter_sampling = iter_sampling_run,
    refresh = refresh_run,
    thin = thin_run,
    adapt_delta = adapt_delta_gp,
    max_treedepth = max_treedepth_gp,
    cmdstan_seed = 201000L + prior_id_run
  )
  
  check_true(
    inherits(mcmc_gp, "BDCT_MCMC"),
    paste0("GP-PH MCMC object has incorrect class for prior ", prior_id_run)
  )
  
  post_gp <- BDCT_POST_ph(
    object = mcmc_gp,
    truth_lookup = NULL,
    tau = tau,
    delta_star = delta_star
  )
  
  missing_gp <- setdiff(
    c(required_common_output, required_eta_output),
    names(post_gp)
  )
  
  check_true(
    length(missing_gp) == 0L,
    paste0(
      "GP-PH output is missing columns: ",
      paste(missing_gp, collapse = ", ")
    )
  )
  
  check_between(
    post_gp$Pi1,
    0,
    1,
    paste0("GP-PH Pi1 is invalid for prior ", prior_id_run)
  )
  
  check_between(
    post_gp$Pi2_cf,
    0,
    1,
    paste0("GP-PH Pi2_cf is invalid for prior ", prior_id_run)
  )
  
  check_true(
    is.finite(post_gp$post_eta_PH_trt_mean_1) &&
      is.finite(post_gp$post_eta_PH_trt_mean_2),
    paste0("GP-PH treatment eta is not finite for prior ", prior_id_run)
  )
  
  # --------------------------------------------------------------------------
  # Independent GP
  # --------------------------------------------------------------------------
  
  mcmc_indep <- BDCT_MCMC_indep(
    dat = dat_one,
    rep_id = 1L,
    mod_single_gp = mod_single_gp,
    tau = tau,
    delta_star = delta_star,
    B0_fun = B0_fun,
    B0_median_pfs = B0_median_pfs,
    a0_gp = a0_gp,
    grid_width_gp = grid_width_gp,
    chains = chains_run,
    parallel_chains = parallel_chains_run,
    iter_warmup = iter_warmup_run,
    iter_sampling = iter_sampling_run,
    refresh = refresh_run,
    thin = thin_run,
    adapt_delta = adapt_delta_indep,
    max_treedepth = max_treedepth_indep,
    cmdstan_seed = 301000L + prior_id_run
  )
  
  check_true(
    inherits(mcmc_indep, "BDCT_MCMC_indep"),
    paste0("Independent GP MCMC object has incorrect class for prior ", prior_id_run)
  )
  
  post_indep <- BDCT_POST_indep(
    object = mcmc_indep,
    truth_lookup = NULL,
    tau = tau,
    delta_star = delta_star
  )
  
  missing_indep <- setdiff(
    c(required_common_output, required_eta_output),
    names(post_indep)
  )
  
  check_true(
    length(missing_indep) == 0L,
    paste0(
      "Independent GP output is missing columns: ",
      paste(missing_indep, collapse = ", ")
    )
  )
  
  check_between(
    post_indep$Pi1,
    0,
    1,
    paste0("Independent GP Pi1 is invalid for prior ", prior_id_run)
  )
  
  check_between(
    post_indep$Pi2_cf,
    0,
    1,
    paste0("Independent GP Pi2_cf is invalid for prior ", prior_id_run)
  )
  
  check_true(
    is.na(post_indep$post_eta_PH_trt_mean_1) &&
      is.na(post_indep$post_eta_PH_trt_mean_2),
    "Independent GP PH treatment eta columns should be NA."
  )
  
  # Old coefficient-output names must be absent.
  old_output_names_gp <- grep(
    "post_(beta|gamma)_PH",
    names(post_gp),
    value = TRUE
  )
  
  old_output_names_indep <- grep(
    "post_(beta|gamma)_PH",
    names(post_indep),
    value = TRUE
  )
  
  check_true(
    length(old_output_names_gp) == 0L,
    paste0(
      "Old GP-PH coefficient output names remain: ",
      paste(old_output_names_gp, collapse = ", ")
    )
  )
  
  check_true(
    length(old_output_names_indep) == 0L,
    paste0(
      "Old Independent-GP coefficient output names remain: ",
      paste(old_output_names_indep, collapse = ", ")
    )
  )
  
  expected_spec <- get_design_prior_spec(prior_id_run)
  
  for (obj_name in c("post_gp", "post_indep")) {
    obj <- get(obj_name)
    
    check_equal(
      obj$prior_id,
      prior_id_run,
      paste0(obj_name, " prior_id is incorrect.")
    )
    
    check_equal(
      obj$scenario_id,
      prior_id_run,
      paste0(obj_name, " scenario_id is incorrect.")
    )
    
    check_equal(
      obj$omega_D,
      expected_spec$omega_D[[1]],
      paste0(obj_name, " omega_D is incorrect.")
    )
    
    check_equal(
      obj$H_D,
      expected_spec$H_D[[1]],
      paste0(obj_name, " H_D is incorrect.")
    )
  }
  
  mcmc_rows[[ii]] <- data.frame(
    prior_id = prior_id_run,
    design_prior_type = post_gp$design_prior_type,
    profile = post_gp$profile,
    GP_Pi1 = post_gp$Pi1,
    GP_Pi2 = post_gp$Pi2_cf,
    Indep_Pi1 = post_indep$Pi1,
    Indep_Pi2 = post_indep$Pi2_cf,
    GP_max_Rhat_1 = post_gp$max_Rhat_1,
    GP_max_Rhat_2 = post_gp$max_Rhat_2,
    Indep_max_Rhat_1 = post_indep$max_Rhat_1,
    Indep_max_Rhat_2 = post_indep$max_Rhat_2,
    stringsAsFactors = FALSE
  )
  
  rm(dat_one, mcmc_gp, post_gp, mcmc_indep, post_indep)
  invisible(gc())
}

mcmc_check <- do.call(rbind, mcmc_rows)

cat("\nShort MCMC smoke-test summary:\n")
print(mcmc_check, row.names = FALSE)

cat(
  "\n[PASS] GP-PH and Independent GP completed for one representative trial ",
  "from each of the four design priors.\n",
  sep = ""
)


# ==============================================================================
# 10. THRESHOLD-GRID POST-PROCESSING ON SYNTHETIC DATA
# ==============================================================================

cat("\n")
cat("------------------------------------------------------------\n")
cat("[10] CHECKING THRESHOLD-GRID POST-PROCESSING\n")
cat("------------------------------------------------------------\n")

synthetic_post <- data.frame(
  rep = 1:4,
  prior_id = rep(3L, 4),
  scenario_id = rep(3L, 4),
  prior_label = rep("benefit_positive_truncnorm_sd010", 4),
  design_prior_type = rep("benefit_positive_truncnorm_sd010", 4),
  omega_D = rep(1, 4),
  profile = rep("constant", 4),
  H_D = rep(1L, 4),
  K1 = rep(18L, 4),
  K2 = rep(30L, 4),
  Pi1 = c(0.996, 0.980, 0.997, 0.800),
  Pi2_cf = c(0.999, 0.990, 0.900, 0.995),
  Pi_beneficial_med_1 = c(0.996, 0.970, 0.998, 0.700),
  Pi_beneficial_med_2 = c(0.999, 0.990, 0.850, 0.996),
  bias_rep_1 = c(0.1, -0.1, 0.2, -0.2),
  bias_rep_2 = c(0.05, -0.05, 0.10, -0.10),
  sqerr_rep_2 = c(0.0025, 0.0025, 0.01, 0.01),
  cover_rep_2 = c(1, 1, 1, 0),
  bias_med_rep_1 = c(0.2, -0.2, 0.3, -0.3),
  bias_med_rep_2 = c(0.1, -0.1, 0.2, -0.2),
  sqerr_med_rep_2 = c(0.01, 0.01, 0.04, 0.04),
  cover_med_rep_2 = c(1, 1, 1, 0),
  post_width_1 = rep(2, 4),
  post_width_2 = rep(1.5, 4),
  post_width_med_1 = rep(4, 4),
  post_width_med_2 = rep(3, 4),
  na_rate_med_1 = rep(0, 4),
  na_rate_med_2 = rep(0, 4),
  stringsAsFactors = FALSE
)

threshold_grid_test <- data.frame(
  nu1_E = c(0.995, 0.990),
  nu2_E = c(0.985, 0.980),
  stringsAsFactors = FALSE
)

synthetic_grid <- apply_threshold_grid(
  res_prob = synthetic_post,
  threshold_grid = threshold_grid_test
)

check_true(
  nrow(synthetic_grid) == nrow(synthetic_post) * nrow(threshold_grid_test),
  "apply_threshold_grid() returned an incorrect number of rows."
)

required_threshold_columns <- c(
  "E1_success",
  "Final_success",
  "Terminal_final_success",
  "S_R",
  "E1_success_med",
  "Final_success_med",
  "Terminal_final_success_med",
  "S_R_med",
  "nu1_E",
  "nu2_E"
)

check_true(
  all(required_threshold_columns %in% names(synthetic_grid)),
  paste0(
    "Threshold-grid output is missing columns: ",
    paste(
      setdiff(required_threshold_columns, names(synthetic_grid)),
      collapse = ", "
    )
  )
)

synthetic_oc <- summarize_oc_by_threshold(synthetic_grid)

check_true(
  nrow(synthetic_oc) == nrow(threshold_grid_test),
  "summarize_oc_by_threshold() returned an incorrect number of rows."
)

check_true(
  all(is.finite(synthetic_oc$BP)),
  "Synthetic benefit setting should produce finite BP values."
)

check_true(
  all(is.na(synthetic_oc$FPR)),
  "Synthetic benefit setting should have FPR = NA."
)

cat("\nSynthetic threshold-grid OC check:\n")
print(
  synthetic_oc[
    ,
    c(
      "prior_id",
      "profile",
      "nu1_E",
      "nu2_E",
      "PES",
      "PFS",
      "BP",
      "Interim_bias",
      "Final_bias_cf",
      "Final_RMSE_cf",
      "Final_coverage_cf"
    )
  ],
  row.names = FALSE
)

cat("\n[PASS] Threshold-grid post-processing functions work correctly.\n")


# ==============================================================================
# 11. DPPS / PPV / NPV IDENTITIES
# ==============================================================================

cat("\n")
cat("------------------------------------------------------------\n")
cat("[11] CHECKING DPPS / PPV / NPV IDENTITIES\n")
cat("------------------------------------------------------------\n")

BP_test <- 0.80
FPR_test <- 0.025
omega_test <- c(0.2, 0.4, 0.5, 0.6, 0.8)

predictive_test <- predictive_ocs_from_components(
  BP = BP_test,
  FPR = FPR_test,
  omega_D = omega_test
)

expected_DPPS <-
  omega_test * BP_test +
  (1 - omega_test) * FPR_test

expected_PPV <-
  omega_test * BP_test /
  expected_DPPS

expected_NPV <-
  (1 - omega_test) * (1 - FPR_test) /
  (
    (1 - omega_test) * (1 - FPR_test) +
      omega_test * (1 - BP_test)
  )

check_equal(
  predictive_test$DPPS,
  expected_DPPS,
  "DPPS identity is incorrect.",
  tolerance = 1e-12
)

check_equal(
  predictive_test$PPV,
  expected_PPV,
  "PPV identity is incorrect.",
  tolerance = 1e-12
)

check_equal(
  predictive_test$NPV,
  expected_NPV,
  "NPV identity is incorrect.",
  tolerance = 1e-12
)

cat("\nPredictive-OC check:\n")
print(predictive_test, row.names = FALSE)

cat("\n[PASS] DPPS/PPV/NPV post-processing identities are correct.\n")


# ==============================================================================
# 12. SCHOENFELD FORMULA
# ==============================================================================

cat("\n")
cat("------------------------------------------------------------\n")
cat("[12] SCHOENFELD FORMULA\n")
cat("------------------------------------------------------------\n")

calculate_events <- function(
    hr,
    alpha = 0.025,
    power = 0.80,
    one_sided = TRUE
) {
  if (one_sided) {
    z_alpha <- qnorm(1 - alpha)
  } else {
    z_alpha <- qnorm(1 - alpha / 2)
  }
  
  z_beta <- qnorm(power)
  
  d <-
    4 * (z_alpha + z_beta)^2 /
    log(hr)^2
  
  ceiling(d)
}

hr_target <- 0.75
alpha_one_sided <- 0.025

n_events_80 <- calculate_events(
  hr = hr_target,
  alpha = alpha_one_sided,
  power = 0.80,
  one_sided = TRUE
)

n_events_90 <- calculate_events(
  hr = hr_target,
  alpha = alpha_one_sided,
  power = 0.90,
  one_sided = TRUE
)

check_equal(
  n_events_80,
  380L,
  "Schoenfeld 80% event calculation no longer equals 380."
)

check_equal(
  n_events_90,
  508L,
  "Schoenfeld 90% event calculation no longer equals 508."
)

N_80 <- ceiling(n_events_80 / (1 - target_censor))
N_90 <- ceiling(n_events_90 / (1 - target_censor))

cat(
  sprintf(
    paste0(
      "One-sided HR %.3f, alpha %.3f, power 0.80 ",
      "-> required events = %d, approximate N at %.0f%% censoring = %d\n"
    ),
    hr_target,
    alpha_one_sided,
    n_events_80,
    100 * target_censor,
    N_80
  )
)

cat(
  sprintf(
    paste0(
      "One-sided HR %.3f, alpha %.3f, power 0.90 ",
      "-> required events = %d, approximate N at %.0f%% censoring = %d\n"
    ),
    hr_target,
    alpha_one_sided,
    n_events_90,
    100 * target_censor,
    N_90
  )
)

cat("[PASS] Schoenfeld calculations completed.\n")


# ==============================================================================
# FINAL STATUS
# ==============================================================================

cat("\n")
cat("============================================================\n")
cat("FINAL CHECK\n")
cat("============================================================\n")

cat("\nALL FOUR-PRIOR BDCT SMOKE TESTS COMPLETED SUCCESSFULLY.\n")

cat(
  "\nChecked:\n",
  "  - BDCT project/function naming\n",
  "  - no old BGSD project names in BDCT_r.R\n",
  "  - no obsolete uniform-prior machinery\n",
  "  - prior_id = 1, ..., 4\n",
  "  - Prior 1: exact-null point mass\n",
  "  - Prior 2: TN(0, 0.30^2; -0.25, 0) independently by interval\n",
  "  - Prior 3: positive TN with SD = 0.10\n",
  "  - Prior 4: positive TN with SD = 0.30\n",
  "  - common normalized temporal profile shapes\n",
  "  - Prior-3 and Prior-4 calibrated kappas\n",
  "  - Prior-3 and Prior-4 location vectors\n",
  "  - calibration-report contents\n",
  "  - direct Prior-3 E_eta[Delta_RMST(24)] = 2 calibration\n",
  "  - direct Prior-4 E_eta[Delta_RMST(24)] = 2 calibration\n",
  "  - independent interval-specific design-prior draws\n",
  "  - all eight DGP evaluation settings\n",
  "  - event-driven interim/final data\n",
  "  - RMST truth classification\n",
  "  - GP-PH Stan compilation\n",
  "  - Independent GP Stan compilation\n",
  "  - GP-PH short MCMC/post-processing for prior_id 1-4\n",
  "  - Independent GP short MCMC/post-processing for prior_id 1-4\n",
  "  - BDCT MCMC object classes\n",
  "  - post_eta_PH_* output naming\n",
  "  - threshold-grid post-processing\n",
  "  - DPPS/PPV/NPV identities\n",
  "  - Schoenfeld event anchors: 380 and 508\n",
  sep = ""
)

cat(
  "\nNOTE: The short MCMC settings above are for execution checking only. ",
  "Do not use their Rhat, ESS, posterior probabilities, bias, RMSE, coverage, ",
  "or other summaries as production operating-characteristic results.\n",
  sep = ""
)
