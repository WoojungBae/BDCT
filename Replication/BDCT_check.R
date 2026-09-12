# ==============================================================================
# BDCT profile-specific threshold calibration -- FPR-first selection
#
# FINAL PRIOR NUMBERING
#   Prior 1: exact-null point mass, eta_j = 0
#   Prior 2: continuous no-benefit prior, eta_j ~ TN(0, 0.30^2; -0.25, 0)
#   Prior 3: continuous positive-benefit prior,
#            eta_j ~ TN(kappa_s * c_{s,j}, 0.10^2; 0, Inf)
#   Prior 4: NOT INCLUDED because simulation results are not yet available.
#
# EXISTING RESULT FILES ARE REUSED. NO MCMC IS RERUN.
#
# Legacy saved IDs are accepted and remapped:
#   old 3 -> new 1
#   old 2 -> new 2
#   old 7 -> new 3
#
# THRESHOLD GRID
#   nu1_E, nu2_E in {0.950, 0.951, ..., 0.995}
#   nu1_E > nu2_E
#   46 threshold values -> 1,035 ordered threshold pairs
#
# FPR-FIRST PROFILE-SPECIFIC CALIBRATION
#
# For each Estimand x model x N:
#   1. Evaluate exact-null FPR for all threshold pairs.
#   2. Restrict to FPR <= 2.5%.
#   3. If exact-null FPR = 2.5% is attainable, use 2.5%.
#      Otherwise use the largest attainable FPR below 2.5% (for example, 2.4%).
#   4. Keep only threshold pairs attaining that selected target FPR.
#
# Then, separately for each benefit profile:
#   5. Among those target-FPR pairs, maximize profile-specific BP.
#   6. Break exact BP ties by the prespecified threshold-grid order.
#
# Thus FPR calibration is primary and profile-specific BP is secondary.
#
# Prior 2 is evaluated secondarily at the selected thresholds.
#
# OUTPUT
#   One compact standalone HTML report is written to Results/.../HTML.
#   Only essential checking/result tables are shown directly.
#   The full 24,840-row calibration surface is NOT embedded in the HTML.
# ==============================================================================

suppressPackageStartupMessages({
  library(dplyr); library(purrr); library(tidyr); library(tibble); library(knitr); library(htmltools)
})

# ==============================================================================
# 1. Settings
# ==============================================================================

results_dir_candidates <- c(
  "/Users/woojung/Documents/Rproject/BDCT/SimulationStudy/Results/ResultsThin1",
  "/Users/woojung/Documents/Rproject/BDCT/SimulationStudy/Results",
  "/Users/woojung/Documents/Rproject/BDCT/Results/ResultsThin1",
  "/Users/woojung/Documents/Rproject/BDCT/Results",
  "/home/woojung.bae/BDCT/Results_thin1",
  "/home/woojung.bae/BDCT/Results"
)

existing_results_dirs <- results_dir_candidates[dir.exists(results_dir_candidates)]
results_dir <- if (length(existing_results_dirs) > 0L) existing_results_dirs[1L] else results_dir_candidates[1L]

if (!dir.exists(results_dir)) {
  stop("Results directory was not found. Checked: ", paste(results_dir_candidates, collapse = "; "))
}

model_prefixes <- c("Original GP-PH" = "GP_POST_delta0", "Independent GP" = "IndepGP_POST_delta0")
N_keep <- c(700L, 1000L)
benefit_profiles <- c("constant", "increasing", "waning")
all_profile_levels <- c(benefit_profiles, "none")

exact_null_prior_id <- 1L
no_benefit_prior_id <- 2L
benefit_prior_id <- 3L

design_prior_specs <- tibble::tribble(
  ~prior_id, ~prior_label,                       ~prior_component, ~prior_family, ~omega_D, ~profile,     ~output_tag,                                  ~saved_prior_ids,
  1L,        "exact_null_point_mass",            "no_benefit",     "point_mass",  0.00,     "none",       "exact_null",                                  "3,1",
  2L,        "no_benefit_continuous",            "no_benefit",     "continuous",  0.00,     "none",       "no_benefit_continuous",                       "2",
  3L,        "benefit_positive_truncnorm_sd010", "benefit",        "continuous",  1.00,     "constant",   "benefit_positive_truncnorm_constant",         "7,3",
  3L,        "benefit_positive_truncnorm_sd010", "benefit",        "continuous",  1.00,     "increasing", "benefit_positive_truncnorm_increasing",       "7,3",
  3L,        "benefit_positive_truncnorm_sd010", "benefit",        "continuous",  1.00,     "waning",     "benefit_positive_truncnorm_waning",           "7,3"
)

n_rep_expected <- 1000L
expected_rep_ids <- seq_len(n_rep_expected)

rmst_interim_prob_col <- "Pi1"
rmst_final_prob_col <- "Pi2_cf"
median_interim_prob_col <- "Pi_beneficial_med_1"
median_final_prob_col <- "Pi_beneficial_med_2"
rmst_true_col <- "true_delta_rmst_rep"
median_true_col <- "true_delta_med_rep"

fpr_target_pct <- 2.5
numeric_tolerance <- 1e-10

# ==============================================================================
# 2. Threshold grid
# ==============================================================================

threshold_values <- round(seq(0.950, 0.995, by = 0.001), 3)

threshold_grid <- tidyr::expand_grid(nu1_E = threshold_values, nu2_E = threshold_values) %>%
  dplyr::filter(nu1_E > nu2_E) %>%
  dplyr::arrange(nu1_E, nu2_E) %>%
  dplyr::mutate(grid_order = dplyr::row_number())

expected_threshold_pairs <- choose(length(threshold_values), 2)

if (length(threshold_values) != 46L) stop("Threshold values should contain exactly 46 values.")
if (nrow(threshold_grid) != expected_threshold_pairs) {
  stop("Threshold grid should contain exactly ", expected_threshold_pairs, " ordered pairs.")
}

# ==============================================================================
# 3. Utility functions
# ==============================================================================

safe_mean <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) == 0L) return(NA_real_)
  mean(x)
}

safe_min <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) == 0L) return(NA_real_)
  min(x)
}

safe_max <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) == 0L) return(NA_real_)
  max(x)
}

safe_rmse <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) == 0L) return(NA_real_)
  sqrt(mean(x^2))
}

pct <- function(x) 100 * x

assert_columns <- function(dat, columns, context) {
  missing_columns <- setdiff(columns, names(dat))
  if (length(missing_columns) > 0L) {
    stop(context, " is missing columns: ", paste(missing_columns, collapse = ", "))
  }
  invisible(TRUE)
}

parse_saved_prior_ids <- function(x) {
  x <- unlist(x, recursive = TRUE, use.names = FALSE)
  x <- unlist(strsplit(as.character(x), ",", fixed = TRUE), use.names = FALSE)
  x <- trimws(x)
  out <- suppressWarnings(as.integer(x))
  out <- unique(out[is.finite(out)])
  if (length(out) == 0L) stop("Could not parse `saved_prior_ids`.")
  out
}

collapse_ids <- function(x, max_show = 30L) {
  x <- sort(unique(as.integer(x)))
  x <- x[is.finite(x)]
  if (length(x) == 0L) return("")
  shown <- head(x, max_show)
  suffix <- if (length(x) > max_show) paste0(" ... (+", length(x) - max_show, " more)") else ""
  paste0(paste(shown, collapse = ", "), suffix)
}

# ==============================================================================
# 4. Expected files
# ==============================================================================

expected_file_grid <- tidyr::crossing(model = names(model_prefixes), N = N_keep) %>%
  dplyr::mutate(prefix = unname(model_prefixes[model])) %>%
  tidyr::crossing(design_prior_specs) %>%
  dplyr::mutate(
    expected_file = paste0(prefix, "_", output_tag, "_N", N, ".txt"),
    expected_path = file.path(results_dir, expected_file)
  )

file_inventory <- expected_file_grid %>%
  dplyr::transmute(model, N, prior_id, profile, expected_file, file_found = file.exists(expected_path))

missing_files <- file_inventory %>% dplyr::filter(!file_found)

if (nrow(missing_files) > 0L) {
  print(missing_files, n = Inf)
  stop("At least one required Prior 1-3 result file is missing. Prior 4 files are intentionally not required.")
}

# ==============================================================================
# 5. Load existing Priors 1-3 and remap legacy IDs
# ==============================================================================

read_one_result_file <- function(
    file, model_name, N_expected, prior_id_final, prior_label_final, prior_component_final,
    prior_family_final, omega_D_final, profile_final, allowed_saved_prior_ids
) {
  allowed_saved_prior_ids <- parse_saved_prior_ids(allowed_saved_prior_ids)
  
  df <- read.table(
    file = file, header = TRUE, sep = "\t", fill = TRUE, quote = "", comment.char = "",
    check.names = FALSE, stringsAsFactors = FALSE
  )
  
  if (nrow(df) == 0L) stop("Result file is empty: ", basename(file))
  
  replicate_column <- if ("rep" %in% names(df)) {
    "rep"
  } else if ("run_ID" %in% names(df)) {
    "run_ID"
  } else {
    stop("Neither `rep` nor `run_ID` was found in: ", basename(file))
  }
  
  df <- df[as.character(df[[replicate_column]]) != replicate_column, , drop = FALSE]
  df[] <- lapply(df, function(x) type.convert(x, as.is = TRUE))
  df$rep_id <- suppressWarnings(as.integer(df[[replicate_column]]))
  df <- df[is.finite(df$rep_id), , drop = FALSE]
  
  if (nrow(df) == 0L) stop("No finite replicate rows remained in: ", basename(file))
  
  df$prior_id_saved <- if ("prior_id" %in% names(df)) suppressWarnings(as.integer(df$prior_id)) else NA_integer_
  
  bad_saved_prior_id <- !is.na(df$prior_id_saved) & !(df$prior_id_saved %in% allowed_saved_prior_ids)
  
  if (any(bad_saved_prior_id)) {
    stop(
      "Unexpected saved prior_id in ", basename(file),
      ". Observed: ", paste(sort(unique(df$prior_id_saved[bad_saved_prior_id])), collapse = ", "),
      "; allowed: ", paste(allowed_saved_prior_ids, collapse = ", ")
    )
  }
  
  df$model <- model_name
  df$N <- as.integer(N_expected)
  df$prior_id <- as.integer(prior_id_final)
  df$prior_label <- prior_label_final
  df$design_prior_type <- prior_label_final
  df$prior_component <- prior_component_final
  df$prior_family <- prior_family_final
  df$omega_D <- as.numeric(omega_D_final)
  df$profile <- profile_final
  df$source_file <- basename(file)
  df
}

dat_raw <- purrr::pmap_dfr(
  expected_file_grid,
  function(
    model, N, prefix, prior_id, prior_label, prior_component, prior_family, omega_D,
    profile, output_tag, saved_prior_ids, expected_file, expected_path
  ) {
    read_one_result_file(
      file = expected_path, model_name = model, N_expected = N, prior_id_final = prior_id,
      prior_label_final = prior_label, prior_component_final = prior_component,
      prior_family_final = prior_family, omega_D_final = omega_D, profile_final = profile,
      allowed_saved_prior_ids = saved_prior_ids
    )
  }
)

# ==============================================================================
# 6. Required columns and replicate checks
# ==============================================================================

required_columns <- c(
  "rep_id", "model", "N", "prior_id", "prior_label", "profile",
  rmst_interim_prob_col, rmst_final_prob_col, median_interim_prob_col, median_final_prob_col,
  rmst_true_col, median_true_col, "events_interim", "events_final",
  "post_mean_1", "post_q025_2", "post_q975_2", "post_mean_2",
  "post_mean_med_1", "post_q025_med_2", "post_q975_med_2", "post_mean_med_2"
)

assert_columns(dat_raw, required_columns, "Loaded simulation results")

duplicate_check <- dat_raw %>%
  dplyr::count(model, N, prior_id, profile, rep_id, name = "n_rows") %>%
  dplyr::filter(n_rows > 1L)

dat <- dat_raw %>%
  dplyr::group_by(model, N, prior_id, profile, rep_id) %>%
  dplyr::slice_tail(n = 1L) %>%
  dplyr::ungroup()

replicate_id_check <- expected_file_grid %>%
  dplyr::select(model, N, prior_id, profile) %>%
  dplyr::distinct() %>%
  dplyr::left_join(
    dat %>%
      dplyr::group_by(model, N, prior_id, profile) %>%
      dplyr::summarise(
        observed_ids = list(sort(unique(rep_id))),
        n_rep = dplyr::n_distinct(rep_id),
        .groups = "drop"
      ),
    by = c("model", "N", "prior_id", "profile")
  ) %>%
  dplyr::mutate(
    n_rep = tidyr::replace_na(n_rep, 0L),
    observed_ids = purrr::map(observed_ids, function(x) if (is.null(x)) integer(0) else as.integer(x)),
    missing_ids = purrr::map(observed_ids, function(x) setdiff(expected_rep_ids, x)),
    extra_ids = purrr::map(observed_ids, function(x) setdiff(x, expected_rep_ids)),
    n_missing = purrr::map_int(missing_ids, length),
    n_extra = purrr::map_int(extra_ids, length),
    Complete = n_rep == n_rep_expected & n_missing == 0L & n_extra == 0L,
    missing_rep_ids = purrr::map_chr(missing_ids, collapse_ids),
    extra_rep_ids = purrr::map_chr(extra_ids, collapse_ids)
  ) %>%
  dplyr::arrange(model, N, prior_id, factor(profile, levels = all_profile_levels))

# ==============================================================================
# 7. Estimand-specific mapping
# ==============================================================================

get_estimand_columns <- function(estimand) {
  if (identical(estimand, "RMST")) {
    return(list(
      pi1 = rmst_interim_prob_col, pi2 = rmst_final_prob_col, true = rmst_true_col,
      mean1 = "post_mean_1", mean2 = "post_mean_2", q025_2 = "post_q025_2", q975_2 = "post_q975_2"
    ))
  }
  
  if (identical(estimand, "Median")) {
    return(list(
      pi1 = median_interim_prob_col, pi2 = median_final_prob_col, true = median_true_col,
      mean1 = "post_mean_med_1", mean2 = "post_mean_med_2",
      q025_2 = "post_q025_med_2", q975_2 = "post_q975_med_2"
    ))
  }
  
  stop("Unknown estimand: ", estimand)
}

# ==============================================================================
# 8. Component OCs for one threshold pair
# ==============================================================================

make_component_oc_one <- function(dat, estimand, nu1_E, nu2_E) {
  cols <- get_estimand_columns(estimand)
  
  dat %>%
    dplyr::mutate(
      pi1 = .data[[cols$pi1]],
      pi2 = .data[[cols$pi2]],
      I = pi1 >= nu1_E,
      F_cf = pi2 >= nu2_E,
      F_terminal = !I & F_cf,
      S = I | F_terminal,
      terminal_events = dplyr::if_else(I, as.numeric(events_interim), as.numeric(events_final), missing = NA_real_)
    ) %>%
    dplyr::group_by(model, N, prior_id, prior_label, profile) %>%
    dplyr::summarise(
      n_rep = dplyr::n(),
      n_rep_OC = sum(!is.na(S)),
      ESP = pct(safe_mean(I)),
      FSP = pct(safe_mean(F_terminal)),
      Success = pct(safe_mean(S)),
      ENE = safe_mean(terminal_events),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      BP = dplyr::if_else(prior_id == benefit_prior_id, Success, NA_real_),
      FPR = dplyr::if_else(prior_id %in% c(exact_null_prior_id, no_benefit_prior_id), Success, NA_real_),
      Estimand = estimand,
      nu1_E = nu1_E,
      nu2_E = nu2_E,
      success_identity_error = abs(Success - (ESP + FSP)) > 1e-8
    ) %>%
    dplyr::select(
      Estimand, model, N, prior_id, prior_label, profile, nu1_E, nu2_E,
      n_rep, n_rep_OC, ESP, FSP, ENE, BP, FPR, Success, success_identity_error
    )
}

make_component_oc_grid <- function(dat, estimand) {
  purrr::pmap_dfr(
    threshold_grid %>% dplyr::select(nu1_E, nu2_E),
    function(nu1_E, nu2_E) make_component_oc_one(dat, estimand, nu1_E, nu2_E)
  )
}

cat("Evaluating ", nrow(threshold_grid), " threshold pairs for RMST...\n", sep = "")
oc_rmst <- make_component_oc_grid(dat, "RMST")

cat("Evaluating ", nrow(threshold_grid), " threshold pairs for Median...\n", sep = "")
oc_median <- make_component_oc_grid(dat, "Median")

oc_all <- dplyr::bind_rows(oc_rmst, oc_median)

if (any(oc_all$success_identity_error, na.rm = TRUE)) {
  stop("At least one OC row violates Success = ESP + FSP.")
}

# ==============================================================================
# 9. Exact-null FPR surface -- Prior 1
# ==============================================================================

exact_null_calibration <- oc_all %>%
  dplyr::filter(prior_id == exact_null_prior_id) %>%
  dplyr::group_by(Estimand, model, N, nu1_E, nu2_E) %>%
  dplyr::summarise(
    ESP_exact = safe_mean(ESP),
    FSP_exact = safe_mean(FSP),
    FPR_exact = safe_mean(FPR),
    ENE_exact = safe_mean(ENE),
    .groups = "drop"
  ) %>%
  dplyr::mutate(Pass_FPR = FPR_exact <= fpr_target_pct)

# ==============================================================================
# 10. Identify exact-null FPR boundary first
# ==============================================================================

fpr_boundary_summary <- exact_null_calibration %>%
  dplyr::filter(Pass_FPR, is.finite(FPR_exact)) %>%
  dplyr::group_by(Estimand, model, N) %>%
  dplyr::summarise(
    exact_2_5_available = any(abs(FPR_exact - fpr_target_pct) <= numeric_tolerance),
    target_FPR = max(FPR_exact, na.rm = TRUE),
    n_pairs_at_target = sum(abs(FPR_exact - max(FPR_exact, na.rm = TRUE)) <= numeric_tolerance),
    distance_below_2_5 = fpr_target_pct - max(FPR_exact, na.rm = TRUE),
    fallback_used = max(FPR_exact, na.rm = TRUE) < fpr_target_pct - numeric_tolerance,
    .groups = "drop"
  ) %>%
  dplyr::arrange(
    factor(Estimand, levels = c("RMST", "Median")),
    factor(model, levels = names(model_prefixes)),
    N
  )

expected_fpr_settings <- tidyr::expand_grid(
  Estimand = c("RMST", "Median"),
  model = names(model_prefixes),
  N = N_keep
)

missing_fpr_settings <- expected_fpr_settings %>%
  dplyr::anti_join(fpr_boundary_summary, by = c("Estimand", "model", "N"))

if (nrow(missing_fpr_settings) > 0L) {
  stop(
    "No threshold pair satisfied exact-null FPR <= ", fpr_target_pct, "% for: ",
    paste0(
      missing_fpr_settings$Estimand, " / ", missing_fpr_settings$model,
      " / N=", missing_fpr_settings$N,
      collapse = "; "
    )
  )
}

if (any(fpr_boundary_summary$exact_2_5_available &
        abs(fpr_boundary_summary$target_FPR - fpr_target_pct) > numeric_tolerance)) {
  stop("An exact 2.5% FPR pair exists but target_FPR was not set to 2.5%.")
}

if (any(fpr_boundary_summary$fallback_used &
        fpr_boundary_summary$target_FPR >= fpr_target_pct - numeric_tolerance)) {
  stop("Fallback flag is inconsistent with target_FPR.")
}

target_fpr_pairs <- exact_null_calibration %>%
  dplyr::inner_join(fpr_boundary_summary, by = c("Estimand", "model", "N")) %>%
  dplyr::filter(abs(FPR_exact - target_FPR) <= numeric_tolerance) %>%
  dplyr::left_join(threshold_grid, by = c("nu1_E", "nu2_E")) %>%
  dplyr::arrange(
    factor(Estimand, levels = c("RMST", "Median")),
    factor(model, levels = names(model_prefixes)),
    N, grid_order
  )

target_fpr_pair_check <- target_fpr_pairs %>%
  dplyr::count(Estimand, model, N, name = "n_target_pairs") %>%
  dplyr::left_join(fpr_boundary_summary, by = c("Estimand", "model", "N"))

if (any(target_fpr_pair_check$n_target_pairs != target_fpr_pair_check$n_pairs_at_target)) {
  stop("Target-FPR threshold-pair counts are inconsistent.")
}

# ==============================================================================
# 11. Profile-specific benefit BP surface -- Prior 3
# ==============================================================================

benefit_power <- oc_all %>%
  dplyr::filter(prior_id == benefit_prior_id) %>%
  dplyr::select(Estimand, model, N, profile, nu1_E, nu2_E, BP, ESP, FSP, ENE)

benefit_profile_check <- benefit_power %>%
  dplyr::distinct(Estimand, model, N, profile) %>%
  dplyr::count(Estimand, model, N, name = "n_profiles")

if (any(benefit_profile_check$n_profiles != length(benefit_profiles))) {
  stop("Benefit calibration must contain exactly three Prior-3 profiles for every estimand x model x N setting.")
}

# ==============================================================================
# 12. Restrict to target-FPR pairs and select profile-specific thresholds
# ==============================================================================

profile_target_candidates <- benefit_power %>%
  dplyr::inner_join(
    target_fpr_pairs %>%
      dplyr::select(
        Estimand, model, N, nu1_E, nu2_E, grid_order,
        FPR_exact, ESP_exact, FSP_exact, ENE_exact, target_FPR
      ),
    by = c("Estimand", "model", "N", "nu1_E", "nu2_E")
  ) %>%
  dplyr::mutate(profile = factor(profile, levels = benefit_profiles)) %>%
  dplyr::group_by(Estimand, model, N, profile) %>%
  dplyr::arrange(dplyr::desc(BP), grid_order, .by_group = TRUE) %>%
  dplyr::mutate(candidate_rank = dplyr::row_number()) %>%
  dplyr::ungroup()

selected_thresholds <- profile_target_candidates %>%
  dplyr::filter(candidate_rank == 1L) %>%
  dplyr::select(
    Estimand, model, N, profile, nu1_E, nu2_E,
    FPR_exact, BP, ESP, FSP, ENE, grid_order
  ) %>%
  dplyr::arrange(
    factor(Estimand, levels = c("RMST", "Median")),
    factor(model, levels = names(model_prefixes)),
    N, factor(profile, levels = benefit_profiles)
  )

expected_calibration_settings <- tidyr::expand_grid(
  Estimand = c("RMST", "Median"),
  model = names(model_prefixes),
  N = N_keep,
  profile = benefit_profiles
)

missing_selected_settings <- expected_calibration_settings %>%
  dplyr::anti_join(selected_thresholds, by = c("Estimand", "model", "N", "profile"))

if (nrow(missing_selected_settings) > 0L) {
  stop(
    "No profile-specific threshold was selected for: ",
    paste0(
      missing_selected_settings$Estimand, " / ", missing_selected_settings$model,
      " / N=", missing_selected_settings$N, " / ", missing_selected_settings$profile,
      collapse = "; "
    )
  )
}

if (nrow(selected_thresholds) != 24L) {
  stop("Expected 24 selected profile-specific threshold pairs, but found ", nrow(selected_thresholds), ".")
}

# ==============================================================================
# 13. Tie summary within fixed target-FPR set
# ==============================================================================

tie_summary <- profile_target_candidates %>%
  dplyr::group_by(Estimand, model, N, profile) %>%
  dplyr::mutate(best_BP = max(BP, na.rm = TRUE)) %>%
  dplyr::filter(abs(BP - best_BP) <= numeric_tolerance) %>%
  dplyr::summarise(
    selected_FPR = dplyr::first(FPR_exact),
    best_BP = dplyr::first(best_BP),
    n_BP_ties = as.integer(dplyr::n()),
    selected_nu1_E = nu1_E[which.min(grid_order)],
    selected_nu2_E = nu2_E[which.min(grid_order)],
    .groups = "drop"
  ) %>%
  dplyr::arrange(
    factor(Estimand, levels = c("RMST", "Median")),
    factor(model, levels = names(model_prefixes)),
    N, factor(profile, levels = benefit_profiles)
  )

# ==============================================================================
# 14. Prior-2 no-benefit FPR at selected thresholds
# ==============================================================================

no_benefit_surface <- oc_all %>%
  dplyr::filter(prior_id == no_benefit_prior_id) %>%
  dplyr::select(
    Estimand, model, N, nu1_E, nu2_E,
    ESP_no_benefit = ESP,
    FSP_no_benefit = FSP,
    FPR_no_benefit = FPR,
    ENE_no_benefit = ENE
  )

selected_null_comparison <- selected_thresholds %>%
  dplyr::left_join(
    no_benefit_surface,
    by = c("Estimand", "model", "N", "nu1_E", "nu2_E")
  ) %>%
  dplyr::select(
    Estimand, model, N, profile, nu1_E, nu2_E,
    FPR_exact, FPR_no_benefit, ESP_no_benefit, FSP_no_benefit, ENE_no_benefit
  ) %>%
  dplyr::arrange(
    factor(Estimand, levels = c("RMST", "Median")),
    factor(model, levels = names(model_prefixes)),
    N, factor(profile, levels = benefit_profiles)
  )

# ==============================================================================
# 15. Estimation summaries for Prior 3
# ==============================================================================

make_estimation_summary <- function(dat, estimand) {
  cols <- get_estimand_columns(estimand)
  
  dat %>%
    dplyr::filter(prior_id == benefit_prior_id) %>%
    dplyr::mutate(
      true_value = .data[[cols$true]],
      interim_estimate = .data[[cols$mean1]],
      final_estimate = .data[[cols$mean2]],
      final_lower = .data[[cols$q025_2]],
      final_upper = .data[[cols$q975_2]],
      interim_error = interim_estimate - true_value,
      final_error = final_estimate - true_value,
      final_cover = true_value >= final_lower & true_value <= final_upper
    ) %>%
    dplyr::group_by(model, N, profile) %>%
    dplyr::summarise(
      `Interim bias` = safe_mean(interim_error),
      `Final bias` = safe_mean(final_error),
      `Final RMSE` = safe_rmse(final_error),
      `Final coverage` = pct(safe_mean(final_cover)),
      .groups = "drop"
    ) %>%
    dplyr::mutate(Estimand = estimand) %>%
    dplyr::select(Estimand, model, N, profile, `Interim bias`, `Final bias`, `Final RMSE`, `Final coverage`)
}

estimation_summary <- dplyr::bind_rows(
  make_estimation_summary(dat, "RMST"),
  make_estimation_summary(dat, "Median")
)

selected_benefit_full <- selected_thresholds %>%
  dplyr::left_join(estimation_summary, by = c("Estimand", "model", "N", "profile")) %>%
  dplyr::arrange(
    factor(Estimand, levels = c("RMST", "Median")),
    factor(model, levels = names(model_prefixes)),
    N, factor(profile, levels = benefit_profiles)
  )

# ==============================================================================
# 16. Compact checking summaries
# ==============================================================================

replicate_summary <- replicate_id_check %>%
  dplyr::group_by(model, N) %>%
  dplyr::summarise(
    all_complete = all(Complete),
    min_reps = min(n_rep),
    max_reps = max(n_rep),
    total_missing_ids = sum(n_missing),
    total_extra_ids = sum(n_extra),
    .groups = "drop"
  )

duplicate_summary <- tibble::tibble(
  Check = "Duplicated model x N x prior x profile x replicate rows",
  Count = nrow(duplicate_check)
)

report_summary <- tibble::tibble(
  Item = c(
    "Threshold range",
    "Threshold increment",
    "Threshold values",
    "Ordered threshold pairs",
    "Primary FPR rule",
    "Secondary BP rule",
    "Profile-specific settings"
  ),
  Value = c(
    "0.950 to 0.995",
    "0.001",
    as.character(length(threshold_values)),
    as.character(nrow(threshold_grid)),
    "Use 2.5% if attainable; otherwise use the largest attainable FPR below 2.5%",
    "Within target-FPR pairs, maximize BP separately by profile",
    "24 = 2 estimands x 2 models x 2 N x 3 profiles"
  )
)

# ==============================================================================
# 17. Console summary
# ==============================================================================

cat("\n============================================================\n")
cat("BDCT PROFILE-SPECIFIC THRESHOLD CALIBRATION -- FPR FIRST\n")
cat("============================================================\n")
cat("Results directory: ", results_dir, "\n", sep = "")
cat("Threshold pairs: ", nrow(threshold_grid), "\n", sep = "")

cat("\n==================== FPR BOUNDARY ====================\n")
print(fpr_boundary_summary, n = Inf)
cat("Fallback settings (target FPR below 2.5%): ", sum(fpr_boundary_summary$fallback_used), "\n", sep = "")

cat("\n==================== SELECTED THRESHOLDS ====================\n")
print(selected_thresholds, n = Inf)

cat("\n==================== TIES AT TARGET FPR ====================\n")
print(tie_summary, n = Inf)

cat("\n==================== EXACT NULL VS NO-BENEFIT REGION ====================\n")
print(selected_null_comparison, n = Inf)

cat("\n==================== SELECTED OCs + ESTIMATION ====================\n")
print(selected_benefit_full, n = Inf)

# ==============================================================================
# 18. Compact standalone HTML report
# ==============================================================================

html_dir <- file.path(results_dir, "HTML")
if (!dir.exists(html_dir)) dir.create(html_dir, recursive = TRUE)
html_file <- file.path(html_dir, "BDCT_profile_specific_threshold_calibration_FPR_target.html")

format_html_data <- function(dat) {
  dat <- as.data.frame(dat, stringsAsFactors = FALSE)
  if (nrow(dat) == 0L) return(data.frame(Message = "No rows.", check.names = FALSE))
  
  for (j in seq_along(dat)) {
    if (is.list(dat[[j]])) {
      dat[[j]] <- vapply(
        dat[[j]],
        function(x) if (length(x) == 0L || all(is.na(x))) "" else paste(x, collapse = ", "),
        character(1)
      )
    }
  }
  
  threshold_cols <- intersect(c("nu1_E", "nu2_E", "selected_nu1_E", "selected_nu2_E"), names(dat))
  percentage_cols <- grep("FPR|BP|ESP|FSP|coverage|distance_below", names(dat), ignore.case = TRUE, value = TRUE)
  ene_cols <- grep("ENE", names(dat), ignore.case = TRUE, value = TRUE)
  estimation_cols <- grep("bias|RMSE", names(dat), ignore.case = TRUE, value = TRUE)
  
  for (nm in threshold_cols) {
    if (is.numeric(dat[[nm]])) dat[[nm]] <- ifelse(is.na(dat[[nm]]), "", sprintf("%.3f", dat[[nm]]))
  }
  
  for (nm in setdiff(percentage_cols, threshold_cols)) {
    if (is.numeric(dat[[nm]])) dat[[nm]] <- ifelse(is.na(dat[[nm]]), "", sprintf("%.1f", dat[[nm]]))
  }
  
  for (nm in ene_cols) {
    if (is.numeric(dat[[nm]])) dat[[nm]] <- ifelse(is.na(dat[[nm]]), "", sprintf("%.1f", dat[[nm]]))
  }
  
  for (nm in estimation_cols) {
    if (is.numeric(dat[[nm]])) dat[[nm]] <- ifelse(is.na(dat[[nm]]), "", sprintf("%.3f", dat[[nm]]))
  }
  
  integer_like <- c(
    "N", "n_pairs_at_target", "n_BP_ties", "min_reps", "max_reps",
    "total_missing_ids", "total_extra_ids", "Count"
  )
  
  for (nm in intersect(integer_like, names(dat))) {
    if (is.numeric(dat[[nm]])) dat[[nm]] <- ifelse(is.na(dat[[nm]]), "", sprintf("%.0f", dat[[nm]]))
  }
  
  remaining_numeric <- names(dat)[vapply(dat, is.numeric, logical(1))]
  for (nm in remaining_numeric) {
    dat[[nm]] <- ifelse(is.na(dat[[nm]]), "", sprintf("%.3f", dat[[nm]]))
  }
  
  dat
}

html_table <- function(title, dat, note = NULL) {
  dat_display <- format_html_data(dat)
  
  table_html <- knitr::kable(
    dat_display,
    format = "html",
    escape = TRUE,
    row.names = FALSE,
    table.attr = 'class="bdct-table"'
  )
  
  htmltools::tags$section(
    class = "table-section",
    htmltools::tags$h2(title),
    if (!is.null(note)) htmltools::tags$p(class = "note", note),
    htmltools::tags$div(class = "table-wrap", htmltools::HTML(table_html))
  )
}

html_doc <- htmltools::tags$html(
  htmltools::tags$head(
    htmltools::tags$meta(charset = "utf-8"),
    htmltools::tags$meta(name = "viewport", content = "width=device-width, initial-scale=1"),
    htmltools::tags$title("BDCT FPR-first profile-specific threshold calibration"),
    htmltools::tags$style(htmltools::HTML("
      body { font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Arial, sans-serif; margin: 28px; color: #1f2937; background: #fff; }
      h1 { margin-bottom: 6px; }
      h2 { margin-top: 32px; margin-bottom: 10px; font-size: 1.22rem; }
      .subtitle { color: #4b5563; margin-top: 0; margin-bottom: 24px; }
      .note { color: #4b5563; margin-top: -3px; margin-bottom: 10px; }
      .table-section { margin-bottom: 34px; }
      .table-wrap { overflow-x: auto; width: 100%; }
      table.bdct-table { border-collapse: collapse; font-size: 13px; white-space: nowrap; width: auto; min-width: 55%; }
      table.bdct-table th { background: #f3f4f6; font-weight: 600; }
      table.bdct-table th, table.bdct-table td { border: 1px solid #d1d5db; padding: 6px 9px; text-align: right; }
      table.bdct-table th:first-child, table.bdct-table td:first-child { text-align: left; }
      table.bdct-table tr:nth-child(even) td { background: #fafafa; }
      .footer { color: #6b7280; margin-top: 36px; font-size: 12px; }
    "))
  ),
  htmltools::tags$body(
    htmltools::tags$h1("BDCT FPR-first profile-specific threshold calibration"),
    htmltools::tags$p(
      class = "subtitle",
      "Exact-null FPR is calibrated first; profile-specific BP is optimized only among threshold pairs attaining the selected FPR boundary."
    ),
    
    html_table("1. Calibration specification", report_summary),
    
    html_table(
      "2. Data checks",
      dplyr::bind_rows(
        duplicate_summary,
        tibble::tibble(
          Check = "All expected replicate sets complete",
          Count = ifelse(all(replicate_summary$all_complete), 1L, 0L)
        )
      )
    ),
    
    html_table("3. Replicate summary", replicate_summary),
    
    html_table(
      "4. Exact-null FPR boundary",
      fpr_boundary_summary,
      "If exact 2.5% is attainable, target_FPR is 2.5%. Otherwise the largest attainable value below 2.5% is used automatically. exact_2_5_available and fallback_used show which case applies."
    ),
    
    html_table(
      "5. Selected profile-specific thresholds",
      selected_thresholds,
      "Within the target-FPR threshold set, BP is maximized separately for constant, increasing, and waning profiles."
    ),
    
    html_table(
      "6. BP ties at the target FPR",
      tie_summary,
      "n_BP_ties counts threshold pairs at the selected FPR that attain the same maximum BP. Grid order resolves exact BP ties."
    ),
    
    html_table(
      "7. Exact-null versus continuous no-benefit FPR",
      selected_null_comparison
    ),
    
    html_table(
      "8. Selected decision and estimation OCs",
      selected_benefit_full
    ),
    
    htmltools::tags$p(
      class = "footer",
      paste0("Generated from existing simulation results. No MCMC was rerun. HTML file: ", html_file)
    )
  )
)

htmltools::save_html(html_doc, file = html_file, background = "white")

if (!file.exists(html_file)) stop("HTML report was not written: ", html_file)

cat("\n============================================================\n")
cat("CALIBRATION COMPLETE\n")
cat("============================================================\n")
cat("Selected profile-specific threshold pairs: ", nrow(selected_thresholds), "\n", sep = "")
cat("HTML report: ", html_file, "\n", sep = "")
cat("The full 24,840-row calibration surface is not embedded in the HTML.\n")
