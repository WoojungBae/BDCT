# ==============================================================================
# BDCT full checked summary and manuscript tables/figures
# Manuscript terminology: efficacy/no-efficacy and design prior.
# Legacy saved-file identifiers containing "benefit"/"no_benefit" are retained
# internally only so existing simulation result files remain compatible.
# Two profile-specific calibration strategies; Priors 1-3 available; no MCMC rerun.
#
# Threshold grid: 0.950-0.995 by 0.001, nu1_E > nu2_E (1,035 pairs).
#
# Strategy 1: FPR-only / FPR-first calibration
#   For each Estimand x model x N, use exact FPR 2.5% when attainable;
#   otherwise use the largest attainable exact-null FPR below 2.5%.
#   Within that FPR set, maximize BP separately by efficacy profile.
#   Conclusion OC order: FPR, BP, PPV, NPV.
#
# Strategy 2: PPV + FPR calibration at omega = 0.5
#   Retain candidate rules satisfying exact-null FPR <= 2.5% and
#   PPV(0.5) >= 97.5%, then maximize BP separately by efficacy profile.
#   BP >= 90% is used to assess full design adequacy.
#   Conclusion OC order: PPV, FPR, BP, NPV.
#
# Both strategies are calibrated for RMST and mPFS; manuscript tables split the two estimands.
#   Candidate setting: Estimand, Model, N, Profile
#   Calibrated thresholds: nu1, nu2
#   Success by look: SP1, SP2
#   Estimation performance: interim bias, final bias, RMSE, CP
#   PPV/NPV are evaluated at omega = 0.5.
#   SP1/SP2 and estimation summaries are evaluated under the corresponding efficacy distribution.
# Main Tables 2 and 3 are not generated.
# Main Table 1: RMST only.
# Supplement S1: RMST exact-null vs no-efficacy-region FPR for FPR-only calibration.
# Supplement S2: mPFS FPR-only decision/estimation OCs.
# Supplement S3: mPFS PPV+FPR decision/estimation OCs.
# Figure 1 only; Figure 2 is not generated.
# ==============================================================================
suppressPackageStartupMessages({ library(dplyr); library(purrr); library(tidyr); library(tibble); library(ggplot2); library(scales); library(knitr); library(kableExtra) })
# ==============================================================================
# 1. Settings
# ==============================================================================
results_dir_candidates <- c( "/Users/woojung/Documents/Rproject/BDCT/SimulationStudy/Results/ResultsThin1",
                             "/Users/woojung/Documents/Rproject/BDCT/SimulationStudy/Results",
                             "/Users/woojung/Documents/Rproject/BDCT/Results/ResultsThin1", "/Users/woojung/Documents/Rproject/BDCT/Results",
                             "/home/woojung.bae/BDCT/Results_thin1", "/home/woojung.bae/BDCT/Results" )
existing_results_dirs <- results_dir_candidates[dir.exists(results_dir_candidates)]
results_dir <- if (length(existing_results_dirs) > 0L) {
  existing_results_dirs[1L]
} else {
  results_dir_candidates[1L]
}
if (!dir.exists(results_dir)) {
  stop( "Results directory was not found. Checked: ", paste(results_dir_candidates, collapse = "; ") )
}
# Prespecified design-prior efficacy mixture used for PPV calibration and manuscript predictive summaries.
ppv_calibration_omega <- 0.5

# Candidate-level PPV calibration target.
ppv_calibration_target <- 0.975

# Keep the manuscript reference design-prior mixture identical to the calibration mixture.
predictive_reference_omega <- ppv_calibration_omega
# ==============================================================================
# FIGURE 1. Predictive reliability at FPR = 0.025
#
# Reliability requirements:
#   PPV >= 97.5%
#   NPV >= 90%
#
# The vertical reference line is omega = 0.5.
# PPV and NPV at omega = 0.5 are shown beside the y-axis in each BP panel.
# ==============================================================================
predictive_figure_BP <- c(0.80, 0.90)
predictive_figure_FPR <- 0.025
predictive_figure_omega <- seq(0, 1, by = 0.001)

predictive_figure_data <- tidyr::expand_grid(
  omega_D = predictive_figure_omega,
  BP = predictive_figure_BP
) %>%
  dplyr::mutate(
    FPR = predictive_figure_FPR,
    PPV_denominator =
      omega_D * BP +
      (1 - omega_D) * FPR,
    NPV_denominator =
      (1 - omega_D) * (1 - FPR) +
      omega_D * (1 - BP),
    PPV = dplyr::if_else(
      PPV_denominator > 0,
      (omega_D * BP) / PPV_denominator,
      NA_real_
    ),
    NPV = dplyr::if_else(
      NPV_denominator > 0,
      ((1 - omega_D) * (1 - FPR)) /
        NPV_denominator,
      NA_real_
    ),
    BP_panel = factor(
      BP,
      levels = predictive_figure_BP,
      labels = c("BP = 0.80", "BP = 0.90")
    )
  ) %>%
  dplyr::select(
    omega_D,
    BP,
    BP_panel,
    PPV,
    NPV
  ) %>%
  tidyr::pivot_longer(
    cols = c(PPV, NPV),
    names_to = "Predictive_value",
    values_to = "Value"
  ) %>%
  dplyr::mutate(
    Predictive_value = factor(
      Predictive_value,
      levels = c("PPV", "NPV")
    )
  )

predictive_reference_points <- tibble::tibble(
  BP = predictive_figure_BP,
  omega_D = predictive_reference_omega
) %>%
  dplyr::mutate(
    FPR = predictive_figure_FPR,
    PPV =
      (omega_D * BP) /
      (
        omega_D * BP +
          (1 - omega_D) * FPR
      ),
    NPV =
      ((1 - omega_D) * (1 - FPR)) /
      (
        (1 - omega_D) * (1 - FPR) +
          omega_D * (1 - BP)
      ),
    BP_panel = factor(
      BP,
      levels = predictive_figure_BP,
      labels = c("BP = 0.80", "BP = 0.90")
    )
  ) %>%
  dplyr::select(
    BP,
    BP_panel,
    omega_D,
    PPV,
    NPV
  ) %>%
  tidyr::pivot_longer(
    cols = c(PPV, NPV),
    names_to = "Predictive_value",
    values_to = "Value"
  ) %>%
  dplyr::mutate(
    Predictive_value = factor(
      Predictive_value,
      levels = c("PPV", "NPV")
    ),
    label = scales::percent(
      Value,
      accuracy = 0.1
    )
  )

# Expected values at omega = 0.5:
#   BP = 0.80: PPV = 97.0%, NPV = 83.0%
#   BP = 0.90: PPV = 97.3%, NPV = 90.7%
predictive_reference_check <- predictive_reference_points %>%
  dplyr::mutate(
    Predictive_value = as.character(Predictive_value)
  ) %>%
  dplyr::select(
    BP,
    Predictive_value,
    Value
  ) %>%
  tidyr::pivot_wider(
    names_from = Predictive_value,
    values_from = Value
  ) %>%
  dplyr::arrange(BP)

expected_reference_PPV <- c(
  0.969696969696970,
  0.972972972972973
)
expected_reference_NPV <- c(
  0.829787234042553,
  0.906976744186047
)

if (nrow(predictive_reference_check) != 2L ||
    any(
      abs(
        predictive_reference_check$PPV -
        expected_reference_PPV
      ) > 1e-10
    ) ||
    any(
      abs(
        predictive_reference_check$NPV -
        expected_reference_NPV
      ) > 1e-10
    )) {
  stop(
    "Figure 1 reference values at omega = 0.5 ",
    "do not match the expected PPV/NPV values."
  )
}

predictive_values_plot <- ggplot2::ggplot(
  predictive_figure_data,
  ggplot2::aes(
    x = omega_D,
    y = Value,
    color = Predictive_value,
    linetype = Predictive_value
  )
) +
  ggplot2::geom_hline(
    yintercept = c(0.90, 0.975),
    color = "grey72",
    linewidth = 0.35,
    linetype = "dotted"
  ) +
  ggplot2::geom_vline(
    xintercept = predictive_reference_omega,
    color = "grey25",
    linewidth = 0.50,
    linetype = "twodash"
  ) +
  ggplot2::geom_segment(
    data = predictive_reference_points,
    ggplot2::aes(
      x = 0.035,
      xend = omega_D,
      y = Value,
      yend = Value,
      color = Predictive_value
    ),
    inherit.aes = FALSE,
    linewidth = 0.35,
    linetype = "dotted",
    show.legend = FALSE
  ) +
  ggplot2::geom_segment(
    data = predictive_reference_points,
    ggplot2::aes(
      x = 0,
      xend = 0.018,
      y = Value,
      yend = Value,
      color = Predictive_value
    ),
    inherit.aes = FALSE,
    linewidth = 0.75,
    show.legend = FALSE
  ) +
  ggplot2::geom_line(
    linewidth = 1.05,
    lineend = "round"
  ) +
  ggplot2::geom_point(
    data = predictive_reference_points,
    ggplot2::aes(
      x = omega_D,
      y = Value,
      color = Predictive_value
    ),
    inherit.aes = FALSE,
    size = 2.5,
    show.legend = FALSE
  ) +
  ggplot2::geom_label(
    data = predictive_reference_points,
    ggplot2::aes(
      x = 0.025,
      y = Value,
      label = label,
      color = Predictive_value
    ),
    inherit.aes = FALSE,
    hjust = 0,
    size = 2.9,
    label.size = 0,
    fill = "white",
    label.padding = grid::unit(0.08, "lines"),
    show.legend = FALSE
  ) +
  ggplot2::facet_wrap(
    ggplot2::vars(BP_panel),
    nrow = 1
  ) +
  ggplot2::scale_color_manual(
    values = c(
      "PPV" = "#1F77B4",
      "NPV" = "#D55E00"
    ),
    breaks = c("PPV", "NPV"),
    name = NULL
  ) +
  ggplot2::scale_linetype_manual(
    values = c(
      "PPV" = "solid",
      "NPV" = "dashed"
    ),
    breaks = c("PPV", "NPV"),
    name = NULL
  ) +
  ggplot2::scale_x_continuous(
    name = expression(
      "Design-prior probability of efficacy, " * omega[D]
    ),
    breaks = sort(
      unique(
        c(
          seq(0, 1, by = 0.2),
          predictive_reference_omega
        )
      )
    ),
    labels = scales::number_format(
      accuracy = 0.1
    ),
    expand = ggplot2::expansion(
      mult = c(0, 0)
    )
  ) +
  ggplot2::scale_y_continuous(
    name = "Predictive value",
    breaks = c(
      0.80,
      0.85,
      0.90,
      0.95,
      0.975,
      1.00
    ),
    labels = scales::percent_format(
      accuracy = 0.1
    ),
    expand = ggplot2::expansion(
      mult = c(0, 0)
    )
  ) +
  ggplot2::coord_cartesian(
    xlim = c(0, 1),
    ylim = c(0.80, 1.00),
    clip = "on"
  ) +
  ggplot2::theme_bw(
    base_size = 12
  ) +
  ggplot2::theme(
    strip.background = ggplot2::element_rect(
      fill = "grey92",
      color = "grey55",
      linewidth = 0.5
    ),
    strip.text = ggplot2::element_text(
      size = 11.5,
      face = "bold"
    ),
    axis.title.x = ggplot2::element_text(
      size = 12,
      margin = ggplot2::margin(t = 10)
    ),
    axis.title.y = ggplot2::element_text(
      size = 12,
      margin = ggplot2::margin(r = 10)
    ),
    axis.text = ggplot2::element_text(
      size = 10
    ),
    panel.grid = ggplot2::element_blank(),
    panel.spacing = grid::unit(
      1.0,
      "lines"
    ),
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.text = ggplot2::element_text(
      size = 10.5
    ),
    legend.key.width = grid::unit(
      1.5,
      "cm"
    ),
    plot.margin = ggplot2::margin(
      10,
      12,
      10,
      10
    )
  )

figure_output_dir <- file.path(results_dir, "LaTeX")
if (!dir.exists(figure_output_dir)) {
  dir.create(figure_output_dir, recursive = TRUE)
}

predictive_figure_pdf <- file.path(
  figure_output_dir,
  "BDCT_predictive_values.pdf"
)
predictive_figure_tex <- file.path(
  figure_output_dir,
  "figureBDCT_predictive_values.tex"
)

ggplot2::ggsave(
  filename = predictive_figure_pdf,
  plot = predictive_values_plot,
  width = 8.8,
  height = 4.8,
  units = "in",
  device = "pdf"
)

predictive_figure_caption <- paste0(
  "PPV and NPV as functions of the design-prior probability of efficacy, ",
  "$\\omega_{\\mathcal{D}}$, at $\\mathrm{FPR}=0.025$ and ",
  "$\\mathrm{BP}=0.80$ or $0.90$. ",
  "Horizontal reference lines denote PPV $97.5\\%$ and NPV $90\\%$. ",
  "The vertical reference line denotes $\\omega=0.5$. ",
  "At $\\omega=0.5$, the corresponding PPV and NPV values are marked ",
  "and labeled at the left edge of each panel."
)

predictive_figure_tex_lines <- c(
  "\\begin{figure}[!tbp]",
  "\\centering",
  "\\includegraphics[width=\\textwidth]{image/BDCT_predictive_values.pdf}",
  paste0(
    "\\caption{\\label{figureBDCT_predictive_values}",
    predictive_figure_caption,
    "}"
  ),
  "\\end{figure}"
)

writeLines(
  predictive_figure_tex_lines,
  con = predictive_figure_tex,
  useBytes = TRUE
)

if (!file.exists(predictive_figure_pdf)) {
  stop("Predictive-values PDF figure was not written: ", predictive_figure_pdf)
}
if (!file.exists(predictive_figure_tex)) {
  stop("Predictive-values LaTeX figure file was not written: ", predictive_figure_tex)
}

# ==============================================================================
# TABLE WORKFLOW STARTS HERE
# ==============================================================================
# ==============================================================================
model_prefixes <- c("Original GP-PH" = "GP_POST_delta0", "Independent GP" = "IndepGP_POST_delta0")
N_keep <- c(700L, 1000L)
benefit_profile_levels <- c("constant", "increasing", "waning")
all_profile_levels <- c(benefit_profile_levels, "none")
# ------------------------------------------------------------------------------
# Current BDCT Priors 1-3
#
# Prior 4 is intentionally excluded until its simulation results are available.
# Existing filenames are retained so Priors 1-3 can be reused without rerunning.
# ------------------------------------------------------------------------------
design_prior_specs <- tibble::tribble(
  ~prior_id, ~design_prior_type,                    ~prior_component, ~prior_family, ~omega_D, ~profile,      ~output_tag,
  1L,        "exact_null_point_mass",               "no_benefit",     "point_mass",  0.00,     "none",        "exact_null",
  2L,        "no_benefit_continuous",               "no_benefit",     "continuous",  0.00,     "none",        "no_benefit_continuous",
  3L,        "benefit_positive_truncnorm_sd010",    "benefit",        "continuous",  1.00,     "constant",    "benefit_positive_truncnorm_constant",
  3L,        "benefit_positive_truncnorm_sd010",    "benefit",        "continuous",  1.00,     "increasing",  "benefit_positive_truncnorm_increasing",
  3L,        "benefit_positive_truncnorm_sd010",    "benefit",        "continuous",  1.00,     "waning",      "benefit_positive_truncnorm_waning"
)
exact_null_prior_id <- 1L
no_benefit_prior_id <- 2L
benefit_prior_id <- 3L
# ------------------------------------------------------------------------------
# Requested design-prior efficacy:no-efficacy ratios
#
# The no-efficacy component for this table is the exact null.
# ------------------------------------------------------------------------------
mixture_grid_base <- tibble::tibble(
  Ratio = c("100:0", "75:25", "50:50", "25:75", "0:100"),
  omega_D = c(1.00, 0.75, 0.50, 0.25, 0.00)
)

mixture_grid <- dplyr::bind_rows(
  mixture_grid_base,
  tibble::tibble(Ratio = "Reference", omega_D = predictive_reference_omega)
) %>%
  dplyr::arrange(dplyr::desc(omega_D)) %>%
  dplyr::distinct(omega_D, .keep_all = TRUE)
# ------------------------------------------------------------------------------
# Expected replicates
# ------------------------------------------------------------------------------
n_rep_expected <- 1000L
expected_rep_ids <- seq_len(n_rep_expected)
# ------------------------------------------------------------------------------
# Event-driven design
# ------------------------------------------------------------------------------
target_censor <- 0.40
tau_expected <- 24
# ------------------------------------------------------------------------------
# Numerical tolerance
# ------------------------------------------------------------------------------
numeric_tolerance <- 1e-10
# ------------------------------------------------------------------------------
# Saved posterior-probability columns
# ------------------------------------------------------------------------------
interim_prob_col <- "Pi1"
final_prob_col <- "Pi2_cf"
rmst_true_col <- "true_delta_rmst_rep"
med_true_col <- "true_delta_med_rep"
# ------------------------------------------------------------------------------
# Calibration constraints
# ------------------------------------------------------------------------------
fpr_target_pct <- 2.5
bp_target_pct <- 90.0
ppv_target_pct <- 100 * ppv_calibration_target
# ------------------------------------------------------------------------------
# MCMC diagnostic cutoffs
# ------------------------------------------------------------------------------
rhat_warn <- 1.05
ess_warn <- 100
median_na_warn_pct <- 5
# ==============================================================================
# ==============================================================================
# 2. Threshold grid
# ==============================================================================
# Fine grid: 0.950, 0.951, ..., 0.995. Only nu1_E > nu2_E is retained.
threshold_values <- round(seq(0.950, 0.995, by = 0.001), 3)
threshold_grid_full <- tidyr::expand_grid(nu1_E = threshold_values, nu2_E = threshold_values) %>%
  dplyr::filter(nu1_E > nu2_E) %>%
  dplyr::arrange(nu1_E, nu2_E)
threshold_grid_selection <- threshold_grid_full
threshold_grid_selection_order <- threshold_grid_selection %>% dplyr::mutate(grid_order = dplyr::row_number())
expected_threshold_pairs <- choose(length(threshold_values), 2)
calibration_index <- tidyr::expand_grid(Estimand = c("RMST", "Median"), model = names(model_prefixes), N = N_keep)
# 3. Utility functions
# ==============================================================================
safe_mean <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) == 0L) return(NA_real_)
  mean(x)
}
safe_sum <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) == 0L) return(NA_real_)
  sum(x)
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
safe_ratio <- function(numerator, denominator) {
  ifelse(is.na(denominator) | abs(denominator) <= numeric_tolerance, NA_real_, numerator / denominator)
}
pct <- function(x) 100 * x
format_two_decimals <- function(x) {
  ifelse(is.na(x), "--", formatC(x, format = "f", digits = 2, big.mark = ""))
}
format_four_decimals <- function(x) {
  ifelse(is.na(x), "--", formatC(x, format = "f", digits = 4, big.mark = ""))
}
format_table_numbers <- function(dat) {
  threshold_columns <- intersect(c("nu1_E", "nu2_E"), names(dat))
  numeric_columns <- names(dat)[vapply(dat, is.numeric, logical(1))]
  two_decimal_columns <- setdiff(numeric_columns, threshold_columns)
  out <- dat
  if (length(two_decimal_columns) > 0L) {
    out <- out %>% dplyr::mutate(dplyr::across(dplyr::all_of(two_decimal_columns), format_two_decimals))
  }
  if (length(threshold_columns) > 0L) {
    out <- out %>% dplyr::mutate(dplyr::across(dplyr::all_of(threshold_columns), format_four_decimals))
  }
  out
}
assert_columns <- function(dat, columns, context) {
  missing_columns <- setdiff(columns, names(dat))
  if (length(missing_columns) > 0L) {
    stop(context, " is missing the following columns: ", paste(missing_columns, collapse = ", "))
  }
  invisible(TRUE)
}
numeric_mismatch <- function(x, y, tolerance = numeric_tolerance) {
  xor(is.na(x), is.na(y)) | (!is.na(x) & !is.na(y) & abs(x - y) > tolerance)
}
character_mismatch <- function(x, y) {
  x <- as.character(x)
  y <- as.character(y)
  xor(is.na(x), is.na(y)) | (!is.na(x) & !is.na(y) & x != y)
}
row_max_na <- function(dat, columns) {
  if (length(columns) == 0L) return(rep(NA_real_, nrow(dat)))
  mat <- as.matrix(dat[, columns, drop = FALSE])
  apply( mat, 1L, function(x) {
    if (all(is.na(x))) return(NA_real_)
    max(x, na.rm = TRUE)
  } )
}
row_min_na <- function(dat, columns) {
  if (length(columns) == 0L) return(rep(NA_real_, nrow(dat)))
  mat <- as.matrix(dat[, columns, drop = FALSE])
  apply( mat, 1L, function(x) {
    if (all(is.na(x))) return(NA_real_)
    min(x, na.rm = TRUE)
  } )
}
row_sum_na <- function(dat, columns) {
  if (length(columns) == 0L) return(rep(NA_real_, nrow(dat)))
  mat <- as.matrix(dat[, columns, drop = FALSE])
  apply( mat, 1L, function(x) {
    if (all(is.na(x))) return(NA_real_)
    sum(x, na.rm = TRUE)
  } )
}
collapse_ids <- function(x, max_show = 30L) {
  x <- sort(unique(as.integer(x)))
  x <- x[is.finite(x)]
  if (length(x) == 0L) return("")
  shown <- head(x, max_show)
  suffix <- if (length(x) > max_show) {
    paste0(" ... (+", length(x) - max_show, " more)")
  } else {
    ""
  }
  paste0(paste(shown, collapse = ", "), suffix)
}
count_true <- function(x) {
  sum(x %in% TRUE, na.rm = TRUE)
}
# ==============================================================================
# ==============================================================================
# 4. Validate settings
# ==============================================================================
validate_settings <- function() {
  if (!dir.exists(results_dir)) stop("Results directory was not found. Checked: ", paste(results_dir_candidates, collapse = "; "))
  if (length(threshold_values) != 46L) stop("The threshold grid must contain exactly 46 candidate values.")
  if (nrow(threshold_grid_full) != expected_threshold_pairs || expected_threshold_pairs != 1035L) {
    stop("The threshold grid must contain exactly 1,035 ordered pairs with nu1_E > nu2_E.")
  }
  if (any(threshold_grid_full$nu1_E <= threshold_grid_full$nu2_E)) stop("Every threshold pair must satisfy nu1_E > nu2_E.")
  if (any(threshold_grid_full$nu1_E <= 0 | threshold_grid_full$nu1_E >= 1 |
          threshold_grid_full$nu2_E <= 0 | threshold_grid_full$nu2_E >= 1)) {
    stop("All thresholds must lie strictly between 0 and 1.")
  }
  if (!all(mixture_grid$omega_D >= 0 & mixture_grid$omega_D <= 1)) stop("All mixture omega_D values must lie in [0, 1].")
  if (ppv_calibration_omega <= 0 || ppv_calibration_omega >= 1) {
    stop("`ppv_calibration_omega` must lie strictly between 0 and 1.")
  }
  if (ppv_calibration_target <= 0 || ppv_calibration_target >= 1) {
    stop("`ppv_calibration_target` must lie strictly between 0 and 1.")
  }
  if (abs(predictive_reference_omega - ppv_calibration_omega) > numeric_tolerance) {
    stop("`predictive_reference_omega` must equal `ppv_calibration_omega` in this workflow.")
  }
  if (!any(abs(mixture_grid$omega_D - predictive_reference_omega) <= numeric_tolerance)) {
    stop("`mixture_grid` must include the prespecified reference mixture omega_D = ",
         format(predictive_reference_omega, trim = TRUE), ".")
  }
  invisible(TRUE)
}
validate_settings()
# 5. Expected output files
# ==============================================================================
expected_file_grid <- tidyr::crossing(model = names(model_prefixes), N = N_keep) %>%
  dplyr::mutate(prefix = unname(model_prefixes[model])) %>% tidyr::crossing(design_prior_specs) %>% dplyr::mutate(
    expected_file = paste0(prefix, "_", output_tag, "_N", N, ".txt"),
    expected_path = file.path(results_dir, expected_file) ) %>% dplyr::select( model, prefix, N, prior_id,
                                                                               design_prior_type, prior_component, prior_family, omega_D, profile, output_tag, expected_file, expected_path )
# ==============================================================================
# 6. File inventory and loading
# ==============================================================================
inspect_expected_file <- function(file, model, N, prior_id, profile) {
  if (!file.exists(file)) {
    return( tibble::tibble( model = model, N = N, prior_id = prior_id, profile = profile, source_file = basename(file),
                            file_found = FALSE, n_rows_raw = 0L, n_rep_raw = 0L, min_rep = NA_integer_, max_rep = NA_integer_ ) )
  }
  df <- read.table( file = file, header = TRUE, sep = "\t", fill = TRUE, quote = "", comment.char = "",
                    check.names = FALSE, stringsAsFactors = FALSE )
  replicate_column <- if ("rep" %in% names(df)) {
    "rep"
  } else if ("run_ID" %in% names(df)) {
    "run_ID"
  } else {
    stop("Neither `rep` nor `run_ID` was found in: ", basename(file))
  }
  rep_id <- suppressWarnings(as.integer(df[[replicate_column]]))
  tibble::tibble( model = model, N = N, prior_id = prior_id, profile = profile, source_file = basename(file),
                  file_found = TRUE, n_rows_raw = nrow(df), n_rep_raw = length(unique(rep_id[is.finite(rep_id)])),
                  min_rep = safe_min(rep_id), max_rep = safe_max(rep_id) )
}
file_inventory <- purrr::pmap_dfr( expected_file_grid, function( model, prefix, N, prior_id, design_prior_type,
                                                                 prior_component, prior_family, omega_D, profile, output_tag, expected_file, expected_path ) {
  inspect_expected_file( file = expected_path, model = model, N = N, prior_id = prior_id, profile = profile )
} )
read_one_expected_file <- function( file, model_name, prefix, N_expected, prior_id_expected, design_prior_type_expected,
                                    prior_component_expected, prior_family_expected, omega_D_expected, profile_expected, output_tag_expected ) {
  if (!file.exists(file)) return(NULL)
  base_file <- basename(file)
  df <- read.table( file = file, header = TRUE, sep = "\t", fill = TRUE, quote = "", comment.char = "",
                    check.names = FALSE, stringsAsFactors = FALSE )
  if (nrow(df) == 0L) {
    warning("Empty result file: ", base_file)
    return(NULL)
  }
  replicate_column <- if ("rep" %in% names(df)) {
    "rep"
  } else if ("run_ID" %in% names(df)) {
    "run_ID"
  } else {
    stop("Neither `rep` nor `run_ID` was found in: ", base_file)
  }
  df <- df[as.character(df[[replicate_column]]) != replicate_column,, drop = FALSE]
  if (nrow(df) == 0L) {
    warning("No data rows remained in: ", base_file)
    return(NULL)
  }
  df[] <- lapply( df, function(x) {
    type.convert(x, as.is = TRUE)
  } )
  prior_id_saved <- if ("prior_id" %in% names(df)) {
    suppressWarnings(as.integer(df$prior_id))
  } else {
    rep(NA_integer_, nrow(df))
  }
  scenario_id_saved <- if ("scenario_id" %in% names(df)) {
    suppressWarnings(as.integer(df$scenario_id))
  } else {
    rep(NA_integer_, nrow(df))
  }
  prior_label_saved <- if ("prior_label" %in% names(df)) {
    as.character(df$prior_label)
  } else {
    rep(NA_character_, nrow(df))
  }
  design_prior_type_saved <- if ("design_prior_type" %in% names(df)) {
    as.character(df$design_prior_type)
  } else {
    rep(NA_character_, nrow(df))
  }
  prior_family_saved <- if ("prior_family" %in% names(df)) {
    as.character(df$prior_family)
  } else {
    rep(NA_character_, nrow(df))
  }
  omega_D_saved <- if ("omega_D" %in% names(df)) {
    suppressWarnings(as.numeric(df$omega_D))
  } else {
    rep(NA_real_, nrow(df))
  }
  profile_saved <- if ("profile" %in% names(df)) {
    as.character(df$profile)
  } else {
    rep(NA_character_, nrow(df))
  }
  N_planned_saved <- if ("N_planned" %in% names(df)) {
    suppressWarnings(as.integer(df$N_planned))
  } else {
    rep(NA_integer_, nrow(df))
  }
  df$rep_id <- suppressWarnings(as.integer(df[[replicate_column]]))
  df %>% dplyr::mutate( model = model_name, file_prefix = prefix, source_file = base_file,
                        source_row = dplyr::row_number(), output_tag = output_tag_expected, prior_id_saved = prior_id_saved,
                        scenario_id_saved = scenario_id_saved, prior_label_saved = prior_label_saved,
                        design_prior_type_saved = design_prior_type_saved, prior_family_saved = prior_family_saved,
                        omega_D_saved = omega_D_saved, profile_saved = profile_saved, N_planned_saved = N_planned_saved,
                        prior_id_expected = prior_id_expected, design_prior_type_expected = design_prior_type_expected,
                        prior_component_expected = prior_component_expected, prior_family_expected = prior_family_expected,
                        omega_D_expected = omega_D_expected, profile_expected = profile_expected, N_expected = N_expected,
                        prior_id = prior_id_expected, design_prior_type = design_prior_type_expected,
                        prior_component = prior_component_expected, prior_family = prior_family_expected, omega_D = omega_D_expected,
                        profile = profile_expected, N = N_expected, rep_id = as.integer(rep_id) ) %>% dplyr::filter(is.finite(rep_id))
}
load_all_results <- function() {
  dat <- purrr::pmap_dfr( expected_file_grid, function( model, prefix, N, prior_id, design_prior_type, prior_component,
                                                        prior_family, omega_D, profile, output_tag, expected_file, expected_path ) {
    read_one_expected_file( file = expected_path, model_name = model, prefix = prefix, N_expected = N,
                            prior_id_expected = prior_id, design_prior_type_expected = design_prior_type,
                            prior_component_expected = prior_component, prior_family_expected = prior_family, omega_D_expected = omega_D,
                            profile_expected = profile, output_tag_expected = output_tag )
  } )
  if (nrow(dat) == 0L) {
    stop("No result rows were loaded from the current BDCT Prior 1-3 file set.")
  }
  dat
}
# ==============================================================================
# 7. Required saved columns
# ==============================================================================
check_required_columns <- function(dat) {
  required_columns <- c( "rep_id", "prior_id", "design_prior_type", "prior_family", "omega_D", "profile", "N",
                         "prior_label", "component", "component_id", "H_D", "N_planned", "K1", "K2", "n_interim", "n_final",
                         "events_interim", "events_final", "censor_rate_interim", "censor_rate_final", "calendar_cutoff_interim",
                         "calendar_cutoff_final", "tau", "delta_star", interim_prob_col, final_prob_col, rmst_true_col, med_true_col,
                         "Pi_beneficial_med_1", "Pi_beneficial_med_2", "post_mean_1", "post_q025_1", "post_q500_1", "post_q975_1",
                         "post_mean_2", "post_q025_2", "post_q500_2", "post_q975_2", "post_mean_med_1", "post_q025_med_1", "post_q500_med_1",
                         "post_q975_med_1", "post_mean_med_2", "post_q025_med_2", "post_q500_med_2", "post_q975_med_2", "na_rate_med_1",
                         "na_rate_med_2", "cover_rep_1", "bias_rep_1", "sqerr_rep_1", "cover_rep_2", "bias_rep_2", "sqerr_rep_2",
                         "post_eta_PH_trt_mean_1", "post_eta_PH_trt_q025_1", "post_eta_PH_trt_q500_1", "post_eta_PH_trt_q975_1",
                         "post_eta_PH_trt_mean_2", "post_eta_PH_trt_q025_2", "post_eta_PH_trt_q500_2", "post_eta_PH_trt_q975_2",
                         "max_Rhat_1", "min_ESS_1", "n_divergent_1", "max_Rhat_2", "min_ESS_2", "n_divergent_2" )
  assert_columns(dat = dat, columns = required_columns, context = "Loaded simulation results")
  invisible(TRUE)
}
# ==============================================================================
# 8. Replicate and duplicate checks
# ==============================================================================
make_duplicate_check <- function(dat_raw) {
  dat_raw %>% dplyr::group_by(model, N, prior_id, profile, rep_id) %>% dplyr::summarise( n_rows = dplyr::n(),
                                                                                         n_source_files = dplyr::n_distinct(source_file), source_files = paste(sort(unique(source_file)), collapse = "; "),
                                                                                         .groups = "drop" ) %>% dplyr::filter(n_rows > 1L) %>%
    dplyr::arrange(model, N, prior_id, factor(profile, levels = all_profile_levels), rep_id)
}
deduplicate_results <- function(dat_raw) {
  dat_raw %>% dplyr::group_by(model, N, prior_id, profile, rep_id) %>% dplyr::slice_tail(n = 1L) %>% dplyr::ungroup()
}
normalize_rep_ids <- function(x) {
  if (is.null(x) || length(x) == 0L || all(is.na(x))) return(integer(0))
  as.integer(x[is.finite(x)])
}
make_replicate_id_check <- function(dat) {
  expected_grid <- expected_file_grid %>% dplyr::select(model, N, prior_id, design_prior_type, profile) %>%
    dplyr::distinct()
  observed <- dat %>% dplyr::group_by(model, N, prior_id, profile) %>% dplyr::summarise(
    observed_ids = list(sort(unique(rep_id))), n_rep = dplyr::n_distinct(rep_id), min_rep = safe_min(rep_id),
    max_rep = safe_max(rep_id), .groups = "drop" )
  expected_grid %>% dplyr::left_join(observed, by = c("model", "N", "prior_id", "profile")) %>% dplyr::mutate(
    observed_ids = purrr::map(observed_ids, normalize_rep_ids), n_rep = tidyr::replace_na(n_rep, 0L),
    missing_ids = purrr::map(observed_ids, function(x) setdiff(expected_rep_ids, x)),
    extra_ids = purrr::map(observed_ids, function(x) setdiff(x, expected_rep_ids)),
    n_missing = purrr::map_int(missing_ids, length), n_extra = purrr::map_int(extra_ids, length),
    missing_rep_ids = purrr::map_chr(missing_ids, collapse_ids),
    extra_rep_ids = purrr::map_chr(extra_ids, collapse_ids), completion_pct = 100 * n_rep / n_rep_expected,
    Complete = n_missing == 0L & n_extra == 0L & n_rep == n_rep_expected ) %>% dplyr::select(
      model, N, prior_id, design_prior_type, profile, n_rep, min_rep, max_rep,
      n_missing, n_extra, completion_pct, Complete, missing_rep_ids, extra_rep_ids ) %>%
    dplyr::arrange(model, N, prior_id, factor(profile, levels = all_profile_levels))
}
# ==============================================================================
# 9. Metadata and row-level integrity checks
# ==============================================================================
make_metadata_check <- function(dat) {
  dat %>% dplyr::mutate(
    # Existing files were generated before final BDCT renumbering:
    #   old 3 -> new 1, old 2 -> new 2, old 7 -> new 3.
    # Files generated after renumbering may already contain the new IDs.
    prior_id_mismatch = !is.na(prior_id_saved) & !dplyr::case_when(
      prior_id_expected == 1L ~ prior_id_saved %in% c(3L, 1L), prior_id_expected == 2L ~ prior_id_saved %in% c(2L),
      prior_id_expected == 3L ~ prior_id_saved %in% c(7L, 3L), TRUE ~ FALSE ), scenario_id_mismatch =
      !is.na(scenario_id_saved) & !dplyr::case_when( prior_id_expected == 1L ~ scenario_id_saved %in% c(3L, 1L),
                                                     prior_id_expected == 2L ~ scenario_id_saved %in% c(2L),
                                                     prior_id_expected == 3L ~ scenario_id_saved %in% c(7L, 3L), TRUE ~ FALSE ), prior_label_mismatch =
      !is.na(prior_label_saved) & !dplyr::case_when( prior_id_expected == 1L ~
                                                       prior_label_saved %in% c("exact_null_point_mass"), prior_id_expected == 2L ~
                                                       prior_label_saved %in% c("no_benefit_continuous"), prior_id_expected == 3L ~ prior_label_saved %in% c(
                                                         "benefit_positive_truncnorm", "benefit_positive_truncnorm_sd010" ), TRUE ~ FALSE ),
    design_prior_type_mismatch = !is.na(design_prior_type_saved) & !dplyr::case_when( prior_id_expected == 1L ~
                                                                                        design_prior_type_saved %in% c("exact_null_point_mass"), prior_id_expected == 2L ~
                                                                                        design_prior_type_saved %in% c("no_benefit_continuous"), prior_id_expected == 3L ~
                                                                                        design_prior_type_saved %in% c( "benefit_positive_truncnorm", "benefit_positive_truncnorm_sd010" ),
                                                                                      TRUE ~ FALSE ), prior_family_mismatch = !is.na(prior_family_saved) & prior_family_saved !=
      prior_family_expected, omega_D_mismatch = !is.na(omega_D_saved) & abs(omega_D_saved - omega_D_expected) >
      numeric_tolerance, profile_mismatch = !is.na(profile_saved) & profile_saved != profile_expected, N_mismatch =
      !is.na(N_planned_saved) & N_planned_saved != N_expected, any_metadata_mismatch = prior_id_mismatch |
      scenario_id_mismatch | prior_label_mismatch | design_prior_type_mismatch | prior_family_mismatch |
      omega_D_mismatch | profile_mismatch | N_mismatch ) %>% dplyr::filter(any_metadata_mismatch) %>% dplyr::select(
        model, N, prior_id, profile, rep_id, source_file, dplyr::ends_with("_mismatch") ) %>% dplyr::arrange( model, N,
                                                                                                              prior_id, factor(profile, levels = all_profile_levels), rep_id )
}
make_integrity_detail <- function(dat) {
  expected_K2 <- as.integer(round(dat$N * (1 - target_censor)))
  expected_K1 <- as.integer(floor(0.6 * dat$K2))
  expected_N_from_K2 <- as.integer(ceiling(dat$K2 / (1 - target_censor)))
  rmst_cover_1_expected <- as.integer( dat$post_q025_1 <= dat[[rmst_true_col]] & dat[[rmst_true_col]] <= dat$post_q975_1
  )
  rmst_cover_2_expected <- as.integer( dat$post_q025_2 <= dat[[rmst_true_col]] & dat[[rmst_true_col]] <= dat$post_q975_2
  )
  rmst_bias_1_expected <- dat$post_mean_1 - dat[[rmst_true_col]]
  rmst_bias_2_expected <- dat$post_mean_2 - dat[[rmst_true_col]]
  dat %>% dplyr::mutate( design_N_planned_mismatch = N_planned != N, design_K2_mismatch = K2 != expected_K2,
                         design_K1_mismatch = K1 != expected_K1, design_N_from_K2_mismatch = N_planned != expected_N_from_K2,
                         design_n_interim_mismatch = n_interim != N_planned, design_n_final_mismatch = n_final != N_planned,
                         design_interim_events_mismatch = events_interim != K1, design_final_events_too_small = events_final < K2,
                         design_final_cutoff_too_early = calendar_cutoff_final < tau, design_tau_mismatch = abs(tau - tau_expected) >
                           numeric_tolerance, rmst_probability_out_of_range = is.na(.data[[interim_prob_col]]) |
                           is.na(.data[[final_prob_col]]) | .data[[interim_prob_col]] < -numeric_tolerance | .data[[interim_prob_col]] >
                           1 + numeric_tolerance | .data[[final_prob_col]] < -numeric_tolerance | .data[[final_prob_col]] > 1 +
                           numeric_tolerance, median_probability_out_of_range = ( !is.na(Pi_beneficial_med_1) &
                                                                                    (Pi_beneficial_med_1 < -numeric_tolerance | Pi_beneficial_med_1 > 1 + numeric_tolerance) ) | (
                                                                                      !is.na(Pi_beneficial_med_2) &
                                                                                        (Pi_beneficial_med_2 < -numeric_tolerance | Pi_beneficial_med_2 > 1 + numeric_tolerance) ),
                         rmst_interval_order_error_1 = post_q025_1 > post_q500_1 | post_q500_1 > post_q975_1, rmst_interval_order_error_2 =
                           post_q025_2 > post_q500_2 | post_q500_2 > post_q975_2, median_interval_order_error_1 =
                           (!is.na(post_q025_med_1) & !is.na(post_q500_med_1) & post_q025_med_1 > post_q500_med_1) |
                           (!is.na(post_q500_med_1) & !is.na(post_q975_med_1) & post_q500_med_1 > post_q975_med_1),
                         median_interval_order_error_2 =
                           (!is.na(post_q025_med_2) & !is.na(post_q500_med_2) & post_q025_med_2 > post_q500_med_2) |
                           (!is.na(post_q500_med_2) & !is.na(post_q975_med_2) & post_q500_med_2 > post_q975_med_2), H_D_prior_mismatch =
                           H_D != as.integer(prior_id == benefit_prior_id), component_id_prior_mismatch = component_id !=
                           as.integer(prior_id == benefit_prior_id), exact_null_truth_mismatch = prior_id == exact_null_prior_id &
                           abs(.data[[rmst_true_col]]) > 1e-8, benefit_truth_mismatch = prior_id == benefit_prior_id &
                           .data[[rmst_true_col]] <= delta_star, no_benefit_truth_mismatch = prior_id == no_benefit_prior_id &
                           .data[[rmst_true_col]] > delta_star, saved_cover_mismatch_1 =
                           numeric_mismatch(cover_rep_1, rmst_cover_1_expected, tolerance = 0), saved_cover_mismatch_2 =
                           numeric_mismatch(cover_rep_2, rmst_cover_2_expected, tolerance = 0), saved_bias_mismatch_1 =
                           numeric_mismatch(bias_rep_1, rmst_bias_1_expected), saved_bias_mismatch_2 =
                           numeric_mismatch(bias_rep_2, rmst_bias_2_expected), saved_sqerr_mismatch_1 =
                           numeric_mismatch(sqerr_rep_1, rmst_bias_1_expected^2), saved_sqerr_mismatch_2 =
                           numeric_mismatch(sqerr_rep_2, rmst_bias_2_expected^2) ) %>% dplyr::mutate( any_design_error =
                                                                                                        design_N_planned_mismatch | design_K2_mismatch | design_K1_mismatch | design_N_from_K2_mismatch |
                                                                                                        design_n_interim_mismatch | design_n_final_mismatch | design_interim_events_mismatch |
                                                                                                        design_final_events_too_small | design_final_cutoff_too_early | design_tau_mismatch, any_probability_error =
                                                                                                        rmst_probability_out_of_range | median_probability_out_of_range, any_interval_error =
                                                                                                        rmst_interval_order_error_1 | rmst_interval_order_error_2 | median_interval_order_error_1 |
                                                                                                        median_interval_order_error_2, any_truth_error = H_D_prior_mismatch | component_id_prior_mismatch |
                                                                                                        exact_null_truth_mismatch | benefit_truth_mismatch | no_benefit_truth_mismatch, any_saved_summary_error =
                                                                                                        saved_cover_mismatch_1 | saved_cover_mismatch_2 | saved_bias_mismatch_1 | saved_bias_mismatch_2 |
                                                                                                        saved_sqerr_mismatch_1 | saved_sqerr_mismatch_2, any_integrity_error = any_design_error |
                                                                                                        any_probability_error | any_interval_error | any_truth_error | any_saved_summary_error ) %>%
    dplyr::filter(any_integrity_error) %>% dplyr::select( model, N, prior_id, design_prior_type, profile, rep_id,
                                                          source_file, any_design_error, any_probability_error, any_interval_error, any_truth_error,
                                                          any_saved_summary_error, dplyr::starts_with("design_"), dplyr::contains("probability_out_of_range"),
                                                          dplyr::contains("interval_order_error"), dplyr::contains("truth_mismatch"), H_D_prior_mismatch,
                                                          component_id_prior_mismatch, dplyr::starts_with("saved_") ) %>%
    dplyr::arrange(model, N, prior_id, factor(profile, levels = all_profile_levels), rep_id)
}
make_integrity_summary <- function(metadata_detail, integrity_detail) {
  tibble::tibble( Check = c( "Metadata mismatches", "Any design error", "Any probability error",
                             "Any interval-order error", "Any design-prior truth error", "Any saved RMST summary error",
                             "Any row-level integrity error" ), n_flagged_rows = c( nrow(metadata_detail), if (nrow(integrity_detail) > 0L) {
                               count_true(integrity_detail$any_design_error)
                             } else {
                               0L
                             }, if (nrow(integrity_detail) > 0L) {
                               count_true(integrity_detail$any_probability_error)
                             } else {
                               0L
                             }, if (nrow(integrity_detail) > 0L) {
                               count_true(integrity_detail$any_interval_error)
                             } else {
                               0L
                             }, if (nrow(integrity_detail) > 0L) {
                               count_true(integrity_detail$any_truth_error)
                             } else {
                               0L
                             }, if (nrow(integrity_detail) > 0L) {
                               count_true(integrity_detail$any_saved_summary_error)
                             } else {
                               0L
                             }, nrow(integrity_detail) ) )
}
# ==============================================================================
# 10. Cross-model shared-data consistency
# ==============================================================================
make_cross_model_check <- function(dat) {
  key_columns <- c("N", "prior_id", "profile", "rep_id")
  shared_columns <- c( "design_prior_type", "prior_family", "omega_D", "component", "component_id", "H_D",
                       rmst_true_col, med_true_col, "N_planned", "K1", "K2", "n_interim", "n_final", "events_interim", "events_final",
                       "censor_rate_interim", "censor_rate_final", "calendar_cutoff_interim", "calendar_cutoff_final", "tau", "delta_star"
  )
  assert_columns( dat = dat, columns = c(key_columns, shared_columns), context = "Cross-model shared-data check" )
  gp <- dat %>% dplyr::filter(model == "Original GP-PH") %>%
    dplyr::select(dplyr::all_of(c(key_columns, shared_columns)))
  indep <- dat %>% dplyr::filter(model == "Independent GP") %>%
    dplyr::select(dplyr::all_of(c(key_columns, shared_columns)))
  paired <- dplyr::full_join(gp, indep, by = key_columns, suffix = c("_gp", "_indep"))
  paired %>% dplyr::mutate( missing_gp = is.na(N_planned_gp), missing_indep = is.na(N_planned_indep),
                            design_prior_type_mismatch = character_mismatch(design_prior_type_gp, design_prior_type_indep),
                            prior_family_mismatch = character_mismatch(prior_family_gp, prior_family_indep), omega_D_mismatch =
                              numeric_mismatch(omega_D_gp, omega_D_indep), component_mismatch =
                              character_mismatch(component_gp, component_indep), component_id_mismatch =
                              numeric_mismatch(component_id_gp, component_id_indep, tolerance = 0), H_D_mismatch =
                              numeric_mismatch(H_D_gp, H_D_indep, tolerance = 0), true_rmst_mismatch =
                              numeric_mismatch(.data[[paste0(rmst_true_col, "_gp")]], .data[[paste0(rmst_true_col, "_indep")]]),
                            true_median_mismatch =
                              numeric_mismatch(.data[[paste0(med_true_col, "_gp")]], .data[[paste0(med_true_col, "_indep")]]),
                            N_planned_mismatch = numeric_mismatch(N_planned_gp, N_planned_indep, tolerance = 0), K1_mismatch =
                              numeric_mismatch(K1_gp, K1_indep, tolerance = 0), K2_mismatch =
                              numeric_mismatch(K2_gp, K2_indep, tolerance = 0), n_interim_mismatch =
                              numeric_mismatch(n_interim_gp, n_interim_indep, tolerance = 0), n_final_mismatch =
                              numeric_mismatch(n_final_gp, n_final_indep, tolerance = 0), events_interim_mismatch =
                              numeric_mismatch(events_interim_gp, events_interim_indep, tolerance = 0), events_final_mismatch =
                              numeric_mismatch(events_final_gp, events_final_indep, tolerance = 0), censor_rate_interim_mismatch =
                              numeric_mismatch(censor_rate_interim_gp, censor_rate_interim_indep), censor_rate_final_mismatch =
                              numeric_mismatch(censor_rate_final_gp, censor_rate_final_indep), calendar_cutoff_interim_mismatch =
                              numeric_mismatch(calendar_cutoff_interim_gp, calendar_cutoff_interim_indep), calendar_cutoff_final_mismatch =
                              numeric_mismatch(calendar_cutoff_final_gp, calendar_cutoff_final_indep), tau_mismatch =
                              numeric_mismatch(tau_gp, tau_indep), delta_star_mismatch = numeric_mismatch(delta_star_gp, delta_star_indep),
                            any_shared_data_mismatch = missing_gp | missing_indep | design_prior_type_mismatch | prior_family_mismatch |
                              omega_D_mismatch | component_mismatch | component_id_mismatch | H_D_mismatch | true_rmst_mismatch |
                              true_median_mismatch | N_planned_mismatch | K1_mismatch | K2_mismatch | n_interim_mismatch | n_final_mismatch |
                              events_interim_mismatch | events_final_mismatch | censor_rate_interim_mismatch | censor_rate_final_mismatch |
                              calendar_cutoff_interim_mismatch | calendar_cutoff_final_mismatch | tau_mismatch | delta_star_mismatch ) %>%
    dplyr::filter(any_shared_data_mismatch) %>%
    dplyr::select(dplyr::all_of(key_columns), missing_gp, missing_indep, dplyr::ends_with("_mismatch")) %>%
    dplyr::arrange(N, prior_id, factor(profile, levels = all_profile_levels), rep_id)
}
make_cross_model_summary <- function(cross_model_detail, dat) {
  n_unique_keys <- dat %>% dplyr::distinct(N, prior_id, profile, rep_id) %>% nrow()
  tibble::tibble( Check = c("Unique shared-data replicate keys", "Replicate keys with a GP-PH/Independent GP mismatch"),
                  Count = c(n_unique_keys, nrow(cross_model_detail)) )
}
# ==============================================================================
# 11. MCMC diagnostics
# ==============================================================================
make_mcmc_diagnostics <- function(dat) {
  rhat_columns <- intersect(c("max_Rhat_1", "max_Rhat_2"), names(dat))
  ess_columns <- intersect(c("min_ESS_1", "min_ESS_2"), names(dat))
  divergent_columns <- intersect(c("n_divergent_1", "n_divergent_2"), names(dat))
  tmp <- dat %>% dplyr::mutate( max_Rhat_any = row_max_na(., rhat_columns), min_ESS_any = row_min_na(., ess_columns),
                                n_divergent_any = row_sum_na(., divergent_columns), bad_Rhat = max_Rhat_any > rhat_warn, bad_ESS = min_ESS_any <
                                  ess_warn, any_divergence = n_divergent_any > 0 )
  tmp %>% dplyr::group_by(model, N, prior_id, design_prior_type, profile) %>% dplyr::summarise( n_rep = dplyr::n(),
                                                                                                max_Rhat = safe_max(max_Rhat_any), min_ESS = safe_min(min_ESS_any), total_divergent = safe_sum(n_divergent_any),
                                                                                                pct_bad_Rhat = pct(safe_mean(bad_Rhat)), pct_bad_ESS = pct(safe_mean(bad_ESS)),
                                                                                                pct_any_divergence = pct(safe_mean(any_divergence)), .groups = "drop" ) %>%
    dplyr::arrange(model, N, prior_id, factor(profile, levels = all_profile_levels))
}
# ==============================================================================
# 12. Median-survival posterior availability
# ==============================================================================
make_median_availability <- function(dat) {
  dat %>% dplyr::group_by(model, N, prior_id, design_prior_type, profile) %>% dplyr::summarise( n_rep = dplyr::n(),
                                                                                                pct_true_median_missing = pct(safe_mean(is.na(.data[[med_true_col]]))),
                                                                                                pct_Pi_median_1_missing = pct(safe_mean(is.na(Pi_beneficial_med_1))),
                                                                                                pct_Pi_median_2_missing = pct(safe_mean(is.na(Pi_beneficial_med_2))),
                                                                                                mean_draw_NA_rate_1 = pct(safe_mean(na_rate_med_1)), max_draw_NA_rate_1 = pct(safe_max(na_rate_med_1)),
                                                                                                mean_draw_NA_rate_2 = pct(safe_mean(na_rate_med_2)), max_draw_NA_rate_2 = pct(safe_max(na_rate_med_2)),
                                                                                                .groups = "drop" ) %>% dplyr::mutate( Warning = pct_true_median_missing > 0 | pct_Pi_median_1_missing > 0 |
                                                                                                                                        pct_Pi_median_2_missing > 0 | max_draw_NA_rate_1 > median_na_warn_pct | max_draw_NA_rate_2 > median_na_warn_pct
                                                                                                ) %>% dplyr::arrange(model, N, prior_id, factor(profile, levels = all_profile_levels))
}
# ==============================================================================
# 13. Estimand-specific columns
# ==============================================================================
get_estimand_columns <- function(estimand) {
  if (identical(estimand, "RMST")) {
    return( list( pi1 = interim_prob_col, pi2 = final_prob_col, true = rmst_true_col, mean1 = "post_mean_1",
                  mean2 = "post_mean_2", q025_2 = "post_q025_2", q975_2 = "post_q975_2" ) )
  }
  if (identical(estimand, "Median")) {
    return( list( pi1 = "Pi_beneficial_med_1", pi2 = "Pi_beneficial_med_2", true = med_true_col,
                  mean1 = "post_mean_med_1", mean2 = "post_mean_med_2", q025_2 = "post_q025_med_2", q975_2 = "post_q975_med_2" ) )
  }
  stop("Unknown estimand: ", estimand)
}
# ==============================================================================
# 14. Component operating characteristics
# ==============================================================================
make_component_oc_one <- function(dat, estimand, nu1_E, nu2_E) {
  estimand_columns <- get_estimand_columns(estimand)
  
  dat %>%
    dplyr::mutate(
      pi_interim = .data[[estimand_columns$pi1]],
      pi_final = .data[[estimand_columns$pi2]],
      I = pi_interim >= nu1_E,
      F_cf = pi_final >= nu2_E,
      S = I | (!I & F_cf)
    ) %>%
    dplyr::group_by(model, N, prior_id, design_prior_type, prior_family, profile) %>%
    dplyr::summarise(
      n_rep = dplyr::n(),
      n_rep_OC = sum(!is.na(S)),
      ESP = pct(safe_mean(I)),
      FSP = pct(safe_mean(!I & F_cf)),
      Success = pct(safe_mean(S)),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      BP = dplyr::if_else(prior_id == benefit_prior_id, Success, NA_real_),
      FPR = dplyr::if_else(
        prior_id %in% c(no_benefit_prior_id, exact_null_prior_id),
        Success,
        NA_real_
      ),
      Estimand = estimand,
      nu1_E = nu1_E,
      nu2_E = nu2_E,
      success_identity_error = abs(Success - (ESP + FSP)) > 1e-8
    ) %>%
    dplyr::select(
      Estimand, model, N, prior_id, design_prior_type, prior_family, profile,
      nu1_E, nu2_E, n_rep, n_rep_OC, ESP, FSP, BP, FPR, Success,
      success_identity_error
    )
}
make_component_oc_grid <- function(dat, estimand, threshold_grid) {
  purrr::pmap_dfr( threshold_grid, function(nu1_E, nu2_E) {
    make_component_oc_one(dat = dat, estimand = estimand, nu1_E = nu1_E, nu2_E = nu2_E)
  } )
}
# ==============================================================================
# ==============================================================================
# 15A. FPR-only / FPR-first, profile-specific threshold calibration
# ==============================================================================
make_fpr_boundary_summary <- function(exact_null_surface) {
  out <- exact_null_surface %>%
    dplyr::filter(
      is.finite(FPR_exact),
      FPR_exact <= fpr_target_pct + numeric_tolerance
    ) %>%
    dplyr::group_by(Estimand, model, N) %>%
    dplyr::summarise(
      exact_2_5_available =
        any(abs(FPR_exact - fpr_target_pct) <= numeric_tolerance),
      target_FPR = max(FPR_exact, na.rm = TRUE),
      n_pairs_at_target = sum(
        abs(FPR_exact - max(FPR_exact, na.rm = TRUE)) <= numeric_tolerance
      ),
      distance_below_2_5 =
        fpr_target_pct - max(FPR_exact, na.rm = TRUE),
      fallback_used =
        max(FPR_exact, na.rm = TRUE) <
        fpr_target_pct - numeric_tolerance,
      .groups = "drop"
    ) %>%
    dplyr::arrange(
      factor(Estimand, levels = c("RMST", "Median")),
      factor(model, levels = names(model_prefixes)),
      N
    )
  
  missing_settings <- calibration_index %>%
    dplyr::anti_join(
      out,
      by = c("Estimand", "model", "N")
    )
  
  if (nrow(missing_settings) > 0L) {
    stop(
      "No threshold pair satisfied exact-null FPR <= ",
      fpr_target_pct,
      "% for FPR-only calibration: ",
      paste0(
        missing_settings$Estimand, " / ",
        missing_settings$model, " / N=",
        missing_settings$N,
        collapse = "; "
      )
    )
  }
  
  if (any(
    out$exact_2_5_available &
    abs(out$target_FPR - fpr_target_pct) > numeric_tolerance
  )) {
    stop(
      "An exact 2.5% FPR pair exists but the FPR-only target was not set to 2.5%."
    )
  }
  
  out
}

select_profile_specific_thresholds_fpr_only <- function(
    exact_null_surface,
    benefit_power_surface
) {
  fpr_boundary <- make_fpr_boundary_summary(exact_null_surface)
  
  target_fpr_pairs <- exact_null_surface %>%
    dplyr::inner_join(
      fpr_boundary,
      by = c("Estimand", "model", "N")
    ) %>%
    dplyr::filter(
      abs(FPR_exact - target_FPR) <= numeric_tolerance
    ) %>%
    dplyr::arrange(
      factor(Estimand, levels = c("RMST", "Median")),
      factor(model, levels = names(model_prefixes)),
      N,
      grid_order
    )
  
  target_pair_count <- target_fpr_pairs %>%
    dplyr::count(
      Estimand, model, N,
      name = "n_target_pairs"
    ) %>%
    dplyr::left_join(
      fpr_boundary %>%
        dplyr::select(
          Estimand, model, N,
          n_pairs_at_target
        ),
      by = c("Estimand", "model", "N")
    )
  
  if (any(
    target_pair_count$n_target_pairs !=
    target_pair_count$n_pairs_at_target
  )) {
    stop(
      "FPR-only target-FPR threshold-pair counts are inconsistent."
    )
  }
  
  candidates <- benefit_power_surface %>%
    dplyr::inner_join(
      target_fpr_pairs %>%
        dplyr::select(
          Estimand, model, N,
          nu1_E, nu2_E,
          grid_order,
          FPR_exact,
          target_FPR
        ),
      by = c(
        "Estimand", "model", "N",
        "nu1_E", "nu2_E"
      )
    ) %>%
    dplyr::mutate(
      profile = factor(
        profile,
        levels = benefit_profile_levels
      )
    ) %>%
    dplyr::group_by(
      Estimand, model, N, profile
    ) %>%
    dplyr::arrange(
      dplyr::desc(BP),
      grid_order,
      .by_group = TRUE
    ) %>%
    dplyr::mutate(
      candidate_rank = dplyr::row_number(),
      best_BP = dplyr::first(BP),
      n_BP_ties = sum(
        abs(BP - dplyr::first(BP)) <= numeric_tolerance
      )
    ) %>%
    dplyr::ungroup()
  
  selected <- candidates %>%
    dplyr::filter(candidate_rank == 1L) %>%
    dplyr::transmute(
      Estimand,
      model,
      N,
      profile = as.character(profile),
      nu1_E,
      nu2_E,
      FPR_exact,
      BP,
      ESP,
      FSP,
      target_FPR,
      Pass_FPR =
        FPR_exact <= fpr_target_pct + numeric_tolerance,
      Pass_BP =
        BP >= bp_target_pct - numeric_tolerance,
      grid_order,
      n_BP_ties
    ) %>%
    dplyr::arrange(
      factor(Estimand, levels = c("RMST", "Median")),
      factor(model, levels = names(model_prefixes)),
      N,
      factor(profile, levels = benefit_profile_levels)
    )
  
  expected_settings <- tidyr::expand_grid(
    Estimand = c("RMST", "Median"),
    model = names(model_prefixes),
    N = N_keep,
    profile = benefit_profile_levels
  )
  
  missing_selected <- expected_settings %>%
    dplyr::anti_join(
      selected,
      by = c("Estimand", "model", "N", "profile")
    )
  
  if (nrow(missing_selected) > 0L ||
      nrow(selected) != 24L) {
    stop(
      "FPR-only calibration must select exactly 24 threshold pairs."
    )
  }
  
  adequacy <- selected %>%
    dplyr::transmute(
      Estimand,
      model,
      N,
      profile,
      nu1_E,
      nu2_E,
      FPR_exact,
      BP,
      Pass_FPR,
      Pass_BP,
      Adequacy_status = dplyr::if_else(
        Pass_BP,
        "BP target met",
        "BP target not met"
      )
    )
  
  list(
    fpr_boundary = fpr_boundary,
    target_fpr_pairs = target_fpr_pairs,
    candidates = candidates,
    selected = selected,
    adequacy = adequacy
  )
}

make_bp_tie_summary_fpr_only <- function(threshold_candidates) {
  threshold_candidates %>%
    dplyr::filter(
      abs(BP - best_BP) <= numeric_tolerance
    ) %>%
    dplyr::group_by(
      Estimand, model, N, profile
    ) %>%
    dplyr::summarise(
      selected_FPR = dplyr::first(FPR_exact),
      best_BP = dplyr::first(best_BP),
      n_BP_ties = dplyr::n(),
      .groups = "drop"
    )
}

# ==============================================================================
# 15B. PPV + FPR-constrained, profile-specific threshold calibration
# ==============================================================================
make_exact_null_surface <- function(oc_grid) {
  oc_grid %>%
    dplyr::filter(prior_id == exact_null_prior_id) %>%
    dplyr::group_by(Estimand, model, N, nu1_E, nu2_E) %>%
    dplyr::summarise(
      FPR_exact = safe_mean(FPR),
      ESP_exact = safe_mean(ESP),
      FSP_exact = safe_mean(FSP),
      .groups = "drop"
    ) %>%
    dplyr::left_join(
      threshold_grid_selection_order,
      by = c("nu1_E", "nu2_E")
    ) %>%
    dplyr::arrange(
      factor(Estimand, levels = c("RMST", "Median")),
      factor(model, levels = names(model_prefixes)),
      N,
      grid_order
    )
}

make_benefit_power_surface <- function(oc_grid) {
  oc_grid %>%
    dplyr::filter(prior_id == benefit_prior_id) %>%
    dplyr::select(
      Estimand, model, N, profile,
      nu1_E, nu2_E, BP, ESP, FSP
    )
}

make_ppv_calibration_surface <- function(
    exact_null_surface,
    benefit_power_surface
) {
  out <- benefit_power_surface %>%
    dplyr::left_join(
      exact_null_surface %>%
        dplyr::select(
          Estimand, model, N,
          nu1_E, nu2_E,
          FPR_exact, grid_order
        ),
      by = c("Estimand", "model", "N", "nu1_E", "nu2_E")
    ) %>%
    dplyr::mutate(
      BP_fraction = BP / 100,
      FPR_fraction = FPR_exact / 100,
      PPV_calibration_fraction = safe_ratio(
        ppv_calibration_omega * BP_fraction,
        ppv_calibration_omega * BP_fraction +
          (1 - ppv_calibration_omega) * FPR_fraction
      ),
      PPV_calibration = pct(PPV_calibration_fraction),
      
      # Algebraically equivalent candidate-specific FPR upper bound induced
      # by PPV >= ppv_calibration_target:
      #
      # FPR <= [omega/(1-omega)] * BP * [(1-PPV_target)/PPV_target].
      #
      # BP and FPR are both on the percentage scale here, so the same formula
      # can be applied directly to BP.
      FPR_limit_from_PPV =
        (ppv_calibration_omega / (1 - ppv_calibration_omega)) *
        BP *
        ((1 - ppv_calibration_target) / ppv_calibration_target),
      
      Pass_FPR = is.finite(FPR_exact) &
        FPR_exact <= fpr_target_pct + numeric_tolerance,
      Pass_PPV = is.finite(PPV_calibration) &
        PPV_calibration >= ppv_target_pct - numeric_tolerance,
      Pass_PPV_equivalent = is.finite(FPR_exact) &
        is.finite(FPR_limit_from_PPV) &
        FPR_exact <= FPR_limit_from_PPV + numeric_tolerance,
      Pass_BP = is.finite(BP) &
        BP >= bp_target_pct - numeric_tolerance,
      Pass_reliability = Pass_FPR & Pass_PPV,
      Pass_all = Pass_FPR & Pass_PPV & Pass_BP
    ) %>%
    dplyr::mutate(
      profile = factor(profile, levels = benefit_profile_levels)
    ) %>%
    dplyr::arrange(
      factor(Estimand, levels = c("RMST", "Median")),
      factor(model, levels = names(model_prefixes)),
      N,
      profile,
      grid_order
    )
  
  ppv_equivalence_error <- out %>%
    dplyr::filter(
      is.finite(PPV_calibration),
      is.finite(FPR_limit_from_PPV),
      Pass_PPV != Pass_PPV_equivalent
    )
  
  if (nrow(ppv_equivalence_error) > 0L) {
    stop(
      "Direct PPV constraint and its algebraically equivalent FPR bound ",
      "do not agree for ", nrow(ppv_equivalence_error), " candidate rows."
    )
  }
  
  out
}

select_profile_specific_thresholds_ppv_fpr <- function(calibration_surface) {
  expected_settings <- tidyr::expand_grid(
    Estimand = c("RMST", "Median"),
    model = names(model_prefixes),
    N = N_keep,
    profile = benefit_profile_levels
  )
  
  # First impose the two reliability/protection constraints that define the
  # admissible candidate set:
  #   FPR <= 2.5%
  #   PPV(omega) >= 97.5%
  #
  # Then maximize BP. If the maximized BP is >= 90%, the exact three-constraint
  # problem is feasible. If the maximized BP is < 90%, no threshold pair can
  # satisfy all three requirements for that setting, but the best admissible
  # rule is retained for transparent diagnostics and manuscript comparison.
  reliability_candidates <- calibration_surface %>%
    dplyr::filter(Pass_reliability)
  
  reliability_counts <- reliability_candidates %>%
    dplyr::count(
      Estimand, model, N, profile,
      name = "n_reliability_feasible"
    )
  
  missing_reliability <- expected_settings %>%
    dplyr::anti_join(
      reliability_counts,
      by = c("Estimand", "model", "N", "profile")
    )
  
  if (nrow(missing_reliability) > 0L) {
    stop(
      "No threshold pair satisfies both FPR <= ",
      fpr_target_pct,
      "% and PPV(omega=",
      format(ppv_calibration_omega, digits = 6, trim = TRUE),
      ") >= ",
      ppv_target_pct,
      "% for: ",
      paste0(
        missing_reliability$Estimand, " / ",
        missing_reliability$model, " / N=",
        missing_reliability$N, " / ",
        missing_reliability$profile,
        collapse = "; "
      )
    )
  }
  
  ranked_candidates <- reliability_candidates %>%
    dplyr::group_by(Estimand, model, N, profile) %>%
    dplyr::arrange(
      dplyr::desc(BP),
      grid_order,
      .by_group = TRUE
    ) %>%
    dplyr::mutate(
      candidate_rank = dplyr::row_number(),
      best_BP = dplyr::first(BP),
      n_BP_ties = sum(
        abs(BP - dplyr::first(BP)) <= numeric_tolerance
      )
    ) %>%
    dplyr::ungroup()
  
  selected <- ranked_candidates %>%
    dplyr::filter(candidate_rank == 1L) %>%
    dplyr::transmute(
      Estimand,
      model,
      N,
      profile = as.character(profile),
      nu1_E,
      nu2_E,
      FPR_exact,
      BP,
      ESP,
      FSP,
      PPV_calibration,
      FPR_limit_from_PPV,
      Pass_FPR,
      Pass_PPV,
      Pass_BP,
      Pass_all,
      grid_order,
      n_BP_ties
    ) %>%
    dplyr::arrange(
      factor(Estimand, levels = c("RMST", "Median")),
      factor(model, levels = names(model_prefixes)),
      N,
      factor(profile, levels = benefit_profile_levels)
    )
  
  missing_selected <- expected_settings %>%
    dplyr::anti_join(
      selected,
      by = c("Estimand", "model", "N", "profile")
    )
  
  if (nrow(missing_selected) > 0L || nrow(selected) != 24L) {
    stop(
      "Profile-specific calibration must return exactly 24 selected ",
      "reliability-admissible threshold pairs."
    )
  }
  
  feasibility_summary <- selected %>%
    dplyr::transmute(
      Estimand,
      model,
      N,
      profile,
      nu1_E,
      nu2_E,
      BP,
      FPR_exact,
      PPV_calibration,
      FPR_limit_from_PPV,
      Pass_BP,
      Pass_FPR,
      Pass_PPV,
      Pass_all,
      Full_constraint_status = dplyr::if_else(
        Pass_all,
        "Feasible",
        "Infeasible: max BP < target"
      )
    )
  
  list(
    candidates = ranked_candidates,
    selected = selected,
    feasibility = feasibility_summary
  )
}

# 16. Apply each selected profile-specific threshold to all required components
# ==============================================================================
make_selected_component_oc <- function(dat, selected_thresholds) {
  purrr::pmap_dfr(
    selected_thresholds %>% dplyr::select(Estimand, model, N, calibration_profile = profile, nu1_E, nu2_E),
    function(Estimand, model, N, calibration_profile, nu1_E, nu2_E) {
      estimand_i <- Estimand; model_i <- model; N_i <- N; profile_i <- calibration_profile
      tmp <- make_component_oc_one(
        dat = dat %>% dplyr::filter(.data$model == .env$model_i, .data$N == .env$N_i),
        estimand = estimand_i, nu1_E = nu1_E, nu2_E = nu2_E
      )
      tmp %>%
        dplyr::filter(prior_id %in% c(exact_null_prior_id, no_benefit_prior_id) |
                        (prior_id == benefit_prior_id & profile == profile_i)) %>%
        dplyr::mutate(calibration_profile = profile_i)
    }
  )
}

# ==============================================================================
# 17. Calibration tie summary for internal checking only
# ==============================================================================
make_bp_tie_summary_ppv_fpr <- function(threshold_candidates) {
  threshold_candidates %>%
    dplyr::filter(abs(BP - best_BP) <= numeric_tolerance) %>%
    dplyr::group_by(Estimand, model, N, profile) %>%
    dplyr::summarise(
      selected_FPR = dplyr::first(FPR_exact),
      selected_PPV = dplyr::first(PPV_calibration),
      best_BP = dplyr::first(best_BP),
      n_BP_ties = dplyr::n(),
      .groups = "drop"
    )
}
# 18. Estimation summaries
# ==============================================================================
make_estimation_one <- function(dat, estimand) {
  estimand_columns <- get_estimand_columns(estimand)
  dat %>% dplyr::mutate( true_value = .data[[estimand_columns$true]], interim_estimate =
                           .data[[estimand_columns$mean1]], final_estimate = .data[[estimand_columns$mean2]], final_lower =
                           .data[[estimand_columns$q025_2]], final_upper = .data[[estimand_columns$q975_2]], interim_error =
                           interim_estimate - true_value, final_error = final_estimate - true_value, final_cover = true_value >=
                           final_lower & true_value <= final_upper ) %>%
    dplyr::group_by(model, N, prior_id, design_prior_type, prior_family, profile) %>% dplyr::summarise(
      n_rep_total = dplyr::n(), n_rep_interim_estimation = sum(!is.na(true_value) & !is.na(interim_estimate)),
      n_rep_final_estimation = sum(!is.na(true_value) & !is.na(final_estimate)),
      n_rep_coverage = sum(!is.na(final_cover)), `Interim bias` = safe_mean(interim_error),
      `Final bias` = safe_mean(final_error), `Final RMSE` = safe_rmse(final_error),
      `Final coverage` = pct(safe_mean(final_cover)), .groups = "drop" ) %>% dplyr::mutate(Estimand = estimand) %>%
    dplyr::select( Estimand, model, N, prior_id, design_prior_type, prior_family, profile, n_rep_total,
                   n_rep_interim_estimation, n_rep_final_estimation, n_rep_coverage, `Interim bias`, `Final bias`, `Final RMSE`,
                   `Final coverage` )
}
make_all_estimation <- function(dat) {
  dplyr::bind_rows( make_estimation_one(dat = dat, estimand = "RMST"),
                    make_estimation_one(dat = dat, estimand = "Median") )
}
# ==============================================================================
# ==============================================================================
# 19. Selected efficacy table with exact-null FPR and estimation summaries
# ==============================================================================
make_selected_benefit_table <- function(selected_component_oc, estimation_summary) {
  benefit <- selected_component_oc %>%
    dplyr::filter(prior_id == benefit_prior_id) %>%
    dplyr::select(Estimand, model, N, profile, calibration_profile, nu1_E, nu2_E, ESP, FSP, BP)
  exact <- selected_component_oc %>%
    dplyr::filter(prior_id == exact_null_prior_id) %>%
    dplyr::select(Estimand, model, N, calibration_profile, nu1_E, nu2_E, FPR_exact = FPR)
  benefit %>%
    dplyr::left_join(exact, by = c("Estimand", "model", "N", "calibration_profile", "nu1_E", "nu2_E")) %>%
    dplyr::left_join(
      estimation_summary %>% dplyr::filter(prior_id == benefit_prior_id) %>%
        dplyr::select(Estimand, model, N, profile, `Interim bias`, `Final bias`, `Final RMSE`, `Final coverage`),
      by = c("Estimand", "model", "N", "profile")
    ) %>%
    dplyr::mutate(Profile = calibration_profile) %>%
    dplyr::arrange(factor(Estimand, levels = c("RMST", "Median")), factor(model, levels = names(model_prefixes)), N,
                   factor(Profile, levels = benefit_profile_levels))
}

# ==============================================================================
# 20. Predictive reliability at each profile-specific selected rule
# ============================================================================== 
# PPV and NPV are calculated from BP, FPR, and omega_D. No separate predictive
# probability-of-success quantity is defined or reported.
make_predictive_reliability_table <- function(selected_component_oc) {
  benefit <- selected_component_oc %>%
    dplyr::filter(prior_id == benefit_prior_id) %>%
    dplyr::select(
      Estimand, model, N, calibration_profile, nu1_E, nu2_E,
      BP
    )
  
  exact_null <- selected_component_oc %>%
    dplyr::filter(prior_id == exact_null_prior_id) %>%
    dplyr::select(
      Estimand, model, N, calibration_profile, nu1_E, nu2_E, FPR
    )
  
  benefit %>%
    dplyr::left_join(
      exact_null,
      by = c("Estimand", "model", "N", "calibration_profile", "nu1_E", "nu2_E")
    ) %>%
    tidyr::crossing(mixture_grid) %>%
    dplyr::mutate(
      BP_fraction = BP / 100,
      FPR_fraction = FPR / 100,
      positive_conclusion_prob =
        omega_D * BP_fraction +
        (1 - omega_D) * FPR_fraction,
      PPV_fraction = safe_ratio(
        omega_D * BP_fraction,
        positive_conclusion_prob
      ),
      NPV_fraction = safe_ratio(
        (1 - omega_D) * (1 - FPR_fraction),
        (1 - omega_D) * (1 - FPR_fraction) +
          omega_D * (1 - BP_fraction)
      ),
      PPV = pct(PPV_fraction),
      NPV = pct(NPV_fraction)
    ) %>%
    dplyr::arrange(
      factor(Estimand, levels = c("RMST", "Median")),
      factor(model, levels = names(model_prefixes)),
      N,
      factor(calibration_profile, levels = benefit_profile_levels),
      dplyr::desc(omega_D)
    ) %>%
    dplyr::select(
      Estimand,
      Model = model,
      N,
      Profile = calibration_profile,
      nu1_E,
      nu2_E,
      Ratio,
      omega_D,
      BP,
      FPR,
      PPV,
      NPV
    )
}

# ==============================================================================
# 21. Sample-size BP summaries after profile-specific calibration
# ==============================================================================
make_sample_size_summary <- function(selected_component_oc) {
  selected_component_oc %>%
    dplyr::filter(prior_id == benefit_prior_id) %>%
    dplyr::group_by(Estimand, model, N) %>%
    dplyr::summarise(n_profiles = dplyr::n_distinct(calibration_profile), BP_min = safe_min(BP), BP_mean = safe_mean(BP),
                     BP_max = safe_max(BP), FNR_max = 100 - BP_min, ESP_min = safe_min(ESP), ESP_mean = safe_mean(ESP),
                     .groups = "drop") %>%
    dplyr::mutate(BP_target = bp_target_pct, Margin_to_target = BP_min - BP_target, Eligible = BP_min >= BP_target) %>%
    dplyr::arrange(factor(Estimand, levels = c("RMST", "Median")), factor(model, levels = names(model_prefixes)), N)
}

# ==============================================================================
# 22. Exact-null versus continuous no-efficacy-region FPR at selected rules
# ==============================================================================
make_null_comparison_table <- function(selected_component_oc) {
  exact <- selected_component_oc %>%
    dplyr::filter(prior_id == exact_null_prior_id) %>%
    dplyr::select(Estimand, model, N, calibration_profile, nu1_E, nu2_E, Exact_ESP = ESP, Exact_FSP = FSP,
                  Exact_FPR = FPR)
  region <- selected_component_oc %>%
    dplyr::filter(prior_id == no_benefit_prior_id) %>%
    dplyr::select(Estimand, model, N, calibration_profile, nu1_E, nu2_E, Region_ESP = ESP, Region_FSP = FSP,
                  Region_FPR = FPR)
  exact %>%
    dplyr::left_join(region, by = c("Estimand", "model", "N", "calibration_profile", "nu1_E", "nu2_E")) %>%
    dplyr::arrange(factor(Estimand, levels = c("RMST", "Median")), factor(model, levels = names(model_prefixes)), N,
                   factor(calibration_profile, levels = benefit_profile_levels))
}
# ==============================================================================
# 24. Main execution
# ==============================================================================
dat_raw <- load_all_results()
check_required_columns(dat_raw)
duplicate_check <- make_duplicate_check(dat_raw)
dat <- deduplicate_results(dat_raw)
replicate_id_check <- make_replicate_id_check(dat)
metadata_detail <- make_metadata_check(dat)
integrity_detail <- make_integrity_detail(dat)
integrity_summary <- make_integrity_summary(
  metadata_detail = metadata_detail,
  integrity_detail = integrity_detail
)
cross_model_detail <- make_cross_model_check(dat)
cross_model_summary <- make_cross_model_summary(
  cross_model_detail = cross_model_detail,
  dat = dat
)
mcmc_diagnostics <- make_mcmc_diagnostics(dat)
median_availability <- make_median_availability(dat)

cat(
  "Evaluating ",
  nrow(threshold_grid_full),
  " threshold pairs for RMST...\n",
  sep = ""
)
oc_rmst_full <- make_component_oc_grid(
  dat = dat,
  estimand = "RMST",
  threshold_grid = threshold_grid_full
)

cat(
  "Evaluating ",
  nrow(threshold_grid_full),
  " threshold pairs for mPFS...\n",
  sep = ""
)
oc_median_full <- make_component_oc_grid(
  dat = dat,
  estimand = "Median",
  threshold_grid = threshold_grid_full
)

oc_all_full <- dplyr::bind_rows(
  oc_rmst_full,
  oc_median_full
)

if (any(
  oc_all_full$success_identity_error,
  na.rm = TRUE
)) {
  stop(
    "At least one OC row violates Success = ESP + FSP."
  )
}

exact_null_surface <- make_exact_null_surface(oc_all_full)
benefit_power_surface <- make_benefit_power_surface(oc_all_full)
estimation_summary <- make_all_estimation(dat)

# ------------------------------------------------------------------------------
# Strategy 1: FPR-only / FPR-first calibration
# ------------------------------------------------------------------------------
calibration_fpr <- select_profile_specific_thresholds_fpr_only(
  exact_null_surface = exact_null_surface,
  benefit_power_surface = benefit_power_surface
)

fpr_boundary_summary <- calibration_fpr$fpr_boundary
target_fpr_pairs <- calibration_fpr$target_fpr_pairs
threshold_candidates_fpr <- calibration_fpr$candidates
selected_thresholds_fpr <- calibration_fpr$selected
calibration_adequacy_fpr <- calibration_fpr$adequacy
bp_tie_summary_fpr <- make_bp_tie_summary_fpr_only(
  threshold_candidates_fpr
)

selected_component_oc_fpr <- make_selected_component_oc(
  dat = dat,
  selected_thresholds = selected_thresholds_fpr
)
selected_benefit_table_fpr <- make_selected_benefit_table(
  selected_component_oc_fpr,
  estimation_summary
)
null_comparison_table_fpr <- make_null_comparison_table(
  selected_component_oc_fpr
)
predictive_reliability_table_fpr <- make_predictive_reliability_table(
  selected_component_oc_fpr
)
sample_size_summary_fpr <- make_sample_size_summary(
  selected_component_oc_fpr
)
rmst_sample_size_summary_fpr <- sample_size_summary_fpr %>%
  dplyr::filter(Estimand == "RMST") %>%
  dplyr::arrange(
    factor(model, levels = names(model_prefixes)),
    N
  )

# ------------------------------------------------------------------------------
# Strategy 2: PPV + FPR calibration at omega = 0.5
# ------------------------------------------------------------------------------
ppv_calibration_surface <- make_ppv_calibration_surface(
  exact_null_surface = exact_null_surface,
  benefit_power_surface = benefit_power_surface
)

calibration_ppv_fpr <- select_profile_specific_thresholds_ppv_fpr(
  calibration_surface = ppv_calibration_surface
)

threshold_candidates_ppv_fpr <- calibration_ppv_fpr$candidates
selected_thresholds_ppv_fpr <- calibration_ppv_fpr$selected
calibration_feasibility_ppv_fpr <- calibration_ppv_fpr$feasibility
bp_tie_summary_ppv_fpr <- make_bp_tie_summary_ppv_fpr(
  threshold_candidates_ppv_fpr
)

selected_component_oc_ppv_fpr <- make_selected_component_oc(
  dat = dat,
  selected_thresholds = selected_thresholds_ppv_fpr
)
selected_benefit_table_ppv_fpr <- make_selected_benefit_table(
  selected_component_oc_ppv_fpr,
  estimation_summary
)
predictive_reliability_table_ppv_fpr <- make_predictive_reliability_table(
  selected_component_oc_ppv_fpr
)
sample_size_summary_ppv_fpr <- make_sample_size_summary(
  selected_component_oc_ppv_fpr
)
rmst_sample_size_summary_ppv_fpr <- sample_size_summary_ppv_fpr %>%
  dplyr::filter(Estimand == "RMST") %>%
  dplyr::arrange(
    factor(model, levels = names(model_prefixes)),
    N
  )

# ------------------------------------------------------------------------------
# Required row-count checks for BOTH strategies
# ------------------------------------------------------------------------------
for (strategy_name in c("FPR-only", "PPV+FPR")) {
  selected_thresholds_i <- if (strategy_name == "FPR-only") {
    selected_thresholds_fpr
  } else {
    selected_thresholds_ppv_fpr
  }
  
  selected_benefit_i <- if (strategy_name == "FPR-only") {
    selected_benefit_table_fpr
  } else {
    selected_benefit_table_ppv_fpr
  }
  
  if (nrow(selected_thresholds_i) != 24L) {
    stop(
      strategy_name,
      " calibration must return 24 selected threshold rows."
    )
  }
  if (nrow(selected_benefit_i) != 24L) {
    stop(
      strategy_name,
      " calibration must return 24 selected efficacy rows."
    )
  }
}

# The exact-null versus continuous no-efficacy-region comparison is retained
# only for the FPR-only calibration (Supplementary Table S1).
if (nrow(null_comparison_table_fpr) != 24L) {
  stop(
    "FPR-only calibration must return 24 null-comparison rows for Supplementary Table S1."
  )
}

# Backward-compatible aliases point to the FPR-only manuscript-primary result.
selected_thresholds <- selected_thresholds_fpr
selected_component_oc <- selected_component_oc_fpr
selected_benefit_table <- selected_benefit_table_fpr
null_comparison_table <- null_comparison_table_fpr
predictive_reliability_table <- predictive_reliability_table_fpr
sample_size_summary <- sample_size_summary_fpr
rmst_sample_size_summary <- rmst_sample_size_summary_fpr

# ==============================================================================
# 25. Concise validation summary
# ==============================================================================
missing_files <- file_inventory %>%
  dplyr::filter(!file_found)
if (nrow(missing_files) > 0L) {
  warning(
    "Missing expected result files: ",
    paste(missing_files$source_file, collapse = "; ")
  )
}

if (nrow(duplicate_check) > 0L) {
  warning(
    "Duplicated model x N x prior x profile x replicate rows were found. ",
    "Inspect `duplicate_check`."
  )
}

incomplete_settings <- replicate_id_check %>%
  dplyr::filter(!Complete)
if (nrow(incomplete_settings) > 0L) {
  warning(
    "Incomplete replicate-ID settings were found. ",
    "Inspect `incomplete_settings`."
  )
}

if (nrow(metadata_detail) > 0L) {
  warning(
    "Metadata mismatches were found. Inspect `metadata_detail`."
  )
}

if (nrow(integrity_detail) > 0L) {
  warning(
    "Row-level integrity issues were found. Inspect `integrity_detail`."
  )
}

if (nrow(cross_model_detail) > 0L) {
  warning(
    "Cross-model shared-data mismatches were found. ",
    "Inspect `cross_model_detail`."
  )
}

# ==============================================================================
# 26. Manuscript-table display helpers
# ==============================================================================
estimand_display <- function(x) {
  dplyr::recode(as.character(x), RMST = "RMST", Median = "mPFS")
}
profile_display <- function(x) {
  dplyr::recode( as.character(x), constant = "Constant", increasing = "Increasing", waning = "Waning", none = "None" )
}
model_display <- function(x) {
  dplyr::recode( as.character(x), `Original GP-PH` = "GP-PH", `Independent GP` = "I-GP" )
}
suppress_repeated_labels <- function(dat, hierarchy) {
  if (nrow(dat) <= 1L || length(hierarchy) == 0L) return(dat)
  missing_columns <- setdiff(hierarchy, names(dat))
  if (length(missing_columns) > 0L) {
    stop( "Cannot suppress repeated labels. Missing columns: ", paste(missing_columns, collapse = ", ") )
  }
  reference_dat <- dat
  output_dat <- dat
  for (j in seq_along(hierarchy)) {
    current_col <- hierarchy[j]
    key_cols <- hierarchy[seq_len(j)]
    suppress_row <- rep(FALSE, nrow(reference_dat))
    for (i in 2:nrow(reference_dat)) {
      suppress_row[i] <- all( vapply( key_cols, function(k) {
        identical( as.character(reference_dat[[k]][i]), as.character(reference_dat[[k]][i - 1L]) )
      }, logical(1) ) )
    }
    output_dat[[current_col]] <- as.character(output_dat[[current_col]])
    output_dat[[current_col]][suppress_row] <- ""
  }
  output_dat
}
format_K2_N <- function(N) {
  N_numeric <- suppressWarnings(as.numeric(as.character(N)))
  K2_numeric <- as.integer(round(N_numeric * (1 - target_censor)))
  ifelse(
    is.na(N_numeric),
    "--",
    paste0(K2_numeric, " (", as.integer(N_numeric), ")")
  )
}
format_threshold <- function(x) {
  ifelse(is.na(x), "--", sprintf("%.3f", x))
}
format_percent <- function(x) {
  ifelse(is.na(x), "--", sprintf("%.1f", x))
}
format_estimation <- function(x) {
  ifelse(is.na(x), "--", sprintf("%.2f", x))
}
format_coefficient <- function(x) {
  ifelse(is.na(x), "--", sprintf("%.2f", x))
}
format_rhat <- function(x) {
  ifelse(is.na(x), "--", sprintf("%.3f", x))
}
format_integer <- function(x) {
  ifelse(is.na(x), "--", sprintf("%.0f", x))
}
format_divergence_fit_pct <- function(x) {
  dplyr::case_when( is.na(x) ~ "--", x > 0 & x < 0.1 ~ "<0.1", TRUE ~ sprintf("%.1f", x) )
}
# ==============================================================================
# ==============================================================================
# 27. Reference-mixture PPV/NPV for both calibration strategies
# ==============================================================================
make_reference_predictive_table <- function(
    predictive_reliability_table,
    selected_thresholds,
    strategy_name
) {
  out <- predictive_reliability_table %>%
    dplyr::filter(
      abs(omega_D - predictive_reference_omega) <= numeric_tolerance
    ) %>%
    dplyr::transmute(
      Estimand,
      model = Model,
      N,
      Profile,
      nu1_E,
      nu2_E,
      PPV_reference = PPV,
      NPV_reference = NPV
    ) %>%
    dplyr::arrange(
      factor(Estimand, levels = c("RMST", "Median")),
      factor(model, levels = names(model_prefixes)),
      N,
      factor(Profile, levels = benefit_profile_levels)
    )
  
  if (nrow(out) != 24L) {
    stop(
      strategy_name,
      " reference PPV/NPV table must contain exactly 24 rows."
    )
  }
  
  selected_key <- selected_thresholds %>%
    dplyr::transmute(
      Estimand,
      model,
      N,
      Profile = profile,
      nu1_E,
      nu2_E
    )
  
  missing_reference <- selected_key %>%
    dplyr::anti_join(
      out,
      by = c(
        "Estimand", "model", "N", "Profile",
        "nu1_E", "nu2_E"
      )
    )
  
  if (nrow(missing_reference) > 0L) {
    stop(
      strategy_name,
      " reference PPV/NPV could not be matched to every selected rule."
    )
  }
  
  out
}

reference_predictive_table_fpr <- make_reference_predictive_table(
  predictive_reliability_table = predictive_reliability_table_fpr,
  selected_thresholds = selected_thresholds_fpr,
  strategy_name = "FPR-only"
)

reference_predictive_table_ppv_fpr <- make_reference_predictive_table(
  predictive_reliability_table = predictive_reliability_table_ppv_fpr,
  selected_thresholds = selected_thresholds_ppv_fpr,
  strategy_name = "PPV+FPR"
)

# Candidate-level PPV used during PPV+FPR calibration must equal the reported
# PPV at the same omega = 0.5.
selected_ppv_validation <- selected_thresholds_ppv_fpr %>%
  dplyr::select(
    Estimand, model, N,
    Profile = profile,
    nu1_E, nu2_E,
    PPV_calibration
  ) %>%
  dplyr::left_join(
    reference_predictive_table_ppv_fpr %>%
      dplyr::select(
        Estimand, model, N, Profile,
        nu1_E, nu2_E,
        PPV_reference
      ),
    by = c(
      "Estimand", "model", "N", "Profile",
      "nu1_E", "nu2_E"
    )
  ) %>%
  dplyr::mutate(
    PPV_difference =
      PPV_calibration - PPV_reference
  )

if (any(
  abs(selected_ppv_validation$PPV_difference) > 1e-8,
  na.rm = TRUE
)) {
  stop(
    "Candidate-level PPV used in PPV+FPR calibration does not match ",
    "the reported reference PPV at omega = 0.5."
  )
}

# ==============================================================================
# 28. Main Table objects for BOTH strategies
# ==============================================================================

make_main_table_numeric <- function(
    selected_benefit_table,
    reference_predictive_table,
    strategy_name
) {
  out <- selected_benefit_table %>%
    dplyr::left_join(
      reference_predictive_table,
      by = c(
        "Estimand", "model", "N", "Profile",
        "nu1_E", "nu2_E"
      )
    ) %>%
    dplyr::mutate(
      Estimand = factor(
        Estimand,
        levels = c("RMST", "Median")
      ),
      model = factor(
        model,
        levels = names(model_prefixes)
      ),
      Profile = factor(
        Profile,
        levels = benefit_profile_levels
      )
    ) %>%
    dplyr::arrange(
      Estimand,
      model,
      N,
      Profile
    )
  
  if (nrow(out) != 24L) {
    stop(
      strategy_name,
      " Main Table must contain exactly 24 rows."
    )
  }
  
  if (anyNA(out$PPV_reference) ||
      anyNA(out$NPV_reference)) {
    stop(
      strategy_name,
      " Main Table is missing PPV/NPV values."
    )
  }
  
  out
}

Table1_FPR_all_numeric <- make_main_table_numeric(
  selected_benefit_table = selected_benefit_table_fpr,
  reference_predictive_table = reference_predictive_table_fpr,
  strategy_name = "FPR-only"
)

Table1_PPV_FPR_all_numeric <- make_main_table_numeric(
  selected_benefit_table = selected_benefit_table_ppv_fpr,
  reference_predictive_table = reference_predictive_table_ppv_fpr,
  strategy_name = "PPV+FPR"
)

Table1_FPR_numeric <- Table1_FPR_all_numeric %>%
  dplyr::filter(as.character(Estimand) == "RMST")
Table1_PPV_FPR_numeric <- Table1_PPV_FPR_all_numeric %>%
  dplyr::filter(as.character(Estimand) == "RMST")
if (nrow(Table1_FPR_numeric) != 12L || nrow(Table1_PPV_FPR_numeric) != 12L) {
  stop("Each RMST Main Table must contain exactly 12 rows.")
}

# FPR-only table: requested Conclusion OC order = FPR, BP, PPV, NPV.
Table1_FPR <- Table1_FPR_numeric %>%
  dplyr::transmute(
    Estimand = estimand_display(Estimand),
    Model = model_display(model),
    N = format_K2_N(N),
    Profile = profile_display(Profile),
    nu1 = format_threshold(nu1_E),
    nu2 = format_threshold(nu2_E),
    FPR = format_percent(FPR_exact),
    BP = format_percent(BP),
    PPV = format_percent(PPV_reference),
    NPV = format_percent(NPV_reference),
    SP1 = format_percent(ESP),
    SP2 = format_percent(FSP),
    `Int. bias` = format_estimation(`Interim bias`),
    `Fin. bias` = format_estimation(`Final bias`),
    RMSE = format_estimation(`Final RMSE`),
    CP = format_percent(`Final coverage`)
  ) %>%
  suppress_repeated_labels(
    hierarchy = c("Estimand", "Model", "N")
  )

# PPV+FPR table: requested Conclusion OC order = PPV, FPR, BP, NPV.
Table1_PPV_FPR <- Table1_PPV_FPR_numeric %>%
  dplyr::transmute(
    Estimand = estimand_display(Estimand),
    Model = model_display(model),
    N = format_K2_N(N),
    Profile = profile_display(Profile),
    nu1 = format_threshold(nu1_E),
    nu2 = format_threshold(nu2_E),
    PPV = format_percent(PPV_reference),
    FPR = format_percent(FPR_exact),
    BP = format_percent(BP),
    NPV = format_percent(NPV_reference),
    SP1 = format_percent(ESP),
    SP2 = format_percent(FSP),
    `Int. bias` = format_estimation(`Interim bias`),
    `Fin. bias` = format_estimation(`Final bias`),
    RMSE = format_estimation(`Final RMSE`),
    CP = format_percent(`Final coverage`)
  ) %>%
  suppress_repeated_labels(
    hierarchy = c("Estimand", "Model", "N")
  )

# Backward-compatible manuscript-primary aliases use the FPR-only strategy.
Table1_Main_numeric <- Table1_FPR_numeric
Table1_Main <- Table1_FPR

# Main Tables 2 and 3 are intentionally not generated.

# ==============================================================================
# 30. SUPPLEMENTARY TABLE S1: exact-null versus regional FPR
#     FPR-only calibration only
# ==============================================================================
make_null_comparison_display <- function(
    null_comparison_table,
    strategy_name
) {
  numeric_table <- null_comparison_table %>%
    dplyr::filter(Estimand == "RMST") %>%
    dplyr::mutate(
      Estimand = factor(
        Estimand,
        levels = "RMST"
      ),
      model = factor(
        model,
        levels = names(model_prefixes)
      ),
      calibration_profile = factor(
        calibration_profile,
        levels = benefit_profile_levels
      )
    ) %>%
    dplyr::arrange(
      Estimand,
      model,
      N,
      calibration_profile
    )
  
  if (nrow(numeric_table) != 12L) {
    stop(
      strategy_name,
      " Supplementary Table S1 must contain 12 RMST rows."
    )
  }
  
  display_table <- numeric_table %>%
    dplyr::transmute(
      Estimand = estimand_display(Estimand),
      Model = model_display(model),
      N = format_K2_N(N),
      Profile = profile_display(calibration_profile),
      nu1 = format_threshold(nu1_E),
      nu2 = format_threshold(nu2_E),
      FPR = format_percent(Exact_FPR),
      `Regional FPR` = format_percent(Region_FPR)
    ) %>%
    suppress_repeated_labels(
      hierarchy = c("Estimand", "Model", "N")
    )
  
  list(
    numeric = numeric_table,
    display = display_table
  )
}

TableS1_FPR_object <- make_null_comparison_display(
  null_comparison_table = null_comparison_table_fpr,
  strategy_name = "FPR-only"
)
TableS1_FPR_numeric <- TableS1_FPR_object$numeric
TableS1_FPR <- TableS1_FPR_object$display

# Backward-compatible manuscript-primary aliases.
TableS1_Supp_numeric <- TableS1_FPR_numeric
TableS1_Supp <- TableS1_FPR

# ==============================================================================
# ==============================================================================
# 31. SUPPLEMENTARY TABLE S2: mPFS, FPR-only calibration
# ==============================================================================
TableS2_mPFS_numeric <- Table1_FPR_all_numeric %>%
  dplyr::filter(as.character(Estimand) == "Median")
if (nrow(TableS2_mPFS_numeric) != 12L) stop("Supplementary Table S2 must contain exactly 12 mPFS rows.")
TableS2_mPFS <- TableS2_mPFS_numeric %>%
  dplyr::transmute(
    Estimand = estimand_display(Estimand), Model = model_display(model), N = format_K2_N(N),
    Profile = profile_display(Profile), nu1 = format_threshold(nu1_E), nu2 = format_threshold(nu2_E),
    FPR = format_percent(FPR_exact), BP = format_percent(BP), PPV = format_percent(PPV_reference), NPV = format_percent(NPV_reference),
    SP1 = format_percent(ESP), SP2 = format_percent(FSP), `Int. bias` = format_estimation(`Interim bias`),
    `Fin. bias` = format_estimation(`Final bias`), RMSE = format_estimation(`Final RMSE`), CP = format_percent(`Final coverage`)
  ) %>% suppress_repeated_labels(hierarchy = c("Estimand", "Model", "N"))

# ==============================================================================
# 32. SUPPLEMENTARY TABLE S3: mPFS, PPV+FPR calibration
# ==============================================================================
TableS3_mPFS_numeric <- Table1_PPV_FPR_all_numeric %>%
  dplyr::filter(as.character(Estimand) == "Median")
if (nrow(TableS3_mPFS_numeric) != 12L) stop("Supplementary Table S3 must contain exactly 12 mPFS rows.")
TableS3_mPFS <- TableS3_mPFS_numeric %>%
  dplyr::transmute(
    Estimand = estimand_display(Estimand), Model = model_display(model), N = format_K2_N(N),
    Profile = profile_display(Profile), nu1 = format_threshold(nu1_E), nu2 = format_threshold(nu2_E),
    PPV = format_percent(PPV_reference), FPR = format_percent(FPR_exact), BP = format_percent(BP), NPV = format_percent(NPV_reference),
    SP1 = format_percent(ESP), SP2 = format_percent(FSP), `Int. bias` = format_estimation(`Interim bias`),
    `Fin. bias` = format_estimation(`Final bias`), RMSE = format_estimation(`Final RMSE`), CP = format_percent(`Final coverage`)
  ) %>% suppress_repeated_labels(hierarchy = c("Estimand", "Model", "N"))

# ==============================================================================
# 33. ADDITIONAL SUPPLEMENTARY TABLE: MCMC diagnostics
# ==============================================================================
mcmc_rhat_columns <- intersect(
  c("max_Rhat_1", "max_Rhat_2"),
  names(dat)
)
mcmc_ess_columns <- intersect(
  c("min_ESS_1", "min_ESS_2"),
  names(dat)
)
mcmc_divergent_columns <- intersect(
  c("n_divergent_1", "n_divergent_2"),
  names(dat)
)

Table_MCMC_numeric <- dat %>%
  dplyr::mutate(
    max_Rhat_any = row_max_na(., mcmc_rhat_columns),
    min_ESS_any = row_min_na(., mcmc_ess_columns),
    n_divergent_any = row_sum_na(., mcmc_divergent_columns),
    any_divergence = n_divergent_any > 0
  ) %>%
  dplyr::group_by(model, N) %>%
  dplyr::summarise(
    Max_Rhat = safe_max(max_Rhat_any),
    Min_bulk_ESS = safe_min(min_ESS_any),
    Divergences = safe_sum(n_divergent_any),
    Replicates_any_divergence = pct(
      safe_mean(any_divergence)
    ),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    model = factor(
      model,
      levels = names(model_prefixes)
    )
  ) %>%
  dplyr::arrange(model, N)

if (nrow(Table_MCMC_numeric) != 4L) {
  stop("MCMC diagnostics table must contain exactly 4 model x N rows.")
}

Table_MCMC <- Table_MCMC_numeric %>%
  dplyr::transmute(
    Model = model_display(model),
    N = format_K2_N(N),
    `Max Rhat` = ifelse(
      is.na(Max_Rhat),
      "--",
      sprintf("%.3f", Max_Rhat)
    ),
    `Min ESS` = ifelse(
      is.na(Min_bulk_ESS),
      "--",
      sprintf("%.0f", Min_bulk_ESS)
    ),
    Total = ifelse(
      is.na(Divergences),
      "--",
      sprintf("%.0f", Divergences)
    ),
    `Fits (%)` = dplyr::case_when(
      is.na(Replicates_any_divergence) ~ "--",
      Replicates_any_divergence > 0 &
        Replicates_any_divergence < 0.1 ~ "<0.1",
      TRUE ~ sprintf("%.1f", Replicates_any_divergence)
    )
  ) %>%
  suppress_repeated_labels(
    hierarchy = "Model"
  )

# 41. MACHINE-READABLE SUPPLEMENTARY DATA
# ==============================================================================
#
# No machine-readable supplementary CSV files are produced.
# ==============================================================================
# 42. LaTeX output directories and generic writer
# ==============================================================================
latex_dir <- file.path(results_dir, "LaTeX")
if (!dir.exists(latex_dir)) {
  dir.create(latex_dir, recursive = TRUE)
}
main_tex_name <- function(name) {
  paste0("tableBDCT_Main_", name, ".tex")
}
supp_tex_name <- function(name) {
  paste0("tableBDCT_Supp_", name, ".tex")
}

render_latex_table <- function( dat, caption, label = NULL, column_names = NULL, header_above = NULL,
                                scale_down = FALSE, table_position = "!htbp", bold_rows = NULL, font_size = NULL, escape = TRUE,
                                suppress_list_entry = FALSE, latex_after_rows = NULL, table_note = NULL ) {
  if (is.null(dat) || nrow(dat) == 0L) {
    stop("Cannot render an empty LaTeX table.")
  }
  caption_placeholder <- "BDSCAPTIONPLACEHOLDER"
  label_placeholder <- "BDSLABELPLACEHOLDER"
  display_column_names <- if (is.null(column_names)) {
    names(dat)
  } else {
    column_names
  }
  if (length(display_column_names) != ncol(dat)) {
    stop("`column_names` must have one entry per table column.")
  }
  if (!is.null(header_above)) {
    if (sum(as.integer(header_above)) != ncol(dat)) {
      stop( "`header_above` spans must sum to the number of table columns." )
    }
  }
  latex_options <- character(0)
  if (isTRUE(scale_down)) {
    latex_options <- c(latex_options, "scale_down")
  }
  tab <- kableExtra::kbl( dat, format = "latex", booktabs = TRUE, caption = caption_placeholder,
                          label = if (is.null(label)) NULL else label_placeholder, col.names = display_column_names, align = "c",
                          escape = escape, linesep = "" )
  if (!is.null(header_above)) {
    tab <- tab %>% kableExtra::add_header_above( header = header_above, escape = FALSE )
  }
  # kable_styling() cannot take latex_options = character(0).
  if (length(latex_options) > 0L) {
    if (is.null(font_size)) {
      tab <- tab %>% kableExtra::kable_styling( latex_options = latex_options, position = "center" )
    } else {
      tab <- tab %>% kableExtra::kable_styling( latex_options = latex_options, position = "center",
                                                font_size = font_size )
    }
  } else {
    if (is.null(font_size)) {
      tab <- tab %>% kableExtra::kable_styling( position = "center" )
    } else {
      tab <- tab %>% kableExtra::kable_styling( position = "center", font_size = font_size )
    }
  }
  if (!is.null(bold_rows) && length(bold_rows) > 0L) {
    bold_rows <- as.integer(bold_rows)
    bold_rows <- bold_rows[
      is.finite(bold_rows) & bold_rows >= 1L & bold_rows <= nrow(dat) ]
    if (length(bold_rows) > 0L) {
      tab <- kableExtra::row_spec( tab, row = bold_rows, bold = TRUE )
    }
  }
  
  # Optional LaTeX rules after selected body rows.
  # Insert separators directly into the completed LaTeX code. Do not use
  # row_spec(extra_latex_after=...), because that mechanism adds vertical
  # spacing around the injected rule and creates an apparent blank row.
  latex_code <- as.character(tab)
  
  if (!is.null(latex_after_rows) && length(latex_after_rows) > 0L) {
    if (is.null(names(latex_after_rows)) || any(names(latex_after_rows) == "")) {
      stop("`latex_after_rows` must be a named character vector with row numbers as names.")
    }
    
    requested_rows <- suppressWarnings(as.integer(names(latex_after_rows)))
    if (any(!is.finite(requested_rows)) ||
        any(requested_rows < 1L) ||
        any(requested_rows >= nrow(dat))) {
      stop("Invalid row number in `latex_after_rows`.")
    }
    
    latex_lines <- strsplit(latex_code, "\n", fixed = TRUE)[[1]]
    body_start <- which(trimws(latex_lines) == "\\midrule")
    body_end <- which(trimws(latex_lines) == "\\bottomrule")
    
    if (length(body_start) == 0L || length(body_end) == 0L) {
      stop("Could not locate the LaTeX table body boundaries.")
    }
    
    body_start <- body_start[1L]
    body_end <- body_end[body_end > body_start][1L]
    
    candidate_lines <- seq.int(body_start + 1L, body_end - 1L)
    data_line_index <- candidate_lines[
      vapply(
        latex_lines[candidate_lines],
        function(x) {
          grepl("&", x, fixed = TRUE) &&
            endsWith(trimws(x), "\\\\")
        },
        logical(1)
      )
    ]
    
    if (length(data_line_index) != nrow(dat)) {
      stop(
        "Could not uniquely identify all LaTeX body rows: expected ",
        nrow(dat), ", found ", length(data_line_index), "."
      )
    }
    
    rule_by_line <- stats::setNames(
      as.character(latex_after_rows),
      as.character(data_line_index[requested_rows])
    )
    
    rebuilt <- character(0)
    for (i in seq_along(latex_lines)) {
      rebuilt <- c(rebuilt, latex_lines[i])
      key <- as.character(i)
      if (key %in% names(rule_by_line)) {
        rebuilt <- c(rebuilt, unname(rule_by_line[[key]]))
      }
    }
    
    latex_code <- paste(rebuilt, collapse = "\n")
  }
  
  if (!is.null(table_note) && nzchar(table_note)) {
    # Insert literal LaTeX lines immediately before \\end{table}.
    # Using sub()/gsub() replacement here can strip LaTeX backslashes.
    latex_lines <- strsplit(latex_code, "\n", fixed = TRUE)[[1]]
    end_table_index <- which(trimws(latex_lines) == "\\end{table}")
    if (length(end_table_index) != 1L) {
      stop("Could not uniquely locate \\\\end{table} when inserting the table note.")
    }
    
    note_lines <- c(
      "\\par\\vspace{2pt}",
      "\\begin{minipage}{\\linewidth}",
      "\\footnotesize\\raggedright",
      paste0("\\textit{Note.} ", table_note),
      "\\end{minipage}"
    )
    
    latex_lines <- append(
      latex_lines,
      values = note_lines,
      after = end_table_index - 1L
    )
    latex_code <- paste(latex_lines, collapse = "\n")
  }
  
  if (!is.null(label)) {
    latex_code <- gsub( paste0("tab:", label_placeholder), label, latex_code, fixed = TRUE )
  }
  latex_code <- gsub( caption_placeholder, caption, latex_code, fixed = TRUE )
  latex_code <- sub( "\\\\begin\\{table\\}(\\[[^]]*\\])?", paste0("\\\\begin{table}[", table_position, "]"), latex_code,
                     perl = TRUE )
  if (isTRUE(suppress_list_entry)) {
    latex_code <- sub( "\\\\caption\\{", "\\\\caption[]{", latex_code, perl = TRUE )
  }
  
  if (!is.null(table_note) && nzchar(table_note)) {
    required_note_tokens <- c(
      "\\par\\vspace{2pt}",
      "\\begin{minipage}{\\linewidth}",
      "\\footnotesize\\raggedright",
      "\\textit{Note.}",
      "\\end{minipage}"
    )
    missing_note_tokens <- required_note_tokens[
      !vapply(
        required_note_tokens,
        function(token) grepl(token, latex_code, fixed = TRUE),
        logical(1)
      )
    ]
    if (length(missing_note_tokens) > 0L) {
      stop(
        "Table note LaTeX was malformed. Missing: ",
        paste(missing_note_tokens, collapse = ", ")
      )
    }
  }
  
  latex_code
}
write_latex_table <- function( dat, filename, caption, label, column_names = NULL, header_above = NULL,
                               scale_down = FALSE, table_position = "!htbp", bold_rows = NULL, font_size = NULL, escape = TRUE,
                               latex_after_rows = NULL, table_note = NULL ) {
  latex_code <- render_latex_table(
    dat = dat,
    caption = caption,
    label = label,
    column_names = column_names,
    header_above = header_above,
    scale_down = scale_down,
    table_position = table_position,
    bold_rows = bold_rows,
    font_size = font_size,
    escape = escape,
    suppress_list_entry = FALSE,
    latex_after_rows = latex_after_rows,
    table_note = table_note
  )
  output_file <- file.path(latex_dir, filename)
  writeLines( latex_code, con = output_file, useBytes = TRUE )
  invisible(output_file)
}
# ==============================================================================
# ==============================================================================
# 43. Write MAIN manuscript tables for BOTH calibration strategies
# ==============================================================================

main_table_cmidrules <- c(
  "3" = "\\cmidrule(l{3pt}r{3pt}){3-3} \\cmidrule(l{3pt}r{3pt}){4-4} \\cmidrule(l{3pt}r{3pt}){5-6} \\cmidrule(l{3pt}r{3pt}){7-10} \\cmidrule(l{3pt}r{3pt}){11-12} \\cmidrule(l{3pt}r{3pt}){13-16}",
  "6" = "\\cmidrule(l{3pt}r{3pt}){2-2} \\cmidrule(l{3pt}r{3pt}){3-3} \\cmidrule(l{3pt}r{3pt}){4-4} \\cmidrule(l{3pt}r{3pt}){5-6} \\cmidrule(l{3pt}r{3pt}){7-10} \\cmidrule(l{3pt}r{3pt}){11-12} \\cmidrule(l{3pt}r{3pt}){13-16}",
  "9" = "\\cmidrule(l{3pt}r{3pt}){3-3} \\cmidrule(l{3pt}r{3pt}){4-4} \\cmidrule(l{3pt}r{3pt}){5-6} \\cmidrule(l{3pt}r{3pt}){7-10} \\cmidrule(l{3pt}r{3pt}){11-12} \\cmidrule(l{3pt}r{3pt}){13-16}"
)

# Strategy 1: FPR-only.
# Requested Conclusion OC order: FPR, BP, PPV, NPV.
file_Table1_FPR <- write_latex_table(
  dat = Table1_FPR,
  filename = main_tex_name("decision_OC_FPR"),
  caption = paste0(
    "FPR-first calibrated posterior efficacy thresholds, conclusion OCs, ",
    "success by look, and estimation performance for the 12 RMST profile-specific ",
    "candidate designs."
  ),
  label = "tableBDCT_Main_decision_OC_FPR",
  header_above = c(
    "Candidate setting" = 4,
    "Calibrated thresholds" = 2,
    "Conclusion OCs" = 4,
    "Success by look" = 2,
    "Estimation performance" = 4
  ),
  column_names = c(
    "Estimand", "Model", "$K_{2}$ ($N$)", "Profile",
    "$\\nu_1$", "$\\nu_2$",
    "FPR", "BP", "PPV", "NPV",
    "$\\mathrm{SP}_1$", "$\\mathrm{SP}_2$",
    "Int. bias", "Fin. bias", "RMSE", "CP"
  ),
  scale_down = TRUE,
  font_size = 7,
  escape = FALSE,
  table_note = paste0(
    "PPV and NPV are evaluated at $\\omega=0.5$. ",
    "Int.\\ and Fin.\\ denote the interim and counterfactual final analyses. ",
    "Bias and RMSE are in months; probability entries are percentages."
  ),
  latex_after_rows = main_table_cmidrules
)

# Strategy 2: PPV + FPR.
# Requested Conclusion OC order: PPV, FPR, BP, NPV.
file_Table1_PPV_FPR <- write_latex_table(
  dat = Table1_PPV_FPR,
  filename = main_tex_name("decision_OC_PPV_FPR"),
  caption = paste0(
    "PPV- and FPR-constrained posterior efficacy thresholds, conclusion OCs, ",
    "success by look, and estimation performance for the 12 RMST profile-specific ",
    "candidate designs."
  ),
  label = "tableBDCT_Main_decision_OC_PPV_FPR",
  header_above = c(
    "Candidate setting" = 4,
    "Calibrated thresholds" = 2,
    "Conclusion OCs" = 4,
    "Success by look" = 2,
    "Estimation performance" = 4
  ),
  column_names = c(
    "Estimand", "Model", "$K_{2}$ ($N$)", "Profile",
    "$\\nu_1$", "$\\nu_2$",
    "PPV", "FPR", "BP", "NPV",
    "$\\mathrm{SP}_1$", "$\\mathrm{SP}_2$",
    "Int. bias", "Fin. bias", "RMSE", "CP"
  ),
  scale_down = TRUE,
  font_size = 7,
  escape = FALSE,
  table_note = paste0(
    "PPV and NPV are evaluated at $\\omega=0.5$. ",
    "Int.\\ and Fin.\\ denote the interim and counterfactual final analyses. ",
    "Bias and RMSE are in months; probability entries are percentages."
  ),
  latex_after_rows = main_table_cmidrules
)

# Backward-compatible alias for the manuscript-primary FPR-only table.
file_Table1 <- file_Table1_FPR

# Main Tables 2 and 3 are intentionally not generated.

# ============================================================================== 
# 44. Write SUPPLEMENTARY manuscript tables
# ==============================================================================
file_TableS1_FPR <- write_latex_table(
  dat = TableS1_FPR, filename = supp_tex_name("null_comparison"),
  caption = paste0(
    "Exact-null and regional FPR for the 12 RMST profile-specific candidate designs selected by FPR-first calibration. ",
    "FPR is the exact-null value used for threshold calibration, whereas Regional FPR averages false-positive behavior ",
    "over the prespecified continuous no-efficacy distribution. Entries are percentages."
  ),
  label = "tableBDCT_Supp_null_comparison",
  header_above = c("Candidate setting" = 4, "Calibrated thresholds" = 2, "FPR evaluation" = 2),
  column_names = c("Estimand", "Model", "$K_{2}$ ($N$)", "Profile", "$\\nu_1$", "$\\nu_2$", "FPR", "Regional FPR"),
  scale_down = TRUE, font_size = 8, escape = FALSE
)
file_TableS1 <- file_TableS1_FPR

file_TableS2 <- write_latex_table(
  dat = TableS2_mPFS, filename = supp_tex_name("mPFS_decision_OC"),
  caption = paste0(
    "FPR-first calibrated posterior efficacy thresholds, conclusion OCs, success by look, and estimation performance ",
    "for the 12 mPFS profile-specific candidate designs."
  ),
  label = "tableBDCT_Supp_mPFS_decision_OC",
  header_above = c("Candidate setting" = 4, "Calibrated thresholds" = 2, "Conclusion OCs" = 4, "Success by look" = 2, "Estimation performance" = 4),
  column_names = c(
    "Estimand", "Model", "$K_{2}$ ($N$)", "Profile", "$\\nu_1$", "$\\nu_2$",
    "FPR", "BP", "PPV", "NPV", "$\\mathrm{SP}_1$", "$\\mathrm{SP}_2$", "Int. bias", "Fin. bias", "RMSE", "CP"
  ),
  scale_down = TRUE, font_size = 7, escape = FALSE,
  table_note = paste0(
    "PPV and NPV are evaluated at $\\omega=0.5$. ",
    "Int.\\ and Fin.\\ denote the interim and counterfactual final analyses. Bias and RMSE are in months; probability entries are percentages."
  ),
  latex_after_rows = main_table_cmidrules
)

file_TableS3 <- write_latex_table(
  dat = TableS3_mPFS, filename = supp_tex_name("mPFS_decision_OC_PPV_FPR"),
  caption = paste0(
    "PPV- and FPR-constrained posterior efficacy thresholds, conclusion OCs, success by look, and estimation performance ",
    "for the 12 mPFS profile-specific candidate designs."
  ),
  label = "tableBDCT_Supp_mPFS_decision_OC_PPV_FPR",
  header_above = c("Candidate setting" = 4, "Calibrated thresholds" = 2, "Conclusion OCs" = 4, "Success by look" = 2, "Estimation performance" = 4),
  column_names = c(
    "Estimand", "Model", "$K_{2}$ ($N$)", "Profile", "$\\nu_1$", "$\\nu_2$",
    "PPV", "FPR", "BP", "NPV", "$\\mathrm{SP}_1$", "$\\mathrm{SP}_2$", "Int. bias", "Fin. bias", "RMSE", "CP"
  ),
  scale_down = TRUE, font_size = 7, escape = FALSE,
  table_note = paste0(
    "PPV and NPV are evaluated at $\\omega=0.5$. ",
    "Int.\\ and Fin.\\ denote the interim and counterfactual final analyses. Bias and RMSE are in months; probability entries are percentages."
  ),
  latex_after_rows = main_table_cmidrules
)

file_Table_MCMC <- write_latex_table(
  dat = Table_MCMC,
  filename = supp_tex_name("MCMC_diagnostics"),
  caption = paste0(
    "Posterior-sampling diagnostics by analysis model and candidate design. ",
    "Max. $\\widehat{R}$ and Min. ESS are worst-case values across posterior fits. ",
    "Total is the total number of divergent transitions, and Fits is the ",
    "percentage of posterior fits containing at least one divergence."
  ),
  label = "tableBDCT_Supp_MCMC_diagnostics",
  header_above = c(
    " " = 4,
    "Divergences" = 2
  ),
  column_names = c(
    "Model",
    "$K_{2}$ ($N$)",
    "Max. $\\widehat{R}$",
    "Min. ESS",
    "Total",
    "Fits (\\%)"
  ),
  scale_down = FALSE,
  escape = FALSE
)

# 45. Copy generated files into an existing manuscript tree when available
# ==============================================================================
manuscript_root_candidates <- unique(normalizePath(c(getwd(), dirname(results_dir), dirname(dirname(results_dir))), winslash = "/", mustWork = FALSE))
copy_to_existing_subdir <- function(files, subdir) {
  target_dirs <- file.path(manuscript_root_candidates, subdir)
  target_dirs <- unique(target_dirs[dir.exists(target_dirs)])
  if (length(target_dirs) == 0L) return(invisible(character(0)))
  copied <- character(0)
  for (target_dir in target_dirs) {
    for (f in files) {
      target_file <- file.path(target_dir, basename(f))
      ok <- file.copy(from = f, to = target_file, overwrite = TRUE)
      if (isTRUE(ok)) copied <- c(copied, target_file)
    }
  }
  invisible(copied)
}
main_files <- c(
  file_Table1_FPR,
  file_Table1_PPV_FPR
)

supp_files <- c(
  file_TableS1_FPR,
  file_TableS2,
  file_TableS3,
  file_Table_MCMC
)

# Remove obsolete prior supplementary tables and Figure 2 from previous runs.
obsolete_prior_table_names <- c(supp_tex_name("GPPH_treatment_coefficient"))
obsolete_prior_table_latex <- file.path(latex_dir, obsolete_prior_table_names)
obsolete_prior_table_latex <- obsolete_prior_table_latex[file.exists(obsolete_prior_table_latex)]
if (length(obsolete_prior_table_latex) > 0L) unlink(obsolete_prior_table_latex, force = TRUE)
obsolete_prior_table_targets <- unique(unlist(lapply(manuscript_root_candidates, function(root) file.path(root, "table", obsolete_prior_table_names)), use.names = FALSE))
obsolete_prior_table_targets <- obsolete_prior_table_targets[file.exists(obsolete_prior_table_targets)]
if (length(obsolete_prior_table_targets) > 0L) unlink(obsolete_prior_table_targets, force = TRUE)

obsolete_figure2_names <- c("figureBDCT_fixed_PPV_FPR_NPV.tex", "BDCT_fixed_PPV_FPR_NPV.pdf")
obsolete_figure2_latex <- file.path(latex_dir, obsolete_figure2_names)
obsolete_figure2_latex <- obsolete_figure2_latex[file.exists(obsolete_figure2_latex)]
if (length(obsolete_figure2_latex) > 0L) unlink(obsolete_figure2_latex, force = TRUE)
obsolete_figure2_targets <- unique(c(
  file.path(manuscript_root_candidates, "figure", "figureBDCT_fixed_PPV_FPR_NPV.tex"),
  file.path(manuscript_root_candidates, "image", "BDCT_fixed_PPV_FPR_NPV.pdf")
))
obsolete_figure2_targets <- obsolete_figure2_targets[file.exists(obsolete_figure2_targets)]
if (length(obsolete_figure2_targets) > 0L) unlink(obsolete_figure2_targets, force = TRUE)

# Remove obsolete standalone calibration-method TeX snippets from prior runs.
obsolete_calibration_method_names <- c(
  "BDCT_calibration_rule_FPR_only.tex",
  "BDCT_calibration_constraints_PPV_FPR.tex"
)

obsolete_calibration_method_latex <- file.path(
  latex_dir,
  obsolete_calibration_method_names
)
obsolete_calibration_method_latex <- obsolete_calibration_method_latex[
  file.exists(obsolete_calibration_method_latex)
]
if (length(obsolete_calibration_method_latex) > 0L) {
  unlink(obsolete_calibration_method_latex, force = TRUE)
}

obsolete_calibration_method_targets <- unique(
  unlist(
    lapply(
      manuscript_root_candidates,
      function(root) {
        file.path(
          root,
          "table",
          obsolete_calibration_method_names
        )
      }
    ),
    use.names = FALSE
  )
)
obsolete_calibration_method_targets <- obsolete_calibration_method_targets[
  file.exists(obsolete_calibration_method_targets)
]
if (length(obsolete_calibration_method_targets) > 0L) {
  unlink(obsolete_calibration_method_targets, force = TRUE)
}

# Remove the excluded PPV+FPR null-comparison supplement so a stale file
# from a previous run cannot be included accidentally.
obsolete_ppv_fpr_supp_name <- supp_tex_name("null_comparison_PPV_FPR")

obsolete_ppv_fpr_supp_latex <- file.path(
  latex_dir,
  obsolete_ppv_fpr_supp_name
)
if (file.exists(obsolete_ppv_fpr_supp_latex)) {
  unlink(obsolete_ppv_fpr_supp_latex, force = TRUE)
}

obsolete_ppv_fpr_supp_targets <- unique(
  file.path(
    manuscript_root_candidates,
    "table",
    obsolete_ppv_fpr_supp_name
  )
)
obsolete_ppv_fpr_supp_targets <- obsolete_ppv_fpr_supp_targets[
  file.exists(obsolete_ppv_fpr_supp_targets)
]
if (length(obsolete_ppv_fpr_supp_targets) > 0L) {
  unlink(obsolete_ppv_fpr_supp_targets, force = TRUE)
}

# Remove the old FPR-only Main Table filename from prior runs.
obsolete_old_fpr_main_name <- main_tex_name("decision_OC")

obsolete_old_fpr_main_latex <- file.path(
  latex_dir,
  obsolete_old_fpr_main_name
)
if (file.exists(obsolete_old_fpr_main_latex)) {
  unlink(obsolete_old_fpr_main_latex, force = TRUE)
}

obsolete_old_fpr_main_targets <- unique(
  file.path(
    manuscript_root_candidates,
    "table",
    obsolete_old_fpr_main_name
  )
)
obsolete_old_fpr_main_targets <- obsolete_old_fpr_main_targets[
  file.exists(obsolete_old_fpr_main_targets)
]
if (length(obsolete_old_fpr_main_targets) > 0L) {
  unlink(obsolete_old_fpr_main_targets, force = TRUE)
}

# Remove obsolete Main Tables 2 and 3 so stale files cannot be included accidentally.
obsolete_main_table_names <- c(
  main_tex_name("estimation_OC"),
  main_tex_name("RMST_predictive_OC")
)

obsolete_main_latex <- file.path(latex_dir, obsolete_main_table_names)
obsolete_main_latex <- obsolete_main_latex[file.exists(obsolete_main_latex)]
if (length(obsolete_main_latex) > 0L) unlink(obsolete_main_latex, force = TRUE)

obsolete_main_targets <- unlist(
  lapply(
    obsolete_main_table_names,
    function(x) file.path(manuscript_root_candidates, "table", x)
  ),
  use.names = FALSE
)
obsolete_main_targets <- unique(
  obsolete_main_targets[file.exists(obsolete_main_targets)]
)
if (length(obsolete_main_targets) > 0L) {
  unlink(obsolete_main_targets, force = TRUE)
}

all_written_table_files <- c(main_files, supp_files)
copied_table_files <- copy_to_existing_subdir(
  all_written_table_files,
  "table"
)

all_figure_tex_files <- c(predictive_figure_tex)
all_figure_pdf_files <- c(predictive_figure_pdf)
copied_figure_tex <- copy_to_existing_subdir(all_figure_tex_files, "figure")
copied_figure_pdf <- copy_to_existing_subdir(all_figure_pdf_files, "image")

# ==============================================================================
# 46. Final validation and console summary
# ==============================================================================

excluded_ppv_fpr_supp_file <- file.path(
  latex_dir,
  supp_tex_name("null_comparison_PPV_FPR")
)
if (file.exists(excluded_ppv_fpr_supp_file)) {
  stop(
    "Excluded PPV+FPR null-comparison supplementary table still exists: ",
    excluded_ppv_fpr_supp_file
  )
}

required_supplement_files <- c(
  file_TableS1_FPR,
  file_Table_MCMC
)
if (any(!file.exists(required_supplement_files))) {
  stop(
    "Required supplementary table files were not written: ",
    paste(
      basename(
        required_supplement_files[
          !file.exists(required_supplement_files)
        ]
      ),
      collapse = ", "
    )
  )
}

missing_written_table_files <- all_written_table_files[
  !file.exists(all_written_table_files)
]
if (length(missing_written_table_files) > 0L) {
  stop(
    "The following expected LaTeX table files were not written: ",
    paste(
      basename(missing_written_table_files),
      collapse = ", "
    )
  )
}

# ------------------------------------------------------------------------------
# Validate Main Table 1 note and OC order for FPR-only calibration.
# ------------------------------------------------------------------------------
fpr_table_tex_check <- readLines(
  file_Table1_FPR,
  warn = FALSE,
  encoding = "UTF-8"
)

required_fpr_note_tokens <- c(
  "\\par\\vspace{2pt}",
  "\\begin{minipage}{\\linewidth}",
  "\\footnotesize\\raggedright",
  "\\textit{Note.} PPV and NPV are evaluated at $\\omega=0.5$.",
  "\\end{minipage}"
)
missing_fpr_note_tokens <- required_fpr_note_tokens[
  !vapply(
    required_fpr_note_tokens,
    function(token) {
      any(grepl(
        token,
        fpr_table_tex_check,
        fixed = TRUE
      ))
    },
    logical(1)
  )
]
if (length(missing_fpr_note_tokens) > 0L) {
  stop(
    "FPR-only Main Table note was not written correctly. Missing: ",
    paste(
      missing_fpr_note_tokens,
      collapse = " | "
    )
  )
}

if (!any(grepl(
  "FPR & BP & PPV & NPV",
  fpr_table_tex_check,
  fixed = TRUE
))) {
  stop(
    "FPR-only Main Table Conclusion OC order must be FPR, BP, PPV, NPV."
  )
}

if (!any(grepl(
  "$K_{2}$ ($N$)",
  fpr_table_tex_check,
  fixed = TRUE
))) {
  stop(
    "FPR-only Main Table candidate-size header must be $K_{2}$ ($N$)."
  )
}

# ------------------------------------------------------------------------------
# Validate Main Table 1 note and OC order for PPV+FPR calibration.
# ------------------------------------------------------------------------------
ppv_fpr_table_tex_check <- readLines(
  file_Table1_PPV_FPR,
  warn = FALSE,
  encoding = "UTF-8"
)

required_ppv_fpr_note_tokens <- c(
  "\\par\\vspace{2pt}",
  "\\begin{minipage}{\\linewidth}",
  "\\footnotesize\\raggedright",
  "\\textit{Note.} PPV and NPV are evaluated at $\\omega=0.5$.",
  "\\end{minipage}"
)
missing_ppv_fpr_note_tokens <- required_ppv_fpr_note_tokens[
  !vapply(
    required_ppv_fpr_note_tokens,
    function(token) {
      any(grepl(
        token,
        ppv_fpr_table_tex_check,
        fixed = TRUE
      ))
    },
    logical(1)
  )
]
if (length(missing_ppv_fpr_note_tokens) > 0L) {
  stop(
    "PPV+FPR Main Table note was not written correctly. Missing: ",
    paste(
      missing_ppv_fpr_note_tokens,
      collapse = " | "
    )
  )
}

if (!any(grepl(
  "PPV & FPR & BP & NPV",
  ppv_fpr_table_tex_check,
  fixed = TRUE
))) {
  stop(
    "PPV+FPR Main Table Conclusion OC order must be PPV, FPR, BP, NPV."
  )
}

if (!any(grepl(
  "$K_{2}$ ($N$)",
  ppv_fpr_table_tex_check,
  fixed = TRUE
))) {
  stop(
    "PPV+FPR Main Table candidate-size header must be $K_{2}$ ($N$)."
  )
}

# ------------------------------------------------------------------------------
# Validate RMST/mPFS manuscript split.
# ------------------------------------------------------------------------------
for (f in c(file_Table1_FPR, file_Table1_PPV_FPR, file_TableS1_FPR)) {
  lines <- readLines(f, warn = FALSE, encoding = "UTF-8")
  if (any(grepl("Median", lines, fixed = TRUE)) || any(grepl("mPFS", lines, fixed = TRUE))) {
    stop("RMST-only table contains an unexpected mPFS/Median label: ", basename(f))
  }
}
for (f in c(file_TableS2, file_TableS3)) {
  lines <- readLines(f, warn = FALSE, encoding = "UTF-8")
  if (!any(grepl("mPFS", lines, fixed = TRUE))) stop("mPFS table lacks the mPFS label: ", basename(f))
  if (any(grepl("Median", lines, fixed = TRUE))) stop("mPFS table still contains Median: ", basename(f))
}

# ------------------------------------------------------------------------------
# Figure-file validation.
# ------------------------------------------------------------------------------
missing_written_figure_files <- c(
  all_figure_tex_files,
  all_figure_pdf_files
)
missing_written_figure_files <- missing_written_figure_files[
  !file.exists(missing_written_figure_files)
]
if (length(missing_written_figure_files) > 0L) {
  stop(
    "The following expected figure files were not written: ",
    paste(
      basename(missing_written_figure_files),
      collapse = ", "
    )
  )
}

cat("\n============================================================\n")
cat("BDCT MANUSCRIPT OUTPUT GENERATED\n")
cat("============================================================\n")
cat("Results directory: ", results_dir, "\n", sep = "")
cat(
  "FPR-only -- RMST: ", sum(calibration_adequacy_fpr$Pass_BP[calibration_adequacy_fpr$Estimand == "RMST"]),
  "/12 meet BP >= 90%; mPFS: ", sum(calibration_adequacy_fpr$Pass_BP[calibration_adequacy_fpr$Estimand == "Median"]),
  "/12 meet BP >= 90%.\n", sep = ""
)
cat(
  "PPV+FPR -- RMST: ", sum(calibration_feasibility_ppv_fpr$Pass_all[calibration_feasibility_ppv_fpr$Estimand == "RMST"]),
  "/12 satisfy all constraints; mPFS: ", sum(calibration_feasibility_ppv_fpr$Pass_all[calibration_feasibility_ppv_fpr$Estimand == "Median"]),
  "/12 satisfy all constraints.\n", sep = ""
)
cat("\nMain tables (RMST only):\n")
for (f in main_files) cat("  ", basename(f), "\n", sep = "")
cat("\nSupplementary tables:\n")
for (f in supp_files) cat("  ", basename(f), "\n", sep = "")
cat("\nFigure 1 only:\n  ", basename(predictive_figure_tex), "\n  ", basename(predictive_figure_pdf), "\n", sep = "")
cat("\nOC order:\n  FPR-only: FPR, BP, PPV, NPV\n  PPV+FPR:  PPV, FPR, BP, NPV\n", sep = "")

# ==============================================================================
# OPTIONAL FIGURE 2
# Set `make_figure2 <- TRUE` only when Figure 2 is needed.
# By default this block does not generate any Figure 2 files.
# ==============================================================================
make_figure2 <- FALSE

if (isTRUE(make_figure2)) {
  figure2_PPV <- ppv_calibration_target
  figure2_BP <- c(0.80, 0.90)
  figure2_omega <- seq(0.001, 0.999, by = 0.001)
  figure2_target_FPR <- 0.025
  figure2_target_NPV <- 0.90
  
  figure2_FPR_from_PPV_BP <- function(omega_D, PPV, BP) {
    omega_D * BP * (1 - PPV) / (PPV * (1 - omega_D))
  }
  
  figure2_NPV_from_BP_FPR <- function(omega_D, BP, FPR) {
    ((1 - omega_D) * (1 - FPR)) /
      ((1 - omega_D) * (1 - FPR) + omega_D * (1 - BP))
  }
  
  figure2_omega_at_NPV <- function(NPV_target, PPV, BP) {
    f <- function(omega_D) {
      FPR <- figure2_FPR_from_PPV_BP(
        omega_D = omega_D,
        PPV = PPV,
        BP = BP
      )
      NPV <- figure2_NPV_from_BP_FPR(
        omega_D = omega_D,
        BP = BP,
        FPR = FPR
      )
      NPV - NPV_target
    }
    stats::uniroot(f, interval = c(0.001, 0.999))$root
  }
  
  figure2_data <- tidyr::expand_grid(
    omega_D = figure2_omega,
    BP = figure2_BP
  ) %>%
    dplyr::mutate(
      PPV = figure2_PPV,
      FPR = figure2_FPR_from_PPV_BP(
        omega_D = omega_D,
        PPV = PPV,
        BP = BP
      ),
      NPV = figure2_NPV_from_BP_FPR(
        omega_D = omega_D,
        BP = BP,
        FPR = FPR
      ),
      valid =
        is.finite(FPR) &
        is.finite(NPV) &
        FPR >= 0 &
        FPR <= 1 &
        NPV >= 0 &
        NPV <= 1,
      FPR = dplyr::if_else(valid, FPR, NA_real_),
      NPV = dplyr::if_else(valid, NPV, NA_real_),
      BP_panel = factor(
        BP,
        levels = figure2_BP,
        labels = c("BP = 0.80", "BP = 0.90")
      )
    ) %>%
    dplyr::select(omega_D, BP_panel, FPR, NPV) %>%
    tidyr::pivot_longer(
      cols = c(FPR, NPV),
      names_to = "Quantity",
      values_to = "Value"
    ) %>%
    dplyr::mutate(
      Quantity = factor(Quantity, levels = c("FPR", "NPV"))
    )
  
  figure2_reference_values <- tibble::tibble(BP = figure2_BP) %>%
    dplyr::mutate(
      omega_NPV090 = purrr::map_dbl(
        BP,
        ~ figure2_omega_at_NPV(
          NPV_target = figure2_target_NPV,
          PPV = figure2_PPV,
          BP = .x
        )
      ),
      BP_panel = factor(
        BP,
        levels = figure2_BP,
        labels = c("BP = 0.80", "BP = 0.90")
      )
    )
  
  figure2_reference_points <- figure2_reference_values %>%
    dplyr::transmute(
      BP_panel,
      omega_D = omega_NPV090,
      Quantity = factor("NPV", levels = c("FPR", "NPV")),
      Value = figure2_target_NPV,
      label = sprintf("%.2f", omega_D)
    )
  
  figure2_equal_reference <- tibble::tibble(
    BP = figure2_BP,
    omega_D = predictive_reference_omega
  ) %>%
    dplyr::mutate(
      FPR = figure2_FPR_from_PPV_BP(
        omega_D = omega_D,
        PPV = figure2_PPV,
        BP = BP
      ),
      NPV = figure2_NPV_from_BP_FPR(
        omega_D = omega_D,
        BP = BP,
        FPR = FPR
      ),
      BP_panel = factor(
        BP,
        levels = figure2_BP,
        labels = c("BP = 0.80", "BP = 0.90")
      )
    )
  
  figure2_equal_reference_check <- figure2_equal_reference %>%
    dplyr::arrange(BP)
  
  if (nrow(figure2_equal_reference_check) != 2L ||
      any(abs(figure2_equal_reference_check$BP - c(0.80, 0.90)) > 1e-12)) {
    stop("Figure 2 reference-mixture validation could not identify both BP values.")
  }
  
  expected_figure2_FPR <- c(
    0.0205128205128205,
    0.0230769230769231
  )
  expected_figure2_NPV <- c(
    0.830434782608696,
    0.907142857142857
  )
  
  if (any(
    abs(
      figure2_equal_reference_check$FPR -
      expected_figure2_FPR
    ) > 1e-10
  ) ||
  any(
    abs(
      figure2_equal_reference_check$NPV -
      expected_figure2_NPV
    ) > 1e-10
  )) {
    stop(
      "Figure 2 reference-mixture values do not match ",
      "the values stated in the caption."
    )
  }
  
  figure2_equal_reference_points <- figure2_equal_reference %>%
    dplyr::select(
      BP_panel,
      omega_D,
      FPR,
      NPV
    ) %>%
    tidyr::pivot_longer(
      cols = c(FPR, NPV),
      names_to = "Quantity",
      values_to = "Value"
    ) %>%
    dplyr::mutate(
      Quantity = factor(
        Quantity,
        levels = c("FPR", "NPV")
      )
    )
  
  figure2_plot <- ggplot2::ggplot(
    figure2_data,
    ggplot2::aes(
      x = omega_D,
      y = Value,
      color = Quantity,
      linetype = Quantity
    )
  ) +
    ggplot2::geom_hline(
      yintercept = figure2_target_FPR,
      color = "grey35",
      linewidth = 0.45,
      linetype = "dotted"
    ) +
    ggplot2::geom_vline(
      xintercept = predictive_reference_omega,
      color = "grey25",
      linewidth = 0.50,
      linetype = "twodash"
    ) +
    ggplot2::geom_vline(
      data = figure2_reference_points,
      ggplot2::aes(xintercept = omega_D),
      inherit.aes = FALSE,
      color = "grey55",
      linewidth = 0.40,
      linetype = "dotted"
    ) +
    ggplot2::geom_line(
      linewidth = 1.05,
      lineend = "round",
      na.rm = TRUE
    ) +
    ggplot2::geom_point(
      data = figure2_reference_points,
      ggplot2::aes(
        x = omega_D,
        y = Value,
        color = Quantity
      ),
      inherit.aes = FALSE,
      shape = 15,
      size = 2.6
    ) +
    ggplot2::geom_text(
      data = figure2_reference_points,
      ggplot2::aes(
        x = omega_D,
        y = 0.925,
        label = label
      ),
      inherit.aes = FALSE,
      size = 3.1
    ) +
    ggplot2::geom_point(
      data = figure2_equal_reference_points,
      ggplot2::aes(
        x = omega_D,
        y = Value,
        color = Quantity
      ),
      inherit.aes = FALSE,
      shape = 16,
      size = 2.4
    ) +
    ggplot2::facet_wrap(
      ggplot2::vars(BP_panel),
      nrow = 1
    ) +
    ggplot2::scale_color_manual(
      values = c(
        "FPR" = "#D55E00",
        "NPV" = "#1F77B4"
      ),
      breaks = c("FPR", "NPV"),
      name = NULL
    ) +
    ggplot2::scale_linetype_manual(
      values = c(
        "FPR" = "dashed",
        "NPV" = "solid"
      ),
      breaks = c("FPR", "NPV"),
      name = NULL
    ) +
    ggplot2::scale_x_continuous(
      name = expression(
        "Design-prior probability of efficacy, " * omega[D]
      ),
      breaks = sort(
        unique(
          c(
            seq(0, 1, by = 0.2),
            predictive_reference_omega
          )
        )
      ),
      labels = scales::number_format(accuracy = 0.1),
      expand = ggplot2::expansion(mult = c(0, 0))
    ) +
    ggplot2::scale_y_continuous(
      name = "Probability",
      breaks = seq(0, 1, by = 0.1),
      labels = scales::percent_format(accuracy = 1),
      expand = ggplot2::expansion(mult = c(0, 0))
    ) +
    ggplot2::coord_cartesian(
      xlim = c(0, 1),
      ylim = c(0, 1)
    ) +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::theme(
      strip.background = ggplot2::element_rect(
        fill = "grey92",
        color = "grey55",
        linewidth = 0.5
      ),
      strip.text = ggplot2::element_text(
        size = 11.5,
        face = "bold"
      ),
      axis.title.x = ggplot2::element_text(
        size = 12,
        margin = ggplot2::margin(t = 10)
      ),
      axis.title.y = ggplot2::element_text(
        size = 12,
        margin = ggplot2::margin(r = 10)
      ),
      axis.text = ggplot2::element_text(size = 10),
      panel.grid = ggplot2::element_blank(),
      panel.spacing = grid::unit(1.0, "lines"),
      legend.position = "bottom",
      legend.direction = "horizontal",
      legend.text = ggplot2::element_text(size = 10.5),
      legend.key.width = grid::unit(1.6, "cm"),
      plot.margin = ggplot2::margin(10, 12, 10, 10)
    )
  
  figure2_pdf <- file.path(
    latex_dir,
    "BDCT_fixed_PPV_FPR_NPV.pdf"
  )
  figure2_tex <- file.path(
    latex_dir,
    "figureBDCT_fixed_PPV_FPR_NPV.tex"
  )
  
  ggplot2::ggsave(
    filename = figure2_pdf,
    plot = figure2_plot,
    width = 8.8,
    height = 4.8,
    units = "in",
    device = "pdf"
  )
  
  figure2_caption <- paste0(
    "FPR and NPV implied by fixing PPV at $97.5\\%$, shown as functions of ",
    "the design-prior probability of efficacy, $\\omega_{\\mathcal{D}}$, ",
    "for $\\mathrm{BP}=0.80$ and $0.90$. ",
    "Markers indicate the values of $\\omega_{\\mathcal{D}}$ at which NPV ",
    "reaches $90\\%$. ",
    "The vertical reference line denotes ",
    "$\\omega=0.5$, and the horizontal reference ",
    "line denotes $\\mathrm{FPR}=0.025$. ",
    "At the reference mixture, fixed PPV $97.5\\%$ implies FPR and NPV of ",
    "approximately $2.1\\%$ and $83.0\\%$ at BP $0.80$ and $2.3\\%$ and ",
    "$90.7\\%$ at BP $0.90$. ",
    "This relationship underlies the additional PPV constraint in joint ",
    "FPR--PPV calibration."
  )
  
  figure2_tex_lines <- c(
    "\\begin{figure}[!tbp]",
    "\\centering",
    "\\includegraphics[width=\\textwidth]{image/BDCT_fixed_PPV_FPR_NPV.pdf}",
    paste0(
      "\\caption{\\label{figureBDCT_fixed_PPV_FPR_NPV}",
      figure2_caption,
      "}"
    ),
    "\\end{figure}"
  )
  
  writeLines(
    figure2_tex_lines,
    con = figure2_tex,
    useBytes = TRUE
  )
  
  if (!file.exists(figure2_pdf)) {
    stop(
      "Optional Figure 2 PDF was not written: ",
      figure2_pdf
    )
  }
  if (!file.exists(figure2_tex)) {
    stop(
      "Optional Figure 2 LaTeX file was not written: ",
      figure2_tex
    )
  }
  
  copy_to_existing_subdir(
    figure2_tex,
    "figure"
  )
  copy_to_existing_subdir(
    figure2_pdf,
    "image"
  )
  
  cat(
    "\nOptional Figure 2 generated:\n",
    "  ", basename(figure2_tex), "\n",
    "  ", basename(figure2_pdf), "\n",
    sep = ""
  )
}