# ---- 1. input ---------------------------------------------------------

data_path  <- "data.xlsx"
output_dir <- "tte_results"


matching_ratio <- 1                
caliper_width  <- 0.20             
random_seed    <- 34590347

# Variable definitions
# composite: 30-day mortality or 30-day re-presentation (0/1)
# admitted: hospital admission (0/1)
# sex: female = 0, male = 1
# troponin: elevated troponin (0/1)
# htn_med: antihypertensive medication use (0/1)
# charlson_bin: 0 = 0, 1 = 1–2, 2 = 3–4, 3 = ≥5
# arrhythmia: arrhythmia or conduction disorder (0/1)
# meds_n_bin: 0 = 0, 1 = 1–5, 2 = 6–10, 3 = ≥11 medications
# cad: coronary artery disease (0/1)
# age: years
# sbp: systolic blood pressure
# hct: hematocrit, %
# creatinine: serum creatinine


dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
set.seed(random_seed)

# ---- 2. load & clean data ------------------------------------------------

dat_raw <- as.data.frame(
  readxl::read_excel(
    data_path,
    col_types = "text",
    .name_repair = "unique_quiet"
  ),
  stringsAsFactors = FALSE,
  check.names = FALSE
)

analysis_variables <- c(
  "composite", "1yr_mortality", "admitted",
  "sex", "troponin", "htn_med", "charlson_bin",
  "arrhythmia", "meds_n_bin", "cad",
  "age", "sbp", "hct", "creatinine"
)
required_variables <- c("risk", analysis_variables)


dat_raw$source_row_id <- seq_len(nrow(dat_raw))


for (v in required_variables) {
  dat_raw[[v]] <- suppressWarnings(as.numeric(as.character(dat_raw[[v]])))
}

# filter low risk
rows_before_risk_filter <- nrow(dat_raw)
dat_risk0 <- dat_raw[!is.na(dat_raw$risk) & dat_raw$risk == 0, , drop = FALSE]

check_allowed_values <- function(x, allowed, variable_name) {
  observed <- unique(x[!is.na(x)])
  invalid <- setdiff(as.character(observed), as.character(allowed))
  if (length(invalid) > 0) {
    stop(
      variable_name, " invalid ",
      paste(invalid, collapse = ", "),
      ". Allowed values: ", paste(allowed, collapse = ", "), "."
    )
  }
}

binary_variables <- c(
  "composite", "1yr_mortality", "admitted", "sex", "troponin",
  "htn_med", "arrhythmia", "cad"
)

for (v in binary_variables) {
  check_allowed_values(dat_risk0[[v]], c(0, 1), v)
}
check_allowed_values(dat_risk0$charlson_bin, 0:3, "charlson_bin")
check_allowed_values(dat_risk0$meds_n_bin, 0:3, "meds_n_bin")



# ---- 3. cohort and covariates --------------------

dat$admitted      <- as.integer(dat$admitted)
dat$composite     <- as.integer(dat$composite)
dat$`1yr_mortality` <- as.integer(dat$`1yr_mortality`)
dat$sex           <- factor(dat$sex, levels = c(0, 1), labels = c("Female", "Male"))
dat$troponin      <- factor(dat$troponin, levels = c(0, 1))
dat$htn_med       <- factor(dat$htn_med, levels = c(0, 1))
dat$charlson_bin  <- factor(dat$charlson_bin, levels = 0:3)
dat$arrhythmia    <- factor(dat$arrhythmia, levels = c(0, 1))
dat$meds_n_bin    <- factor(dat$meds_n_bin, levels = 0:3)
dat$cad           <- factor(dat$cad, levels = c(0, 1))
dat$log2_creatinine <- log2(dat$creatinine)

# Propensity score model
ps_formula <- admitted ~
  sex + troponin + htn_med + charlson_bin + arrhythmia +
  meds_n_bin + cad + age + sbp + hct + log2_creatinine

# ---- 4. Estimate propensity scores and match -----------------------

# Calculate propensity scores
ps_model <- glm(
  ps_formula,
  data = dat,
  family = binomial(link = "logit")
)
dat$propensity_score <- predict(ps_model, type = "response")
dat$logit_ps <- predict(ps_model, type = "link")

match_object <- MatchIt::matchit(
  formula  = ps_formula,
  data     = dat,
  method   = "nearest",
  distance = dat$logit_ps,
  estimand = "ATT",
  ratio    = matching_ratio,
  replace  = FALSE,
  caliper  = caliper_width,
  std.caliper = TRUE,
  m.order  = "random"
)

matched <- MatchIt::match.data(
  match_object,
  data = dat,
  drop.unmatched = TRUE
)


write.csv(
  matched,
  file.path(output_dir, "matched_dataset.csv"),
  row.names = FALSE
)

capture.output(
  summary(match_object, standardize = TRUE),
  file = file.path(output_dir, "matching_summary.txt")
)

# ---- 5. Covariate balance and propensity-score overlap -----------------

balance <- cobalt::bal.tab(
  match_object,
  un = TRUE,
  m.threshold = 0.10,
  binary = "std"
)

write.csv(
  data.frame(variable = rownames(balance$Balance), balance$Balance, row.names = NULL),
  file.path(output_dir, "covariate_balance.csv"),
  row.names = FALSE
)

png(
  file.path(output_dir, "love_plot.png"),
  width = 1800, height = 1400, res = 180
)
print(
  cobalt::love.plot(
    match_object,
    stats = "mean.diffs",
    abs = TRUE,
    thresholds = c(m = 0.10),
    var.order = "unadjusted",
    binary = "std",
    line = FALSE,
    colors = c("#B2182B", "#2166AC"),
    shapes = c("triangle", "circle")
  )
)
dev.off()

overlap_plot <- ggplot2::ggplot(
  matched,
  ggplot2::aes(
    x = distance,
    weight = weights,
    fill = factor(admitted)
  )
) +
  ggplot2::geom_density(alpha = 0.40) +
  ggplot2::scale_fill_manual(
    values = c("0" = "#2166AC", "1" = "#B2182B"),
    labels = c("Not admitted", "Admitted"),
    name = "Group"
  ) +
  ggplot2::labs(
    x = "Logit of the estimated propensity score",
    y = "Weighted density",
    title = "Propensity-score overlap in the matched sample"
  ) +
  ggplot2::theme_minimal(base_size = 12)

ggplot2::ggsave(
  file.path(output_dir, "propensity_score_overlap.png"),
  plot = overlap_plot,
  width = 7, height = 5, dpi = 300
)

# ---- 6. complete 1:1 matched pairs -----------------------------------

# The same matched cohort is used for the primary and secondary outcomes.
pair_sizes <- table(matched$subclass)
pair_treatment_counts <- tapply(
  matched$admitted,
  matched$subclass,
  function(x) sum(x == 1)
)

valid_pair_ids <- names(pair_sizes)[
  pair_sizes == 2 & pair_treatment_counts[names(pair_sizes)] == 1
]

if (length(valid_pair_ids) == 0) {
  stop("No complete 1:1 matched pairs are available for outcome analysis.")
}

# ---- 7. descriptives ----------------------------

weighted_mean <- function(x, w) sum(x * w) / sum(w)

analyze_matched_outcome <- function(outcome_name) {
  outcome_values <- matched[[outcome_name]]

  group_summary <- do.call(
    rbind,
    lapply(split(seq_len(nrow(matched)), matched$admitted), function(idx) {
      data.frame(
        outcome = outcome_name,
        admitted = unique(matched$admitted[idx]),
        n_rows = length(idx),
        effective_weight = sum(matched$weights[idx]),
        events = sum(outcome_values[idx]),
        weighted_risk = weighted_mean(outcome_values[idx], matched$weights[idx])
      )
    })
  )

  matched_pairs <- matched[
    as.character(matched$subclass) %in% valid_pair_ids,
    c("subclass", "admitted", outcome_name)
  ]
  treated <- matched_pairs[[outcome_name]][matched_pairs$admitted == 1]
  names(treated) <- as.character(
    matched_pairs$subclass[matched_pairs$admitted == 1]
  )
  control <- matched_pairs[[outcome_name]][matched_pairs$admitted == 0]
  names(control) <- as.character(
    matched_pairs$subclass[matched_pairs$admitted == 0]
  )

  pair_ids <- sort(intersect(names(treated), names(control)))
  pair_data <- data.frame(
    outcome = outcome_name,
    subclass = pair_ids,
    treated_outcome = as.integer(treated[pair_ids]),
    control_outcome = as.integer(control[pair_ids]),
    row.names = NULL
  )
  pair_data$pair_difference <-
    pair_data$treated_outcome - pair_data$control_outcome


  n_pairs <- nrow(pair_data)
  rd_estimate <- mean(pair_data$pair_difference)
  rd_se <- stats::sd(pair_data$pair_difference) / sqrt(n_pairs)
  rd_critical <- stats::qt(0.975, df = n_pairs - 1)
  rd_conf_low <- rd_estimate - rd_critical * rd_se
  rd_conf_high <- rd_estimate + rd_critical * rd_se
  if (is.na(rd_se) || rd_se == 0) {
    rd_p <- if (rd_estimate == 0) 1 else 0
  } else {
    rd_p <- 2 * stats::pt(
      abs(rd_estimate / rd_se),
      df = n_pairs - 1,
      lower.tail = FALSE
    )
  }


  b <- sum(pair_data$treated_outcome == 1 & pair_data$control_outcome == 0)
  c <- sum(pair_data$treated_outcome == 0 & pair_data$control_outcome == 1)

  if ((b + c) == 0) {
    matched_or_estimate <- NA_real_
    matched_or_conf_low <- NA_real_
    matched_or_conf_high <- NA_real_
    matched_or_p <- NA_real_
    warning(
      "There are no discordant matched pairs for '", outcome_name,
      "'; the matched odds ratio is undefined."
    )
  } else {

    b_for_or <- if (b == 0 || c == 0) b + 0.5 else b
    c_for_or <- if (b == 0 || c == 0) c + 0.5 else c
    matched_or_estimate <- b_for_or / c_for_or
    log_or_se <- sqrt(1 / b_for_or + 1 / c_for_or)
    matched_or_conf_low <- exp(
      log(matched_or_estimate) - stats::qnorm(0.975) * log_or_se
    )
    matched_or_conf_high <- exp(
      log(matched_or_estimate) + stats::qnorm(0.975) * log_or_se
    )
    matched_or_p <- stats::binom.test(
      x = b, n = b + c, p = 0.5, alternative = "two.sided"
    )$p.value
  }

  effect_results <- rbind(
    data.frame(
      outcome = outcome_name,
      estimand = "ATT",
      effect_measure = "Risk difference",
      estimate = rd_estimate,
      conf_low = rd_conf_low,
      conf_high = rd_conf_high,
      p_value = rd_p,
      discordant_treated_event = NA_integer_,
      discordant_control_event = NA_integer_
    ),
    data.frame(
      outcome = outcome_name,
      estimand = "ATT",
      effect_measure = "Matched odds ratio",
      estimate = matched_or_estimate,
      conf_low = matched_or_conf_low,
      conf_high = matched_or_conf_high,
      p_value = matched_or_p,
      discordant_treated_event = b,
      discordant_control_event = c
    )
  )

  list(
    group_summary = group_summary,
    pair_data = pair_data,
    effect_results = effect_results
  )
}

outcomes <- c("composite", "1yr_mortality")
outcome_analyses <- lapply(outcomes, analyze_matched_outcome)

group_summary <- do.call(
  rbind, lapply(outcome_analyses, function(x) x$group_summary)
)
pair_data <- do.call(
  rbind, lapply(outcome_analyses, function(x) x$pair_data)
)
effect_results <- do.call(
  rbind, lapply(outcome_analyses, function(x) x$effect_results)
)
rownames(group_summary) <- NULL
rownames(pair_data) <- NULL
rownames(effect_results) <- NULL

write.csv(
  group_summary,
  file.path(output_dir, "matched_outcome_summary.csv"),
  row.names = FALSE
)
write.csv(
  pair_data,
  file.path(output_dir, "matched_pair_outcomes.csv"),
  row.names = FALSE
)
write.csv(
  effect_results,
  file.path(output_dir, "treatment_effect_estimates.csv"),
  row.names = FALSE
)

# ---- 8. Print  ------------------------------------------------------

cat("\nMATCHED SAMPLE\n")
print(group_summary, row.names = FALSE)

cat("\nATT TREATMENT-EFFECT ESTIMATES\n")
print(effect_results, row.names = FALSE, digits = 4)

cat("\nFiles saved to:", normalizePath(output_dir), "\n")




sessionInfo()
