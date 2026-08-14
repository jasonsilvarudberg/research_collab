library(readxl)
library(writexl)
library(dplyr)

file_path <- "data.xlsx"
df <- read_excel(file_path)

# 1. Calculate LASTRAIN quartile cut points
# Filter out NA values for tail trimming and quartile calculation
lastrain_valid <- df$LASTRAIN[!is.na(df$LASTRAIN)]

# exclude extreme 0.5% tails
tail_low  <- quantile(lastrain_valid, 0.005)
tail_high <- quantile(lastrain_valid, 0.995)

trimmed_lastrain <- lastrain_valid[lastrain_valid >= tail_low & lastrain_valid <= tail_high]

# calculate quartile cut points on trimmed cohort
cuts <- quantile(trimmed_lastrain, probs = c(0.25, 0.50, 0.75))

# assign full cohort into Q1-Q4 based on cut points
df <- df %>%
  mutate(
    LASTRAIN_Q = case_when(
      is.na(LASTRAIN) ~ NA_character_,
      LASTRAIN <= cuts[1] ~ "Q1",
      LASTRAIN <= cuts[2] ~ "Q2",
      LASTRAIN <= cuts[3] ~ "Q3",
      TRUE ~ "Q4"
    )
  )

# variables
categorical_vars <- c(
  "SEX_V1", "EDUC_V3", "BKGRD_V1", "SMOKER_V3", "ALCOHOL_V3",
  "DIABETES_V3", "HYPERTENSION_NHANES_V3", "HTN_MED_V3", "STATIN_V3",
  "STROKE", "ANYCMB"
)

welch_vars <- c(
  "AGE_V3", "BMI_V3", "HEART_RATE", "SBP", "LDL_V3", "HDL_V3",
  "GFRSCR_V3", "LAVMIN_INDEX", "LVEDM_INDEX", "LVEF", "LVGLS",
  "LOG2_CEREBRUMV2_WMH", "FA", "GCS", "STD_SEVLT_SUM",
  "STD_SEVLT_RECALL", "STD_WF", "STD_DSS"
)

kruskal_vars <- c(
  "TG_V3", "log_prontbnp"
)

# 2. Output
results <- list()

# First Row: N
n_overall <- nrow(df)
n_q1 <- sum(df$LASTRAIN_Q == "Q1", na.rm = TRUE)
n_q2 <- sum(df$LASTRAIN_Q == "Q2", na.rm = TRUE)
n_q3 <- sum(df$LASTRAIN_Q == "Q3", na.rm = TRUE)
n_q4 <- sum(df$LASTRAIN_Q == "Q4", na.rm = TRUE)

results[[length(results) + 1]] <- data.frame(
  Characteristic = "N",
  Overall = as.character(n_overall),
  Q1 = as.character(n_q1),
  Q2 = as.character(n_q2),
  Q3 = as.character(n_q3),
  Q4 = as.character(n_q4),
  `P value` = "",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

# --- categorical variables ---
for (var in categorical_vars) {
  if (!var %in% colnames(df)) next
  
  # chi-square Tests
  sub_df <- df %>% filter(!is.na(.data[[var]]), !is.na(LASTRAIN_Q))
  p_val <- tryCatch({
    as.character(chisq.test(table(sub_df[[var]], sub_df$LASTRAIN_Q))$p.value)
  }, error = function(e) "")
  
  # header row
  results[[length(results) + 1]] <- data.frame(
    Characteristic = var,
    Overall = "", Q1 = "", Q2 = "", Q3 = "", Q4 = "",
    `P value` = p_val,
    check.names = FALSE, stringsAsFactors = FALSE
  )
  
  # category levels
  levels_val <- sort(unique(na.omit(df[[var]])))
  for (lvl in levels_val) {
    # Overall
    total_valid <- sum(!is.na(df[[var]]))
    n_tot <- sum(df[[var]] == lvl, na.rm = TRUE)
    pct_tot <- ifelse(total_valid > 0, (n_tot / total_valid) * 100, 0)
    str_tot <- paste0(n_tot, " (", pct_tot, "%)")
    
    # quartiles
    q_strs <- c()
    for (q in c("Q1", "Q2", "Q3", "Q4")) {
      q_sub <- df %>% filter(LASTRAIN_Q == q)
      q_valid <- sum(!is.na(q_sub[[var]]))
      n_q <- sum(q_sub[[var]] == lvl, na.rm = TRUE)
      pct_q <- ifelse(q_valid > 0, (n_q / q_valid) * 100, 0)
      q_strs <- c(q_strs, paste0(n_q, " (", pct_q, "%)"))
    }
    
    results[[length(results) + 1]] <- data.frame(
      Characteristic = paste0("  ", lvl),
      Overall = str_tot,
      Q1 = q_strs[1],
      Q2 = q_strs[2],
      Q3 = q_strs[3],
      Q4 = q_strs[4],
      `P value` = "",
      check.names = FALSE, stringsAsFactors = FALSE
    )
  }
}

# --- welch anova ---
for (var in welch_vars) {
  if (!var %in% colnames(df)) next
  
  # Welch ANOVA (excluding NAs)
  sub_df <- df %>% filter(!is.na(.data[[var]]), !is.na(LASTRAIN_Q))
  p_val <- tryCatch({
    as.character(oneway.test(as.formula(paste(var, "~ LASTRAIN_Q")), data = sub_df, var.equal = FALSE)$p.value)
  }, error = function(e) "")
  
  # Overall Mean (SD)
  x_tot <- na.omit(df[[var]])
  str_tot <- paste0(mean(x_tot), " (", sd(x_tot), ")")
  
  # Quartile Mean (SD)
  q_strs <- c()
  for (q in c("Q1", "Q2", "Q3", "Q4")) {
    x_q <- na.omit(df[[var]][df$LASTRAIN_Q == q])
    if (length(x_q) > 0) {
      q_strs <- c(q_strs, paste0(mean(x_q), " (", sd(x_q), ")"))
    } else {
      q_strs <- c(q_strs, "N/A")
    }
  }
  
  results[[length(results) + 1]] <- data.frame(
    Characteristic = paste0(var, ", Mean (SD)"),
    Overall = str_tot,
    Q1 = q_strs[1], Q2 = q_strs[2], Q3 = q_strs[3], Q4 = q_strs[4],
    `P value` = p_val,
    check.names = FALSE, stringsAsFactors = FALSE
  )
}

# --- Kruskal-wallis ---
for (var in kruskal_vars) {
  if (!var %in% colnames(df)) next
  
  # Kruskal-Wallis Test
  sub_df <- df %>% filter(!is.na(.data[[var]]), !is.na(LASTRAIN_Q))
  p_val <- tryCatch({
    as.character(kruskal.test(as.formula(paste(var, "~ LASTRAIN_Q")), data = sub_df)$p.value)
  }, error = function(e) "")
  
  # Overall Median [Q1, Q3]
  x_tot <- na.omit(df[[var]])
  q_tot <- quantile(x_tot, probs = c(0.25, 0.5, 0.75))
  str_tot <- paste0(q_tot[2], " [", q_tot[1], ", ", q_tot[3], "]")
  
  # Quartile Median [Q1, Q3]
  q_strs <- c()
  for (q in c("Q1", "Q2", "Q3", "Q4")) {
    x_q <- na.omit(df[[var]][df$LASTRAIN_Q == q])
    if (length(x_q) > 0) {
      q_sub_vals <- quantile(x_q, probs = c(0.25, 0.5, 0.75))
      q_strs <- c(q_strs, paste0(q_sub_vals[2], " [", q_sub_vals[1], ", ", q_sub_vals[3], "]"))
    } else {
      q_strs <- c(q_strs, "N/A")
    }
  }
  
  results[[length(results) + 1]] <- data.frame(
    Characteristic = paste0(var, ", Median [Q1, Q3]"),
    Overall = str_tot,
    Q1 = q_strs[1], Q2 = q_strs[2], Q3 = q_strs[3], Q4 = q_strs[4],
    `P value` = p_val,
    check.names = FALSE, stringsAsFactors = FALSE
  )
}

final_table <- bind_rows(results)

output_path <- "table1.xlsx"
write_xlsx(final_table, path = output_path)
