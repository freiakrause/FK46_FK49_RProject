###############################################################################
# FK49_Metabolomics_2_Statistics.R
#
# Loads preprocessed RDS files (one per dataset), runs per-metabolite
# statistics, and saves one combined CSV + RDS.
#
# Two analysis modes (set `analysis` at top):
#   "CDHFD" - FK49, both sexes, model value ~ Treatment * Sex
#             censored metabolites: cen2way (NADA2) with interact = TRUE
#   "ND"    - BH, females only, model value ~ Treatment
#             censored metabolites: cens1way (NADA2)
#
# Uncensored metabolites: lm on log2-transformed values with emmeans
#   (pairwise Treatment contrast, Cohen's d, sex-stratified contrasts).
# Censored metabolites: cen2way/cens1way on raw LOD-imputed values with
#   LOG = TRUE; internal survreg reproduced for GMR and CI.
#
# FDR (BH) applied within each dataset separately.
###############################################################################

rm(list = ls())
gc()

library(tidyverse)
library(emmeans)
library(NADA2)
library(survival)
source("FK49_Definitions.R")

# ============================================================
# CONFIGURATION  -- change analysis to "CDHFD" or "ND"
# ============================================================

analysis <- "CDHFD"

rawdata_pwd <- PATHS$metabolomics$rawdata
output_pwd  <- PATHS$metabolomics$output

analysis_folder <- if (analysis == "CDHFD") "CDHFD" else "ND"
analysis_pwd <- file.path(output_pwd, analysis_folder)
if (!dir.exists(analysis_pwd)) dir.create(analysis_pwd, recursive = TRUE)

# Analysis-specific settings
if (analysis == "CDHFD") {
  diet_filter  <- "CDHFD13"
  expid_filter <- "FK49"
  has_sex      <- TRUE
} else if (analysis == "ND") {
  diet_filter  <- "ND"
  expid_filter <- "BH"
  has_sex      <- FALSE
} else {
  stop('analysis must be "CDHFD" or "ND"')
}

# ============================================================
# LOAD PREPROCESSED RDS
# ============================================================

pos <- readRDS(file.path(rawdata_pwd, "FK49_metabolome_positive_processed.rds"))
neg <- readRDS(file.path(rawdata_pwd, "FK49_metabolome_negative_processed.rds"))
tar <- readRDS(file.path(rawdata_pwd, "FK49_metabolome_targeted_processed.rds"))

datasets <- list(positive = pos, negative = neg, targeted = tar)

# ============================================================
# HELPER: cen2way path (CDHFD, two factors)
# ============================================================

run_cen2way <- function(value_raw, cens, Treatment, Sex) {
  # --- cen2way for p-values ---
  cen_result <- capture.output(
    suppressWarnings(
      cen2way(
        ifelse(cens, value_raw * 2, value_raw),
        cens,
        Treatment,
        Sex,
        LOG = TRUE,
        interact = TRUE
      )
    )
  )

  p_lines <- cen_result[grep("Treatment|Sex|interaction", cen_result)]
  p_values <- as.numeric(sub(".*\\s([0-9]+\\.[0-9]+)$", "\\1", p_lines))

  # p_values[3] = Treatment, [4] = Sex, [5] = interaction
  p_Treatment      <- p_values[3]
  p_Sex            <- p_values[4]
  p_Treatment_Sex  <- p_values[5]

  # --- Reproduce internal survreg for coefficients and GMR ---
  y1 <- ifelse(cens, value_raw * 2, value_raw)
  y2 <- cens
  fac1 <- factor(Treatment)
  fac2 <- factor(Sex)
  keep <- complete.cases(y1, y2, fac1, fac2)
  y1 <- y1[keep]; y2 <- y2[keep]; fac1 <- fac1[keep]; fac2 <- fac2[keep]

  levels_1 <- levels(fac1)  # Ctrl, TAM
  levels_2 <- levels(fac2)  # female, male

  # Effect coding: first level = +1, last level = -1
  e <- ifelse(fac1 == levels_1[1], 1, -1)
  d <- ifelse(fac2 == levels_2[1], 1, -1)
  int <- e * d

  # Censored log-normal model (same as cen2way internal)
  lnvar     <- log(y1)
  fconst    <- max(lnvar)
  flip.log  <- fconst + 1 - lnvar
  detect    <- !y2
  logCensData <- survival::Surv(flip.log, detect, type = "right")
  cen_model <- survival::survreg(logCensData ~ e + d + int, dist = "gaussian")

  # Reverse the flipping (same as cen2way)
  beta <- coef(cen_model)
  beta <- -beta
  beta[1] <- fconst + 1 + beta[1]

  b0    <- beta["(Intercept)"]
  b_tr  <- beta["e"]
  b_sex <- beta["d"]
  b_int <- beta["int"]

  # Marginal geometric means for Treatment (averaged over Sex)
  # Ctrl = +1 (first level), TAM = -1 (second level)
  # female = +1 (first level), male = -1 (second level)
  log_GM_Ctrl_Female <- b0 + b_tr + b_sex + b_int
  log_GM_Ctrl_Male   <- b0 + b_tr - b_sex - b_int
  log_GM_TAM_Female   <- b0 - b_tr + b_sex - b_int
  log_GM_TAM_Male     <- b0 - b_tr - b_sex + b_int

  log_GM_Ctrl <- mean(c(log_GM_Ctrl_Female, log_GM_Ctrl_Male))
  log_GM_TAM  <- mean(c(log_GM_TAM_Female,  log_GM_TAM_Male))

  GM_ctrl <- exp(log_GM_Ctrl)
  GM_tam  <- exp(log_GM_TAM)

  # GMR = TAM / Ctrl
  log_GMR <- -2 * b_tr
  GMR     <- exp(log_GMR)

  se_b_tr     <- sqrt(vcov(cen_model)["e", "e"])
  se_log_GMR  <- 2 * se_b_tr
  CI_low      <- exp(log_GMR - 1.96 * se_log_GMR)
  CI_high     <- exp(log_GMR + 1.96 * se_log_GMR)

  list(
    p_Treatment     = p_Treatment,
    p_Sex           = p_Sex,
    p_Treatment_Sex = p_Treatment_Sex,
    mean_Ctrl       = GM_ctrl,
    mean_TAM        = GM_tam,
    effect_size     = GMR,
    effect_CI_low   = CI_low,
    effect_CI_high  = CI_high,
    effect_size_type = "GMR",
    fold_change     = GMR,
    log2FC          = log2(GMR)
  )
}

# ============================================================
# HELPER: cens1way path (ND, single factor)
# ============================================================

run_cens1way <- function(value_raw, cens, Treatment) {
  # --- cens1way for p-value ---
  # Same censoring adjustment as cen2way: censored values doubled
  # (imputed value was min/2, doubling gives back min = censoring limit)
  cen_vals <- ifelse(cens, value_raw * 2, value_raw)
  cen_result <- capture.output(
    suppressWarnings(
      cens1way(cen_vals, cens, Treatment, LOG = TRUE)
    )
  )

  # Parse p-value from output
  p_lines <- cen_result[grep("p-value|Treatment", cen_result, ignore.case = TRUE)]
  p_values <- as.numeric(sub(".*\\s([0-9]+\\.[0-9]+)$", "\\1", p_lines))
  p_Treatment <- if (length(p_values) > 0 && !is.na(p_values[length(p_values)])) {
    p_values[length(p_values)]
  } else NA_real_

  # --- Reproduce survreg for GMR ---
  y1 <- ifelse(cens, value_raw * 2, value_raw)
  y2 <- cens
  fac1 <- factor(Treatment)
  keep <- complete.cases(y1, y2, fac1)
  y1 <- y1[keep]; y2 <- y2[keep]; fac1 <- fac1[keep]

  levels_1 <- levels(fac1)  # Ctrl, TAM
  e <- ifelse(fac1 == levels_1[1], 1, -1)

  lnvar     <- log(y1)
  fconst    <- max(lnvar)
  flip.log  <- fconst + 1 - lnvar
  detect    <- !y2
  logCensData <- survival::Surv(flip.log, detect, type = "right")
  cen_model <- survival::survreg(logCensData ~ e, dist = "gaussian")

  beta <- coef(cen_model)
  beta <- -beta
  beta[1] <- fconst + 1 + beta[1]

  b0   <- beta["(Intercept)"]
  b_tr <- beta["e"]

  log_GM_Ctrl <- b0 + b_tr
  log_GM_TAM  <- b0 - b_tr

  GM_ctrl <- exp(log_GM_Ctrl)
  GM_tam  <- exp(log_GM_TAM)

  log_GMR <- -2 * b_tr
  GMR     <- exp(log_GMR)

  se_b_tr    <- sqrt(vcov(cen_model)["e", "e"])
  se_log_GMR <- 2 * se_b_tr
  CI_low     <- exp(log_GMR - 1.96 * se_log_GMR)
  CI_high    <- exp(log_GMR + 1.96 * se_log_GMR)

  # Fallback p-value via LRT if parsing failed
  if (is.na(p_Treatment)) {
    cen_null <- survival::survreg(logCensData ~ 1, dist = "gaussian")
    lr_stat  <- as.numeric(2 * (logLik(cen_model) - logLik(cen_null)))
    p_Treatment <- pchisq(lr_stat, df = 1, lower.tail = FALSE)
  }

  list(
    p_Treatment      = p_Treatment,
    mean_Ctrl        = GM_ctrl,
    mean_TAM         = GM_tam,
    effect_size      = GMR,
    effect_CI_low    = CI_low,
    effect_CI_high   = CI_high,
    effect_size_type = "GMR",
    fold_change      = GMR,
    log2FC           = log2(GMR)
  )
}

# ============================================================
# HELPER: LM path (uncensored)
# ============================================================

run_lm <- function(value_log, Treatment, Sex, use_sex, Batch = NULL, use_batch = FALSE) {
  if (use_sex && use_batch) {
    df_metab <- data.frame(value_log = value_log, Treatment = Treatment,
                           Sex = Sex, Batch = factor(Batch))
    fit <- tryCatch(lm(value_log ~ Treatment * Sex + Batch, data = df_metab),
                    error = function(e) NULL)
  } else if (use_sex) {
    df_metab <- data.frame(value_log = value_log, Treatment = Treatment, Sex = Sex)
    fit <- tryCatch(lm(value_log ~ Treatment * Sex, data = df_metab),
                    error = function(e) NULL)
  } else if (use_batch) {
    df_metab <- data.frame(value_log = value_log, Treatment = Treatment,
                           Batch = factor(Batch))
    fit <- tryCatch(lm(value_log ~ Treatment + Batch, data = df_metab),
                    error = function(e) NULL)
  } else {
    df_metab <- data.frame(value_log = value_log, Treatment = Treatment)
    fit <- tryCatch(lm(value_log ~ Treatment, data = df_metab),
                    error = function(e) NULL)
  }

  if (is.null(fit)) {
    return(list(method = "model_failed"))
  }

  # ANOVA p-values
  anova_result <- tryCatch(as.data.frame(anova(fit)), error = function(e) NULL)
  p_Treatment <- p_Sex <- p_Treatment_Sex <- p_Batch <- NA_real_

  if (!is.null(anova_result)) {
    p_Treatment <- anova_result$`Pr(>F)`[rownames(anova_result) == "Treatment"]
    if (use_batch) {
      p_Batch <- anova_result$`Pr(>F)`[rownames(anova_result) == "Batch"]
    }
    if (use_sex) {
      p_Sex           <- anova_result$`Pr(>F)`[rownames(anova_result) == "Sex"]
      p_Treatment_Sex <- anova_result$`Pr(>F)`[rownames(anova_result) == "Treatment:Sex"]
    }
  }

  # Means and SDs (on log2 scale)
  mean_Ctrl <- mean(value_log[Treatment == "Ctrl"], na.rm = TRUE)
  sd_Ctrl   <- sd(value_log[Treatment == "Ctrl"], na.rm = TRUE)
  mean_TAM  <- mean(value_log[Treatment == "TAM"],  na.rm = TRUE)
  sd_TAM    <- sd(value_log[Treatment == "TAM"],  na.rm = TRUE)

  log2FC    <- mean_TAM - mean_Ctrl
  fold_change <- 2^log2FC

  # emmeans pairwise Treatment contrast
  emm_treatment <- tryCatch(emmeans(fit, ~ Treatment), error = function(e) NULL)
  con_treatment <- if (!is.null(emm_treatment))
    tryCatch(as.data.frame(contrast(emm_treatment, method = "pairwise",
                                    adjust = "none", infer = TRUE)),
             error = function(e) NULL) else NULL

  if (!is.null(con_treatment) && nrow(con_treatment) > 0) {
    p_Treatment <- con_treatment$p.value[1]
  }

  # Cohen's d via eff_size (negated so positive = TAM higher)
  effect_size <- effect_CI_low <- effect_CI_high <- NA_real_
  effect_size_type <- NA_character_

  if (!is.null(emm_treatment)) {
    eff <- tryCatch(eff_size(emm_treatment, sigma = sigma(fit),
                             edf = df.residual(fit)), error = function(e) NULL)
    if (!is.null(eff)) {
      eff_summary <- as.data.frame(summary(eff))
      effect_size     <- -eff_summary$effect.size[1]
      effect_CI_low   <- -eff_summary$lower.CL[1]
      effect_CI_high  <- -eff_summary$upper.CL[1]
      effect_size_type <- "Cohen's d"
    }
  }

  # Sex-stratified contrasts (only when use_sex)
  p_Treatment_female <- treatment_effect_female <- NA_real_
  CI_low_female <- CI_high_female <- NA_real_
  p_Treatment_male <- treatment_effect_male <- NA_real_
  CI_low_male <- CI_high_male <- NA_real_

  if (use_sex) {
    emm_sex <- tryCatch(emmeans(fit, ~ Treatment | Sex), error = function(e) NULL)
    con_sex <- if (!is.null(emm_sex))
      tryCatch(as.data.frame(contrast(emm_sex, method = "pairwise",
                                      adjust = "none", infer = TRUE)),
               error = function(e) NULL) else NULL

    if (!is.null(con_sex) && nrow(con_sex) > 0) {
      cond_female <- con_sex$Sex == "female"
      cond_male   <- con_sex$Sex == "male"

      if (any(cond_female)) {
        p_Treatment_female     <- con_sex$p.value[cond_female][1]
        treatment_effect_female <- con_sex$estimate[cond_female][1]
        CI_low_female          <- con_sex$lower.CL[cond_female][1]
        CI_high_female         <- con_sex$upper.CL[cond_female][1]
      }
      if (any(cond_male)) {
        p_Treatment_male     <- con_sex$p.value[cond_male][1]
        treatment_effect_male <- con_sex$estimate[cond_male][1]
        CI_low_male          <- con_sex$lower.CL[cond_male][1]
        CI_high_male         <- con_sex$upper.CL[cond_male][1]
      }
    }
  }

  method <- if (use_sex) "linear_model_Treatment_x_Sex" else "linear_model_Treatment"

  list(
    method                = method,
    p_Treatment           = p_Treatment,
    p_Sex                 = p_Sex,
    p_Treatment_Sex       = p_Treatment_Sex,
    p_Batch               = p_Batch,
    p_Treatment_female    = p_Treatment_female,
    treatment_effect_female = treatment_effect_female,
    CI_low_female         = CI_low_female,
    CI_high_female        = CI_high_female,
    p_Treatment_male      = p_Treatment_male,
    treatment_effect_male = treatment_effect_male,
    CI_low_male           = CI_low_male,
    CI_high_male          = CI_high_male,
    mean_Ctrl             = mean_Ctrl,
    sd_Ctrl               = sd_Ctrl,
    mean_TAM              = mean_TAM,
    sd_TAM                = sd_TAM,
    log2FC                = log2FC,
    fold_change           = fold_change,
    effect_size           = effect_size,
    effect_CI_low         = effect_CI_low,
    effect_CI_high        = effect_CI_high,
    effect_size_type      = effect_size_type
  )
}

# ============================================================
# MAIN LOOP
# ============================================================

all_results <- list()

for (ds_name in names(datasets)) {
  ds <- datasets[[ds_name]]
  cat("\n\n========== Dataset:", ds_name, "==========\n")

  # Filter to analysis subset
  idx <- ds$metadata$Diet == diet_filter & ds$metadata$ExpID == expid_filter
  metadata   <- ds$metadata[idx, , drop = FALSE]
  raw_values <- ds$raw_values[idx, , drop = FALSE]
  log_values <- ds$log_values[idx, , drop = FALSE]
  censored   <- ds$censored[idx, , drop = FALSE]

  metab_names <- colnames(raw_values)
  results <- list()

  # Batch covariate: include only when > 1 level present
  Batch      <- factor(metadata$Batch)
  use_batch  <- n_distinct(Batch) > 1

  for (metab in metab_names) {
    cat(".")

    value_raw <- raw_values[, metab]
    value_log <- log_values[, metab]
    cens      <- censored[, metab]

    Treatment <- factor(metadata$Treatment, levels = c("Ctrl", "TAM"))
    Sex       <- factor(metadata$Sex, levels = c("female", "male"))

    n_obs      <- length(value_raw)
    n_animals  <- n_distinct(metadata$Animal)
    n_ctrl     <- sum(Treatment == "Ctrl", na.rm = TRUE)
    n_tam      <- sum(Treatment == "TAM",  na.rm = TRUE)
    n_female   <- sum(Sex == "female", na.rm = TRUE)
    n_male     <- sum(Sex == "male",   na.rm = TRUE)
    n_censored <- sum(cens, na.rm = TRUE)

    # Initialise result row
    result <- data.frame(
      Dataset = ds_name, Metabolite = metab,
      n = n_obs, n_Ctrl = n_ctrl, n_TAM = n_tam,
      n_female = n_female, n_male = n_male, n_censored = n_censored,
      method = NA_character_,
      p_Treatment = NA_real_, p_Sex = NA_real_, p_Treatment_Sex = NA_real_,
      p_Batch = NA_real_,
      p_Treatment_female = NA_real_, treatment_effect_female = NA_real_,
      CI_low_female = NA_real_, CI_high_female = NA_real_,
      p_Treatment_male = NA_real_, treatment_effect_male = NA_real_,
      CI_low_male = NA_real_, CI_high_male = NA_real_,
      mean_Ctrl = NA_real_, sd_Ctrl = NA_real_,
      mean_TAM = NA_real_, sd_TAM = NA_real_,
      log2FC = NA_real_, fold_change = NA_real_,
      effect_size = NA_real_, effect_CI_low = NA_real_,
      effect_CI_high = NA_real_, effect_size_type = NA_character_,
      stringsAsFactors = FALSE
    )

    # Check sufficient data
    if (n_obs < 4 || n_animals < 2 || n_distinct(Treatment) < 2) {
      result$method <- "insufficient_data"
      results[[metab]] <- result
      next
    }

    # Determine whether Sex can be used (need >= 2 levels)
    use_sex <- has_sex && n_distinct(Sex) >= 2

    if (n_censored > 0) {
      # --- Censored path ---
      # Check sufficient uncensored data
      uncensored_vals <- value_raw[!cens]
      sufficient <- length(uncensored_vals) >= 4

      if (sufficient) {
        if (use_sex) {
          cen_out <- tryCatch(run_cen2way(value_raw, cens, Treatment, Sex),
                              error = function(e) NULL)
          if (!is.null(cen_out)) {
            result$method           <- "cen2way"
            result$p_Treatment      <- cen_out$p_Treatment
            result$p_Sex            <- cen_out$p_Sex
            result$p_Treatment_Sex  <- cen_out$p_Treatment_Sex
            result$mean_Ctrl        <- cen_out$mean_Ctrl
            result$mean_TAM         <- cen_out$mean_TAM
            result$effect_size      <- cen_out$effect_size
            result$effect_CI_low    <- cen_out$effect_CI_low
            result$effect_CI_high   <- cen_out$effect_CI_high
            result$effect_size_type <- cen_out$effect_size_type
            result$fold_change      <- cen_out$fold_change
            result$log2FC           <- cen_out$log2FC
          } else {
            result$method <- "cen2way_failed"
          }
        } else {
          cen_out <- tryCatch(run_cens1way(value_raw, cens, Treatment),
                              error = function(e) NULL)
          if (!is.null(cen_out)) {
            result$method           <- "cens1way"
            result$p_Treatment      <- cen_out$p_Treatment
            result$mean_Ctrl        <- cen_out$mean_Ctrl
            result$mean_TAM         <- cen_out$mean_TAM
            result$effect_size      <- cen_out$effect_size
            result$effect_CI_low    <- cen_out$effect_CI_low
            result$effect_CI_high   <- cen_out$effect_CI_high
            result$effect_size_type <- cen_out$effect_size_type
            result$fold_change      <- cen_out$fold_change
            result$log2FC           <- cen_out$log2FC
          } else {
            result$method <- "cens1way_failed"
          }
        }
      } else {
        result$method <- "insufficient_uncensored"
      }
    } else {
      # --- Uncensored path (LM) ---
      lm_out <- run_lm(value_log, Treatment, Sex, use_sex, Batch, use_batch)
      if (lm_out$method == "model_failed") {
        result$method <- "model_failed"
      } else {
        result$method                <- lm_out$method
        result$p_Treatment           <- lm_out$p_Treatment
        result$p_Sex                 <- lm_out$p_Sex
        result$p_Treatment_Sex       <- lm_out$p_Treatment_Sex
        result$p_Batch               <- lm_out$p_Batch
        result$p_Treatment_female    <- lm_out$p_Treatment_female
        result$treatment_effect_female <- lm_out$treatment_effect_female
        result$CI_low_female         <- lm_out$CI_low_female
        result$CI_high_female        <- lm_out$CI_high_female
        result$p_Treatment_male      <- lm_out$p_Treatment_male
        result$treatment_effect_male <- lm_out$treatment_effect_male
        result$CI_low_male           <- lm_out$CI_low_male
        result$CI_high_male          <- lm_out$CI_high_male
        result$mean_Ctrl             <- lm_out$mean_Ctrl
        result$sd_Ctrl               <- lm_out$sd_Ctrl
        result$mean_TAM              <- lm_out$mean_TAM
        result$sd_TAM                <- lm_out$sd_TAM
        result$log2FC                <- lm_out$log2FC
        result$fold_change           <- lm_out$fold_change
        result$effect_size           <- lm_out$effect_size
        result$effect_CI_low         <- lm_out$effect_CI_low
        result$effect_CI_high        <- lm_out$effect_CI_high
        result$effect_size_type      <- lm_out$effect_size_type
      }
    }

    results[[metab]] <- result
  }

  cat("\n")

  # FDR within dataset
  results_df <- bind_rows(results)
  results_df <- results_df %>%
    mutate(
      adj_p_Treatment        = p.adjust(p_Treatment,        method = "fdr"),
      adj_p_Sex              = p.adjust(p_Sex,              method = "fdr"),
      adj_p_Treatment_Sex    = p.adjust(p_Treatment_Sex,    method = "fdr"),
      adj_p_Batch            = p.adjust(p_Batch,            method = "fdr"),
      adj_p_Treatment_female = p.adjust(p_Treatment_female, method = "fdr"),
      adj_p_Treatment_male   = p.adjust(p_Treatment_male,   method = "fdr")
    )

  all_results[[ds_name]] <- results_df
}

# ============================================================
# COMBINE AND ANNOTATE
# ============================================================

final_results <- bind_rows(all_results)

final_results <- final_results %>%
  mutate(
    significant = !is.na(adj_p_Treatment) & adj_p_Treatment < 0.05 &
                  abs(log2FC) > 0.5,
    trend = !is.na(adj_p_Treatment) & adj_p_Treatment < 0.1 &
            abs(log2FC) > 0.3,
    direction = case_when(
      significant & log2FC > 0  ~ "UP",
      significant & log2FC < 0  ~ "DOWN",
      TRUE                      ~ "NS"
    )
  )

# ============================================================
# SAVE
# ============================================================

csv_path <- file.path(analysis_pwd, "FK49_metabolome_statistics.csv")
rds_path <- file.path(analysis_pwd, "FK49_metabolome_statistics.rds")

write.csv2(final_results, file = csv_path, row.names = FALSE)
saveRDS(final_results, file = rds_path)

cat("\n=== Statistics complete (", analysis, ") ===\n")
cat("CSV written to:", csv_path, "\n")
cat("RDS written to:", rds_path, "\n")

for (ds_name in names(all_results)) {
  df <- all_results[[ds_name]]
  cat("  ", ds_name, ":", sum(df$method == "linear_model_Treatment_x_Sex" |
                              df$method == "linear_model_Treatment"),
      "LM,", sum(df$method == "cen2way" | df$method == "cens1way"),
      "censored,", sum(df$method == "insufficient_data" |
                       df$method == "insufficient_uncensored" |
                       df$method == "model_failed" |
                       df$method == "cen2way_failed" |
                       df$method == "cens1way_failed"), "skipped\n")
}
