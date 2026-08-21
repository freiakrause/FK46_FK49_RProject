gc()
rm(list = ls())

library(tidyverse)
library(lme4)
library(lmerTest)
library(emmeans)
library(NADA2)
source("FK49_Definitions.R")

ExpID <- "FK49"
output_pwd <- PATHS$BA$output
BA_filtered <- readRDS(file = paste0(dirname(dirname(output_pwd)), "/01_RawData/FK49_BA_preprocessed.rds"))
#BA_filtered <- readRDS(file = paste0(dirname(dirname(output_pwd)), "/01_RawData/FK49_BAfiltered_preprocessed.rds"))
BA_to_test <- PARAMETERS$BA$BA_sort

results <- list()

# ============================================================
# LOOP OVER BILE ACIDS
# ============================================================

for (ba in BA_to_test) {
  cat("\n\n", ba, "\n")
  censored_col <- paste0(ba, "_censored")
  direction_col <- paste0(ba, "_direction")
  df_ba <- BA_filtered %>%
    mutate(censored = as.logical(.data[[censored_col]]),
           direction = as.character(.data[[direction_col]]),
           cens_logical = censored) %>%
    select(Animal, Timepoint, Treatment, Sex, value = all_of(ba), censored, direction, cens_logical) %>%
    filter(!is.na(value), Timepoint == "11") %>%
    mutate(Animal = factor(Animal),
           Treatment = factor(Treatment, levels = c("Ctrl", "TAM")),
           Sex = factor(Sex, levels = c("male", "female")))
  
  n_obs <- nrow(df_ba)
  n_animals <- n_distinct(df_ba$Animal)
  n_ctrl <- sum(df_ba$Treatment == "Ctrl", na.rm = TRUE)
  n_tam <- sum(df_ba$Treatment == "TAM", na.rm = TRUE)
  n_male <- sum(df_ba$Sex == "male", na.rm = TRUE)
  n_female <- sum(df_ba$Sex == "female", na.rm = TRUE)
  n_censored <- sum(df_ba$cens_logical, na.rm = TRUE)
  n_below_LOD <- sum(df_ba$direction == "<", na.rm = TRUE)
  n_above_ULOQ <- sum(df_ba$direction == ">", na.rm = TRUE)
  
  result <- data.frame(
    BA = ba, n = n_obs, n_animals = n_animals, n_ctrl = n_ctrl, n_tam = n_tam, n_male = n_male, n_female = n_female,
    n_censored = n_censored, n_below_LOD = n_below_LOD, n_above_ULOQ = n_above_ULOQ,
    method = NA_character_,
    p_Treatment = NA_real_, p_Sex = NA_real_, p_Treatment_Sex = NA_real_,
    p_Treatment_female = NA_real_, treatment_effect_female = NA_real_, CI_low_Treatment_female = NA_real_, CI_high_Treatment_female = NA_real_,
    p_Treatment_male = NA_real_, treatment_effect_male = NA_real_, CI_low_Treatment_male = NA_real_, CI_high_Treatment_male = NA_real_,
    mean_Ctrl = NA_real_, sd_Ctrl = NA_real_, mean_TAM = NA_real_, sd_TAM = NA_real_,
    mean_male = NA_real_, sd_male = NA_real_, mean_female = NA_real_, sd_female = NA_real_,
    log2FC = NA_real_, fold_change = NA_real_,
    effect_size = NA_real_, effect_CI_low = NA_real_, effect_CI_high = NA_real_,
    effect_size_type = NA_character_, effect_size_calc = NA_character_,
    stringsAsFactors = FALSE
  )
  
  if (n_obs < 4 || n_animals < 2 || n_distinct(df_ba$Treatment) < 2 || n_distinct(df_ba$Sex) < 2) {
    result$method <- "insufficient_data"
    results[[ba]] <- result
    cat("Insufficient data\n")
    next
  }
  
  if (n_censored > 0) {
    if (n_above_ULOQ > 0) {
      result$method <- "censored_data_with_ULOQ_not_analyzed"
      results[[ba]] <- result
      cat("Contains >ULOQ values; cen2way not performed\n")
      next
    }
    
    uncensored <- df_ba %>% filter(!cens_logical)
    sufficient <- nrow(uncensored) >= 4 &&
      all(uncensored %>% group_by(Treatment) %>%
            summarise(n_dist = n_distinct(value), .groups = "drop") %>%
            pull(n_dist) >= 2)
    
    if (sufficient) {
      # Censored two-way model --------------------------------------------
      
      cen_result <- capture.output(
        suppressWarnings(
          with(
            df_ba,
            cen2way(
              ifelse(cens_logical, value * 2, value),
              cens_logical,
              Treatment,
              Sex,
              LOG = TRUE,
              interact = TRUE
            )
          )
        )
      )
      
      p_lines <- cen_result[grep("Treatment|Sex|interaction", cen_result)]
      print(p_lines)
      
      p_values <- as.numeric(
        sub(".*\\s([0-9]+\\.[0-9]+)$", "\\1", p_lines)
      )
      
      result$method <- "cen2way"
      result$p_Treatment <- p_values[3]
      result$p_Sex <- p_values[4]
      result$p_Treatment_Sex <- p_values[5]
      
      
      # ------------------------------------------------------------
      # Estimated geometric means for Treatment
      # marginal over Sex
      # ------------------------------------------------------------
      
      # Same input used by cen2way:
      y1 <- ifelse(
        df_ba$cens_logical,
        df_ba$value * 2,
        df_ba$value
      )
      
      y2 <- df_ba$cens_logical
      
      # Same factor ordering as used by the model
      fac1 <- factor(df_ba$Treatment)
      fac2 <- factor(df_ba$Sex)
      
      # Remove incomplete observations exactly as cen2way does
      keep <- complete.cases(y1, y2, fac1, fac2)
      
      y1 <- y1[keep]
      y2 <- y2[keep]
      fac1 <- fac1[keep]
      fac2 <- fac2[keep]
      
      # ------------------------------------------------------------
      # Reproduce cen2way's effect coding
      # ------------------------------------------------------------
      
      levels_1 <- levels(fac1)
      levels_2 <- levels(fac2)
      
      # Treatment coding:
      # first level = +1
      # last level  = -1
      e <- ifelse(fac1 == levels_1[1], 1, -1)
      
      # Sex coding:
      # first level = +1
      # last level  = -1
      d <- ifelse(fac2 == levels_2[1], 1, -1)
      
      # Treatment x Sex interaction
      int <- e * d
      
      # ------------------------------------------------------------
      # Same censored log-normal model as cen2way
      # ------------------------------------------------------------
      
      lnvar <- log(y1)
      
      fconst <- max(lnvar)
      
      flip.log <- fconst + 1 - lnvar
      
      detect <- !y2
      
      logCensData <- survival::Surv(
        flip.log,
        detect,
        type = "right"
      )
      
      cen_model <- survival::survreg(
        logCensData ~ e + d + int,
        dist = "gaussian"
      )
      
      # ------------------------------------------------------------
      # Convert coefficients back to the original log scale
      # exactly as cen2way does
      # ------------------------------------------------------------
      
      beta <- coef(cen_model)
      
      # cen2way reverses the flipping
      beta <- -beta
      # Correct intercept
      beta[1] <- fconst + 1 + beta[1]
      names(beta)
      # ------------------------------------------------------------
      # Marginal geometric means for Treatment
      # averaged over Sex
      # ------------------------------------------------------------
      b0   <- beta["(Intercept)"]
      b_tr <- beta["e"]
      b_sex <- beta["d"]
      b_int <- beta["int"]
      
      # Expected log means for the four Treatment x Sex combinations
      log_GM_Ctrl_Female <- b0 + b_tr + b_sex + b_int
      log_GM_Ctrl_Male   <- b0 + b_tr - b_sex - b_int
      log_GM_TAM_Female  <- b0 - b_tr + b_sex - b_int
      log_GM_TAM_Male    <- b0 - b_tr - b_sex + b_int
      
      # Marginal means on the log scale
      log_GM_Ctrl <- mean(c(log_GM_Ctrl_Female, log_GM_Ctrl_Male) )
      
      log_GM_TAM <- mean(c(log_GM_TAM_Female, log_GM_TAM_Male))
      
      # Back-transform -> geometric means
      GM_ctrl <- exp(log_GM_Ctrl)
      GM_tam  <- exp(log_GM_TAM)
      # Treatment fold change / geometric mean ratio
      GMR <- GM_tam / GM_ctrl
      # Store results
      result$mean_Ctrl <- GM_ctrl
      result$mean_TAM <- GM_tam
      result$effect_size <- GMR
      result$effect_size_type <- "GMR"
      result$effect_size_calc <-"censored log-normal two-way model; marginal Treatment GMR TAM/Ctrl"
    } else {
      result$method <- "cen2way_not_performed"
    }
    
  } else {
    fit <- tryCatch(lm(value ~ Treatment * Sex, data = df_ba), error = function(e) NULL)
    if (is.null(fit)) {
      result$method <- "model_failed"
      results[[ba]] <- result
      cat("Model failed\n")
      next
    }
    
    result$method <- "linear_model_Treatment_x_Sex"
    anova_result <- tryCatch(as.data.frame(anova(fit)), error = function(e) NULL)
    if (!is.null(anova_result)) {
      result$p_Treatment <- anova_result$`Pr(>F)`[rownames(anova_result) == "Treatment"]
      result$p_Sex <- anova_result$`Pr(>F)`[rownames(anova_result) == "Sex"]
      result$p_Treatment_Sex <- anova_result$`Pr(>F)`[rownames(anova_result) == "Treatment:Sex"]
    }
    
    result$mean_Ctrl <- mean(df_ba$value[df_ba$Treatment == "Ctrl"], na.rm = TRUE)
    result$sd_Ctrl <- sd(df_ba$value[df_ba$Treatment == "Ctrl"], na.rm = TRUE)
    result$mean_TAM <- mean(df_ba$value[df_ba$Treatment == "TAM"], na.rm = TRUE)
    result$sd_TAM <- sd(df_ba$value[df_ba$Treatment == "TAM"], na.rm = TRUE)
    result$mean_male <- mean(df_ba$value[df_ba$Sex == "male"], na.rm = TRUE)
    result$sd_male <- sd(df_ba$value[df_ba$Sex == "male"], na.rm = TRUE)
    result$mean_female <- mean(df_ba$value[df_ba$Sex == "female"], na.rm = TRUE)
    result$sd_female <- sd(df_ba$value[df_ba$Sex == "female"], na.rm = TRUE)
    
    if (is.finite(result$mean_Ctrl) && is.finite(result$mean_TAM) && result$mean_Ctrl > 0 && result$mean_TAM > 0) {
      result$fold_change <- result$mean_TAM / result$mean_Ctrl
      result$log2FC <- log2(result$fold_change)
    }
    
    emm_treatment <- tryCatch(emmeans(fit, ~ Treatment), error = function(e) NULL)
    con_treatment <- if (!is.null(emm_treatment)) tryCatch(as.data.frame(contrast(emm_treatment, method = "pairwise", adjust = "none", infer = TRUE)), error = function(e) NULL) else NULL
    
    if (!is.null(con_treatment) && nrow(con_treatment) > 0) {
      result$p_Treatment <- con_treatment$p.value[1]
    }
    
    if (!is.null(emm_treatment)) {
      eff <- tryCatch(eff_size(emm_treatment, sigma = sigma(fit), edf = df.residual(fit)), error = function(e) NULL)
      
      if (!is.null(eff)) {
        eff_summary <- as.data.frame(summary(eff))
        result$effect_size <- -eff_summary$effect.size[1]
        result$effect_CI_low <- -eff_summary$lower.CL[1]
        result$effect_CI_high <- -eff_summary$upper.CL[1]
        result$effect_size_type <- "standardized model effect"
        result$effect_size_calc <- "emmeans + eff_size; Cohen's d"
      }
    }
    
    emm_sex <- tryCatch(emmeans(fit, ~ Treatment | Sex), error = function(e) NULL)
    con_sex <- if (!is.null(emm_sex)) tryCatch(as.data.frame(contrast(emm_sex, method = "pairwise", adjust = "none", infer = TRUE)), error = function(e) NULL) else NULL
    
    if (!is.null(con_sex) && nrow(con_sex) > 0) {
      cond_female <- con_sex$Sex == "female"
      cond_male <- con_sex$Sex == "male"
      
      if (any(cond_female)) {
        result$p_Treatment_female <- con_sex$p.value[cond_female][1]
        result$treatment_effect_female <- con_sex$estimate[cond_female][1]
        result$CI_low_Treatment_female <- con_sex$lower.CL[cond_female][1]
        result$CI_high_Treatment_female <- con_sex$upper.CL[cond_female][1]
      }
      
      if (any(cond_male)) {
        result$p_Treatment_male <- con_sex$p.value[cond_male][1]
        result$treatment_effect_male <- con_sex$estimate[cond_male][1]
        result$CI_low_Treatment_male <- con_sex$lower.CL[cond_male][1]
        result$CI_high_Treatment_male <- con_sex$upper.CL[cond_male][1]
      }
    }
    
    cat("\n")
    print(ba)
    # print(summary(fit))
    # print(anova(fit))
    cat("\n")
  }
  
  results[[ba]] <- result
}

# ============================================================
# COMBINE ALL BILE ACIDS
# ============================================================

results_df <- bind_rows(results)

# ============================================================
# FDR ADJUSTMENT
# ============================================================

results_df <- results_df %>% mutate(
  adj_p_Treatment = p.adjust(p_Treatment, method = "fdr"),
  adj_p_Sex = p.adjust(p_Sex, method = "fdr"),
  adj_p_Treatment_Sex = p.adjust(p_Treatment_Sex, method = "fdr"),
  adj_p_Treatment_female = p.adjust(p_Treatment_female, method = "fdr"),
  adj_p_Treatment_male = p.adjust(p_Treatment_male, method = "fdr")
)

# ============================================================
# RELEVANT COMPARISONS
# ============================================================

results_df <- results_df %>% rowwise() %>% mutate(
  Relevant_comparisons = {
    statements <- c(
      if (!is.na(adj_p_Treatment) && adj_p_Treatment < 0.05) "Significant Treatment effect at TP11",
      if (!is.na(adj_p_Treatment) && adj_p_Treatment >= 0.05) "No significant Treatment effect at TP11",
      if (!is.na(adj_p_Treatment_female) && adj_p_Treatment_female < 0.05) "Significant Treatment effect at TP11 in females",
      if (!is.na(adj_p_Treatment_female) && adj_p_Treatment_female >= 0.05) "No significant Treatment effect at TP11 in females",
      if (!is.na(adj_p_Treatment_male) && adj_p_Treatment_male < 0.05) "Significant Treatment effect at TP11 in males",
      if (!is.na(adj_p_Treatment_male) && adj_p_Treatment_male >= 0.05) "No significant Treatment effect at TP11 in males",
      if (!is.na(adj_p_Treatment_Sex) && adj_p_Treatment_Sex < 0.05) "Treatment effect differs significantly between sexes",
      if (!is.na(adj_p_Treatment_Sex) && adj_p_Treatment_Sex >= 0.05) "Treatment effect does not differ significantly between sexes",
      if (!is.na(adj_p_Sex) && adj_p_Sex < 0.05) "Significant overall Sex effect",
      if (!is.na(adj_p_Sex) && adj_p_Sex >= 0.05) "No significant overall Sex effect"
    )
    statements <- statements[!is.na(statements) & statements != ""]
    if (length(statements) == 0) NA_character_ else paste(statements, collapse = "; ")
  }
) %>% ungroup()

# ============================================================
# WRITE CSV
# ============================================================

csv_path <- file.path(output_pwd, paste0(ExpID, "_BA_TP11_Statistics.csv"))
write.csv2(results_df, file = csv_path, row.names = FALSE)

# ============================================================
# SAVE R OBJECT
# ============================================================

rds_path <- file.path(output_pwd, paste0(ExpID, "_BA_TP11_Statistics.rds"))
saveRDS(results_df, file = rds_path)

cat("\nCSV written to:\n", csv_path, "\n")
cat("\nR object written to:\n", rds_path, "\n")