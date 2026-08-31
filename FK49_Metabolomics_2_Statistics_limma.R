###############################################################################
# FK49_Metabolomics_2_Statistics.R
###############################################################################

rm(list = ls())
gc()

library(tidyverse)
library(emmeans)
library(NADA2)
library(survival)
library(car)
library(limma)
source("FK49_Definitions.R")

analysis <- "CDHFD"

rawdata_pwd <- PATHS$metabolomics$rawdata
output_pwd  <- PATHS$metabolomics$output

analysis_folder <- if (analysis == "CDHFD") "CDHFD" else "ND"
analysis_pwd <- file.path(output_pwd, analysis_folder)
if (!dir.exists(analysis_pwd)) dir.create(analysis_pwd, recursive = TRUE)

if (analysis == "CDHFD") {
  diet_filter <- "CDHFD13"; expid_filter <- "FK49"; has_sex <- TRUE
} else if (analysis == "ND") {
  diet_filter <- "ND"; expid_filter <- "BH"; has_sex <- FALSE
} else stop('analysis must be "CDHFD" or "ND"')

pos <- readRDS(file.path(rawdata_pwd, "FK49_metabolome_positive_processed.rds"))
neg <- readRDS(file.path(rawdata_pwd, "FK49_metabolome_negative_processed.rds"))
tar <- readRDS(file.path(rawdata_pwd, "FK49_metabolome_targeted_processed.rds"))

datasets <- list(positive = pos, negative = neg, targeted = tar)

# ============================================================
# CEN2WAY
# ============================================================

run_cen2way <- function(value_raw, cens, Treatment, Sex) {
  
  cen_result <- capture.output(
    suppressWarnings(
      cen2way(
        ifelse(cens, value_raw * 2, value_raw),
        cens, Treatment, Sex,
        LOG = TRUE, interact = TRUE
      )
    )
  )
  
  p_lines <- cen_result[grep("Treatment|Sex|interaction", cen_result)]
  p_values <- as.numeric(sub(".*\\s([0-9]+\\.[0-9]+)$", "\\1", p_lines))
  
  p_Treatment     <- p_values[3]
  p_Sex           <- p_values[4]
  p_Treatment_Sex <- p_values[5]
  
  y1 <- ifelse(cens, value_raw * 2, value_raw)
  y2 <- cens
  fac1 <- factor(Treatment)
  fac2 <- factor(Sex)
  
  keep <- complete.cases(y1, y2, fac1, fac2)
  y1 <- y1[keep]; y2 <- y2[keep]
  fac1 <- fac1[keep]; fac2 <- fac2[keep]
  
  levels_1 <- levels(fac1)
  levels_2 <- levels(fac2)
  
  e <- ifelse(fac1 == levels_1[1], 1, -1)
  d <- ifelse(fac2 == levels_2[1], 1, -1)
  int <- e * d
  
  lnvar <- log(y1)
  fconst <- max(lnvar)
  flip.log <- fconst + 1 - lnvar
  detect <- !y2
  
  logCensData <- survival::Surv(flip.log, detect, type = "right")
  cen_model <- survival::survreg(
    logCensData ~ e + d + int,
    dist = "gaussian"
  )
  
  beta <- -coef(cen_model)
  beta[1] <- fconst + 1 + beta[1]
  
  b0 <- beta["(Intercept)"]
  b_tr <- beta["e"]
  b_sex <- beta["d"]
  b_int <- beta["int"]
  
  log_GM_Ctrl_Female <- b0 + b_tr + b_sex + b_int
  log_GM_Ctrl_Male   <- b0 + b_tr - b_sex - b_int
  log_GM_TAM_Female  <- b0 - b_tr + b_sex - b_int
  log_GM_TAM_Male    <- b0 - b_tr - b_sex + b_int
  
  log_GM_Ctrl <- mean(c(log_GM_Ctrl_Female, log_GM_Ctrl_Male))
  log_GM_TAM  <- mean(c(log_GM_TAM_Female, log_GM_TAM_Male))
  
  GM_ctrl <- exp(log_GM_Ctrl)
  GM_tam  <- exp(log_GM_TAM)
  
  log_GMR <- -2 * b_tr
  GMR <- exp(log_GMR)
  
  se_log_GMR <- 2 * sqrt(vcov(cen_model)["e", "e"])
  
  list(
    p_Treatment = p_Treatment,
    p_Sex = p_Sex,
    p_Treatment_Sex = p_Treatment_Sex,
    mean_Ctrl = GM_ctrl,
    mean_TAM = GM_tam,
    effect_size = GMR,
    effect_CI_low = exp(log_GMR - 1.96 * se_log_GMR),
    effect_CI_high = exp(log_GMR + 1.96 * se_log_GMR),
    effect_size_type = "GMR",
    fold_change = GMR,
    log2FC = log2(GMR)
  )
}

# ============================================================
# CENS1WAY
# ============================================================

run_cens1way <- function(value_raw, cens, Treatment) {
  
  cen_vals <- ifelse(cens, value_raw * 2, value_raw)
  
  cen_result <- capture.output(
    suppressWarnings(
      cens1way(cen_vals, cens, Treatment, LOG = TRUE)
    )
  )
  
  p_lines <- cen_result[
    grep("p-value|Treatment", cen_result, ignore.case = TRUE)
  ]
  
  p_values <- as.numeric(
    sub(".*\\s([0-9]+\\.[0-9]+)$", "\\1", p_lines)
  )
  
  p_Treatment <- if (
    length(p_values) > 0 &&
    !is.na(p_values[length(p_values)])
  ) p_values[length(p_values)] else NA_real_
  
  y1 <- ifelse(cens, value_raw * 2, value_raw)
  y2 <- cens
  fac1 <- factor(Treatment)
  
  keep <- complete.cases(y1, y2, fac1)
  
  y1 <- y1[keep]; y2 <- y2[keep]; fac1 <- fac1[keep]
  
  levels_1 <- levels(fac1)
  e <- ifelse(fac1 == levels_1[1], 1, -1)
  
  lnvar <- log(y1)
  fconst <- max(lnvar)
  flip.log <- fconst + 1 - lnvar
  detect <- !y2
  
  logCensData <- survival::Surv(flip.log, detect, type = "right")
  
  cen_model <- survival::survreg(
    logCensData ~ e,
    dist = "gaussian"
  )
  
  beta <- -coef(cen_model)
  beta[1] <- fconst + 1 + beta[1]
  
  b0 <- beta["(Intercept)"]
  b_tr <- beta["e"]
  
  log_GM_Ctrl <- b0 + b_tr
  log_GM_TAM <- b0 - b_tr
  
  GMR <- exp(-2 * b_tr)
  se_log_GMR <- 2 * sqrt(vcov(cen_model)["e", "e"])
  
  if (is.na(p_Treatment)) {
    cen_null <- survival::survreg(
      logCensData ~ 1,
      dist = "gaussian"
    )
    
    lr_stat <- as.numeric(
      2 * (logLik(cen_model) - logLik(cen_null))
    )
    
    p_Treatment <- pchisq(
      lr_stat, df = 1, lower.tail = FALSE
    )
  }
  
  list(
    p_Treatment = p_Treatment,
    mean_Ctrl = exp(log_GM_Ctrl),
    mean_TAM = exp(log_GM_TAM),
    effect_size = GMR,
    effect_CI_low = exp(-2 * b_tr - 1.96 * se_log_GMR),
    effect_CI_high = exp(-2 * b_tr + 1.96 * se_log_GMR),
    effect_size_type = "GMR",
    fold_change = GMR,
    log2FC = log2(GMR)
  )
}

# ============================================================
# MAIN LOOP
# ============================================================

all_results <- list()

for (ds_name in names(datasets)) {
  
  ds <- datasets[[ds_name]]
  
  cat("\n\n========== Dataset:", ds_name, "==========\n")
  
  idx <- ds$metadata$Diet == diet_filter &
    ds$metadata$ExpID == expid_filter
  
  metadata   <- ds$metadata[idx, , drop = FALSE]
  raw_values <- ds$raw_values[idx, , drop = FALSE]
  log_values <- ds$log_values[idx, , drop = FALSE]
  censored   <- ds$censored[idx, , drop = FALSE]
  
  metab_names <- colnames(raw_values)
  
  Treatment <- factor(
    metadata$Treatment,
    levels = c("Ctrl", "TAM")
  )
  
  Sex <- factor(
    metadata$Sex,
    levels = c("female", "male")
  )
  
  Batch <- factor(metadata$Batch)
  use_batch <- n_distinct(Batch) > 1
  use_sex <- has_sex && n_distinct(Sex) >= 2
  
  # ----------------------------------------------------------
  # LIMMA: all uncensored metabolites simultaneously
  # ----------------------------------------------------------
  
  uncensored_metabs <- metab_names[
    colSums(censored[, metab_names, drop = FALSE], na.rm = TRUE) == 0
  ]
  
  if (length(uncensored_metabs) > 0) {
    
    if (use_sex && use_batch) {
      design <- model.matrix(
        ~ Treatment * Sex + Batch,
        data = metadata
      )
    } else if (use_sex) {
      design <- model.matrix(
        ~ Treatment * Sex,
        data = metadata
      )
    } else if (use_batch) {
      design <- model.matrix(
        ~ Treatment + Batch,
        data = metadata
      )
    } else {
      design <- model.matrix(
        ~ Treatment,
        data = metadata
      )
    }
    
    protein_matrix <- t(
      as.matrix(log_values[, uncensored_metabs, drop = FALSE])
    )
    
    fit <- lmFit(protein_matrix, design)
    fit <- eBayes(fit)
    
    coef_names <- colnames(design)
    
    treatment_coef <- match("TreatmentTAM", coef_names)
    sex_coef <- match("Sexmale", coef_names)
    interaction_coef <- match("TreatmentTAM:Sexmale", coef_names)
    batch_coef <- grep("^Batch", coef_names)[1]
    
    if (use_sex) {
      
      # Overall Treatment effect = average effect across sexes
      C_overall <- rep(0, ncol(design))
      C_overall[treatment_coef] <- 1
      C_overall[interaction_coef] <- 0.5
      
      fit_overall <- eBayes(
        contrasts.fit(
          fit,
          contrasts = matrix(
            C_overall,
            ncol = 1,
            dimnames = list(coef_names, "Treatment_overall")
          )
        )
      )
      
      # Difference in Treatment effect between sexes
      C_interaction <- rep(0, ncol(design))
      C_interaction[interaction_coef] <- 1
      
      fit_interaction <- eBayes(
        contrasts.fit(
          fit,
          contrasts = matrix(
            C_interaction,
            ncol = 1,
            dimnames = list(coef_names, "Treatment_x_Sex")
          )
        )
      )
      
      # Treatment effect in females
      C_female <- rep(0, ncol(design))
      C_female[treatment_coef] <- 1
      
      fit_female <- eBayes(
        contrasts.fit(
          fit,
          contrasts = matrix(
            C_female,
            ncol = 1,
            dimnames = list(coef_names, "Treatment_female")
          )
        )
      )
      
      # Treatment effect in males
      C_male <- rep(0, ncol(design))
      C_male[treatment_coef] <- 1
      C_male[interaction_coef] <- 1
      
      fit_male <- eBayes(
        contrasts.fit(
          fit,
          contrasts = matrix(
            C_male,
            ncol = 1,
            dimnames = list(coef_names, "Treatment_male")
          )
        )
      )
      
    } else {
      
      C_overall <- rep(0, ncol(design))
      C_overall[treatment_coef] <- 1
      
      fit_overall <- eBayes(
        contrasts.fit(
          fit,
          contrasts = matrix(
            C_overall,
            ncol = 1,
            dimnames = list(coef_names, "Treatment_overall")
          )
        )
      )
    }
  }
  
  # ----------------------------------------------------------
  # RESULT LOOP
  # ----------------------------------------------------------
  
  results <- list()
  
  for (metab in metab_names) {
    
    cat(".")
    
    value_raw <- raw_values[, metab]
    value_log <- log_values[, metab]
    cens <- censored[, metab]
    
    n_obs <- length(value_raw)
    n_animals <- n_distinct(metadata$Animal)
    
    n_ctrl <- sum(Treatment == "Ctrl", na.rm = TRUE)
    n_tam <- sum(Treatment == "TAM", na.rm = TRUE)
    n_female <- sum(Sex == "female", na.rm = TRUE)
    n_male <- sum(Sex == "male", na.rm = TRUE)
    n_censored <- sum(cens, na.rm = TRUE)
    
    result <- data.frame(
      Dataset = ds_name,
      Metabolite = metab,
      n = n_obs,
      n_Ctrl = n_ctrl,
      n_TAM = n_tam,
      n_female = n_female,
      n_male = n_male,
      n_censored = n_censored,
      method = NA_character_,
      p_Treatment = NA_real_,
      p_Sex = NA_real_,
      p_Treatment_Sex = NA_real_,
      p_Batch = NA_real_,
      p_Treatment_female = NA_real_,
      treatment_effect_female = NA_real_,
      CI_low_female = NA_real_,
      CI_high_female = NA_real_,
      p_Treatment_male = NA_real_,
      treatment_effect_male = NA_real_,
      CI_low_male = NA_real_,
      CI_high_male = NA_real_,
      mean_Ctrl = NA_real_,
      sd_Ctrl = NA_real_,
      mean_TAM = NA_real_,
      sd_TAM = NA_real_,
      log2FC = NA_real_,
      fold_change = NA_real_,
      effect_size = NA_real_,
      effect_CI_low = NA_real_,
      effect_CI_high = NA_real_,
      effect_size_type = NA_character_,
      stringsAsFactors = FALSE
    )
    
    if (
      n_obs < 4 ||
      n_animals < 2 ||
      n_distinct(Treatment) < 2
    ) {
      result$method <- "insufficient_data"
      results[[metab]] <- result
      next
    }
    
    # --------------------------------------------------------
    # CENSORED
    # --------------------------------------------------------
    
    if (n_censored > 0) {
      
      uncensored_vals <- value_raw[!cens]
      
      if (length(uncensored_vals) >= 4) {
        
        if (use_sex) {
          
          cen_out <- tryCatch(
            run_cen2way(
              value_raw, cens,
              Treatment, Sex
            ),
            error = function(e) NULL
          )
          
          if (!is.null(cen_out)) {
            
            result$method <- "cen2way"
            result$p_Treatment <- cen_out$p_Treatment
            result$p_Sex <- cen_out$p_Sex
            result$p_Treatment_Sex <- cen_out$p_Treatment_Sex
            result$mean_Ctrl <- cen_out$mean_Ctrl
            result$mean_TAM <- cen_out$mean_TAM
            result$effect_size <- cen_out$effect_size
            result$effect_CI_low <- cen_out$effect_CI_low
            result$effect_CI_high <- cen_out$effect_CI_high
            result$effect_size_type <- cen_out$effect_size_type
            result$fold_change <- cen_out$fold_change
            result$log2FC <- cen_out$log2FC
            
          } else result$method <- "cen2way_failed"
          
        } else {
          
          cen_out <- tryCatch(
            run_cens1way(
              value_raw, cens,
              Treatment
            ),
            error = function(e) NULL
          )
          
          if (!is.null(cen_out)) {
            
            result$method <- "cens1way"
            result$p_Treatment <- cen_out$p_Treatment
            result$mean_Ctrl <- cen_out$mean_Ctrl
            result$mean_TAM <- cen_out$mean_TAM
            result$effect_size <- cen_out$effect_size
            result$effect_CI_low <- cen_out$effect_CI_low
            result$effect_CI_high <- cen_out$effect_CI_high
            result$effect_size_type <- cen_out$effect_size_type
            result$fold_change <- cen_out$fold_change
            result$log2FC <- cen_out$log2FC
            
          } else result$method <- "cens1way_failed"
        }
        
      } else {
        result$method <- "insufficient_uncensored"
      }
      
      # --------------------------------------------------------
      # UNCENSORED = LIMMA
      # --------------------------------------------------------
      
    } else {
      
      result$method <- if (use_sex)
        "limma_Treatment_x_Sex"
      else
        "limma_Treatment"
      
      # Overall Treatment effect
      result$p_Treatment <- fit_overall$p.value[
        metab, "Treatment_overall"
      ]
      
      if (use_sex) {
        
        # Difference in Treatment effect between sexes
        result$p_Treatment_Sex <- fit_interaction$p.value[
          metab, "Treatment_x_Sex"
        ]
        
        # Treatment effect in females
        result$p_Treatment_female <- fit_female$p.value[
          metab, "Treatment_female"
        ]
        
        result$treatment_effect_female <- fit_female$coefficients[
          metab, "Treatment_female"
        ]
        
        result$CI_low_female <-
          result$treatment_effect_female -
          1.96 *
          fit_female$stdev.unscaled[
            metab, "Treatment_female"
          ] *
          fit_female$sigma[metab]
        
        result$CI_high_female <-
          result$treatment_effect_female +
          1.96 *
          fit_female$stdev.unscaled[
            metab, "Treatment_female"
          ] *
          fit_female$sigma[metab]
        
        # Treatment effect in males
        result$p_Treatment_male <- fit_male$p.value[
          metab, "Treatment_male"
        ]
        
        result$treatment_effect_male <- fit_male$coefficients[
          metab, "Treatment_male"
        ]
        
        result$CI_low_male <-
          result$treatment_effect_male -
          1.96 *
          fit_male$stdev.unscaled[
            metab, "Treatment_male"
          ] *
          fit_male$sigma[metab]
        
        result$CI_high_male <-
          result$treatment_effect_male +
          1.96 *
          fit_male$stdev.unscaled[
            metab, "Treatment_male"
          ] *
          fit_male$sigma[metab]
      }
      
      if (use_sex)
        result$p_Sex <- fit$p.value[
          metab, sex_coef
        ]
      
      if (use_batch)
        result$p_Batch <- fit$p.value[
          metab, batch_coef
        ]
      
      result$mean_Ctrl <- mean(
        value_log[Treatment == "Ctrl"],
        na.rm = TRUE
      )
      
      result$sd_Ctrl <- sd(
        value_log[Treatment == "Ctrl"],
        na.rm = TRUE
      )
      
      result$mean_TAM <- mean(
        value_log[Treatment == "TAM"],
        na.rm = TRUE
      )
      
      result$sd_TAM <- sd(
        value_log[Treatment == "TAM"],
        na.rm = TRUE
      )
      
      result$log2FC <-
        result$mean_TAM -
        result$mean_Ctrl
      
      result$fold_change <-
        2^result$log2FC
      
      result$effect_size <-
        (result$mean_TAM - result$mean_Ctrl) /
        sqrt(
          (
            (n_tam - 1) * result$sd_TAM^2 +
              (n_ctrl - 1) * result$sd_Ctrl^2
          ) /
            (n_tam + n_ctrl - 2)
        )
      
      result$effect_size_type <- "Cohen's d"
    }
    
    results[[metab]] <- result
  }
  
  cat("\n")
  
  results_df <- bind_rows(results) %>%
    mutate(
      adj_p_Treatment =
        p.adjust(p_Treatment, method = "fdr"),
      
      adj_p_Sex =
        p.adjust(p_Sex, method = "fdr"),
      
      adj_p_Treatment_Sex =
        p.adjust(p_Treatment_Sex, method = "fdr"),
      
      adj_p_Batch =
        p.adjust(p_Batch, method = "fdr"),
      
      adj_p_Treatment_female =
        p.adjust(p_Treatment_female, method = "fdr"),
      
      adj_p_Treatment_male =
        p.adjust(p_Treatment_male, method = "fdr")
    )
  
  all_results[[ds_name]] <- results_df
}

# ============================================================
# COMBINE
# ============================================================

final_results <- bind_rows(all_results) %>%
  mutate(
    significant =
      !is.na(adj_p_Treatment) &
      adj_p_Treatment < 0.05 &
      abs(log2FC) > 0.5,
    
    trend =
      !is.na(adj_p_Treatment) &
      adj_p_Treatment < 0.1 &
      abs(log2FC) > 0.3,
    
    direction = case_when(
      significant & log2FC > 0 ~ "UP",
      significant & log2FC < 0 ~ "DOWN",
      TRUE ~ "NS"
    )
  )

# ============================================================
# SAVE
# ============================================================

csv_path <- file.path(
  analysis_pwd,
  "FK49_metabolome_statistics.csv"
)

rds_path <- file.path(
  analysis_pwd,
  "FK49_metabolome_statistics.rds"
)

write.csv2(
  final_results,
  file = csv_path,
  row.names = FALSE
)

saveRDS(
  final_results,
  file = rds_path
)

cat(
  "\n=== Statistics complete (",
  analysis,
  ") ===\n"
)

cat(
  "CSV written to:",
  csv_path,
  "\n"
)

cat(
  "RDS written to:",
  rds_path,
  "\n"
)

for (ds_name in names(all_results)) {
  
  df <- all_results[[ds_name]]
  
  cat(
    "  ", ds_name, ":",
    sum(
      df$method == "limma_Treatment_x_Sex" |
        df$method == "limma_Treatment"
    ),
    "limma,",
    sum(
      df$method == "cen2way" |
        df$method == "cens1way"
    ),
    "censored,",
    sum(
      df$method == "insufficient_data" |
        df$method == "insufficient_uncensored" |
        df$method == "cen2way_failed" |
        df$method == "cens1way_failed"
    ),
    "skipped\n"
  )
}