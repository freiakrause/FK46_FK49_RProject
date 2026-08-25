###############################################################################
# FK49_Metabolomics_1_Preprocess.R
#
# Reads metabolomics XLSX files (untargeted: Positive/Negative sheets,
# targeted: Concentrations + SampleData sheets), applies LOD/censored
# handling (BA-style pivot-long pivot-wide pattern), normalises,
# log2-transforms, and saves one RDS per dataset.
#
# RDS structure (list per dataset):
#   metadata    - data.frame: Animal, Sample, Sex, Treatment, Diet, ExpID,
#                 Batch, T_D_S, T_D, T_S
#   raw_values  - numeric matrix: LOD-imputed values (for cen2way / dot plots)
#   log_values  - numeric matrix: log2-transformed (for LM / PCA)
#   norm_values - numeric matrix or NULL: row-sum normalised (untargeted only)
#   censored    - logical matrix: TRUE where original value was <LOD
#   metabolite_abbrev - named vector: RefMet abbreviated names (names = original)
###############################################################################

rm(list = ls())
gc()

library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)
library(readxl)
source("FK49_Definitions.R")

# ============================================================
# CONFIGURATION
# ============================================================

meta_cols <- c("Animal", "Sample", "Sex", "Treatment", "Diet", "ExpID", "Batch",
               "T_D_S", "T_D", "T_S")

rawdata_pwd <- PATHS$metabolomics$rawdata
output_pwd  <- PATHS$metabolomics$output

# ============================================================
# READ AND PROCESS XLSX FILES
# ============================================================

# --- SampleData from targeted XLSX (single source of truth for metadata + Batch) ---
sample_data  <- read_excel(PATHS$metabolomics$targeted_xlsx, sheet = "SampleData",
                           .name_repair = "minimal")
batch_lookup <- sample_data %>% select(Animal, Batch)

# --- Positive (untargeted, from untargeted XLSX) ---
metabolome_positive <- read_excel(PATHS$metabolomics$untargeted_xlsx, sheet = "Positive",
                                  .name_repair = "minimal") %>%
  left_join(batch_lookup, by = "Animal") %>%
  mutate(Sample = SampleNo) %>%
  select(-SampleCode, -SampleNo) %>%
  process_metabolome(meta_cols)

# --- Negative (untargeted, from untargeted XLSX) ---
metabolome_negative <- read_excel(PATHS$metabolomics$untargeted_xlsx, sheet = "Negative",
                                  .name_repair = "minimal") %>%
  left_join(batch_lookup, by = "Animal") %>%
  mutate(Sample = SampleNo) %>%
  select(-`Sample Code`, -SampleNo) %>%
  process_metabolome(meta_cols)

# --- Targeted (HILIC09, from targeted XLSX: join Concentrations + SampleData) ---
metabolome_targeted <- read_excel(PATHS$metabolomics$targeted_xlsx,
                                  sheet = "Concentrations, pM-mg",
                                  .name_repair = "minimal") %>%
  left_join(sample_data, by = "sampleID") %>%
  select(-sampleID, -`Sample Idxx`) %>%
  rename(Sample = SampleNo) %>%
  process_metabolome(meta_cols)

# ============================================================
# HELPER: Extract matrices, normalise, build RDS
# ============================================================

prepare_dataset <- function(processed, meta_cols, method, dataset_name) {
  cat("\n--- Preparing", dataset_name, "(", method, ") ---\n")

  # Identify metabolite names (columns without _censored/_direction/_raw suffix)
  all_cols <- colnames(processed)
  metab_names <- all_cols[!all_cols %in% meta_cols &
                          !grepl("_censored$|_direction$|_raw$", all_cols)]
  cat("Metabolites found:", length(metab_names), "\n")

  # Extract metadata
  metadata <- processed %>% select(all_of(meta_cols)) %>% as.data.frame()

  # Extract raw values (LOD-imputed numeric, pre-normalisation)
  raw_values <- processed %>% select(all_of(metab_names)) %>% as.matrix()
  rownames(raw_values) <- metadata$Animal

  # Extract censored matrix
  cens_cols <- paste0(metab_names, "_censored")
  censored <- processed %>% select(all_of(cens_cols)) %>% as.matrix()
  colnames(censored) <- metab_names
  rownames(censored) <- metadata$Animal

  # Remove zero-variance metabolites (NA variance → treat as zero variance)
  nzv <- apply(raw_values, 2, var, na.rm = TRUE) > 0
  nzv[is.na(nzv)] <- FALSE
  cat("Zero-variance metabolites removed:", sum(!nzv), "\n")
  raw_values <- raw_values[, nzv, drop = FALSE]
  censored   <- censored[, nzv, drop = FALSE]
  metab_names <- metab_names[nzv]

  # Normalisation
  if (method == "untargeted") {
    # Row-sum normalisation then log2
    rs <- rowSums(raw_values, na.rm = TRUE)
    rs[rs == 0] <- NA
    norm_values <- raw_values / rs
    log_values  <- log2(norm_values + 1e-9)
  } else if (method == "targeted") {
    # Log2 with epsilon (no row-sum normalisation)
    eps <- min(raw_values[raw_values > 0], na.rm = TRUE) / 2
    log_values  <- log2(raw_values + eps)
    norm_values <- NULL
  }

  # QC histogram
  raw_long <- as.data.frame(raw_values) %>%
    pivot_longer(everything(), names_to = "Metabolite", values_to = "Value") %>%
    mutate(type = "raw")

  log_long <- as.data.frame(log_values) %>%
    pivot_longer(everything(), names_to = "Metabolite", values_to = "Value") %>%
    mutate(type = "log2-transformed")

  if (!is.null(norm_values)) {
    norm_long <- as.data.frame(norm_values) %>%
      pivot_longer(everything(), names_to = "Metabolite", values_to = "Value") %>%
      mutate(type = "normalized")
    combined <- bind_rows(raw_long, norm_long, log_long) %>%
      mutate(type = factor(type, levels = c("raw", "normalized", "log2-transformed")))
  } else {
    combined <- bind_rows(raw_long, log_long) %>%
      mutate(type = factor(type, levels = c("raw", "log2-transformed")))
  }

  plot_qc <- ggplot(combined, aes(x = Value, fill = type)) +
    geom_histogram(bins = 60, alpha = 0.6, position = "identity") +
    facet_wrap(~ type, scales = "free_x", ncol = 4) +
    scale_fill_manual(values = c("gray", "blue", "violet")) +
    labs(title = paste0("Distribution - ", dataset_name),
         x = "Value", y = "Frequency") +
    theme_minimal() +
    theme(legend.position = "none")

  # Build abbrev lookup for this dataset's metabolites (from RefMet mapping)
  metab_abbrev <- PARAMETERS$metabolomics$metabolite_abbrev[colnames(raw_values)]
  metab_abbrev[is.na(metab_abbrev)] <- colnames(raw_values)[is.na(metab_abbrev)]

  # Build RDS list
  rds_data <- list(
    metadata    = metadata,
    raw_values  = raw_values,
    log_values  = log_values,
    norm_values = norm_values,
    censored    = censored,
    metabolite_abbrev = metab_abbrev
  )

  list(rds = rds_data, qc_plot = plot_qc)
}

# ============================================================
# PROCESS ALL THREE DATASETS
# ============================================================

pos <- prepare_dataset(metabolome_positive, meta_cols, "untargeted", "positive")
neg <- prepare_dataset(metabolome_negative, meta_cols, "untargeted", "negative")
tar <- prepare_dataset(metabolome_targeted, meta_cols, "targeted",  "targeted")

# ============================================================
# SAVE RDS AND QC PLOTS
# ============================================================

saveRDS(pos$rds, file = file.path(rawdata_pwd, "FK49_metabolome_positive_processed.rds"))
saveRDS(neg$rds, file = file.path(rawdata_pwd, "FK49_metabolome_negative_processed.rds"))
saveRDS(tar$rds, file = file.path(rawdata_pwd, "FK49_metabolome_targeted_processed.rds"))

create_output_folders(output_pwd, c("positive", "negative", "targeted"))

ggsave(plot = pos$qc_plot, filename = "FK49_metabolome_positive_QC.png",
       path = file.path(output_pwd, "positive"), width = 12, height = 6, dpi = 300)
ggsave(plot = neg$qc_plot, filename = "FK49_metabolome_negative_QC.png",
       path = file.path(output_pwd, "negative"), width = 12, height = 6, dpi = 300)
ggsave(plot = tar$qc_plot, filename = "FK49_metabolome_targeted_QC.png",
       path = file.path(output_pwd, "targeted"), width = 9, height = 6, dpi = 300)

cat("\n=== Preprocessing complete ===\n")
cat("RDS files saved to:", rawdata_pwd, "\n")
cat("QC plots saved to:", output_pwd, "\n")
cat("  positive:", ncol(pos$rds$raw_values), "metabolites,",
    nrow(pos$rds$raw_values), "samples\n")
cat("  negative:", ncol(neg$rds$raw_values), "metabolites,",
    nrow(neg$rds$raw_values), "samples\n")
cat("  targeted:", ncol(tar$rds$raw_values), "metabolites,",
    nrow(tar$rds$raw_values), "samples\n")
