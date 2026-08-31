###############################################################################
# FK49_Microbiome_1_Preprocessing.R
#
# Reads raw taxa_SV CSV files (OTU/ASV tables from collaboration partner),
# builds a phyloseq object, performs QC, filtering, and saves intermediate
# RDS files for downstream statistics and plotting.
#
# Key decisions:
#   - EtOH renamed to Ctrl (consistent with all other FK49 scripts)
#   - F5 excluded entirely (systematic sequencing failure: 163_F5=149 reads,
#     164_F5=1 read, 160_F5=1033 reads; male TAM at F5 = n=1 after QC)
#   - Prevalence filter: taxa present in >= 5% of samples
#   - Min read filter: samples with < 2000 total reads removed
#
# Inputs:  PATHS$microbiome$input_16S  (folder with *_taxa_SV.csv)
#          PATHS$microbiome$input_meta (metadata CSV)
# Outputs: ps.rds, ps_rel.rds, ps_genus.rds, ps_clr.rds, metadata.rds,
#          Seq_Info.csv, QC plots
###############################################################################

rm(list = ls())
gc()

library(tidyverse)
library(phyloseq)
library(tools)

source("FK49_Definitions.R")

# ============================================================
# CONFIGURATION
# ============================================================
mb_params <- PARAMETERS$microbiome
input_16S  <- PATHS$microbiome$input_16S
input_meta <- PATHS$microbiome$input_meta
output_dir <- PATHS$microbiome$output

# Create output folders before saving
create_output_folders(output_dir, c("Statistics", "Plots", "Plots/Exploratory"))

cat("=== FK49 Microbiome Preprocessing ===\n")

# ============================================================
# 1. READ AND COMBINE TAXA CSV FILES
# ============================================================
files <- list.files(input_16S, pattern = "_taxa_SV\\.csv$", full.names = TRUE)
cat("Found", length(files), "taxa CSV files\n")

df_list <- lapply(files, function(f) {
  read.csv(f, sep = ";", header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
})

# Add SampleID from filename (extract Animal_Feces pattern, e.g. "157_F1")
df_list <- Map(function(df, fname) {
  df$SampleID <- file_path_sans_ext(basename(fname))
  df
}, df_list, files)

all_df <- do.call(rbind, df_list)
colnames(all_df)[1] <- "Sequence"

# Parse SampleID: from "20250310-iAL157-F1-FK491" to "157_F1"
all_df <- all_df %>%
  mutate(SampleID = gsub(".*iAL([0-9]+)-([A-Z0-9]+)-.*", "\\1_\\2", SampleID))

cat("Combined data:", nrow(all_df), "rows,", n_distinct(all_df$SampleID), "samples\n")

# ============================================================
# 2. BUILD COUNT MATRIX
# ============================================================
count_mat <- all_df %>%
  dplyr::select(Sequence, SampleID, abundance) %>%
  pivot_wider(names_from = SampleID, values_from = abundance, values_fill = 0) %>%
  column_to_rownames(var = "Sequence") %>%
  as.matrix()

# Save sequence mapping (assign short OTU IDs)
sequences <- rownames(count_mat)
sequence_ids <- paste0("OTU", seq_along(sequences))
seq_info <- data.frame(seq_ids = sequence_ids, Sequence = sequences,
                       stringsAsFactors = FALSE)
write.csv(seq_info, file.path(output_dir, "Seq_Info.csv"), row.names = FALSE)

# Convert to numeric and reassign OTU IDs as rownames
count_mat <- apply(count_mat, 2, as.numeric)
rownames(count_mat) <- sequence_ids

# ============================================================
# 3. BUILD TAXONOMY MATRIX
# ============================================================
tax_mat_unique <- all_df %>%
  dplyr::select(Sequence, Kingdom:Species) %>%
  distinct(Sequence, .keep_all = TRUE) %>%
  column_to_rownames(var = "Sequence") %>%
  as.matrix()

# Align taxonomy rows to count matrix sequence order
tax_mat_unique <- tax_mat_unique[sequences, ]
rownames(tax_mat_unique) <- sequence_ids

# ============================================================
# 4. READ METADATA
# ============================================================
metadata <- read.csv(input_meta, sep = ";", stringsAsFactors = FALSE) %>%
  mutate(
    # Rename EtOH -> Ctrl (consistent with all FK49 scripts)
    Treatment = ifelse(Treatment == "EtOH", "Ctrl", Treatment),
    # Create derived grouping columns
    T_T_D_S = paste0(Treatment, "_", Feces, "_", Diet, "_", Sex),
    T_T_D   = paste0(Treatment, "_", Feces, "_", Diet),
    T_T     = paste0(Treatment, "_", Feces),
    T_S     = paste0(Treatment, "_", Sex)
  ) %>%
  mutate(
    Animal = factor(Animal),
    Treatment = factor(Treatment, levels = c("Ctrl", "TAM")),
    Sex       = factor(Sex,       levels = c("female", "male")),
    Feces     = factor(Feces,     levels = mb_params$timepoints),
    Diet      = factor(Diet,      levels = c("ND", "ND_TAM", "CDHFD_3", "CDHFD_7", "CDHFD_11")),
    Diet_short = factor(Diet_short, levels = c("ND", "CDHFD"))
  )
rownames(metadata) <- metadata$SampleID

cat("Metadata:", nrow(metadata), "samples\n")

# ============================================================
# 5. CONSTRUCT PHYLOSEQ OBJECT
# ============================================================
ps <- phyloseq(
  otu_table(count_mat, taxa_are_rows = TRUE),
  tax_table(tax_mat_unique),
  sample_data(metadata)
)

# Remove taxa with zero total counts
ps <- prune_taxa(taxa_sums(ps) > 0, ps)
cat("Initial phyloseq:", nsamples(ps), "samples,", ntaxa(ps), "taxa\n")

# ============================================================
# 6. EXCLUDE F5 SAMPLES
# ============================================================
ps <- prune_samples(sample_data(ps)$Feces %in% mb_params$timepoints, ps)
cat("After F5 exclusion:", nsamples(ps), "samples\n")

# ============================================================
# 7. QC: LIBRARY SIZE PLOT
# ============================================================
lib_sizes <- data.frame(
  SampleID = sample_names(ps),
  Reads    = sample_sums(ps),
  Treatment = as.character(sample_data(ps)$Treatment),
  Feces     = as.character(sample_data(ps)$Feces)
)

p_lib <- ggplot(lib_sizes, aes(x = SampleID, y = Reads, fill = Treatment)) +
  geom_col() +
  facet_grid(~ Feces, scales = "free_x", space = "free_x") +
  geom_hline(yintercept = mb_params$min_reads, linetype = "dashed", color = "red") +
  scale_fill_manual(values = Treatment_colors) +
  labs(title = "Library size per sample", x = NULL, y = "Total reads") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 7),
        legend.position = "bottom")
ggsave(file.path(output_dir, "QC_LibrarySize.png"), p_lib, width = 10, height = 5, dpi = 300)

# ============================================================
# 8. QC: PREVALENCE-ABUNDANCE SCATTER
# ============================================================
prevdf <- data.frame(
  Prevalence      = apply(otu_table(ps), 1, function(x) sum(x > 0)),
  TotalAbundance  = taxa_sums(ps),
  tax_table(ps)
)

p_prev <- ggplot(prevdf, aes(TotalAbundance, Prevalence / nsamples(ps), color = Phylum)) +
  geom_hline(yintercept = mb_params$prevalence_threshold, alpha = 0.5, linetype = 2) +
  geom_point(size = 2, alpha = 0.7) +
  scale_x_log10() +
  labs(x = "Total Abundance", y = "Prevalence [Frac. Samples]") +
  facet_wrap(~ Phylum) +
  theme(legend.position = "none")
ggsave(file.path(output_dir, "QC_PrevalenceAbundance.png"), p_prev, width = 10, height = 8, dpi = 300)

# ============================================================
# 9. FILTERING
# ============================================================
# 9a. Remove taxa with NA or empty Phylum
ps <- subset_taxa(ps, !is.na(Phylum) & !Phylum %in% c("", "uncharacterized"))
cat("After Phylum filter:", ntaxa(ps), "taxa\n")

# 9b. Remove samples below min_reads threshold
low_samples <- sample_names(ps)[sample_sums(ps) < mb_params$min_reads]
if (length(low_samples) > 0) {
  cat("Removing low-read samples (<", mb_params$min_reads, "reads):",
      paste(low_samples, collapse = ", "), "\n")
  ps <- prune_samples(sample_sums(ps) >= mb_params$min_reads, ps)
}
cat("After read filter:", nsamples(ps), "samples\n")

# 9c. Prevalence filter: keep taxa present in >= threshold fraction of samples
prev_threshold <- mb_params$prevalence_threshold * nsamples(ps)
keep_taxa <- taxa_names(ps)[apply(otu_table(ps), 1, function(x) sum(x > 0)) >= prev_threshold]
ps <- prune_taxa(keep_taxa, ps)
cat("After prevalence filter (>=", mb_params$prevalence_threshold * 100, "%):",
    ntaxa(ps), "taxa\n")

# ============================================================
# 10. AGGLOMERATE AT GENUS LEVEL
# ============================================================
ps_genus <- tax_glom(ps, taxrank = "Genus", NArm = TRUE)
cat("After genus agglomeration:", ntaxa(ps_genus), "genera\n")

# ============================================================
# 11. TRANSFORMS
# ============================================================
# Relative abundance
ps_rel <- transform_sample_counts(ps, function(x) x / sum(x))
ps_genus_rel <- transform_sample_counts(ps_genus, function(x) x / sum(x))

# CLR transform (compositional) at genus level
# Add pseudocount of 0.5 to handle zeros before CLR
genus_counts <- as(otu_table(ps_genus), "matrix")
if (!taxa_are_rows(ps_genus)) genus_counts <- t(genus_counts)
genus_clr <- apply(genus_counts + 0.5, 2, function(x) log(x) - mean(log(x)))
rownames(genus_clr) <- taxa_names(ps_genus)

# ============================================================
# 12. SAVE OUTPUTS
# ============================================================
saveRDS(ps, file.path(output_dir, "ps.rds"))
saveRDS(ps_rel, file.path(output_dir, "ps_rel.rds"))
saveRDS(ps_genus, file.path(output_dir, "ps_genus.rds"))
saveRDS(ps_genus_rel, file.path(output_dir, "ps_genus_rel.rds"))
saveRDS(metadata, file.path(output_dir, "metadata.rds"))
saveRDS(genus_clr, file.path(output_dir, "genus_clr.rds"))

# Save CLR as phyloseq object too (for Aitchison distance in beta diversity)
ps_clr <- ps_genus
otu_table(ps_clr) <- otu_table(genus_clr, taxa_are_rows = TRUE)
saveRDS(ps_clr, file.path(output_dir, "ps_clr.rds"))

cat("\n=== Preprocessing complete ===\n")
cat("Saved to:", output_dir, "\n")
cat("  ps.rds          - raw counts (filtered, F1-F4)\n")
cat("  ps_rel.rds      - relative abundance (filtered)\n")
cat("  ps_genus.rds    - genus-glommed raw counts\n")
cat("  ps_genus_rel.rds- genus-glommed relative abundance\n")
cat("  ps_clr.rds      - genus-glommed CLR transformed\n")
cat("  genus_clr.rds   - CLR matrix (standalone)\n")
cat("  metadata.rds    - sample metadata\n")
cat("  Seq_Info.csv    - OTU ID to sequence mapping\n")
cat("  QC_LibrarySize.png  - library size barplot\n")
cat("  QC_PrevalenceAbundance.png - prevalence scatter\n")
