# FK49_Microbiome_3_Plots.R -----
# 
#
# Visualization of FK49 16S microbiome data.
# Loads preprocessed phyloseq objects + statistics results.
#
# Plot strategy:
#   - Default: both sexes together (Sex as shape, Treatment as color)
#   - Sex-faceted versions generated ONLY if Sex interaction was significant
#     (flag read from Statistics/sex_interaction_flag.csv)
#   - F5 excluded (no data in preprocessed objects)
#
# Plot categories:
#   Publication (9 plots): QC, composition, alpha, beta, DA, Lactobacillus
#   Exploratory (4 plots): per-phylum violin, Lactobacillales drill-down, absolute counts
#
# Proposed cuts from original ~50 plots (user confirms before final):
#   - Per-Order/Per-Family abundance loops (~15 plots) -> consolidated per-phylum violin
#   - Duplicate absolute+relative barplot pairs (~20 plots) -> one per level
#   - plot_abundance() violin calls duplicating plot_Micro_Bars() -> merged
#   - Lactobacillus-specific plots (~8) -> consolidated into 2-3 panels
#
# Inputs:  ps.rds, ps_rel.rds, ps_genus.rds, ps_genus_rel.rds, ps_clr.rds,
#          metadata.rds, stats CSVs (from Scripts 1 & 2)
# Outputs: PNG plots in PATHS$microbiome$output_plots
###############################################################################
###############################################################################
# FK49_Microbiome_3_Plots.R
###############################################################################

rm(list = ls())
gc()

library(tidyverse)
library(phyloseq)
library(microbiome)
library(vegan)
library(ggrepel)
library(pheatmap)
library(gridExtra)
library(RColorBrewer)
library(emmeans)

source("FK49_Definitions.R")

# CONFIGURATION -----
mb_params <- PARAMETERS$microbiome
output_dir <- PATHS$microbiome$output
stats_dir <- PATHS$microbiome$output_stats
plots_dir <- PATHS$microbiome$output_plots
explor_dir <- file.path(plots_dir, "Exploratory")

if (!dir.exists(plots_dir)) dir.create(plots_dir, recursive = TRUE)
if (!dir.exists(explor_dir)) dir.create(explor_dir, recursive = TRUE)

ps <- readRDS(file.path(output_dir, "ps.rds"))
ps_rel <- readRDS(file.path(output_dir, "ps_rel.rds"))
ps_genus <- readRDS(file.path(output_dir, "ps_genus.rds"))
ps_genus_rel <- readRDS(file.path(output_dir, "ps_genus_rel.rds"))
ps_clr <- readRDS(file.path(output_dir, "ps_clr.rds"))
metadata <- readRDS(file.path(output_dir, "metadata.rds"))

feces_labels <- c(F1 = "-1wks", F2 = "0wks", F3 = "3wks", F4 = "7wks")
feces_labeller <- as_labeller(feces_labels)
feces_order <- names(feces_labels)

sex_flag_file <- file.path(stats_dir, "sex_interaction_flag.csv")
sex_stratify <- FALSE
if (file.exists(sex_flag_file)) {
  sex_flag <- read.csv(sex_flag_file)
  sex_stratify <- sex_flag$sex_interaction_significant[1]
}

cat("Sex stratification for plots:", sex_stratify, "\n")

theme_mb <- function() {
  theme_minimal() +
    theme(
      axis.line = element_line(color = "black", linewidth = 0.5),
      axis.ticks = element_line(color = "black", linewidth = 0.5),
      axis.title = element_text(size = 14, face = "bold"),
      axis.text = element_text(size = 12, face = "bold"),
      plot.title = element_blank(),
      legend.title = element_text(size = 10),
      legend.text = element_text(size = 10),
      panel.grid = element_blank()
    )
}

treat_cols <- Treatment_colors[c("Ctrl", "TAM")]
sex_shapes <- Sex_shape

sig_stars <- function(p) {
  ifelse(is.na(p), "NA",
         ifelse(p < 0.001, "***",
                ifelse(p < 0.01, "**",
                       ifelse(p < 0.05, "*", "ns"))))
}

cat("=== FK49 Microbiome Plots ===\n\n")



# PLOT 1: LIBRARY SIZE --------------
cat("Plot 1: Library size\n")

lib_sizes <- data.frame(
  SampleID = sample_names(ps),
  Reads = sample_sums(ps),
  Treatment = as.character(sample_data(ps)$Treatment),
  Feces = factor(as.character(sample_data(ps)$Feces), levels = feces_order)
)

p1 <- ggplot(lib_sizes, aes(x = SampleID, y = Reads, fill = Treatment)) +
  geom_col() +
  facet_grid(~Feces, scales = "free_x", space = "free_x") +
  scale_fill_manual(values = treat_cols) +
  labs(x = NULL, y = "Total reads") +
  theme_mb() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 7),
        legend.position = "bottom")

ggsave(file.path(plots_dir, "01_QC_LibrarySize.png"),
       p1, width = 10, height = 5, dpi = 300)



# PLOT 2: PHYLUM COMPOSITION -----
cat("Plot 2: Phylum composition\n")

ps_phylum <- tax_glom(ps_rel, "Phylum", NArm = TRUE)

phylum_abund <- psmelt(ps_phylum) %>%
  group_by(Phylum) %>%
  summarise(mean_abund = mean(Abundance), .groups = "drop") %>%
  arrange(desc(mean_abund))

top_phyla <- phylum_abund$Phylum[1:min(8, nrow(phylum_abund))]

ps_phylum_filt <- ps_phylum
tax_table(ps_phylum_filt)[
  !tax_table(ps_phylum_filt)[, "Phylum"] %in% top_phyla, "Phylum"
] <- "Other"

p2_data <- psmelt(ps_phylum_filt) %>%
  mutate(
    Animal = factor(Animal),
    Phylum = factor(Phylum, levels = c(top_phyla, "Other")),
    Feces = factor(Feces, levels = feces_order),
    Treatment = factor(Treatment, levels = c("Ctrl", "TAM"))
  ) %>%
  group_by(SampleID, Animal, Phylum, Treatment, Feces, Sex) %>%
  summarise(Abundance = sum(Abundance, na.rm = TRUE), .groups = "drop") %>%
  arrange(Feces, Treatment, Animal) %>%
  group_by(Feces) %>%
  mutate(x_pos = match(Animal, unique(Animal)) + ifelse(Treatment == "TAM", 0.1, 0)) %>%
  ungroup()

treatment_labels <- p2_data %>%
  distinct(Feces, Treatment, Animal, x_pos) %>%
  group_by(Feces, Treatment) %>%
  summarise(x = mean(x_pos), .groups = "drop")

p2 <- ggplot(p2_data, aes(x = x_pos, y = Abundance, fill = Phylum)) +
  geom_col(width = 1) +
  facet_wrap(~Feces, scales = "free_x", space = "free_x") +
  scale_fill_manual(values=Phylum_colors)+
  scale_x_continuous(
    breaks = treatment_labels$x,
    labels = treatment_labels$Treatment,
    expand = expansion(mult = c(0.01, 0.01))
  ) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(x = NULL, y = "Relative abundance") +
  theme_mb() +
  theme(legend.position = "right")

ggsave(file.path(plots_dir, "02_Composition_Phylum.png"),
       p2, width = 10, height = 6, dpi = 300)



# PLOT 3: TOP-15 GENERA ------

cat("Plot 3: Top-15 Genera composition\n")

genus_abund <- psmelt(ps_genus_rel) %>%
  group_by(Genus) %>%
  summarise(mean_abund = mean(Abundance), .groups = "drop") %>%
  arrange(desc(mean_abund))

top_genera <- genus_abund$Genus[1:min(15, nrow(genus_abund))]

ps_genus_filt <- ps_genus_rel
tax_table(ps_genus_filt)[
  !tax_table(ps_genus_filt)[, "Genus"] %in% top_genera, "Genus"
] <- "Other"

p3_data <- psmelt(ps_genus_filt) %>%
  mutate(
    Animal = factor(Animal),
    Genus = factor(Genus, levels = c(top_genera, "Other")),
    Feces = factor(Feces, levels = feces_order),
    Treatment = factor(Treatment, levels = c("Ctrl", "TAM"))
  ) %>%
  group_by(SampleID, Animal, Genus, Treatment, Feces, Sex) %>%
  summarise(Abundance = sum(Abundance, na.rm = TRUE), .groups = "drop") %>%
  arrange(Feces, Treatment, Animal) %>%
  group_by(Feces) %>%
  mutate(x_pos = match(Animal, unique(Animal)) + ifelse(Treatment == "TAM", 0.1, 0)) %>%
  ungroup()

treatment_labels <- p3_data %>%
  distinct(Feces, Treatment, Animal, x_pos) %>%
  group_by(Feces, Treatment) %>%
  summarise(x = mean(x_pos), .groups = "drop")

p3 <- ggplot(p3_data, aes(x = x_pos, y = Abundance, fill = Genus)) +
  geom_col(width = 1) +
  facet_wrap(~Feces, scales = "free_x", space = "free_x") +
  scale_x_continuous(
    breaks = treatment_labels$x,
    labels = treatment_labels$Treatment,
    expand = expansion(mult = c(0.01, 0.01))
  ) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(x = NULL, y = "Relative abundance") +
  theme_mb() +
  theme(legend.position = "right")

ggsave(file.path(plots_dir, "03_Composition_Top15Genera.png"),
       p3, width = 15, height = 6, dpi = 300)



# PLOT 4: ALPHA DIVERSITY -----
cat("Plot 4: Alpha diversity\n")

alpha_div <- estimate_richness(
  ps,
  measures = c("Observed", "Shannon", "Simpson")
) %>%
  rownames_to_column("SampleID") %>%
  mutate(SampleID = gsub("^X", "", SampleID)) %>%
  left_join(metadata, by = "SampleID") %>%
  mutate(
    Feces = factor(Feces, levels = feces_order),
    Time = as.numeric(dplyr::recode(as.character(Feces), F1 = "-1", F2 = "0", F3 = "3", F4 = "7")),
    Treatment = factor(Treatment, levels = c("Ctrl", "TAM")),
    Sex = factor(Sex, levels = c("female", "male"))
  )

alpha_pwc_file <- file.path(stats_dir, "alpha_diversity_posthoc.csv")
alpha_pwc <- if (file.exists(alpha_pwc_file)) read.csv(alpha_pwc_file) else NULL

for (a in c("Observed", "Shannon", "Simpson")) {
  
  alpha_summary <- alpha_div %>%
    group_by(Time, Treatment) %>%
    summarise(
      mean = mean(.data[[a]], na.rm = TRUE),
      sd = sd(.data[[a]], na.rm = TRUE),
      .groups = "drop"
    )
  
  sig_labels <- if (!is.null(alpha_pwc)) {
    alpha_pwc %>% filter(Measure == a) %>%
      dplyr::select(Feces, significance) %>%
      mutate(Time = as.numeric(dplyr::recode(as.character(Feces), F1 = "-1", F2 = "0", F3 = "3", F4 = "7")))
  } else {
    data.frame(Time = c(-1, 0, 3, 7), significance = rep("", 4))
  }
  
  p4 <- ggplot() +
    geom_jitter(
      data = alpha_div,
      aes(x = Time, y = .data[[a]],
          fill = Treatment, color = Treatment, shape = Sex),
      width = 0.15, size = 2.5, alpha = 0.5
    ) +
    geom_line(
      data = alpha_summary,
      aes(x = Time, y = mean, group = Treatment, color = Treatment),
      linewidth = 1
    ) +
    geom_ribbon(
      data = alpha_summary,
      aes(x = Time, fill = Treatment, group = Treatment,
          ymin = mean - sd, ymax = mean + sd),
      alpha = 0.1, color = NA
    ) +
    geom_text(
      data = sig_labels,
      aes(x = Time,
          y = max(alpha_div[[a]], na.rm = TRUE) * 1.05,
          label = significance),
      size = 3, fontface = "italic"
    ) +
    scale_color_manual(values = treat_cols) +
    scale_fill_manual(values = treat_cols) +
    scale_shape_manual(values = sex_shapes) +
    scale_x_continuous(breaks = c(-1, 0, 3, 7)) +
    labs(x = "Time on CDHFD [wks]", y = a) +
    theme_mb() +
    theme(legend.position = "bottom")
  
  ggsave(
    file.path(plots_dir, paste0("04_AlphaDiversity_", a, ".png")),
    p4, width = 8, height = 5, dpi = 300
  )
  
  if (sex_stratify) {
    p4_sex <- p4 + facet_wrap(~Sex)
    ggsave(
      file.path(plots_dir, paste0("04_AlphaDiversity_", a, "_bySex.png")),
      p4_sex, width = 10, height = 5, dpi = 300
    )
  }
}


# PLOT 5: PCoA -----------

cat("Plot 5: PCoA\n")

plot_pcoa <- function(ps_obj, dist_method, title, filename) {
  
  meta_pcoa <- as(sample_data(ps_obj), "data.frame")
  meta_pcoa$Feces <- factor(meta_pcoa$Feces, levels = feces_order)
  
  if (dist_method == "aitchison") {
    clr_mat <- as(otu_table(ps_clr), "matrix")
    if (!taxa_are_rows(ps_clr)) clr_mat <- t(clr_mat)
    dist_obj <- dist(t(clr_mat), method = "euclidean")
  } else {
    rel_mat <- as(otu_table(ps_genus_rel), "matrix")
    if (!taxa_are_rows(ps_genus_rel)) rel_mat <- t(rel_mat)
    dist_obj <- vegdist(t(rel_mat), method = dist_method)
  }
  
  pcoa_res <- cmdscale(dist_obj, k = 2, eig = TRUE)
  
  pcoa_df <- data.frame(
    PC1 = pcoa_res$points[, 1],
    PC2 = pcoa_res$points[, 2],
    SampleID = rownames(pcoa_res$points)
  ) %>%
    left_join(meta_pcoa, by = "SampleID")
  
  var_explained <- pcoa_res$eig[1:2] / sum(pcoa_res$eig) * 100
  
  gg <- ggplot(
    pcoa_df,
    aes(x = PC1, y = PC2, color = Treatment, shape = Sex)
  ) +
    geom_point(size = 3, alpha = 0.8) +
    stat_ellipse(aes(group = Treatment), level = 0.95, linetype = "dashed") +
    facet_wrap(~Feces, labeller = feces_labeller) +
    scale_color_manual(values = treat_cols) +
    scale_shape_manual(values = sex_shapes) +
    labs(
      title = title,
      x = paste0("PC1 (", round(var_explained[1], 1), "%)"),
      y = paste0("PC2 (", round(var_explained[2], 1), "%)")
    ) +
    theme_bw() +
    theme(text = element_text(family = "Liberation Sans"))
  
  ggsave(file.path(plots_dir, filename), gg, width = 10, height = 8, dpi = 300)
}

plot_pcoa(ps_genus_rel, "bray", "PCoA — Bray-Curtis", "05_PCoA_BrayCurtis.png")
plot_pcoa(ps_genus_rel, "aitchison", "PCoA — Aitchison", "05_PCoA_Aitchison.png")



# PLOT 6: VOLCANO ---------
cat("Plot 6: Volcano plots\n")

ancom_coef_file <- file.path(stats_dir, "ancombc2_model_coefficients.csv")

if (file.exists(ancom_coef_file)) {
  
  ancom_coef <- read.csv(ancom_coef_file)
  volcano_list <- list()
  
  for (f in mb_params$timepoints) {
    
    log2fc_col <- paste0("log2FC_", f)
    q_col <- paste0("q_", f)
    
    if (!(log2fc_col %in% colnames(ancom_coef)) ||
        !(q_col %in% colnames(ancom_coef))) {
      cat("Skipping", f, "— columns not found\n")
      next
    }
    
    vol_df <- ancom_coef %>%
      filter(!is.na(.data[[log2fc_col]]), !is.na(.data[[q_col]])) %>%
      mutate(
        neg_log10q = -log10(.data[[q_col]]),
        significance = sig_stars(.data[[q_col]]),
        label = ifelse(
          .data[[q_col]] < mb_params$fdr_threshold & abs(.data[[log2fc_col]]) > 1,
          as.character(Genus), NA
        )
      )
    
    p_vol <- ggplot(vol_df, aes(x = .data[[log2fc_col]], y = neg_log10q)) +
      geom_point(aes(color = .data[[q_col]] < mb_params$fdr_threshold), size = 2, alpha = 0.7) +
      geom_text_repel(aes(label = label), size = 3, max.overlaps = 15, box.padding = 0.5, na.rm = TRUE) +
      geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
      geom_hline(yintercept = -log10(mb_params$fdr_threshold), linetype = "dashed", color = "grey50") +
      scale_color_manual(
        values = c("TRUE" = "#8B0000", "FALSE" = "grey70"),
        labels = c("TRUE" = "FDR < 0.05", "FALSE" = "ns"), name = NULL
      ) +
      labs(
        title = paste0("Volcano — ", f, " (", feces_labels[f], ")"),
        x = "Log2 Fold Change (TAM vs Ctrl)",
        y = "-log10(q-value)"
      ) +
      theme_bw() +
      theme(text = element_text(family = "Liberation Sans"))
    
    volcano_list[[f]] <- p_vol
    
    ggsave(file.path(plots_dir, paste0("06_Volcano_", f, ".png")),
           p_vol, width = 6, height = 5, dpi = 300)
  }
  
  if (length(volcano_list) > 0) {
    p6_combined <- do.call("grid.arrange", c(volcano_list, ncol = 2))
    ggsave(file.path(plots_dir, "06_Volcanos_combined.png"),
           p6_combined, width = 12, height = 10, dpi = 300)
  }
}


# PLOT 7A: HEATMAP CLR -which genera dominate per sample ------------------
cat("Plot 7: Heatmap\n")

str(tax_table(ps_clr))
tax_table(ps_clr)[, "Genus"]
plot_genera <- c(tax_table(ps_clr)[, "Genus"])

clr_mat <- as(otu_table(ps_clr), "matrix")
if (!taxa_are_rows(ps_clr)) clr_mat <- t(clr_mat)

tax_info <- as.data.frame(tax_table(ps_clr)) %>%
  tibble::rownames_to_column("OTU") %>%
  dplyr::select(OTU, Genus)

clr_df <- as.data.frame(t(clr_mat)) %>%
  tibble::rownames_to_column("SampleID") %>%
  tidyr::pivot_longer(-SampleID, names_to = "OTU", values_to = "CLR") %>%
  dplyr::left_join(tax_info, by = "OTU") %>%
  dplyr::left_join(metadata, by = "SampleID") %>%
  filter(Genus %in% plot_genera)

hm <- clr_df %>%
  dplyr::select(Genus, SampleID, CLR) %>%
  tidyr::pivot_wider(names_from = SampleID, values_from = CLR) %>%
  as.data.frame()

ann <- metadata %>%
  dplyr::select(SampleID, Sex, Treatment, Diet_short, Week) %>%
  dplyr::mutate(Week = as.numeric(Week-1))%>%# Week in metadata is time after tam injection, but experimnt start 1 week earlier
  dplyr::filter(SampleID %in% colnames(hm)) %>%
  dplyr::arrange(match(SampleID, colnames(hm))) %>%
  dplyr::mutate(Treatment = dplyr::if_else(Week == -1, "none", as.character(Treatment)),
                Week = factor(Week, levels = c(-1, 0, 3, 7))) %>%
  as.data.frame()

rownames(ann) <- ann$SampleID
ann$SampleID <- NULL

rownames(hm) <- hm$Genus
hm$Genus <- NULL
hm <- as.matrix(hm)
Diet_colors <- c(
  "CDHFD"   = "#B39BC8",
  "CDHFD13" = "#C8B6D9",
  "ND"      = "#D2A679"
)
Week_colors <- c(
  "-1" = "#E5E7EB",
  "0"  = "#CBD5E1",
  "3"  = "#94A3B8",
  "7"  = "#64748B"
)
ann_colors <- list(
  Treatment = c(
    "none" = "lightgrey",
    "Ctrl" = adjustcolor(Treatment_colors["Ctrl"], alpha.f = 0.8),
    "TAM" = adjustcolor(Treatment_colors["TAM"], alpha.f = 0.8)
  ),
  Diet_short = Diet_colors[c("ND","CDHFD")],
  Sex = Sex_colors,
  Week = Week_colors,
  Phylum=Phylum_colors
)
ann <- ann[, c("Sex" ,"Treatment", "Week", "Diet_short"), drop = FALSE]

ann_row <- as.data.frame(tax_table(ps_clr)) %>%
  tibble::rownames_to_column("OTU") %>%
  dplyr::select(OTU, Genus, Phylum) %>%
  dplyr::filter(Genus %in% rownames(hm)) %>%
  dplyr::distinct(Genus, .keep_all = TRUE) %>%
  tibble::column_to_rownames("Genus") %>%
  dplyr::select(Phylum)

ann_row <- ann_row[rownames(hm), , drop = FALSE]
colnames(ann_row) <- "P."
# Heatmap ohne Row scaling
heatmap1<-pheatmap(
  hm,
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  scale = "none",
  annotation_col = ann,
  annotation_row =ann_row,
  annotation_colors = ann_colors,
  show_colnames = F,
  border_color = NA,
  fontsize = 8,
  main = "Genera CLR"
)
ggsave(file.path(plots_dir, paste0("07A_Heat_CLR_scaled.png")),
       heatmap1, width = 6, height = 6, dpi = 300,bg="white")
# Heatmap mit Row scaling
heatmap2<-pheatmap(
  hm,
  cluster_rows = TRUE,
  cluster_cols = F,
  scale = "row",
  annotation_col = ann,
  annotation_row =ann_row,
  annotation_colors = ann_colors,
  show_colnames = F,
  border_color = NA,
  fontsize = 8,
  main = "Genera CLR-row scaled"
)
ggsave(file.path(plots_dir, paste0("07B_Heat_CLR_scaled.png")),
       heatmap2, width = 6, height = 6, dpi = 300,bg="white")

# Heatmap nur Week 3 und 7
hm_3_7 <- hm[, rownames(ann)[ann$Week %in% c("3", "7")], drop = FALSE]
ann_3_7 <- ann[rownames(ann) %in% colnames(hm_3_7), , drop = FALSE]

heatmap3 <- pheatmap(
  hm_3_7,
  cluster_rows = TRUE,
  cluster_cols = T,
  scale = "none",
  annotation_col = ann_3_7,
  annotation_row = ann_row,
  annotation_colors = ann_colors,
  show_colnames = FALSE,
  border_color = NA,
  fontsize = 8,
  main = "Genera CLR — Week 3 and 7"
)

ggsave(
  file.path(plots_dir, "07C_Heat_CLR_Week3_7.png"),
  heatmap3, width = 5.9, height = 6, dpi = 300, bg = "white"
)
# Heatmap nur Week 3 und 7
hm_3_7 <- hm[, rownames(ann)[ann$Week %in% c("3", "7")], drop = FALSE]
ann_3_7 <- ann[rownames(ann) %in% colnames(hm_3_7), , drop = FALSE]

heatmap4 <- pheatmap(
  hm_3_7,
  cluster_rows = TRUE,
  cluster_cols = T,
  scale = "row",
  annotation_col = ann_3_7,
  annotation_row = ann_row,
  annotation_colors = ann_colors,
  show_colnames = FALSE,
  border_color = NA,
  fontsize = 8,
  main = "Genera CLR — Week 3 and 7"
)

ggsave(
  file.path(plots_dir, "07D_Heat_CLR_Week3_7_scaled.png"),
  heatmap4, width = 5.9, height = 6, dpi = 300, bg = "white"
)
# PLOT 8: SIGNIFICANT GENERA BOXPLOTS — CLR ------
cat("Plot 8: Significant genera boxplots\n")

consensus_file <- file.path(stats_dir, "da_consensus.csv")

if (file.exists(consensus_file)) {
  
  consensus <- read.csv(consensus_file)
  
  sig_genera <- consensus %>%
    filter(consensus %in% c("consensus_significant", "suggestive", "aldex_significant")) %>%
    pull(Genus) %>%
    unique()
}

if (exists("sig_genera") && length(sig_genera) > 0) {
  
  clr_mat <- as(otu_table(ps_clr), "matrix")
  if (!taxa_are_rows(ps_clr)) clr_mat <- t(clr_mat)
  
  tax_info <- as.data.frame(tax_table(ps_clr)) %>%
    rownames_to_column("OTU") %>%
    dplyr::select(OTU, Genus)
  
  clr_df <- as.data.frame(t(clr_mat)) %>%
    rownames_to_column("SampleID") %>%
    pivot_longer(-SampleID, names_to = "OTU", values_to = "CLR") %>%
    left_join(tax_info, by = "OTU") %>%
    left_join(metadata, by = "SampleID") %>%
    filter(Genus %in% sig_genera) %>%
    mutate(
      Feces = factor(Feces, levels = feces_order),
      Treatment = factor(Treatment, levels = c("Ctrl", "TAM")),
      Sex = factor(Sex, levels = c("female", "male"))
    )
  
  boxplot_list <- list()
  
  for (g in sig_genera) {
    
    gen_df <- clr_df %>% filter(Genus == g)
    
    if (nrow(gen_df) == 0) next
    
    p_box <- ggplot(gen_df, aes(x = Feces, y = CLR, fill = Treatment, color = Treatment)) +
      geom_boxplot(
        position = position_dodge(0.8), outlier.shape = NA,
        alpha = 0.3, width = 0.5
      ) +
      geom_point(
        aes(shape = Sex),
        position = position_jitterdodge(jitter.width = 0, dodge.width = 0.8),
        size = 1.8, alpha = 0.7
      ) +
      scale_fill_manual(values = treat_cols) +
      scale_color_manual(values = treat_cols) +
      scale_shape_manual(values = sex_shapes) +
      scale_x_discrete(labels = feces_labels) +
      labs(
        title = g, x = "Timepoint", y = "CLR abundance",
        fill = "Treatment", color = "Treatment", shape = "Sex"
      ) +
      theme_mb() +
      theme(
        legend.position = "bottom",
        plot.title = element_text(face = "italic", size = 12)
      )
    
    boxplot_list[[g]] <- p_box
    
    ggsave(file.path(plots_dir, paste0("08_Boxplot_", g, ".png")),
           p_box, width = 6, height = 5, dpi = 300)
    
    if (sex_stratify) {
      p8_sex <- p_box + facet_wrap(~Sex)
      ggsave(file.path(plots_dir, paste0("08_Boxplot_", g, "_bySex.png")),
             p8_sex, width = 8, height = 5, dpi = 300)
    }
  }
  
  if (length(boxplot_list) > 0) {
    
    n_per_page <- 9
    n_pages <- ceiling(length(boxplot_list) / n_per_page)
    
    for (page in 1:n_pages) {
      
      start_idx <- (page - 1) * n_per_page + 1
      end_idx <- min(page * n_per_page, length(boxplot_list))
      page_plots <- boxplot_list[start_idx:end_idx]
      
      p_combined <- do.call("grid.arrange", c(page_plots, ncol = 3))
      
      ggsave(file.path(plots_dir, paste0("08_Boxplots_page", page, ".png")),
             p_combined, width = 15, height = 10, dpi = 300)
    }
  }
}


# PLOT 9: LACTOBACILLUS --------
cat("Plot 9: Lactobacillus\n")

ps_lacto <- subset_taxa(ps_genus_rel, Genus == "Lactobacillus")

if (ntaxa(ps_lacto) > 0) {
  
  lacto_df <- psmelt(ps_lacto) %>%
    mutate(
      Feces = factor(Feces, levels = feces_order),
      Time = as.numeric(dplyr::recode(as.character(Feces), F1 = "-1", F2 = "0", F3 = "3", F4 = "7")),
      Treatment = factor(Treatment, levels = c("Ctrl", "TAM")),
      Sex = factor(Sex, levels = c("female", "male"))
    ) %>%
    group_by(SampleID, Feces, Time, Treatment, Sex, Animal) %>%
    summarise(Abundance = sum(Abundance), .groups = "drop")
  
  lacto_summary <- lacto_df %>%
    group_by(Time, Treatment) %>%
    summarise(mean = mean(Abundance), sd = sd(Abundance), .groups = "drop")
  
  lacto_pwc_file <- file.path(stats_dir, "lactobacillus_posthoc_pooled.csv")
  lacto_pwc <- if (file.exists(lacto_pwc_file)) read.csv(lacto_pwc_file) else NULL
  
  sig_labels_lacto <- if (!is.null(lacto_pwc)) {
    lacto_pwc %>%
      dplyr::select(Feces, significance) %>%
      mutate(Time = as.numeric(dplyr::recode(as.character(Feces), F1 = "-1", F2 = "0", F3 = "3", F4 = "7")))
  } else {
    data.frame(Time = c(-1, 0, 3, 7), significance = rep("", 4))
  }
  
  p9 <- ggplot() +
    geom_jitter(data = lacto_df, aes(x = Time, y = Abundance, fill = Treatment, color = Treatment, shape = Sex), width = 0.15, size = 2.5, alpha = 0.4) +
    geom_line(data = lacto_summary, aes(x = Time, y = mean, group = Treatment, color = Treatment), linewidth = 1) +
    geom_ribbon(data = lacto_summary, aes(x = Time, fill = Treatment, group = Treatment, ymin = mean - sd, ymax = mean + sd), alpha = 0.1, color = NA) +
    geom_text(data = sig_labels_lacto, aes(x = Time, y = max(lacto_df$Abundance, na.rm = TRUE) * 1.08, label = significance), size = 3.5, fontface = "italic") +
    scale_color_manual(values = treat_cols) +
    scale_fill_manual(values = treat_cols) +
    scale_shape_manual(values = sex_shapes) +
    scale_x_continuous(breaks = c(-1, 0, 3, 7)) +
    labs(x = "Time on CDHFD [wks]", y = "Relative abundance", title = "Lactobacillus") +
    theme_mb() +
    theme(legend.position = "bottom", plot.title = element_text(face = "italic", size = 12))
  
  ggsave(file.path(plots_dir, "09_Lactobacillus_targeted.png"), p9, width = 8, height = 5, dpi = 300)
  
  lacto_male <- lacto_df %>% filter(Sex == "male")
  lacto_male_summary <- lacto_male %>% group_by(Time, Treatment) %>% summarise(mean = mean(Abundance), sd = sd(Abundance), .groups = "drop")
  
  lacto_pwc_male_file <- file.path(stats_dir, "lactobacillus_posthoc_by_sex.csv")
  lacto_pwc_male <- if (file.exists(lacto_pwc_male_file)) read.csv(lacto_pwc_male_file) %>% filter(Sex == "male") else NULL
  
  sig_labels_male <- if (!is.null(lacto_pwc_male)) {
    lacto_pwc_male %>% dplyr::select(Feces, significance) %>% mutate(Time = as.numeric(dplyr::recode(as.character(Feces), F1 = "-1", F2 = "0", F3 = "3", F4 = "7")))
  } else {
    data.frame(Time = c(-1, 0, 3, 7), significance = rep("", 4))
  }
  
  p9_male <- ggplot() +
    geom_jitter(data = lacto_male, aes(x = Time, y = Abundance, fill = Treatment, color = Treatment), width = 0.15, size = 3, alpha = 0.4) +
    geom_line(data = lacto_male_summary, aes(x = Time, y = mean, group = Treatment, color = Treatment), linewidth = 1) +
    geom_ribbon(data = lacto_male_summary, aes(x = Time, fill = Treatment, group = Treatment, ymin = mean - sd, ymax = mean + sd), alpha = 0.1, color = NA) +
    geom_text(data = sig_labels_male, aes(x = Time, y = max(lacto_male$Abundance, na.rm = TRUE) * 1.08, label = significance), size = 3.5, fontface = "italic") +
    scale_color_manual(values = treat_cols) +
    scale_fill_manual(values = treat_cols) +
    scale_x_continuous(breaks = c(-1, 0, 3, 7)) +
    labs(x = "Time on CDHFD [wks]", y = "Relative abundance", title = "Lactobacillus (male)") +
    theme_mb() +
    theme(legend.position = "bottom", plot.title = element_text(face = "italic", size = 12))
  
  ggsave(file.path(plots_dir, "09_Lactobacillus_targeted_male.png"), p9_male, width = 8, height = 5, dpi = 300)
  
  lacto_female <- lacto_df %>% filter(Sex == "female")
  lacto_female_summary <- lacto_female %>% group_by(Time, Treatment) %>% summarise(mean = mean(Abundance), sd = sd(Abundance), .groups = "drop")
  
  lacto_pwc_female <- if (file.exists(lacto_pwc_male_file)) {
    read.csv(lacto_pwc_male_file) %>% filter(Sex == "female")
  } else NULL
  
  sig_labels_female <- if (!is.null(lacto_pwc_female)) {
    lacto_pwc_female %>% dplyr::select(Feces, significance) %>% mutate(Time = as.numeric(dplyr::recode(as.character(Feces), F1 = "-1", F2 = "0", F3 = "3", F4 = "7")))
  } else {
    data.frame(Time = c(-1, 0, 3, 7), significance = rep("", 4))
  }
  
  p9b <- ggplot() +
    geom_jitter(data = lacto_female, aes(x = Time, y = Abundance, fill = Treatment, color = Treatment), width = 0.15, size = 3, alpha = 0.4) +
    geom_line(data = lacto_female_summary, aes(x = Time, y = mean, group = Treatment, color = Treatment), linewidth = 1) +
    geom_ribbon(data = lacto_female_summary, aes(x = Time, fill = Treatment, group = Treatment, ymin = mean - sd, ymax = mean + sd), alpha = 0.1, color = NA) +
    geom_text(data = sig_labels_female, aes(x = Time, y = max(lacto_female$Abundance, na.rm = TRUE) * 1.08, label = significance), size = 3.5, fontface = "italic") +
    scale_color_manual(values = treat_cols) +
    scale_fill_manual(values = treat_cols) +
    scale_x_continuous(breaks = c(-1, 0, 3, 7)) +
    labs(x = "Time on CDHFD [wks]", y = "Relative abundance", title = "Lactobacillus (female)") +
    theme_mb() +
    theme(legend.position = "bottom", plot.title = element_text(face = "italic", size = 12))
  
  ggsave(file.path(plots_dir, "09_Lactobacillus_targeted_female.png"), p9b, width = 8, height = 5, dpi = 300)
}


# EXPLORATORY PLOTS -----
cat("\n--- Exploratory plots ---\n")


# PLOT 10: PER-PHYLUM VIOLIN -----
cat("Plot 10: Per-phylum violin\n")

unique_phyla <- unique(tax_table(ps_rel)[, "Phylum"])
unique_phyla <- unique_phyla[!is.na(unique_phyla)]

for (ph in unique_phyla) {
  
  ps_ph <- subset_taxa(ps_rel, Phylum == ph)
  
  df_ph <- psmelt(ps_ph) %>%
    mutate(
      Feces = factor(Feces, levels = feces_order),
      Treatment = factor(Treatment, levels = c("Ctrl", "TAM")),
      Sex = factor(Sex, levels = c("female", "male"))
    ) %>%
    filter(Abundance > 0)
  
  p10 <- ggplot(df_ph, aes(x = Feces, y = Abundance, color = Treatment, shape = Sex)) +
    geom_violin(fill = NA, aes(group = interaction(Feces, Treatment))) +
    geom_jitter(width = 0.2, size = 1, alpha = 0.3) +
    facet_wrap(~Order, ncol = 3) +
    scale_y_log10() +
    scale_color_manual(values = treat_cols) +
    scale_shape_manual(values = sex_shapes) +
    labs(x = "Timepoint", y = "Relative abundance", title = ph) +
    theme_mb() +
    theme(
      legend.position = "bottom",
      plot.title = element_text(face = "bold", size = 10) )
  
  ggsave(file.path(explor_dir, paste0("10_Violin_Phylum_", ph, ".png")),
         p10, width = 9, height = 6, dpi = 300)
}


# PLOT 11: LACTOBACILLALES DRILL-DOWN -----
cat("Plot 11: Lactobacillales drill-down\n")

ps_lacto_order <- subset_taxa(ps_genus_rel, Order == "Lactobacillales")

if (ntaxa(ps_lacto_order) > 0) {
  
  for (rank in c("Family", "Genus")) {
    
    df_lacto <- psmelt(ps_lacto_order) %>%
      mutate(
        Feces = factor(Feces, levels = feces_order),
        Treatment = factor(Treatment, levels = c("Ctrl", "TAM")),
        Sex = factor(Sex, levels = c("female", "male"))
      ) %>%
      group_by(SampleID, .data[[rank]], Feces, Treatment, Sex) %>%
      summarise(Abundance = sum(Abundance), .groups = "drop") %>%
      filter(Abundance > 0)
    
    p11 <- ggplot(df_lacto, aes(x = Feces, y = Abundance, fill = Treatment)) +
      geom_col(position = "dodge") +
      facet_wrap(~.data[[rank]], scales = "free_y", ncol = 3) +
      scale_fill_manual(values = treat_cols) +
      labs(
        x = "Timepoint", y = "Relative abundance",
        title = paste("Lactobacillales —", rank)
      ) +
      theme_mb() +
      theme(
        legend.position = "bottom",
        plot.title = element_text(face = "bold", size = 10),
        axis.text.x = element_text(angle = 45, hjust = 1)
      )
    
    ggsave(file.path(explor_dir, paste0("11_Lactobacillales_", rank, ".png")),
           p11, width = 9, height = 6, dpi = 300)
  }
}


# PLOT 12: ABSOLUTE COUNTS -----
cat("Plot 12: Absolute counts\n")

for (rank in c("Phylum", "Genus")) {
  
  ps_abs <- tax_glom(ps, rank, NArm = TRUE)
  
  if (rank == "Genus") {
    top_g <- taxa_names(ps_abs)[
      order(taxa_sums(ps_abs), decreasing = TRUE)
    ][1:min(20, ntaxa(ps_abs))]
    
    ps_abs <- prune_taxa(top_g, ps_abs)
  }
  
  df_abs <- psmelt(ps_abs) %>%
    mutate(
      Feces = factor(Feces, levels = feces_order),
      Treatment = factor(Treatment, levels = c("Ctrl", "TAM"))
    )
  
  p12 <- ggplot(df_abs, aes(x = SampleID, y = Abundance, fill = .data[[rank]])) +
    geom_col(width = 1) +
    facet_grid(Treatment ~ Feces, scales = "free_x", space = "free_x") +
    labs(x = NULL, y = "Absolute counts", title = paste("Absolute counts —", rank)) +
    theme_mb() +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      legend.position = "right",
      plot.title = element_text(face = "bold", size = 10)
    )
  
  ggsave(file.path(explor_dir, paste0("12_AbsoluteCounts_", rank, ".png")),
         p12, width = 10, height = 6, dpi = 300)
}

cat("\n=== Plots complete ===\n")
cat("Publication plots saved to:", plots_dir, "\n")
cat("Exploratory plots saved to:", explor_dir, "\n")