###############################################################################
# FK49_Microbiome_3_Plots.R
#
# Visualization of FK49 16S microbiome data.
# Loads preprocessed phyloseq objects + statistical results from Scripts 1 & 2.
#
# Plots:
# 1. Alpha diversity boxplots (per timepoint, Treatment comparison)
# 2. Stacked barplot — Phylum composition (per timepoint × Treatment)
# 3. Stacked barplot — Genus composition (top 15, per timepoint × Treatment)
# 4. Alpha diversity with Sex faceting (if Sex interaction significant)
# 5. PCoA ordination (Bray-Curtis + Aitchison)
# 6. Volcano plots — per-timepoint ANCOM-BC2 absolute treatment effects
# 7. Heatmap — significant genera across timepoints
# 8. Individual boxplots — significant genera (CLR-transformed)
# 9. Lactobacillus targeted plot
#
# Inputs: ps.rds, ps_genus.rds, ps_genus_rel.rds, ps_clr.rds, metadata.rds
# Statistics CSVs from Script 2
# Outputs: PNG figures in PATHS$microbiome$output_plots
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

source("FK49_Definitions.R")

# ============================================================
# CONFIGURATION
# ============================================================
mb_params <- PARAMETERS$microbiome
output_dir <- PATHS$microbiome$output
stats_dir <- PATHS$microbiome$output_stats
plots_dir <- PATHS$microbiome$output_plots

if (!dir.exists(plots_dir)) dir.create(plots_dir, recursive = TRUE)

# Load data
ps <- readRDS(file.path(output_dir, "ps.rds"))
ps_genus <- readRDS(file.path(output_dir, "ps_genus.rds"))
ps_genus_rel <- readRDS(file.path(output_dir, "ps_genus_rel.rds"))
ps_clr <- readRDS(file.path(output_dir, "ps_clr.rds"))
metadata <- readRDS(file.path(output_dir, "metadata.rds"))

# Feces display labels (internal = F1-F4, display = week labels)
feces_labels <- mb_params$feces_labels
feces_labeller <- as_labeller(feces_labels)

# Color palettes (defined in FK49_Definitions.R as top-level variables)
# Treatment_colors, Sex_colors, Sex_shape are sourced above

# Helper: significance stars
sig_stars <- function(p) {
 ifelse(is.na(p), "NA",
 ifelse(p < 0.001, "***",
 ifelse(p < 0.01, "**",
 ifelse(p < 0.05, "*", "ns"))))
}

# Read sex interaction flag
sex_flag <- read.csv(file.path(stats_dir, "sex_interaction_flag.csv"))
sex_interaction <- sex_flag$sex_interaction_significant[1]

cat("=== FK49 Microbiome Plots ===\n")
cat("Sex interaction significant:", sex_interaction, "\n\n")

# ============================================================
# PLOT 1: ALPHA DIVERSITY BOXPLOTS
# ============================================================
cat("--- Plot 1: Alpha Diversity ---\n")

alpha_div <- estimate_richness(ps, measures = mb_params$alpha_measures) %>%
 rownames_to_column(var = "SampleID") %>%
 mutate(SampleID = gsub("^X", "", SampleID)) %>%
 left_join(metadata, by = "SampleID") %>%
 mutate(Feces = factor(Feces, levels = mb_params$feces_order))

alpha_pwc <- read.csv(file.path(stats_dir, "alpha_diversity_posthoc.csv"))

for (a in mb_params$alpha_measures) {
 pwc_sub <- alpha_pwc %>%
 filter(Measure == a) %>%
 mutate(Feces = factor(Feces, levels = mb_params$feces_order))

 p <- ggplot(alpha_div, aes(x = Feces, y = .data[[a]], fill = Treatment)) +
 geom_boxplot(position = position_dodge(0.8), outlier.shape = NA, alpha = 0.7) +
 geom_point(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8),
 aes(shape = Sex), size = 2, alpha = 0.8) +
 geom_text(data = pwc_sub, aes(x = Feces, label = significance, y = max(alpha_div[[a]], na.rm = TRUE) * 1.05),
 vjust = 0, hjust = 0.5, inherit.aes = FALSE, size = 5) +
 scale_fill_manual(values = Treatment_colors) +
 scale_shape_manual(values = Sex_shape) +
 scale_x_discrete(labels = feces_labels) +
 labs(title = a, x = "Timepoint", y = a, fill = "Treatment", shape = "Sex") +
 theme_bw() +
 theme(text = element_text(family = "Liberation Sans"),
 axis.text.x = element_text(angle = 45, hjust = 1))

 ggsave(file.path(plots_dir, paste0("plot1_alpha_", a, ".png")), p,
 width = 6, height = 5, dpi = 300)
}

# ============================================================
# PLOT 2: PHYLUM COMPOSITION (stacked barplot)
# ============================================================
cat("--- Plot 2: Phylum Composition ---\n")

ps_phylum <- tax_glom(ps_genus_rel, taxrank = "Phylum")
phylum_mat <- as(otu_table(ps_phylum), "matrix")
if (!taxa_are_rows(ps_phylum)) phylum_mat <- t(phylum_mat)

phylum_df <- as.data.frame(t(phylum_mat)) %>%
 rownames_to_column("SampleID") %>%
 pivot_longer(-SampleID, names_to = "Phylum", values_to = "Abundance") %>%
 left_join(metadata, by = "SampleID") %>%
 mutate(Feces = factor(Feces, levels = mb_params$feces_order))

# Summarize by group
phylum_summary <- phylum_df %>%
 group_by(Feces, Treatment, Phylum) %>%
 summarise(mean_abundance = mean(Abundance, na.rm = TRUE), .groups = "drop") %>%
 mutate(Phylum = fct_reorder(Phylum, mean_abundance, .desc = TRUE))

# Keep top 8 phyla, rest as "Other"
top_phyla <- phylum_summary %>%
 group_by(Phylum) %>%
 summarise(total = sum(mean_abundance), .groups = "drop") %>%
 arrange(desc(total)) %>%
 slice_head(n = 8) %>%
 pull(Phylum)

phylum_summary <- phylum_summary %>%
 mutate(Phylum_display = ifelse(Phylum %in% top_phyla, as.character(Phylum), "Other"))

# FIX 5: Increased offset from 0.1 to 0.2 for better Ctrl/TAM separation
phylum_summary <- phylum_summary %>%
 mutate(x_pos = as.numeric(Feces) + ifelse(Treatment == "TAM", 0.2, 0))

p2 <- ggplot(phylum_summary, aes(x = x_pos, y = mean_abundance, fill = Phylum_display)) +
 geom_bar(stat = "identity", width = 0.35) +
 scale_x_continuous(
 breaks = as.numeric(factor(mb_params$feces_order)),
 labels = feces_labels
 ) +
 scale_fill_brewer(palette = "Set3") +
 labs(title = "Phylum Composition", x = "Timepoint", y = "Mean Relative Abundance",
 fill = "Phylum") +
 theme_bw() +
 theme(text = element_text(family = "Liberation Sans"),
 axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file.path(plots_dir, "plot2_phylum_composition.png"), p2,
 width = 8, height = 6, dpi = 300)

# ============================================================
# PLOT 3: GENUS COMPOSITION (top 15)
# ============================================================
cat("--- Plot 3: Genus Composition ---\n")

genus_mat <- as(otu_table(ps_genus_rel), "matrix")
if (!taxa_are_rows(ps_genus_rel)) genus_mat <- t(genus_mat)

genus_df <- as.data.frame(t(genus_mat)) %>%
 rownames_to_column("SampleID") %>%
 pivot_longer(-SampleID, names_to = "Genus", values_to = "Abundance") %>%
 left_join(metadata, by = "SampleID") %>%
 mutate(Feces = factor(Feces, levels = mb_params$feces_order))

genus_summary <- genus_df %>%
 group_by(Feces, Treatment, Genus) %>%
 summarise(mean_abundance = mean(Abundance, na.rm = TRUE), .groups = "drop")

top_genera <- genus_summary %>%
 group_by(Genus) %>%
 summarise(total = sum(mean_abundance), .groups = "drop") %>%
 arrange(desc(total)) %>%
 slice_head(n = 15) %>%
 pull(Genus)

genus_summary <- genus_summary %>%
 mutate(Genus_display = ifelse(Genus %in% top_genera, as.character(Genus), "Other"))

# FIX 5: Increased offset from 0.1 to 0.2 for better Ctrl/TAM separation
genus_summary <- genus_summary %>%
 mutate(x_pos = as.numeric(Feces) + ifelse(Treatment == "TAM", 0.2, 0))

p3 <- ggplot(genus_summary, aes(x = x_pos, y = mean_abundance, fill = Genus_display)) +
 geom_bar(stat = "identity", width = 0.35) +
 scale_x_continuous(
 breaks = as.numeric(factor(mb_params$feces_order)),
 labels = feces_labels
 ) +
 scale_fill_brewer(palette = "Set3") +
 labs(title = "Genus Composition (Top 15)", x = "Timepoint",
 y = "Mean Relative Abundance", fill = "Genus") +
 theme_bw() +
 theme(text = element_text(family = "Liberation Sans"),
 axis.text.x = element_text(angle = 45, hjust = 1),
 legend.position = "right",
 legend.text = element_text(size = 7))

ggsave(file.path(plots_dir, "plot3_genus_composition.png"), p3,
 width = 10, height = 6, dpi = 300)

# ============================================================
# PLOT 4: ALPHA DIVERSITY WITH SEX FACETING
# ============================================================
cat("--- Plot 4: Alpha Diversity (Sex-faceted) ---\n")

if (sex_interaction) {
 alpha_stratified <- read.csv(file.path(stats_dir, "alpha_diversity_stratified.csv"))

 for (a in mb_params$alpha_measures) {
 # FIX 3: Use sex-stratified p-values instead of pooled p-values
 pwc_strat <- alpha_stratified %>%
 filter(Measure == a) %>%
 mutate(Feces = factor(Feces, levels = mb_params$feces_order))

 p4 <- ggplot(alpha_div, aes(x = Feces, y = .data[[a]], fill = Treatment)) +
 geom_boxplot(position = position_dodge(0.8), outlier.shape = NA, alpha = 0.7) +
 geom_point(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8),
 size = 1.5, alpha = 0.6) +
 geom_text(data = pwc_strat, aes(x = Feces, label = significance,
 y = max(alpha_div[[a]], na.rm = TRUE) * 1.05),
 vjust = 0, hjust = 0.5, inherit.aes = FALSE, size = 4) +
 facet_wrap(~Sex, labeller = as_labeller(Sex_colors)) +
 scale_fill_manual(values = Treatment_colors) +
 scale_x_discrete(labels = feces_labels) +
 labs(title = paste(a, "(by Sex)"), x = "Timepoint", y = a, fill = "Treatment") +
 theme_bw() +
 theme(text = element_text(family = "Liberation Sans"),
 axis.text.x = element_text(angle = 45, hjust = 1))

 ggsave(file.path(plots_dir, paste0("plot4_alpha_sex_", a, ".png")), p4,
 width = 8, height = 5, dpi = 300)
 }
} else {
 cat(" Sex interaction not significant — skipping sex-faceted alpha plots\n")
}

# ============================================================
# PLOT 5: PCoA ORDINATION
# ============================================================
cat("--- Plot 5: PCoA Ordination ---\n")

plot_pcoa <- function(ps_obj, dist_method, title, filename) {
 meta_pcoa <- as(sample_data(ps_obj), "data.frame")
 meta_pcoa$Feces <- factor(meta_pcoa$Feces, levels = mb_params$feces_order)

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

 p <- ggplot(pcoa_df, aes(x = PC1, y = PC2, color = Treatment, shape = Sex)) +
 geom_point(size = 3, alpha = 0.8) +
 stat_ellipse(aes(group = Treatment), level = 0.95, linetype = "dashed") +
 facet_wrap(~Feces, labeller = feces_labeller) +
 scale_color_manual(values = Treatment_colors) +
 scale_shape_manual(values = Sex_shape) +
 labs(title = title,
 x = paste0("PC1 (", round(var_explained[1], 1), "%)"),
 y = paste0("PC2 (", round(var_explained[2], 1), "%)")) +
 theme_bw() +
 theme(text = element_text(family = "Liberation Sans"))

 ggsave(file.path(plots_dir, filename), p, width = 10, height = 8, dpi = 300)
}

plot_pcoa(ps_genus_rel, "bray", "PCoA — Bray-Curtis", "plot5_pcoa_bray.png")
plot_pcoa(ps_genus_rel, "aitchison", "PCoA — Aitchison", "plot5_pcoa_aitchison.png")

# ============================================================
# PLOT 6: VOLCANO PLOTS (per-timepoint ANCOM-BC2 absolute effects)
# ============================================================
cat("--- Plot 6: Volcano Plots ---\n")

# FIX 1: Read per-timepoint ANCOM-BC2 absolute effects (not interaction terms)
# ancombc2_model_coefficients.csv now contains log2FC_F1, q_F1, log2FC_F2, q_F2, etc.
ancom_coef_file <- file.path(stats_dir, "ancombc2_model_coefficients.csv")
if (!file.exists(ancom_coef_file)) {
 cat(" ancombc2_model_coefficients.csv not found — skipping volcano plots\n")
 cat(" (ANCOMBC not available or per-timepoint section not run)\n")
} else {
ancom_coef <- read.csv(ancom_coef_file)

volcano_list <- list()

for (f in mb_params$timepoints) {
 log2fc_col <- paste0("log2FC_", f)
 q_col <- paste0("q_", f)

 if (!(log2fc_col %in% colnames(ancom_coef)) || !(q_col %in% colnames(ancom_coef))) {
 cat(" Skipping", f, "— columns not found\n")
 next
 }

 vol_df <- ancom_coef %>%
 filter(!is.na(.data[[log2fc_col]]), !is.na(.data[[q_col]])) %>%
 mutate(
 neg_log10q = -log10(.data[[q_col]]),
 significance = sig_stars(.data[[q_col]]),
 label = ifelse(.data[[q_col]] < fdr_thresh & abs(.data[[log2fc_col]]) > 1,
 as.character(Genus), NA)
 )

 p_vol <- ggplot(vol_df, aes(x = .data[[log2fc_col]], y = neg_log10q)) +
 geom_point(aes(color = .data[[q_col]] < fdr_thresh), size = 2, alpha = 0.7) +
 geom_text_repel(aes(label = label), size = 3, max.overlaps = 15,
 box.padding = 0.5, na.rm = TRUE) +
 geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
 geom_hline(yintercept = -log10(fdr_thresh), linetype = "dashed", color = "grey50") +
 scale_color_manual(values = c("TRUE" = "#8B0000", "FALSE" = "grey70"),
 labels = c("TRUE" = "FDR < 0.05", "FALSE" = "ns"), name = NULL) +
 labs(title = paste0("Volcano — ", f, " (", feces_labels[f], ")"),
 x = "Log2 Fold Change (TAM vs Ctrl)",
 y = "-log10(q-value)") +
 theme_bw() +
 theme(text = element_text(family = "Liberation Sans"))

 volcano_list[[f]] <- p_vol
 ggsave(file.path(plots_dir, paste0("plot6_volcano_", f, ".png")), p_vol,
 width = 6, height = 5, dpi = 300)
}

if (length(volcano_list) > 0) {
 p6_combined <- do.call("grid.arrange", c(volcano_list, ncol = 2))
 ggsave(file.path(plots_dir, "plot6_volcanos_combined.png"), p6_combined,
 width = 12, height = 10, dpi = 300)
}
} # end else (ancombc2_model_coefficients.csv exists)

# ============================================================
# PLOT 7: HEATMAP — significant genera across timepoints
# ============================================================
cat("--- Plot 7: Heatmap ---\n")

consensus <- read.csv(file.path(stats_dir, "da_consensus.csv"))

# FIX 2: Filter to significant genera only (not all genera)
sig_genera <- consensus %>%
 filter(consensus %in% c("consensus_significant", "suggestive", "aldex_significant")) %>%
 pull(Genus) %>%
 unique()

if (length(sig_genera) > 0) {
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
 mutate(Feces = factor(Feces, levels = mb_params$feces_order))

 # Helper function for heatmap panels
 make_heatmap <- function(data, genera, title) {
 data_sub <- data %>%
 filter(Genus %in% genera) %>%
 group_by(Genus, Feces, Treatment) %>%
 summarise(mean_CLR = mean(CLR, na.rm = TRUE), .groups = "drop") %>%
 mutate(Feces = factor(Feces, levels = mb_params$feces_order))

 if (nrow(data_sub) == 0) return(NULL)

 heat_mat <- data_sub %>%
 pivot_wider(names_from = c(Feces, Treatment), values_from = mean_CLR) %>%
 column_to_rownames("Genus") %>%
 as.matrix()

 # Set column names to display labels
 colnames(heat_mat) <- gsub("^(F[1-4])_", "\\1: ", colnames(heat_mat))
 # Replace F1-F4 with week labels in column names
 for (f in names(feces_labels)) {
 colnames(heat_mat) <- gsub(paste0("^", f, ":"), paste0(feces_labels[f], " - "), colnames(heat_mat))
 }

 p <- pheatmap(heat_mat,
 cluster_rows = TRUE, cluster_cols = FALSE,
 scale = "row",
 main = title,
 color = colorRampPalette(rev(brewer.pal(11, "RdBu")))(100),
 fontsize = 8,
 border_color = NA)
 return(p)
 }

 # 3-panel: All, CDHFD (F3F4), ND (F1F2)
 p7_all <- make_heatmap(clr_df, sig_genera, "All Timepoints")
 p7_cdhfd <- make_heatmap(
 clr_df %>% filter(Feces %in% c("F3", "F4")),
 sig_genera, "CDHFD (3wks + 7wks)")
 p7_nd <- make_heatmap(
 clr_df %>% filter(Feces %in% c("F1", "F2")),
 sig_genera, "ND (-1wks + 0wks)")

 # Save individual panels
 if (!is.null(p7_all)) {
 png(file.path(plots_dir, "plot7_heatmap_all.png"), width = 800, height = 600)
 print(p7_all)
 dev.off()
 }
 if (!is.null(p7_cdhfd)) {
 png(file.path(plots_dir, "plot7_heatmap_cdhfd.png"), width = 800, height = 600)
 print(p7_cdhfd)
 dev.off()
 }
 if (!is.null(p7_nd)) {
 png(file.path(plots_dir, "plot7_heatmap_nd.png"), width = 800, height = 600)
 print(p7_nd)
 dev.off()
 }
} else {
 cat(" No significant genera for heatmap\n")
}

# ============================================================
# PLOT 8: INDIVIDUAL BOXPLOTS — significant genera (CLR)
# ============================================================
cat("--- Plot 8: Individual Boxplots ---\n")

# FIX 4: Remove arbitrary sig_genera[1] limit; fix jitter width from 0.00 to 0.15
if (length(sig_genera) > 0) {
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
 mutate(Feces = factor(Feces, levels = mb_params$feces_order))

 boxplot_list <- list()

 for (g in sig_genera) {
 gen_df <- clr_df %>% filter(Genus == g)
 if (nrow(gen_df) == 0) next

 p_box <- ggplot(gen_df, aes(x = Feces, y = CLR, fill = Treatment)) +
 geom_boxplot(position = position_dodge(0.8), outlier.shape = NA, alpha = 0.7) +
 # FIX 4: width = 0.15 instead of 0.00 for visible jitter
 geom_point(position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.8),
 aes(shape = Sex), size = 1.5, alpha = 0.6) +
 scale_fill_manual(values = Treatment_colors) +
 scale_shape_manual(values = Sex_shape) +
 scale_x_discrete(labels = feces_labels) +
 labs(title = g, x = "Timepoint", y = "CLR Abundance",
 fill = "Treatment", shape = "Sex") +
 theme_bw() +
 theme(text = element_text(family = "Liberation Sans"),
 axis.text.x = element_text(angle = 45, hjust = 1),
 plot.title = element_text(face = "italic"))

 boxplot_list[[g]] <- p_box
 ggsave(file.path(plots_dir, paste0("plot8_boxplot_", g, ".png")), p_box,
 width = 6, height = 5, dpi = 300)
 }

 # Combined panel (up to 9 per page)
 if (length(boxplot_list) > 0) {
 n_per_page <- 9
 n_pages <- ceiling(length(boxplot_list) / n_per_page)
 for (page in 1:n_pages) {
 start_idx <- (page - 1) * n_per_page + 1
 end_idx <- min(page * n_per_page, length(boxplot_list))
 page_plots <- boxplot_list[start_idx:end_idx]
 p_combined <- do.call("grid.arrange", c(page_plots, ncol = 3))
 ggsave(file.path(plots_dir, paste0("plot8_boxplots_page", page, ".png")),
 p_combined, width = 15, height = 10, dpi = 300)
 }
 }
} else {
 cat(" No significant genera for individual boxplots\n")
}

# ============================================================
# PLOT 9: LACTOBACILLUS TARGETED PLOT
# ============================================================
cat("--- Plot 9: Lactobacillus ---\n")

lacto_anova <- read.csv(file.path(stats_dir, "lactobacillus_anova.csv"))
lacto_pwc <- read.csv(file.path(stats_dir, "lactobacillus_posthoc_pooled.csv"))
lacto_pwc_sex <- read.csv(file.path(stats_dir, "lactobacillus_posthoc_by_sex.csv"))

ps_lacto <- subset_taxa(ps_genus, Genus == "Lactobacillus")
if (ntaxa(ps_lacto) > 0) {
 lacto_counts <- as(otu_table(ps_lacto), "matrix")
 if (!taxa_are_rows(ps_lacto)) lacto_counts <- t(lacto_counts)
 lacto_total <- colSums(lacto_counts)
 lacto_clr <- log(lacto_total + 0.5) - mean(log(lacto_total + 0.5))

 lacto_df <- data.frame(
 SampleID = names(lacto_total),
 Lacto_counts = lacto_total,
 Lacto_CLR = lacto_clr
 ) %>%
 left_join(metadata, by = "SampleID") %>%
 mutate(Feces = factor(Feces, levels = mb_params$feces_order))

 lacto_pwc_plot <- lacto_pwc %>%
 mutate(Feces = factor(Feces, levels = mb_params$feces_order))

 # Pooled plot
 p9a <- ggplot(lacto_df, aes(x = Feces, y = Lacto_CLR, fill = Treatment)) +
 geom_boxplot(position = position_dodge(0.8), outlier.shape = NA, alpha = 0.7) +
 geom_point(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8),
 aes(shape = Sex), size = 2, alpha = 0.8) +
 geom_text(data = lacto_pwc_plot,
 aes(x = Feces, label = significance, y = max(lacto_df$Lacto_CLR, na.rm = TRUE) * 1.05),
 vjust = 0, hjust = 0.5, inherit.aes = FALSE, size = 5) +
 scale_fill_manual(values = Treatment_colors) +
 scale_shape_manual(values = Sex_shape) +
 scale_x_discrete(labels = feces_labels) +
 labs(title = "Lactobacillus (Pooled)", x = "Timepoint",
 y = "CLR Abundance", fill = "Treatment", shape = "Sex") +
 theme_bw() +
 theme(text = element_text(family = "Liberation Sans"),
 axis.text.x = element_text(angle = 45, hjust = 1),
 plot.title = element_text(face = "bold"))

 ggsave(file.path(plots_dir, "plot9_lactobacillus_pooled.png"), p9a,
 width = 6, height = 5, dpi = 300)

 # Sex-faceted plot
 lacto_pwc_sex_plot <- lacto_pwc_sex %>%
 mutate(Feces = factor(Feces, levels = mb_params$feces_order))

 p9b <- ggplot(lacto_df, aes(x = Feces, y = Lacto_CLR, fill = Treatment)) +
 geom_boxplot(position = position_dodge(0.8), outlier.shape = NA, alpha = 0.7) +
 geom_point(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8),
 size = 1.5, alpha = 0.6) +
 geom_text(data = lacto_pwc_sex_plot,
 aes(x = Feces, label = significance, y = max(lacto_df$Lacto_CLR, na.rm = TRUE) * 1.05),
 vjust = 0, hjust = 0.5, inherit.aes = FALSE, size = 4) +
 facet_wrap(~Sex, labeller = as_labeller(Sex_colors)) +
 scale_fill_manual(values = Treatment_colors) +
 scale_x_discrete(labels = feces_labels) +
 labs(title = "Lactobacillus (by Sex)", x = "Timepoint",
 y = "CLR Abundance", fill = "Treatment") +
 theme_bw() +
 theme(text = element_text(family = "Liberation Sans"),
 axis.text.x = element_text(angle = 45, hjust = 1),
 plot.title = element_text(face = "bold"))

 ggsave(file.path(plots_dir, "plot9_lactobacillus_sex.png"), p9b,
 width = 8, height = 5, dpi = 300)

 # Relative abundance plot
 lacto_rel <- lacto_df %>%
 mutate(Lacto_rel = Lacto_counts / sample_sums(ps)[SampleID])

 p9c <- ggplot(lacto_rel, aes(x = Feces, y = Lacto_rel, fill = Treatment)) +
 geom_boxplot(position = position_dodge(0.8), outlier.shape = NA, alpha = 0.7) +
 geom_point(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8),
 aes(shape = Sex), size = 2, alpha = 0.8) +
 scale_fill_manual(values = Treatment_colors) +
 scale_shape_manual(values = Sex_shape) +
 scale_x_discrete(labels = feces_labels) +
 scale_y_continuous(labels = scales::percent_format()) +
 labs(title = "Lactobacillus Relative Abundance", x = "Timepoint",
 y = "Relative Abundance", fill = "Treatment", shape = "Sex") +
 theme_bw() +
 theme(text = element_text(family = "Liberation Sans"),
 axis.text.x = element_text(angle = 45, hjust = 1),
 plot.title = element_text(face = "bold"))

 ggsave(file.path(plots_dir, "plot9_lactobacillus_relabun.png"), p9c,
 width = 6, height = 5, dpi = 300)

} else {
 cat(" No Lactobacillus taxa found\n")
}

cat("\n=== Plots complete ===\n")
cat("Figures saved to:", plots_dir, "\n")
