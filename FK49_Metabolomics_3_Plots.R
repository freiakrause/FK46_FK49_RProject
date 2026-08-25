###############################################################################
# FK49_Metabolomics_3_Plots.R
#
# Loads preprocessed RDS files + statistics CSV and generates all plots:
#   1. PCA (autoplot, Treatment colour, Sex shape)
#   2. Volcano plot
#   3. Per-metabolite dot plots (top N by adj_p_Treatment)
#   4. Conditional sex-stratified plots (when adj_p_Treatment_Sex < 0.05)
#   5. Sex-specific dot plots (when adj_p_Treatment_female or male < 0.05)
#   6. Z-score scaled heatmap with censoring visualisation
#   7. Effect size plots (Cohen's d for LM, GMR for cen2way)
#
# Set `analysis` at top to "CDHFD" or "ND" (must match statistics script).
###############################################################################

rm(list = ls())
gc()

library(tidyverse)
library(pheatmap)
library(ggfortify)
library(factoextra)
library(ggrepel)
library(ggnewscale)
library(patchwork)
source("FK49_Definitions.R")

# ============================================================
# CONFIGURATION
# ============================================================

analysis <- "CDHFD"
top_n    <- 12

rawdata_pwd <- PATHS$metabolomics$rawdata
output_pwd  <- PATHS$metabolomics$output

analysis_folder <- if (analysis == "CDHFD") "CDHFD" else "ND"
analysis_pwd <- file.path(output_pwd, analysis_folder)

has_sex <- (analysis == "CDHFD")

# Analysis-specific filter
if (analysis == "CDHFD") {
  diet_filter <- "CDHFD13"; expid_filter <- "FK49"
} else {
  diet_filter <- "ND";      expid_filter <- "BH"
}

# ============================================================
# LOAD DATA
# ============================================================

pos <- readRDS(file.path(rawdata_pwd, "FK49_metabolome_positive_processed.rds"))
neg <- readRDS(file.path(rawdata_pwd, "FK49_metabolome_negative_processed.rds"))
tar <- readRDS(file.path(rawdata_pwd, "FK49_metabolome_targeted_processed.rds"))

datasets <- list(positive = pos, negative = neg, targeted = tar)

stats <- read.csv2(file.path(analysis_pwd, "FK49_metabolome_statistics.csv"),
                   stringsAsFactors = FALSE)

# ============================================================
# HELPER: Look up abbreviated metabolite names from RDS
# ============================================================
get_abbrev <- function(metab_names, abbrev_map) {
  if (is.null(abbrev_map)) return(metab_names)
  abbrev <- abbrev_map[metab_names]
  abbrev[is.na(abbrev)] <- metab_names[is.na(abbrev)]
  unname(abbrev)
}

# ============================================================
# HELPER: Filter RDS to analysis subset
# ============================================================

filter_subset <- function(rds, diet, expid) {
  idx <- rds$metadata$Diet == diet & rds$metadata$ExpID == expid
  list(
    metadata   = rds$metadata[idx, , drop = FALSE],
    raw_values = rds$raw_values[idx, , drop = FALSE],
    log_values = rds$log_values[idx, , drop = FALSE],
    norm_values = if (is.null(rds$norm_values)) NULL else rds$norm_values[idx, , drop = FALSE],
    censored   = rds$censored[idx, , drop = FALSE]
  )
}

# ============================================================
# HELPER: Per-metabolite dot plot (adapted from BA do_BA)
# ============================================================

do_metab_dotplot <- function(metadata, values_vec, cens_vec, metab_name,
                             sex = "both", p_value = NULL) {
  d <- data.frame(
    Treatment = factor(metadata$Treatment, levels = c("Ctrl", "TAM")),
    Sex       = factor(metadata$Sex,       levels = c("female", "male")),
    value     = values_vec,
    Censoring = factor(cens_vec, levels = c(FALSE, TRUE),
                       labels = c("Not censored", "Censored"))
  )

  if (sex != "both") d <- d %>% filter(Sex == sex)

  y_vals  <- d$value
  y_max   <- max(y_vals, na.rm = TRUE)
  y_min   <- min(y_vals, na.rm = TRUE)
  y_range <- y_max - y_min
  if (!is.finite(y_range) || y_range == 0) y_range <- max(abs(y_max) * 0.1, 1)
  y_pos   <- y_max + 0.15 * y_range

  p1 <- ggplot(d, aes(x = Treatment, y = value)) +
    stat_summary(fun = mean, geom = "bar", aes(fill = Treatment, color = Treatment),
                 alpha = 0.5, width = 0.75) +
    scale_color_manual(values = Treatment_colors[c("Ctrl", "TAM")]) +
    guides(color = "none") +
    ggnewscale::new_scale_color() +
    stat_summary(fun.data = mean_sdl, fun.args = list(mult = 1),
                 geom = "errorbar", width = 0.2, color = "black") +
    scale_fill_manual(name = "Treatment",
                      values = Treatment_colors[c("Ctrl", "TAM")]) +
    ggnewscale::new_scale_fill() +
    ggnewscale::new_scale_color() +
    geom_point(aes(shape = Sex, color = Censoring), fill = "lightgrey",
               alpha = 0.8, size = 5.3,
               position = position_jitter(width = 0.15, height = 0), stroke = 1.8) +
    scale_shape_manual(name = "Sex", values = Sex_shape) +
    scale_color_manual(name = "Censoring",
                       values = c("Not censored" = "black", "Censored" = "blue")) +
    scale_y_continuous(name = metab_name, expand = expansion(mult = c(0.05, 0.20))) +
    theme_minimal() +
    theme(legend.position = "bottom",
          legend.box = "vertical", legend.box.just = "left",
          legend.title = element_text(size = 11, face = "bold"),
          legend.text  = element_text(size = 9),
          axis.line  = element_line(color = "black", linewidth = 0.5),
          axis.ticks = element_line(color = "black", linewidth = 0.5),
          axis.title = element_text(size = 20, face = "bold"),
          axis.title.x = element_blank(),
          axis.text  = element_text(size = 19, face = "bold"),
          panel.grid = element_blank()) +
    guides(shape = guide_legend(title = "Sex", order = 2, nrow = 1, byrow = TRUE),
           color = guide_legend(title = "Censoring", order = 2, nrow = 1, byrow = TRUE),
           fill  = guide_legend(title = "Group", order = 3, nrow = 1, byrow = TRUE))

  if (!is.null(p_value) && is.finite(p_value))
    p1 <- p1 + annotate("text", x = 1.5, y = y_pos,
                        label = paste0("adj p = ", format.pval(p_value, digits = 3)),
                        size = 5.5, fontface = "italic")

  list(plot = p1, y_pos = y_pos)
}

# ============================================================
# PLOT LOOP: per dataset
# ============================================================

for (ds_name in names(datasets)) {
  cat("\n\n========== Plotting:", ds_name, "==========\n")

  ds <- filter_subset(datasets[[ds_name]], diet_filter, expid_filter)
  ds_stats <- stats %>% filter(Dataset == ds_name)

  out_dir <- file.path(analysis_pwd, ds_name)
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

  # Values for dot plots: norm_values for untargeted, raw_values for targeted
  is_untargeted <- !is.null(ds$norm_values)
  dot_values <- if (is_untargeted) ds$norm_values else ds$raw_values
  y_unit <- if (is_untargeted) "" else " (pmol/mg)"

  metab_names <- colnames(ds$log_values)

  # ----------------------------------------------------------
  # 1. PCA
  # ----------------------------------------------------------
  cat("  PCA...\n")

  pca_mat <- ds$log_values
  # Autoscale within subset
  pca_mat_scaled <- scale(pca_mat)

  # Drop samples with NA (prcomp cannot handle NA)
  complete_rows <- complete.cases(pca_mat_scaled)
  if (sum(complete_rows) >= 4) {
    pca_data <- pca_mat_scaled[complete_rows, , drop = FALSE]
    pca_meta <- ds$metadata[complete_rows, , drop = FALSE]

    data.pca <- prcomp(pca_data)

    p_pca <- autoplot(data.pca, data = pca_meta, x = 1, y = 2,
                      size = 3, fill = "Treatment", color = "Treatment",
                      shape = if (has_sex) "Sex" else NULL,
                      frame = TRUE, frame.type = "norm", frame.level = 0.95) +
      theme_bw() + theme_classic() +
      ggtitle(paste0("PCA - ", ds_name, " (", analysis, ")")) +
      scale_color_manual(values = Treatment_colors[c("Ctrl", "TAM")]) +
      scale_fill_manual(values  = Treatment_colors[c("Ctrl", "TAM")]) +
      {if (has_sex) scale_shape_manual(values = Sex_shape)} +
      theme(text = element_text(size = 20))

    ggsave(p_pca, filename = "PCA.png",
           path = out_dir, dpi = 300, width = 9, height = 7)

    # PCA diagnostics
    p_scree <- fviz_eig(data.pca, addlabels = TRUE)
    ggsave(p_scree, filename = "PCA_scree.png",
           path = out_dir, dpi = 300, width = 4, height = 2.5)

    p_var <- fviz_pca_var(data.pca, col.var = "cos2",
                          gradient.cols = c("black", "orange", "green"),
                          repel = TRUE)
    ggsave(p_var, filename = "PCA_variables.png",
           path = out_dir, dpi = 300, width = 5, height = 4)

    # Top 10 PC1/PC2 contributors
    loadings <- data.pca$rotation[, 1:2]
    top_pc1 <- names(sort(abs(loadings[, 1]), decreasing = TRUE))[1:min(10, nrow(loadings))]
    top_pc2 <- names(sort(abs(loadings[, 2]), decreasing = TRUE))[1:min(10, nrow(loadings))]
    top_contrib <- data.frame(PC1 = top_pc1, PC2 = top_pc2)
    write.csv(top_contrib,
              file = file.path(out_dir, "PCA_top10_loadings.csv"),
              row.names = FALSE)
  }

  # ----------------------------------------------------------
  # 2. Volcano plot
  # ----------------------------------------------------------
  cat("  Volcano...\n")

  volcano_data <- ds_stats %>%
    filter(is.finite(log2FC), is.finite(adj_p_Treatment)) %>%
    mutate(
      negLog10FDR = -log10(adj_p_Treatment),
      direction = case_when(
        adj_p_Treatment < 0.05 & log2FC > 0 ~ "TAM higher",
        adj_p_Treatment < 0.05 & log2FC < 0 ~ "TAM lower",
        TRUE ~ "Not significant"
      ),
      abbrev = get_abbrev(Metabolite, ds$metabolite_abbrev)
    )

  p_volcano <- ggplot(volcano_data, aes(x = log2FC, y = negLog10FDR)) +
    geom_point(aes(fill = direction), alpha = 0.5, size = 3, stroke = 0.5,
               shape = 21, color = "black",
               position = position_jitter(width = 0.08)) +
    scale_fill_manual(values = c("TAM lower" = "blue",
                                 "Not significant" = "grey60",
                                 "TAM higher" = "firebrick")) +
    geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed", color = "grey80") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey80") +
    labs(title = paste0("Volcano - ", ds_name, " (", analysis, ")"),
         x = expression(paste("log"[2], " FC (TAM / Ctrl)")),
         y = expression(paste("-log"[10], "(adj.p.value)"))) +
    theme_classic() +
    theme(panel.grid = element_line(color = "grey90", linewidth = 0.1)) +
    geom_text_repel(data = volcano_data %>% filter(adj_p_Treatment < 0.05),
                    aes(label = abbrev), size = 3.5, max.overlaps = 25)

  ggsave(p_volcano, filename = "volcano.png",
         path = out_dir, dpi = 300, width = 6, height = 9)

  # ----------------------------------------------------------
  # 3. Per-metabolite dot plots (top N)
  # ----------------------------------------------------------
  cat("  Dot plots (top", top_n, ")...\n")

  top_metabs <- ds_stats %>%
    filter(!is.na(adj_p_Treatment)) %>%
    arrange(adj_p_Treatment, desc(abs(log2FC))) %>%
    head(top_n) %>%
    pull(Metabolite)

  for (metab in top_metabs) {
    if (!(metab %in% colnames(dot_values))) next
    stat_row <- ds_stats %>% filter(Metabolite == metab)
    p_adj <- stat_row$adj_p_Treatment[1]

    res <- do_metab_dotplot(ds$metadata, dot_values[, metab],
                            ds$censored[, metab],
                            paste0(get_abbrev(metab, ds$metabolite_abbrev), y_unit),
                            sex = "both", p_value = p_adj)
    fname <- substr(gsub("[^[:alnum:]_]", "_", get_abbrev(metab, ds$metabolite_abbrev)), 1, 30)
    ggsave(res$plot, filename = paste0(fname, ".png"),
           path = out_dir, width = 4, height = 11, dpi = 300)
  }

  # ----------------------------------------------------------
  # 4. Conditional sex-stratified plots (CDHFD only)
  # ----------------------------------------------------------
  if (has_sex) {
    cat("  Sex-stratified plots...\n")

    for (metab in metab_names) {
      stat_row <- ds_stats %>% filter(Metabolite == metab)
      if (nrow(stat_row) == 0) next

      p_adj_Treat_Sex <- stat_row$adj_p_Treatment_Sex[1]
      p_adj_female    <- stat_row$adj_p_Treatment_female[1]
      p_adj_male      <- stat_row$adj_p_Treatment_male[1]

      # Interaction significant: facet by Sex
      if (!is.na(p_adj_Treat_Sex) && p_adj_Treat_Sex < 0.05 &&
          metab %in% colnames(dot_values)) {

        cens_col <- paste0(metab, "_censored")
        d_sex <- data.frame(
          Treatment = factor(ds$metadata$Treatment, levels = c("Ctrl", "TAM")),
          Sex       = factor(ds$metadata$Sex,       levels = c("female", "male")),
          value     = dot_values[, metab],
          Censoring = factor(ds$censored[, metab], levels = c(FALSE, TRUE),
                             labels = c("Not censored", "Censored"))
        )

        y_max <- max(d_sex$value, na.rm = TRUE)
        y_min <- min(d_sex$value, na.rm = TRUE)
        y_range <- max(y_max - y_min, abs(y_max) * 0.1, 1)

        p_sex <- ggplot(d_sex, aes(x = Treatment, y = value)) +
          stat_summary(fun = mean, geom = "bar", aes(fill = Treatment),
                       color = "black", alpha = 0.5, width = 0.75) +
          stat_summary(fun.data = mean_sdl, fun.args = list(mult = 1),
                       geom = "errorbar", width = 0.2, color = "black") +
          facet_wrap(~Sex) +
          scale_fill_manual(name = "Treatment",
                            values = Treatment_colors[c("Ctrl", "TAM")]) +
          ggnewscale::new_scale_fill() +
          ggnewscale::new_scale_color() +
          geom_point(aes(shape = Sex, color = Censoring), alpha = 0.8,
                     size = 4.5, position = position_jitter(width = 0.15),
                     stroke = 1.5) +
          scale_color_manual(name = "Censoring",
                             values = c("Not censored" = "black", "Censored" = "blue")) +
          scale_shape_manual(name = "Sex", values = Sex_shape) +
          scale_y_continuous(name = paste0(get_abbrev(metab, ds$metabolite_abbrev), y_unit), expand = expansion(mult = c(0.05, 0.20))) +
          theme_minimal() +
          theme(legend.position = "bottom",
                axis.line  = element_line(color = "black", linewidth = 0.5),
                axis.ticks = element_line(color = "black", linewidth = 0.5),
                axis.title = element_text(size = 20, face = "bold"),
                axis.title.x = element_blank(),
                axis.text   = element_text(size = 17, face = "bold"),
                strip.text  = element_text(size = 15, face = "bold"),
                panel.grid  = element_blank())

        annot_df <- data.frame(
          Sex   = factor(c("female", "male"), levels = c("female", "male")),
          x     = 1.5,
          y     = y_max + 0.15 * (y_max-y_min),
          label = c(paste0("adj p = ", format.pval(p_adj_female, digits = 3)),
                    paste0("adj p = ", format.pval(p_adj_male, digits = 3)))
        )
        p_sex <- p_sex +
          geom_text(data = annot_df, aes(x = x, y = y, label = label),
                    size = 4.5, fontface = "italic", inherit.aes = FALSE)

        fname <- substr(gsub("[^[:alnum:]_]", "_", get_abbrev(metab, ds$metabolite_abbrev)), 1, 30)
        ggsave(p_sex, filename = paste0(fname, "_bySex.png"),
               path = out_dir, width = 8, height = 11, dpi = 300)
      }

      # Sex-specific significant: individual plots
      sex_plots <- list()
      if (!is.na(p_adj_female) && p_adj_female < 0.05 &&
          metab %in% colnames(dot_values)) {
        res_f <- do_metab_dotplot(ds$metadata, dot_values[, metab],
                                  ds$censored[, metab],
                                  paste0(get_abbrev(metab, ds$metabolite_abbrev), y_unit),
                                  sex = "female", p_value = p_adj_female)
        sex_plots$female <- res_f$plot + ggtitle("Female")
      }
      if (!is.na(p_adj_male) && p_adj_male < 0.05 &&
          metab %in% colnames(dot_values)) {
        res_m <- do_metab_dotplot(ds$metadata, dot_values[, metab],
                                  ds$censored[, metab],
                                  paste0(get_abbrev(metab, ds$metabolite_abbrev), y_unit),
                                  sex = "male", p_value = p_adj_male)
        sex_plots$male <- res_m$plot + ggtitle("Male")
      }
      if (length(sex_plots) > 0) {
        p_combined <- wrap_plots(sex_plots)
        fname <- substr(gsub("[^[:alnum:]_]", "_", get_abbrev(metab, ds$metabolite_abbrev)), 1, 30)
        ggsave(p_combined, filename = paste0(fname, "_sexSpecific.png"),
               path = out_dir, width = 4 * length(sex_plots), height = 11, dpi = 300,
               limitsize = FALSE)
      }
    }
  }

  # ----------------------------------------------------------
  # 5. Z-score scaled heatmap with censoring
  # ----------------------------------------------------------
  cat("  Heatmap...\n")

  # Select metabolites: significant + trend, fallback to top N
  heat_metabs <- ds_stats %>%
    filter(significant == TRUE | trend == TRUE) %>%
    pull(Metabolite)

  if (length(heat_metabs) < 3) {
    heat_metabs <- top_metabs
  }
  heat_metabs <- intersect(heat_metabs, colnames(ds$log_values))

  if (length(heat_metabs) >= 2) {
    d_mat <- ds$log_values[, heat_metabs, drop = FALSE]
    rownames(d_mat) <- ds$metadata$Animal
    d_mat <- t(d_mat)  # metabolites as rows, samples as columns

    # Column annotation
    ann <- data.frame(
      Treatment = factor(ds$metadata$Treatment, levels = c("Ctrl", "TAM")),
      Sex       = factor(ds$metadata$Sex,       levels = c("female", "male"))
    )
    rownames(ann) <- ds$metadata$Animal
    ann <- ann[colnames(d_mat), , drop = FALSE]

    # Z-score per metabolite (row)
    d_scaled <- t(scale(t(d_mat)))

    # Sort rows by effect_size_type then adj_p_Treatment
    row_order <- ds_stats %>%
      filter(Metabolite %in% heat_metabs) %>%
      mutate(effect_size_type = ifelse(is.na(effect_size_type), "zzz", effect_size_type)) %>%
      arrange(effect_size_type, adj_p_Treatment) %>%
      pull(Metabolite)
    row_order <- intersect(row_order, rownames(d_scaled))

    # Sort columns by Treatment, Sex
    col_order <- order(ann$Treatment, ann$Sex)

    d_scaled <- d_scaled[row_order, col_order, drop = FALSE]
    ann      <- ann[col_order, , drop = FALSE]

    # Censoring overlay
    cens_mat <- t(ds$censored[, heat_metabs, drop = FALSE])
    rownames(cens_mat) <- heat_metabs
    colnames(cens_mat) <- ds$metadata$Animal
    cens_mat <- cens_mat[row_order, col_order, drop = FALSE]

    below_detection <- cens_mat == TRUE
    d_scaled[below_detection] <- -4

    # Row annotation: adj p categories
    row_annot <- ds_stats %>%
      filter(Metabolite %in% heat_metabs) %>%
      select(Metabolite, adj_p_Treatment) %>%
      mutate(adj.p = cut(adj_p_Treatment,
                         breaks = c(-Inf, 0.0001, 0.001, 0.01, 0.05, Inf),
                         labels = c("<0.0001", "0.0001-0.001", "0.001-0.01",
                                    "0.01-0.05", ">0.05"),
                         include.lowest = TRUE)) %>%
      select(Metabolite, adj.p) %>%
      tibble::column_to_rownames("Metabolite")
    row_annot$adj.p <- as.character(row_annot$adj.p)
    row_annot$adj.p[is.na(row_annot$adj.p)] <- "NA"
    row_annot$adj.p <- factor(row_annot$adj.p,
                              levels = c("<0.0001", "0.0001-0.001", "0.001-0.01",
                                         "0.01-0.05", ">0.05", "NA"))
    row_annot <- row_annot[row_order, , drop = FALSE]

    heatmap_colors <- c("#FFF5CC",
                        colorRampPalette(c("#FFE699", "orange", "red"))(98),
                        "darkred")
    heatmap_breaks <- seq(-4, 4, length.out = length(heatmap_colors) + 1)

    p_color_list <- c(
      "<0.0001" = "darkgreen", "0.0001-0.001" = "forestgreen",
      "0.001-0.01" = "green3", "0.01-0.05" = "green1",
      ">0.05" = "white", "NA" = "grey80"
    )

    ann_colors_list <- list(
      Treatment = Treatment_colors[c("Ctrl", "TAM")],
      Sex       = Sex_colors,
      adj.p     = p_color_list
    )

    # Ordered heatmap (no clustering)
    p_heat <- pheatmap(
      d_scaled, cluster_rows = FALSE, cluster_cols = FALSE,
      fontsize = 10, fontsize_row = 7, fontsize_col = 8,
      annotation_col = ann, annotation_row = row_annot,
      color = heatmap_colors, breaks = heatmap_breaks,
      annotation_colors = ann_colors_list,
      border_color = "black", cellwidth = 10, cellheight = 7,
      angle_col = 90, show_colnames = FALSE,
      labels_row = get_abbrev(rownames(d_scaled), ds$metabolite_abbrev),
      legend_breaks = c(-4, -2, 0, 2, 4),
      legend_labels = c("Below LOD", "-2", "0", "2", "Above ULOQ")
    )
    ggsave(p_heat, filename = "heatmap_ordered.png",
           path = out_dir, width = 12, height = 5, dpi = 300, bg = "white",
           limitsize = FALSE)

    # Clustered heatmap
    p_heat_clust <- pheatmap(
      d_scaled, cluster_rows = TRUE, cluster_cols = TRUE,
      fontsize = 10, fontsize_row = 7, fontsize_col = 8,
      annotation_col = ann, annotation_row = row_annot,
      color = heatmap_colors, breaks = heatmap_breaks,
      annotation_colors = ann_colors_list,
      border_color = "black", cellwidth = 10, cellheight = 7,
      angle_col = 90, show_colnames = FALSE,
      labels_row = get_abbrev(rownames(d_scaled), ds$metabolite_abbrev),
      legend_breaks = c(-4, -2, 0, 2, 4),
      legend_labels = c("Below LOD", "-2", "0", "2", "Above ULOQ")
    )
    ggsave(p_heat_clust, filename = "heatmap_clustered.png",
           path = out_dir, width = 12, height = 5, dpi = 300, bg = "white",
           limitsize = FALSE)
  }

  # ----------------------------------------------------------
  # 6. Effect size plots
  # ----------------------------------------------------------
  cat("  Effect size plots...\n")
  effect_metabs <- ds_stats %>%
    filter(significant == TRUE | trend == TRUE) %>%
    pull(Metabolite)
  stats_eff <- ds_stats %>%
    filter(!is.na(effect_size), is.finite(effect_size),Metabolite %in% effect_metabs)

  # Cohen's d (LM metabolites)
  stats_lm <- stats_eff %>% filter(effect_size_type == "Cohen's d")
  if (nrow(stats_lm) > 0) {
    stats_lm <- stats_lm %>%
      mutate(abbrev = get_abbrev(Metabolite, ds$metabolite_abbrev),
             abbrev = factor(abbrev, levels = unique(abbrev[order(effect_size)])))
    p_lm <- ggplot(stats_lm, aes(x = effect_size, y = abbrev)) +
      geom_vline(xintercept = 0, linetype = "dashed") +
      geom_errorbar(aes(xmin = effect_CI_low, xmax = effect_CI_high), height = 0) +
      geom_point(size = 2) +
      labs(x = "Cohen's d (TAM - Ctrl)", y = NULL) +
      theme_classic()
    ggsave(p_lm, filename = "effectsize_LM.png",
           path = out_dir, dpi = 300, width = 5,
           height = max(3, 0.25 * nrow(stats_lm) + 1),
           limitsize = FALSE)
  }

  # GMR (censored metabolites)
  stats_gmr <- stats_eff %>% filter(effect_size_type == "GMR")
  if (nrow(stats_gmr) > 0) {
    stats_gmr <- stats_gmr %>%
      mutate(abbrev = get_abbrev(Metabolite, ds$metabolite_abbrev),
             abbrev = factor(abbrev, levels = unique(abbrev[order(effect_size)])))
    p_gmr <- ggplot(stats_gmr, aes(x = effect_size, y = abbrev)) +
      geom_vline(xintercept = 1, linetype = "dashed") +
      geom_errorbar(aes(xmin = effect_CI_low, xmax = effect_CI_high), height = 0) +
      geom_point(size = 2) +
      scale_x_log10() +
      labs(x = "Geometric Mean Ratio (TAM / Ctrl)", y = NULL) +
      theme_classic()
    ggsave(p_gmr, filename = "effectsize_GMR.png",
           path = out_dir, dpi = 300, width = 5,
           height = max(3, 0.25 * nrow(stats_gmr) + 1),
           limitsize = FALSE)
  }
}

cat("\n=== Plotting complete (", analysis, ") ===\n")
cat("Plots saved to:", analysis_pwd, "\n")
