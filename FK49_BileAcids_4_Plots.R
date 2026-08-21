rm(list=ls())
gc()

library(tidyverse)
library(pheatmap)
library(ggfortify)
library(factoextra)
library(patchwork)
library(ggnewscale)
source("FK49_Definitions.R")
### hier später weiter achen. statisik wurde auf censoring angepasst.
### 
###  bei exigo prüfen ob ich statt cen2 means cen2way implementiren kann und sollte.
###  volcano plot prüfen, warum censored werte nicht dargestellt werden.
###  bei exigo prüfen ob p werte richtig gezogen, ob volcano da ist und ob censored value im volcano erscheinen
###  
###  # Read Raw Inputdata ------------------------------------------------------

ExpID <- "FK49"
output_pwd <- PATHS$BA$output
#d1 <- readRDS(file = paste0(dirname(dirname(output_pwd)), "/01_RawData/FK49_BAfiltered_preprocessed.rds"))
d1 <-  readRDS(file = paste0(dirname(dirname(output_pwd)), "/01_RawData/FK49_BA_preprocessed.rds"))%>%
filter(Timepoint == "11") 
param_list <- PARAMETERS$BA$BA_sort
stats <- read.csv2(file.path(output_pwd, paste0(ExpID, "_BA_TP11_Statistics.csv")))
head(d1)
d1$GCDCA

# Dotplot Bile Acids ------------------------------------------------------

do_BA <- function(inputdata, value, sex = "both", y_title, path_images, p_value = NULL) {
  censored_col <- paste0(value, "_censored")

  d <- inputdata %>%
    filter(Timepoint == "11", !is.na(.data[[value]])) %>%
    mutate(Treatment = factor(Treatment, levels = c("Ctrl", "TAM")),
           Sex = factor(Sex, levels = c("female", "male")),
           Censoring = factor( .data[[censored_col]],levels = c(FALSE, TRUE),labels = c("Not censored", "Censored")))
  print(value)
  print(summary(d[[value]]))
  print(sort(unique(d[[value]])))
  if (sex != "both") d <- d %>% filter(Sex == sex)
  
  y_vals <- d[[value]]
  y_max <- max(y_vals, na.rm = TRUE)
  y_min <- min(y_vals, na.rm = TRUE)
  y_range <- y_max - y_min
  if (!is.finite(y_range) || y_range == 0) y_range <- max(abs(y_max) * 0.1, 1)
  y_pos <- y_max + 0.15 * y_range
  
  p1 <- ggplot(d, aes(x = Treatment, y = .data[[value]])) +
    stat_summary(fun = mean, geom = "bar", aes(fill = Treatment, color = Treatment), alpha = 0.5, width = 0.75) +
    scale_color_manual(values = Treatment_colors[c("Ctrl","TAM")], labels = c("Ctrl" = "Ctrl", "TAM" = "TAM")) +
    guides(color = "none") +
    ggnewscale::new_scale_color() +
    stat_summary(fun.data = mean_sdl, fun.args = list(mult = 1), geom = "errorbar", width = 0.2, color = "black") +
    scale_fill_manual(name = "Treatment", values = Treatment_colors[c("Ctrl","TAM")], labels = c("Ctrl" = "Ctrl", "TAM" = "TAM")) +
    ggnewscale::new_scale_fill() +
    ggnewscale::new_scale_color() +
    geom_point(aes(shape = Sex,color =Censoring ), fill = "lightgrey", alpha = 0.8, size = 5.3, position = position_jitter(width = 0.15, height = 0),stroke = 1.8) + #position = position_jitter(width = 0.15), 
    scale_shape_manual(name = "Sex", values = Sex_shape) +
    scale_color_manual(name = "Censoring", values = c("Not censored" = "black","Censored"="blue")) +
    scale_y_continuous(name = y_title, expand = expansion(mult = c(0.05, 0.20))) +
    theme_minimal() +
    theme(legend.position = "bottom",
          legend.box = "vertical",
          legend.box.just = "left",
          legend.title = element_text(size = 11, face = "bold"),
          legend.text = element_text(size = 9),
          axis.line = element_line(color = "black", linewidth = 0.5),
          axis.ticks = element_line(color = "black", linewidth = 0.5),
          axis.title = element_text(size = 20, face = "bold"),
          axis.title.x = element_blank(),
          axis.text = element_text(size = 19, face = "bold"),
          panel.grid = element_blank()) +
    guides(shape = guide_legend(title = "Sex", order = 2, nrow = 1, byrow = TRUE),
           color = guide_legend(title = "Censoring", order = 2, nrow = 1, byrow = TRUE),
           fill = guide_legend(title = "Group", order = 3, nrow = 1, byrow = TRUE))
  
  if (!is.null(p_value) && is.finite(p_value))
    p1 <- p1 + annotate("text", x = 1.5, y = y_pos, label = paste0("adj p = ", format.pval(p_value, digits = 3)), size = 5.5, fontface = "italic")
  
  return(list(plot_raw = p1, y_pos = y_pos))
}


# Generate Bile Acid Treatment Plots --------------------------------------

plots <- lapply(param_list, function(p) {
  stat_row <- stats %>% filter(BA == p)
  if (nrow(stat_row) == 0) return(NULL)
  
  p_adj_treatment <- stat_row$adj_p_Treatment[1]
  p_adj_sex <- stat_row$adj_p_Treatment_Sex[1]
  p_adj_female <- stat_row$adj_p_Treatment_female[1]
  p_adj_male <- stat_row$adj_p_Treatment_male[1]
  
  # ------------------------------------------------------------
  # Overall Treatment effect
  # ------------------------------------------------------------
  
  res <- do_BA(d1, p, sex = "both", y_title = p, path_images = output_pwd)
  p_overall <- res$plot_raw + annotate("text", x = 1.5, y = res$y_pos, label = paste0("adj p = ", format.pval(p_adj_treatment, digits = 3)), size = 5.5, fontface = "italic")
  
  fname_val <- gsub("[^[:alnum:]_]", "_", p)
  ggsave(filename = paste0(ExpID, "_BA_", fname_val, "_Treatment.png"), plot = p_overall, path = output_pwd, width = 4, height = 11, dpi = 300)
  
  # ------------------------------------------------------------
  # Treatment effect by sex if Treatment x Sex is significant
  # ------------------------------------------------------------
  
  if (!is.na(p_adj_sex) && p_adj_sex < 0.05) {
    d_sex <- d1 %>%
      filter(Timepoint == "11", !is.na(.data[[p]])) %>%
      mutate(Treatment = factor(Treatment, levels = c("Ctrl", "TAM")),
             Sex = factor(Sex, levels = c("female", "male")), 
             Censoring = factor(.data[[censored_col]],levels = c(FALSE, TRUE),labels = c("Not censored", "Censored")))
    
    
    p_sex <- ggplot(d_sex, aes(x = Treatment, y = .data[[p]])) +
      stat_summary(fun = mean, geom = "bar", aes(fill = Treatment), color="black",alpha = 0.5, width = 0.75) +
      stat_summary(fun.data = mean_sdl, fun.args = list(mult = 1), geom = "errorbar", width = 0.2, color = "black") +
      facet_wrap(~Sex) +
      scale_fill_manual(name = "Treatment", values = Treatment_colors[c("Ctrl","TAM")], labels = c("Ctrl" = "Ctrl", "TAM" = "TAM")) +
      ggnewscale::new_scale_fill() +
      ggnewscale::new_scale_color() +
      geom_point(aes(shape = Sex, fill = Sex,color = .data[[censored_col]]), alpha = 0.8, size = 4.5,position = position_jitter(width = 0.15, height = 0), stroke = 1.5) + # position = position_jitter(width = 0.15),
      scale_color_manual(name = "Censoring", values = c("Not censored"="black","Censored"="blue")) +
      scale_shape_manual(name = "Sex", values = Sex_shape) +
      scale_fill_manual(name = "Sex", values = Sex_colors) +
      scale_y_continuous(name = p, expand = expansion(mult = c(0.05, 0.20))) +
      theme_minimal() +
      theme(legend.position = "bottom",
            axis.line = element_line(color = "black", linewidth = 0.5),
            axis.ticks = element_line(color = "black", linewidth = 0.5),
            axis.title = element_text(size = 20, face = "bold"),
            axis.title.x = element_blank(),
            axis.text = element_text(size = 17, face = "bold"),
            strip.text = element_text(size = 15, face = "bold"),
            panel.grid = element_blank())
    
    y_max <- max(d_sex[[p]], na.rm = TRUE)
    y_min <- min(d_sex[[p]], na.rm = TRUE)
    y_range <- max(y_max - y_min, abs(y_max) * 0.1, 1)
    
    p_sex <- p_sex +
      annotate("text", x = 1.5, y = y_max + 0.15 * y_range, label = paste0("adj p = ", format.pval(p_adj_female, digits = 3)), size = 4.5, fontface = "italic", data = d_sex %>% filter(Sex == "female")) +
      annotate("text", x = 1.5, y = y_max + 0.15 * y_range, label = paste0("adj p = ", format.pval(p_adj_male, digits = 3)), size = 4.5, fontface = "italic", data = d_sex %>% filter(Sex == "male"))
    
    ggsave(filename = paste0(ExpID, "_BA_", fname_val, "_Treatment_by_Sex.png"), plot = p_sex, path = output_pwd, width = 8, height = 11, dpi = 300)
  }
  
  # ------------------------------------------------------------
  # Sex-specific Treatment effects
  # ------------------------------------------------------------
  
  sex_plots <- list()
  
  if (!is.na(p_adj_female) && p_adj_female < 0.05) {
    res_f <- do_BA(d1, p, sex = "female", y_title = p, path_images = output_pwd, p_value = p_adj_female)
    sex_plots$female <- res_f$plot_raw + ggtitle("Female")
  }
  
  if (!is.na(p_adj_male) && p_adj_male < 0.05) {
    res_m <- do_BA(d1, p, sex = "male", y_title = p, path_images = output_pwd, p_value = p_adj_male)
    sex_plots$male <- res_m$plot_raw + ggtitle("Male")
  }
  
  if (length(sex_plots) > 0) {
    p_sex_specific <- wrap_plots(sex_plots)
    ggsave(filename = paste0(ExpID, "_BA_", fname_val, "_Treatment_SexSpecific.png"), plot = p_sex_specific, path = output_pwd, width = 4, height = 11, dpi = 300)
  }
  
  p_overall
})


# Generate Matrix for Heatmaps --------------------------------------------
# Bile Acid heatmaps -------------------------------------------------------

d_mat <- d1 %>%
  filter(Timepoint == "11") %>%
  select(all_of(param_list)) %>%
  as.matrix()

rownames(d_mat) <- d1 %>% 
  filter(Timepoint == "11") %>% 
  pull(Animal)

d_mat <- t(d_mat)
ann <- d1 %>%
  filter(Timepoint == "11") %>%
  select(Treatment, Sex) %>%
  as.data.frame()

rownames(ann) <- d1 %>%
  filter(Timepoint == "11") %>%
  pull(Animal)

ann <- ann[colnames(d_mat), , drop = FALSE]

dim(d_mat)
head(d_mat[, 1:5])


# Correlation heatmap based on bile acids ---------------------------------
#d_mat is oriented so taht animals are columns and animal are going to be correlated against each other
cor_spearman <- cor(d_mat,method = "spearman",use = "pairwise.complete.obs")

round(cor_spearman, 2)
summary(cor_spearman[upper.tri(cor_spearman)])

p_cor <- pheatmap(
  cor_spearman,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  fontsize = 10,
  fontsize_main = 7,
  fontsize_row = 9,
  fontsize_col = 9,
  fontsize_number = 8,
  annotation_col = ann,
  color = colorRampPalette(c("blue", "white", "red"))(100),
  breaks = seq(-1, 1, length.out = 101),
  annotation_colors = list( Treatment = Treatment_colors[c("Ctrl", "TAM")],
    Sex = Sex_colors),
  border_color = "black",
  cellwidth = 10,
  cellheight = 10,
  angle_col = 90,
  display_numbers = FALSE,
  show_rownames = FALSE,
  show_colnames = FALSE,
  legend_breaks = c(-1, -0.5, 0, 0.5, 1),
  legend_labels = c("-1", "-0.5", "0", "0.5", "1"))

ggsave(
  filename = paste0(ExpID, "_BA_Correlation_Animals_ranks.png"),
  plot = p_cor,
  path = output_pwd,
  width = 5.5,
  height = 5,
  dpi = 300,
  bg = "white"
)

p_cor
# Correlation heatmap between bile acids -------------------------------
ba_sd <- apply(d_mat, 1, sd, na.rm = TRUE)

d_mat_cor <- d_mat[ba_sd > 0, , drop = FALSE] # bila acids that have 0 varianze * which are constant (due to omy censoring) are left out

cor_spearman <- cor(         #d_mat_cor is transposed so that bile acids are columns and bile acids are being correlated
  t(d_mat_cor),
  method = "spearman",
  use = "pairwise.complete.obs"
)

sum(is.na(cor_spearman))
p_cor <- pheatmap(
  cor_spearman,
  cluster_rows = F,
  cluster_cols = F,
  fontsize = 10,
  fontsize_main = 7,
  fontsize_row = 9,
  fontsize_col = 9,
  fontsize_number = 8,
  color = colorRampPalette(c("blue", "white", "red"))(100),
  breaks = seq(-1, 1, length.out = 101),
  border_color = "black",
  cellwidth = 10,
  cellheight = 10,
  angle_col = 90,
  display_numbers = FALSE,
  show_rownames = TRUE,
  show_colnames = TRUE
)

ggsave(
  filename = paste0(ExpID, "_BA_Correlation_Parameters_ranks.png"),
  plot = p_cor,
  path = output_pwd,
  width = 8,
  height = 5,
  dpi = 300,
  bg = "white"
)


# Correlation heatmap based on bile acids - Pearson ----------------------

d_mat_pearson <- d1 %>%
  filter(Timepoint == "11") %>%
  select(all_of(param_list)) %>%
  mutate(across(all_of(param_list), log10)) %>% # peasron needs values with linear dependency, hope taht log10 makes them linear
  as.matrix()

rownames(d_mat_pearson) <- d1 %>% 
  filter(Timepoint == "11") %>% 
  pull(Animal)

d_mat_pearson <- t(d_mat_pearson)

dim(d_mat_pearson)
head(d_mat_pearson[, 1:5])

summary(as.vector(as.matrix(d_mat_pearson)))

apply(d_mat_pearson, 2, sd, na.rm = TRUE)

apply(d_mat_pearson,  2,function(x) length(unique(na.omit(x))))

cor_pearson <- cor( d_mat_pearson,method = "pearson", use = "pairwise.complete.obs") #d_mat is oriented so taht animals are columns and animal are going to be correlated against each other


round(cor_pearson, 2)
summary(cor_pearson[upper.tri(cor_pearson)])

p_cor <- pheatmap(
  cor_pearson,
  cluster_rows = F,
  cluster_cols = F,
  fontsize = 10,
  fontsize_main = 7,
  fontsize_row = 9,
  fontsize_col = 9,
  fontsize_number = 8,
  color = colorRampPalette(c("blue", "white", "red"))(100),
  breaks = seq(-1, 1, length.out = 101),
  border_color = "black",
  cellwidth = 10,
  cellheight = 10,
  angle_col = 90,
  annotation_col = ann,
  annotation_colors = list( Treatment = Treatment_colors[c("Ctrl", "TAM")],
                            Sex = Sex_colors),
  display_numbers = FALSE,
  show_rownames = TRUE,
  show_colnames = TRUE
)

ggsave(
  filename = paste0(ExpID, "_BA_Correlation_Animals_values.png"),
  plot = p_cor,
  path = output_pwd,
  width = 8,
  height = 5,
  dpi = 300,
  bg = "white"
)

ba_sd <- apply(d_mat_pearson, 1, sd, na.rm = TRUE)
d_mat_cor <- d_mat_pearson[ba_sd > 0, , drop = FALSE]
d_mat_person_values <- t(d_mat_cor) # d_mat is transposed so that bile acids are columns and bile acids are being correlated
cor_pearson_values <- cor(d_mat_person_values,method = "pearson", use = "pairwise.complete.obs")

round(cor_pearson_values, 2)
summary(cor_pearson_values[upper.tri(cor_pearson_values)])

p_cor <- pheatmap(
  cor_pearson_values,
  cluster_rows = F,
  cluster_cols = F,
  fontsize = 10,
  fontsize_main = 7,
  fontsize_row = 9,
  fontsize_col = 9,
  fontsize_number = 8,
  color = colorRampPalette(c("blue", "white", "red"))(100),
  breaks = seq(-1, 1, length.out = 101),
  border_color = "black",
  cellwidth = 10,
  cellheight = 10,
  angle_col = 90,
  display_numbers = FALSE,
  show_rownames = TRUE,
  show_colnames = TRUE
)

ggsave(
  filename = paste0(ExpID, "_BA_Correlation_Parameters_values.png"),
  plot = p_cor,
  path = output_pwd,
  width = 8,
  height = 5,
  dpi = 300,
  bg = "white"
)


# Z score scaled heatmap - Bile Acids ------------------------------------

d_scaled <- t(scale(t(d_mat)))


# Effect size type wie bei Exigo
stats <- stats %>%
  mutate(
    effect_size_type = factor(
      effect_size_type,
      levels = c("standardized model effect", "NA")
    )
  )


# Sort parameters: first by effect size type,
# then by adjusted Treatment p-value
p_order <- stats %>%
  arrange(
    is.na(effect_size_type),
    effect_size_type,
    adj_p_Treatment
  ) %>%
  pull(BA)

p_order <- intersect(p_order, rownames(d_scaled))


ann$Treatment <- factor(
  ann$Treatment,
  levels = c("Ctrl", "TAM")
)

ann$Sex <- factor(
  ann$Sex,
  levels = c("female", "male")
)

column_order <- order(
  ann$Treatment,
  ann$Sex
)


# Row annotation: adjusted p-value categories
row_annot <- stats %>%
  select(BA, adj_p_Treatment) %>%
  mutate(
    adj.p = cut(
      adj_p_Treatment,
      breaks = c(
        -Inf,
        0.0001,
        0.001,
        0.01,
        0.05,
        Inf
      ),
      labels = c(
        "<0.0001",
        "0.0001–0.001",
        "0.001–0.01",
        "0.01–0.05",
        ">0.05"
      ),
      include.lowest = TRUE
    )
  ) %>%
  select(BA, adj.p) %>%
  tibble::column_to_rownames("BA")

row_annot$adj.p <- as.character(row_annot$adj.p)

row_annot$adj.p[is.na(row_annot$adj.p)] <- "NA"

row_annot$adj.p <- factor(
  row_annot$adj.p,
  levels = c(
    "<0.0001",
    "0.0001–0.001",
    "0.001–0.01",
    "0.01–0.05",
    ">0.05",
    "NA"
  )
)


# Apply column and row order
d_scaled <- d_scaled[
  ,
  column_order,
  drop = FALSE
]

ann <- ann[
  column_order,
  ,
  drop = FALSE
]

d_scaled <- d_scaled[
  p_order,
  ,
  drop = FALSE
]

row_annot <- row_annot[
  p_order,
  ,
  drop = FALSE
]


# Censoring information for bile acids -------------------------------

cens_cols <- paste0(param_list, "_censored")
direction_cols <- paste0(param_list, "_direction")

cens_mat <- d1 %>%
  filter(Timepoint == "11") %>%
  select(all_of(cens_cols)) %>%
  as.matrix()

rownames(cens_mat) <- d1 %>%
  filter(Timepoint == "11") %>%
  pull(Animal)

cens_mat <- t(cens_mat)
rownames(cens_mat) <- param_list


direction_mat <- d1 %>%
  filter(Timepoint == "11") %>%
  select(all_of(direction_cols)) %>%
  as.matrix()

rownames(direction_mat) <- d1 %>%
  filter(Timepoint == "11") %>%
  pull(Animal)

direction_mat <- t(direction_mat)
rownames(direction_mat) <- param_list


# Reorder censoring information in exactly the same way
# as the scaled data
cens_mat <- cens_mat[
  p_order,
  column_order,
  drop = FALSE
]

direction_mat <- direction_mat[
  p_order,
  column_order,
  drop = FALSE
]


# Identify censored values
below_detection <- cens_mat == TRUE & direction_mat == "<"
above_detection <- cens_mat == TRUE & direction_mat == ">"


# Replace censored values ONLY for visualization
d_scaled[below_detection] <- -4
d_scaled[above_detection] <- 4


# Heatmap colour scale
heatmap_colors <- c(
  "#FFF5CC",
  colorRampPalette(
    c("#FFE699", "orange", "red")
  )(98),
  "darkred"
)

heatmap_breaks <- seq(
  -4,
  4,
  length.out = length(heatmap_colors) + 1
)


# p-value annotation colours
p_color_list <- c(
  "<0.0001" = "darkgreen",
  "0.0001–0.001" = "forestgreen",
  "0.001–0.01" = "green3",
  "0.01–0.05" = "green1",
  ">0.05" = "white",
  "NA" = "grey80"
)


# Z score scaled Bile Acid heatmap -------------------------------------

d_s <- pheatmap(
  d_scaled,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  fontsize = 10,
  fontsize_main = 7,
  fontsize_row = 9,
  fontsize_col = 9,
  fontsize_number = 8,
  annotation_col = ann,
  annotation_row = row_annot,
  color = heatmap_colors,
  breaks = heatmap_breaks,
  annotation_colors = list(
    Treatment = Treatment_colors[c("Ctrl", "TAM")],
    Sex = Sex_colors,
    adj.p = p_color_list
  ),
  border_color = "black",
  cellwidth = 10,
  cellheight = 10,
  angle_col = 90,
  gaps_row = 0,
  display_numbers = FALSE,
  number_format = "%.3f",
  show_colnames = FALSE,
  legend_breaks = c(-4, -2, 0, 2, 4),
  legend_labels = c(
    "Below LOD",
    "-2",
    "0",
    "2",
    "Above ULOQ"
  )
)

ggsave(
  filename = paste0(ExpID, "_BA_scaled.png"),
  plot = d_s,
  path = output_pwd,
  width = 6,
  height = 5,
  dpi = 300,
  bg = "white"
)
# Effect size plot --------------------------------------------------------

stats_eff_size <- stats %>%
  filter(!is.na(effect_size)) %>%
  mutate(BA = factor(BA, levels = rev(p_order)))

stats_lm <- stats_eff_size %>%
  filter(effect_size_type == "standardized model effect")

p_lm <- ggplot(stats_lm, aes(x = effect_size, y = BA)) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_errorbar(aes(xmin = effect_CI_low, xmax = effect_CI_high), height = 0) +
  geom_point(size = 2) +
  scale_y_discrete(limits = rev(p_order)) +
  labs(x = "Standardized model-based effect TAM - Ctrl", y = NULL) +
  theme_classic()

ggsave(p_lm, file = paste0(ExpID, "_BA_Effect_size_LM.png"), dpi = 300, width = 5, height = 3, path = output_pwd)


# PCA ----------------------------------------------------------------------

pca_data <- d1 %>%
  filter(Timepoint == "11") %>%
  select(Animal, Sex, Treatment, all_of(param_list),-GCDCA,-GDCA) %>%
  drop_na() 

numerical_d <- pca_data %>%
  select(-Animal, -Sex, -Treatment)

scaled_d <- scale(numerical_d)
data.pca <- prcomp(scaled_d)

summary(data.pca)
data.pca$rotation[, 1:2]

f1 <- fviz_eig(data.pca, addlabels = TRUE)
f2 <- fviz_pca_var(data.pca, col.var = "black")
f3 <- fviz_cos2(data.pca, choice = "var", axes = 1:2)
f4 <- fviz_pca_var(data.pca, col.var = "cos2", gradient.cols = c("black", "orange", "green"), repel = TRUE)
f5 <- fviz_contrib(data.pca, choice = "var", axes = 1, top = 15, sort.val = c("desc"))

f6 <- autoplot(data.pca, data = pca_data, x = 1, y = 2, size = 3, fill = "Treatment", color = "Treatment", shape = "Sex") +
  theme_bw() +
  theme_classic() +
  ggtitle("Principal Component Analysis") +
  scale_color_manual(values = Treatment_colors[c("Ctrl","TAM")]) +
  scale_fill_manual(values = Treatment_colors[c("Ctrl","TAM")]) +
  scale_shape_manual(values = Sex_shape) +
  theme(text = element_text(size = 20))

ggsave(f1, path = output_pwd, file = paste0(ExpID, "_BA_PCA_ScreePlot.png"), dpi = 300, width = 4, height = 2.5)
ggsave(f2, path = output_pwd, file = paste0(ExpID, "_BA_PCA_VariablePlot.png"), dpi = 300, width = 5, height = 4)
ggsave(f3, path = output_pwd, file = paste0(ExpID, "_BA_PCA_Cos2.png"), dpi = 300, width = 5, height = 4)
ggsave(f4, path = output_pwd, file = paste0(ExpID, "_BA_PCA_Variable_Cos2.png"), dpi = 300, width = 5, height = 4)
ggsave(f5, path = output_pwd, file = paste0(ExpID, "_BA_PCA_Variable_Contribution_PC1.png"), dpi = 300, width = 5, height = 4)
ggsave(f6, path = output_pwd, file = paste0(ExpID, "_BA_PCA_Scores.png"), dpi = 300, width = 9, height = 7)


# Save plotting data ------------------------------------------------------

saveRDS(
  list(
    data = d1 %>% filter(Timepoint == "11"),
    statistics = stats,
    PCA = data.pca
  ),
  file = file.path(output_pwd, paste0(ExpID, "_BA_TP11_PlottingData.rds"))
)

cat("\nBile acid plots, heatmaps, effect size plot and PCA saved to:\n", output_pwd, "\n")


# Volcano plot for Bile Acids -----
# Volcano plot based on the Treatment comparison at TP11.
# Fold change is calculated directly from the observed group means.
# log2FC = log2(mean TAM / mean Ctrl).
# Significance is based on the FDR-adjusted Treatment p-value.

volcano_data <- stats %>%
  mutate(
    FC = mean_TAM / mean_Ctrl,
    log2FC = log2(FC),
    negLog10FDR = -log10(adj_p_Treatment),
    direction = case_when(
      adj_p_Treatment < 0.05 & log2FC < 0 ~ "TAM lower",
      adj_p_Treatment < 0.05 & log2FC > 0 ~ "TAM higher",
      TRUE ~ "Not significant"
    ),
    significant = !is.na(adj_p_Treatment) & adj_p_Treatment < 0.05,
    Metabolite = BA
  ) %>%
  filter(
    is.finite(log2FC),
    is.finite(negLog10FDR)
  )

p_volcano <- ggplot(
  volcano_data,
  aes(x = log2FC, y = negLog10FDR)
) +
  geom_point(
    aes(fill = direction),
    alpha = 0.5,
    size = 3,
    stroke = 0.5,
    position = position_jitter(width = 0.08),
    shape = 21,
    color = "black"
  ) +
  scale_fill_manual(
    values = c(
      "TAM lower" = "blue",
      "Not significant" = "grey60",
      "TAM higher" = "firebrick"
    )
  ) +
  geom_vline(
    xintercept = c(-0.5, 0.5),
    linetype = "dashed",
    color = "grey80"
  ) +
  geom_hline(
    yintercept = -log10(0.05),
    linetype = "dashed",
    color = "grey80"
  ) +
  labs(
    title = paste0("Volcano plot - Treatment at TP11"),
    x = expression(paste("log"[2], " FC (TAM / Ctrl)")),
    y = expression(paste("-log"[10], "(adj.p.value)"))
  ) +
  theme_classic() +
  theme(
    panel.grid = element_line(
      color = "grey90",
      linewidth = 0.1
    )
  ) +
  ggrepel::geom_text_repel(
    data = volcano_data %>%
      filter(significant == TRUE),
    aes(label = Metabolite),
    size = 4,
    max.overlaps = 20
  ) +
  coord_cartesian(
    xlim = c(-4, 4),
    ylim = c(0, 3)
  )

print(p_volcano)

ggsave(
  plot = p_volcano,
  filename = paste0(ExpID, "_BA_volcano.png"),
  width = 6,
  height = 9,
  dpi = 300,
  path = output_pwd
)

