###############################################################################
# FK49_Proteomics_3_Plots.R
#
# Plotting of FK49 proteomics statistical results.
###############################################################################

rm(list = ls())
gc()

library(dplyr)
library(tidyr)
library(ggplot2)
library(ggrepel)
library(pheatmap)
library(stringr)
library(tibble)
library(RColorBrewer)
#library(VennDiagram)
source("FK49_Definitions.R")

proteom_output_pwd <- PATHS$proteomics$output

# Load data -------------------------------------------------------------------

Proteins <- read.csv2(file.path(proteom_output_pwd,  "Statistics/FK49_Proteomics_Statistics.csv" ))

protein_matrix <- readRDS(file.path(proteom_output_pwd,"Data/FK49_Proteomics_protein_matrix.rds" ))
meta <- readRDS(file.path(proteom_output_pwd,"Data/FK49_Proteomics_metadata.rds"))
# load calculateion of Correaltion animals
# load calculattion of Crorealtion proteins
# Volcano ---------------------------------------------------------------------
p_volcano <- ggplot(Proteins, aes(x = logFC_Treatment,y = -log10(pValue_Treatment))) +
  geom_point( aes(fill = Direction),alpha = 0.5,size = 3, stroke = 0.5, shape = 21, color = "black" ) +
  scale_fill_manual(values = c("blue","grey60","firebrick") ) +
  geom_vline( xintercept = c(-1, 1),linetype = "dashed",color = "grey80") +
  geom_hline(yintercept = -log10(0.05),linetype = "dashed", color = "grey80") +
  labs( title = "Volcano plot - Treatment",
    x = expression(paste("FC [", log[2], "]")),
    y = expression(paste("-log"[10], "(p.value)"))) +
  theme_classic() +
  theme(panel.grid = element_line(color = "grey90",linewidth = 0.1) ) +
  geom_text_repel(data = Proteins %>%filter(Direction %in% c("Up", "Down")),
                  aes(label = Name),size = 2.5, max.overlaps = 25) +
  coord_cartesian(xlim = c(-4.5, 4.5),ylim = c(0, 15) )

p_volcano

ggsave( plot = p_volcano, filename = "/Plots/01_Prots_a_volcano.png",
  width = 9,height = 9,  dpi = 300,  path = proteom_output_pwd)


# Heatmaps --------------------------------------------------------------------
## Heatmap significant up -----------------------------------------------------

values_for_heatmap_up <- as.data.frame(
  Proteins %>%
    filter( adj_pvalue_Treatment < 0.05 & logFC_Treatment > 1 ) %>%
    dplyr::select( Name, Genes,contains(c("EtOH", "TAM"))))

rownames(values_for_heatmap_up) <- values_for_heatmap_up$Name
values_for_heatmap_up$Name <- NULL

ann_col <- data.frame(Sample = colnames(values_for_heatmap_up %>% dplyr::select(-Genes))) %>%
  separate(  Sample,into = c("Sex", "Treatment", "Replicate"),sep = "_" ) %>%
  dplyr::select(-Replicate) %>%
  mutate( Sex = case_when(Sex == "F" ~ "female",Sex == "M" ~ "male" ),
    Treatment = case_when(Treatment == "EtOH" ~ "Ctrl", TRUE ~ Treatment),
    Sex = factor( Sex,levels = c("female", "male")),Treatment = factor(Treatment,levels = c("Ctrl", "TAM")) )

values_for_heatmap_up$Genes <- NULL
rownames(ann_col) <- colnames(values_for_heatmap_up)

ph <- pheatmap(
  values_for_heatmap_up,
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  show_colnames = FALSE,
  scale = "row",
  annotation_col = ann_col,
  annotation_colors = list( Sex = Sex_colors,
                            Treatment = Treatment_colors[c("Ctrl", "TAM")])
  )

ggsave( ph, file = "/Plots/02_Heat_sig_up.png",  path = proteom_output_pwd,
  width = 12,  height = 1 + nrow(values_for_heatmap_up) / 5,
  dpi = 300,  bg = "white",  limitsize = FALSE)


## Heatmap significant down ---------------------------------------------------
values_for_heatmap_down <- as.data.frame(
  Proteins %>%filter( adj_pvalue_Treatment < 0.05 & logFC_Treatment < -1) %>%
    dplyr::select( Name, Genes,contains(c("EtOH", "TAM"))))

rownames(values_for_heatmap_down) <- values_for_heatmap_down$Name
values_for_heatmap_down$Name <- NULL

ann_col <- data.frame(Sample = colnames(values_for_heatmap_down %>% dplyr::select(-Genes))) %>%
  separate( Sample, into = c("Sex", "Treatment", "Replicate"),sep = "_") %>%
  dplyr::select(-Replicate) %>%
  mutate(Sex = case_when( Sex == "F" ~ "female",Sex == "M" ~ "male"),
    Treatment = case_when(Treatment == "EtOH" ~ "Ctrl", TRUE ~ Treatment),
    Sex = factor( Sex,levels = c("female", "male")),
    Treatment = factor(Treatment,levels = c("Ctrl", "TAM"))
  )

values_for_heatmap_down$Genes <- NULL
rownames(ann_col) <- colnames(values_for_heatmap_down)

ph <- pheatmap(
  values_for_heatmap_down,
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  scale = "row",
  show_colnames = FALSE,
  annotation_col = ann_col,
  annotation_colors = list( Sex = Sex_colors,
    Treatment = Treatment_colors[c("Ctrl", "TAM")])
)

ggsave(
  ph,
  file = "/Plots/02_Heat_sig_down.png",
  path = proteom_output_pwd,
  width = 12,
  height = 1 + nrow(values_for_heatmap_down) / 5,
  dpi = 300,
  bg = "white",
  limitsize = FALSE
)


## Top 50 up + Top 50 down ----------------------------------------------------

values_for_heatmap_top100 <- as.data.frame(
  Proteins %>%
    filter(adj_pvalue_Treatment < 0.05) %>%
    arrange(logFC_Treatment) %>%
    slice_head(n = 50) %>%
    bind_rows(
      Proteins %>%
        filter(adj_pvalue_Treatment < 0.05) %>%
        arrange(desc(logFC_Treatment)) %>%
        slice_head(n = 50)) %>%
    dplyr::select( Name,  Genes, contains(c("EtOH", "TAM"))))

rownames(values_for_heatmap_top100) <-values_for_heatmap_top100$Name

values_for_heatmap_top100$Name <- NULL

ann_col <- data.frame(Sample = colnames(values_for_heatmap_top100 %>% dplyr::select(-Genes))) %>%
  separate(  Sample,into = c("Sex", "Treatment", "Replicate"),sep = "_" ) %>%
  dplyr::select(-Replicate) %>%
  mutate( Sex = case_when(Sex == "F" ~ "female",Sex == "M" ~ "male" ),
    Treatment = case_when(Treatment == "EtOH" ~ "Ctrl", TRUE ~ Treatment),
    Sex = factor(Sex,levels = c("female", "male") ),
    Treatment = factor( Treatment,levels = c("Ctrl", "TAM")))

values_for_heatmap_top100$Genes <- NULL
rownames(ann_col) <- colnames(values_for_heatmap_top100)

ph <- pheatmap(
  values_for_heatmap_top100,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  scale = "row",
  show_colnames = FALSE,
  annotation_col = ann_col,
  annotation_colors = list(
    Sex = Sex_colors,
    Treatment = Treatment_colors[c("Ctrl", "TAM")]
  )
)

ggsave( ph, file = "/Plots/02_Heat_top50up_top50down.png", path = proteom_output_pwd,
  width = 12, height = 2 + nrow(values_for_heatmap_top100) / 8,
  dpi = 300,bg = "white", limitsize = FALSE)


# PCA -------------------------------------------------------------------------
## PCA all proteins -----------------------------------------------------------
pca <- prcomp(t(protein_matrix),center = TRUE,scale. = TRUE)

pca_df <- as.data.frame(pca$x) %>%
  rownames_to_column("Sample") %>%
  separate( Sample, into = c("Sex", "Treatment", "Replicate"),sep = "_" ) %>%
  mutate(Sex = case_when(Sex == "F" ~ "female",Sex == "M" ~ "male"),
    Treatment = case_when(Treatment == "EtOH" ~ "Ctrl",TRUE ~ Treatment),
    Sex = factor(Sex,levels = c("female", "male")),
    Treatment = factor(Treatment,levels = c("Ctrl", "TAM")))

pca_plot <- ggplot( pca_df, aes(x = PC1,y = PC2,shape = Sex,color = Treatment,fill = Treatment)) +
  geom_point( size = 4,  alpha = 0.5) + stat_ellipse(  aes(group = Treatment,fill = Treatment ),
    type = "norm",level = 0.95,  geom = "polygon",alpha = 0.05, linewidth = 0) +
  stat_ellipse( aes(group = Treatment,color = Treatment),
    type = "norm",level = 0.95,  alpha = 0.2,  linewidth = 0.8) +
  scale_color_manual( values = Treatment_colors[c("Ctrl", "TAM")] ) +
  scale_fill_manual(values = Treatment_colors[c("Ctrl", "TAM")] ) +
  scale_shape_manual( values = Sex_shape) +
  theme_classic() +
  labs( x = paste0("PC1 (", round(100 * summary(pca)$importance[2, 1], 1),"%)"),
        y = paste0("PC2 (", round(100 * summary(pca)$importance[2, 2], 1),"%)"))

pca_plot

ggsave( pca_plot, file = "/Plots/03_PCA_Proteins_all.png", path = proteom_output_pwd,
  width = 6.5,height = 6, dpi = 300,bg = "white", limitsize = FALSE)

## PCA significant proteins ---------------------------------------------------

pca_proteins <- Proteins %>%
  filter( adj_pvalue_Treatment < 0.05,  abs(logFC_Treatment) > 1) %>%
  pull(Name)

protein_matrix_sig <- protein_matrix[pca_proteins, ]

pca_sig <- prcomp(t(protein_matrix_sig),center = TRUE,scale. = TRUE)

pca_sig_df <- as.data.frame(pca_sig$x) %>%
  rownames_to_column("Sample") %>%
  separate(  Sample,into = c("Sex", "Treatment", "Replicate"),sep = "_" ) %>%
  mutate( Sex = case_when(Sex == "F" ~ "female",Sex == "M" ~ "male" ),
    Treatment = case_when(Treatment == "EtOH" ~ "Ctrl", TRUE ~ Treatment),
    Sex = factor(Sex,levels = c("female", "male")),
    Treatment = factor(Treatment,levels = c("Ctrl", "TAM"))
  )

pca_sig_plot <- ggplot( pca_sig_df,  aes(x = PC1,y = PC2,shape = Sex,color = Treatment,fill = Treatment)) +
  geom_point( size = 4, alpha = 0.5 ) +
  stat_ellipse( aes(group = Treatment,color = Treatment),
    type = "norm",level = 0.95,  geom = "polygon",  alpha = 0.05,linewidth = 0) +
  stat_ellipse( aes(group = Treatment,color = Treatment),  type = "norm", level = 0.95,
    alpha = 0.2,linewidth = 0.8 ) +
  scale_color_manual( values = Treatment_colors[c("Ctrl", "TAM")] ) +
  scale_fill_manual(values = Treatment_colors[c("Ctrl", "TAM")] ) +
  scale_shape_manual( values = Sex_shape) +
  theme_classic() +
  labs(x = paste0("PC1 (", round(100 * summary(pca_sig)$importance[2, 1], 1),"%)"),
       y = paste0("PC2 (", round(100 * summary(pca_sig)$importance[2, 2], 1),"%)"))

pca_sig_plot

ggsave( pca_sig_plot, file = "/Plots/03_PCA_Proteins_sig.png", path = proteom_output_pwd,
  width = 6.5,height = 6, dpi = 300,bg = "white", limitsize = FALSE)
# Correaltions PLots HEatmaps -----
## Correaltion between ANimals -----
# ggsave(Cor_Animals_spearman, file= file.path(proteom_output_pwd,"Plots/04_Correlation_Animals_Spearman.png, width = 8, height = 8, dpi = 300, bg= "white)
# ggsave(Cor_Animals_pearson, file= file.path(proteom_output_pwd,"Plots/04_Correlation_Animals_Pearson.png, width = 8, height = 8, dpi = 300, bg= "white)

# Correaltion between proteins -----
# ggsave(Cor_Proteins_spearman, file= file.path(proteom_output_pwd,"Plots/04_Correlation_Proteins_Spearman.png, width = 8, height = 8, dpi = 300, bg= "white)
# ggsave(Cor_Proteins_pearson, file= file.path(proteom_output_pwd,"Plots/04_Correlation_Proteins_Pearson.png, width = 8, height = 8, dpi = 300, bg= "white)

# VennDagramm -----
# Venn Proteins Overall in TAM EtOH
# not possible i only have filtered proteins that seem to be present in both sets
# Exploratory heatmaps --------------------------------------------------------

biological_oxidations <- c(
  "Nnmt", "Cyp4a10", "Gsta", "Cyp2a4", "Gstt2",
  "Ugt1a9", "Cyp2c29", "Cyp3a11", "Gstt1", "Fmo2",
  "Ugp2", "Mgst3", "Sult1c2", "Sult1b1", "Adh4"
) 
metabolism <- c(
  "Adh4", "Abcc3", "Acot1", "Ttr", "Alb", "C3", "Ca3",
  "Fabp2", "Cox7a1", "Apoa4", "Fabp5", "Hsd3b5",
  "Acsl6", "Lpgat1", "Tymp", "Slco1a4", "Acot2", "Mvk",
  "Fga", "Orm2", "Orm1", "Gas6", "Itih3", "Fgb", "Fgg"
)

phaseII_conjugation_of_compounds <- c(
  "Nnmt", "Gsta", "Ugp2", "Ugt1a9", "Gstt2",
  "Gstt1", "Mgst3", "Sult1c2", "Sult1b1"
)

metabolism_of_lipids <- c(
  "Fabp2", "Cox7a1", "Apoa4", "Fabp5", "Hsd3b5",
  "Acsl6", "Lpgat1", "Acot2", "Mvk"
)

platelet_degranulation <- c(
  "Fga", "Fgb", "Fgg", "Orm1", "Orm2", "Gas6", "Itih3"
)

response_to_elevated_platelet_cytosolic_Ca2 <- c(
  "Fga", "Fgb", "Fgg", "Gas6", "Itih3"
)

drug_ADME <- c(
  "Cyp4a10", "Cyp2a4", "Ugt1a9", "Cyp2c29",
  "Cyp3a11", "Fmo2", "Mgst3", "Sult1c2", "Adh4"
)

ciprofloxacin_ADME <- c(
  "Fga", "Fgb", "Fgg", "Gas6"
)

GRB2_SOS_integrin_MAPK <- c(
  "Itih3", "Gas6", "Fga"
)

p130Cas_integrin_MAPK <- c(
  "Itih3", "Gas6", "Fga"
)


# Bile acid heatmap -----------------------------------------------------------

Uptake <- c(
  "Slc10a1", "Slco1a1", "Slco1b2", "Slco2b1"
)

Secretion <- c(
  "Abcb11", "Abcc2", "Abcg5", "Abcg8"
)

Synthesis <- c(
  "Cyp7a1", "Cyp8b1", "Cyp27a1", "Cyp7b1",
  "Hsd3b7", "Akr1d1", "Baat"
)

TF <- c(
  "Nr1h4", "Nr0b2", "Rxra", "Hnf4a", "Ppara",
  "Ppard", "Pparg"
)

Conjugation <- c(
  "Slc27a5", "Baat", "Sult2a1", "Sult2a2",
  "Ugt1a1", "Ugt1a6", "Ugt2b1"
)

values_for_heatmap <- as.data.frame(
  Proteins %>%
    dplyr::select(
      Name,
      Genes,
      contains(c("Etoh", "Tam"))
    ) %>%
    filter(
      Genes %in% c(
        Uptake,
        Secretion,
        Synthesis,
        TF,
        Conjugation
      )
    )
)

rownames(values_for_heatmap) <- values_for_heatmap$Name

row_category <- case_when(
  values_for_heatmap$Genes %in% Uptake ~ "Uptake",
  values_for_heatmap$Genes %in% Secretion ~ "Secretion",
  values_for_heatmap$Genes %in% Synthesis ~ "Synthesis",
  values_for_heatmap$Genes %in% TF ~ "TF",
  values_for_heatmap$Genes %in% Conjugation ~ "Conjugation"
)

row_category <- factor(
  row_category,
  levels = c(
    "Synthesis",
    "Secretion",
    "Conjugation",
    "Uptake",
    "TF"
  )
)

ann <- data.frame(
  Sample = colnames(values_for_heatmap)
) %>%
  separate(
    Sample,
    into = c("Sex", "Treatment", "Replicate"),
    sep = "_"
  ) %>%
  dplyr::select(-Replicate) %>%
  mutate(
    Sex = case_when(
      Sex == "F" ~ "female",
      Sex == "M" ~ "male"
    )
  )

annotation_row <- data.frame(
  Uptake = values_for_heatmap$Genes %in% Uptake,
  Secretion = values_for_heatmap$Genes %in% Secretion,
  Synthesis = values_for_heatmap$Genes %in% Synthesis,
  TF = values_for_heatmap$Genes %in% TF,
  Conjugation = values_for_heatmap$Genes %in% Conjugation,
  row.names = values_for_heatmap$Genes
)

annotation_row[] <- lapply(
  annotation_row,
  factor,
  levels = c(FALSE, TRUE)
)

annotation_colors <- list(
  Sex = Sex_colors,
  Treatment = Treatment_colors,
  Diet = Diet_colors
)

annotation_colors$Uptake <-
  c("FALSE" = "grey95", "TRUE" = "#8dd3c7")

annotation_colors$Secretion <-
  c("FALSE" = "grey95", "TRUE" = "#ffffb3")

annotation_colors$Synthesis <-
  c("FALSE" = "grey90", "TRUE" = "#bebada")

annotation_colors$TF <-
  c("FALSE" = "grey90", "TRUE" = "forestgreen")

annotation_colors$Conjugation <-
  c("FALSE" = "grey90", "TRUE" = "brown3")

rownames(annotation_row) <- rownames(values_for_heatmap)
rownames(ann) <- colnames(values_for_heatmap)

ord <- order(row_category)

values_for_heatmap <- values_for_heatmap[ord, ]
annotation_row <- annotation_row[ord, ]

values_for_heatmap <- values_for_heatmap %>%
  dplyr::select(-Name, -Genes)

gaps_row <- cumsum(
  table(row_category[ord])
)[-length(table(row_category))]

p_heat <- pheatmap::pheatmap(
  values_for_heatmap,
  scale = "row",
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  annotation_col = ann,
  annotation_row = annotation_row,
  annotation_names_row = FALSE,
  annotation_colors = annotation_colors,
  show_colnames = FALSE,
  treeheight_row = 0,
  treeheight_col = 5,
  gaps_row = gaps_row,
  main = "Bile Acid Metabolism"
)

print(p_heat)

ggsave(
  plot = p_heat,
  filename = "/Plots/Heatmap_BA_clustered.png",
  limitsize = FALSE,
  width = 9,
  height = nrow(values_for_heatmap) / 7 + 3,
  dpi = 500,
  bg = "white",
  path = proteom_output_pwd
)


# PPAR heatmap ----------------------------------------------------------------

PPARa_FAoxidation <- c(
  "Ppara", "Acox1", "Acot1", "Acot2", "Acsl1", "Acadm",
  "Acadl", "Cpt1a", "Cyp4a10", "Cyp4a14", "Ehhadh",
  "Fabp1", "Fgf21", "Hmgcs2", "Pdk4"
)

PPARg_lipogenesis <- c(
  "Pparg", "Cd36", "Fabp4", "Plin2", "Lpl",
  "Scd1", "Fasn", "Adipoq"
)

PPAR_Coactivators <- c(
  "Rxra", "Ppargc1a", "Ppargc1b",
  "Nr1h4", "Srebf1", "Hnf4a"
)

values_for_heatmap <- as.data.frame(
  Proteins %>%
    dplyr::select(
      Name,
      Genes,
      contains(c("Etoh", "Tam"))
    ) %>%
    filter(
      Genes %in% c(
        PPARa_FAoxidation,
        PPARg_lipogenesis,
        PPAR_Coactivators
      )
    )
)

rownames(values_for_heatmap) <- values_for_heatmap$Name

row_category <- case_when(
  values_for_heatmap$Genes %in% PPARa_FAoxidation ~
    "PPARa_FAoxidation",
  values_for_heatmap$Genes %in% PPAR_Coactivators ~
    "PPAR_Coactivators",
  values_for_heatmap$Genes %in% PPARg_lipogenesis ~
    "PPARg_lipogenesis"
)

row_category <- factor(
  row_category,
  levels = c(
    "PPARa_FAoxidation",
    "PPAR_Coactivators",
    "PPARg_lipogenesis"
  )
)

ann <- data.frame(
  Sample = colnames(values_for_heatmap)
) %>%
  separate(
    Sample,
    into = c("Sex", "Treatment", "Replicate"),
    sep = "_"
  ) %>%
  dplyr::select(-Replicate) %>%
  mutate(
    Sex = case_when(
      Sex == "F" ~ "female",
      Sex == "M" ~ "male"
    )
  )

annotation_row <- data.frame(
  PPARa_FAoxidation =
    values_for_heatmap$Genes %in% PPARa_FAoxidation,
  PPARg_lipogenesis =
    values_for_heatmap$Genes %in% PPARg_lipogenesis,
  PPAR_Coactivators =
    values_for_heatmap$Genes %in% PPAR_Coactivators,
  row.names = values_for_heatmap$Genes
)

annotation_row[] <- lapply(
  annotation_row,
  factor,
  levels = c(FALSE, TRUE)
)

annotation_colors <- list(
  Sex = Sex_colors,
  Treatment = Treatment_colors,
  Diet = Diet_colors
)

annotation_colors$PPARg_lipogenesis <-
  c("FALSE" = "grey95", "TRUE" = "#8dd3c7")

annotation_colors$PPARa_FAoxidation <-
  c("FALSE" = "grey95", "TRUE" = "#ffffb3")

annotation_colors$PPAR_Coactivators <-
  c("FALSE" = "grey95", "TRUE" = "orange2")

rownames(annotation_row) <- rownames(values_for_heatmap)
rownames(ann) <- colnames(values_for_heatmap)

ord <- order(row_category)

values_for_heatmap <- values_for_heatmap[ord, ]
annotation_row <- annotation_row[ord, ]

values_for_heatmap <- values_for_heatmap %>%
  dplyr::select(-Name, -Genes)

gaps_row <- cumsum(
  table(row_category[ord])
)[-length(table(row_category))]

p_heat <- pheatmap::pheatmap(
  values_for_heatmap,
  scale = "row",
  cluster_rows = FALSE,
  cluster_cols = TRUE,
  annotation_col = ann,
  annotation_row = annotation_row,
  annotation_names_row = FALSE,
  annotation_colors = annotation_colors,
  show_colnames = FALSE,
  treeheight_row = 0,
  treeheight_col = 5,
  gaps_row = gaps_row,
  main = "PPARs"
)

print(p_heat)

ggsave(
  plot = p_heat,
  filename = "/Plots/Heatmap_PPAR_clustered.png",
  limitsize = FALSE,
  width = 9,
  height = nrow(values_for_heatmap) / 7 + 3,
  dpi = 500,
  bg = "white",
  path = proteom_output_pwd
)


# Metabolic pathways ----------------------------------------------------------

Metabolic_Pathways <- list(
  "Purine metabolism" = c(
    "Ppat", "Gart", "Pfas", "Paics", "Adsl", "Atic",
    "Impdh1", "Impdh2", "Gmps", "Ampd1", "Ampd2", "Ampd3",
    "Nt5e", "Ada", "Xdh", "Hprt", "Aprt", "Adk", "Ak1", "Ak2",
    "Nme1", "Nme2", "Pnp", "Entpd1", "Entpd2", "Entpd5",
    "Gda", "Adenosine kinase", "Itpa", "Enpp1"
  )
)

Metabolic_Pathways <- unique(
  unlist(Metabolic_Pathways)
)

values_for_heatmap <- as.data.frame(
  Proteins %>%
    dplyr::select(
      Name,
      Genes,
      contains(c("Etoh", "Tam"))
    ) %>%
    filter(
      Genes %in% Metabolic_Pathways
    )
)

rownames(values_for_heatmap) <- values_for_heatmap$Name

annotation_row <- data.frame(
  row.names = values_for_heatmap$Genes
)

annotation_colors <- list(
  Sex = Sex_colors,
  Treatment = Treatment_colors,
  Diet = Diet_colors
)

ann <- data.frame(
  Sample = colnames(values_for_heatmap)
) %>%
  separate(
    Sample,
    into = c("Sex", "Treatment", "Replicate"),
    sep = "_"
  ) %>%
  dplyr::select(-Replicate) %>%
  mutate(
    Sex = case_when(
      Sex == "F" ~ "female",
      Sex == "M" ~ "male"
    )
  )

values_for_heatmap <- values_for_heatmap %>%
  dplyr::select(-Name, -Genes)

p_metabolism <- pheatmap::pheatmap(
  values_for_heatmap,
  scale = "row",
  cluster_rows = FALSE,
  cluster_cols = TRUE,
  annotation_col = ann,
  annotation_names_row = FALSE,
  annotation_colors = annotation_colors,
  show_colnames = FALSE,
  treeheight_row = 0,
  treeheight_col = 5,
  main = "Metabolic pathways"
)

print(p_metabolism)

ggsave(
  plot = p_metabolism,
  filename = "/Plots/Heatmap_Metabolic_pathways_clustered.png",
  limitsize = FALSE,
  width = 10,
  height = nrow(values_for_heatmap) / 7 + 3,
  dpi = 500,
  bg = "white",
  path = proteom_output_pwd
)