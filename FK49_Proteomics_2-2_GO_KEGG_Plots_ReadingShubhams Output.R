rm(list = ls())
gc()

library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(enrichplot)
library(circlize)
library(httr)
library(jsonlite)
library(org.Mm.eg.db)
library(pheatmap)
library(ReactomeContentService4R)
library(AnnotationDbi)
source("FK49_Definitions.R")

proteom_input_pwd  <- PATHS$proteomics$input
proteom_output_pwd <- PATHS$proteomics$output
stats_pwd <- file.path(proteom_output_pwd, "Statistics")

# Load protein statistics
Proteins <- read.csv2(file.path(proteom_output_pwd, "FK49_Proteomics_statistics.csv"))

# Significant proteins
sig_proteins <- Proteins %>%
  filter(adj_pvalue_Treatment < 0.05, abs(logFC_Treatment) > 1)

up_proteins <- sig_proteins %>%
  filter(logFC_Treatment > 1) %>%
  pull(Genes) %>%
  unique()

down_proteins <- sig_proteins %>%
  filter(logFC_Treatment < -1) %>%
  pull(Genes) %>%
  unique()

all_proteins <- sig_proteins %>%
  pull(Genes) %>%
  unique()

background <- Proteins %>%
  pull(Genes) %>%
  unique()

ego_up <- read.csv2(file.path(stats_pwd, "GO_ALL_up.csv"))
ego_up_simplified <- read.csv2(file.path(stats_pwd, "GO_ALL_up_simplified.csv"))

ego_down <- read.csv2(file.path(stats_pwd, "GO_ALL_down.csv"))
ego_down_simplified <- read.csv2(file.path(stats_pwd, "GO_ALL_down_simplified.csv"))

ego_all <- read.csv2(file.path(stats_pwd, "GO_ALL_significant.csv"))
ego_all_simplified <- read.csv2(file.path(stats_pwd, "GO_ALL_significant_simplified.csv"))

kegg_up <- read.csv2(file.path(stats_pwd, "KEGG_up.csv"))
kegg_down <- read.csv2(file.path(stats_pwd, "KEGG_down.csv"))
kegg_all <- read.csv2(file.path(stats_pwd, "KEGG_ALL_significant.csv"))

reactome_up <- read.csv2(file.path(stats_pwd, "Reactome_up.csv"))
reactome_down <- read.csv2(file.path(stats_pwd, "Reactome_down.csv"))
reactome_all <- read.csv2(file.path(stats_pwd, "Reactome_ALL_significant.csv"))



# Chord plot function ---------------------------------------------------------
save_chord <- function(enrich, filename, n = 10) {
  
  if (is.null(enrich) || nrow(as.data.frame(enrich)) == 0) return(NULL)
  
  df <- as.data.frame(enrich) %>%
    filter(p.adjust < 0.05) %>%
    arrange(p.adjust) %>%
    slice_head(n = n)
  
  if (nrow(df) == 0) return(NULL)
  
  links <- df %>%
    dplyr::select(Description, geneID) %>%
    tidyr::separate_rows(geneID, sep = "/") %>%
    filter(!is.na(geneID), geneID != "")
  
  links <- unique(links)
  
  terms <- unique(links$Description)
  proteins <- unique(links$geneID)
  
  grid.col <- c(
    setNames(grDevices::hcl.colors(length(terms), "Set 3"), terms),
    setNames(rep("grey70", length(proteins)), proteins)
  )
  
  png(
    file.path(proteom_output_pwd, filename),
    width = 12, height = 12, units = "in", res = 300
  )
  
  chordDiagram(
    links,
    grid.col = grid.col,
    grid.border = NA,
    transparency = 0,
    directional = 0,
    direction.type = "diffHeight",
    diffHeight = mm_h(2),
    link.sort = FALSE,
    link.decreasing = TRUE,
    link.target.prop = TRUE,
    link.border = NA,
    link.lwd = 1,
    link.lty = 1,
    link.visible = TRUE,
    link.overlap = FALSE,
    scale = FALSE,
    big.gap = 10,
    small.gap = 1,
    symmetric = FALSE,
    keep.diagonal = FALSE,
    self.link = 2,
    preAllocateTracks = list(track.height = mm_h(35)),
    annotationTrack = "grid",
    annotationTrackHeight = mm_h(c(1.5, 1, 1)),
    link.auto = TRUE,
    reduce = 1e-5
  )
  
  circos.trackPlotRegion(
    track.index = 2,
    panel.fun = function(x, y) {
      
      sector.name <- get.cell.meta.data("sector.index") %>%
        stringr::str_wrap(width = 25)
      
      circos.text(
        CELL_META$xcenter,
        mean(CELL_META$ylim)+mm_y(1),
        sector.name,
        facing = "clockwise",
        niceFacing = TRUE,
        adj = c(0, 0.5),
        cex = 1
      )
    },
    bg.border = NA
  )
  
  dev.off()
  
  circos.clear()
}


# GO chord plots --------------------------------------------------------------
save_chord(ego_up,   "Chord_GO_up.png")
save_chord(ego_down, "Chord_GO_down.png")
save_chord(ego_all,  "Chord_GO_ALL_significant.png")


# KEGG chord plots ------------------------------------------------------------
save_chord(kegg_up,   "Chord_KEGG_up.png")
save_chord(kegg_down, "Chord_KEGG_down.png")
save_chord(kegg_all,  "Chord_KEGG_ALL_significant.png")


# Reactome chord plots --------------------------------------------------------
save_chord(reactome_up,   "Chord_Reactome_up.png")
save_chord(reactome_down, "Chord_Reactome_down.png")
save_chord(reactome_all,  "Chord_Reactome_ALL_significant.png")


# Reactome proteins -----------------------------------------------------------
reactome_df <- as.data.frame(reactome_all) %>%
  filter(p.adjust < 0.05) %>%
  arrange(p.adjust) %>%
  slice_head(n = 10) %>%
  select(Description, geneID) %>%
  separate_rows(geneID, sep = "/")

# Protein abundances in long format ------------------------------------------
protein_long <- Proteins %>%
  select(Genes,adj_pvalue_Treatment, starts_with("F_"), starts_with("M_")) %>%
  pivot_longer(c(-Genes, -adj_pvalue_Treatment),  
                 names_to = "Sample",values_to = "ProteinValue") %>%
  separate(  Sample,into = c("Sex", "Treatment", "Replicate"),sep = "_" ) %>%
  select(-Replicate) %>%
  mutate( Treatment = case_when(Treatment == "EtOH" ~ "Ctrl", TRUE ~Treatment),
          Treatment = factor(Treatment, levels=c("Ctrl","TAM")),
          Sex = case_when(Sex == "F" ~ "female", Sex == "M" ~ "male"),
          Sex= factor(Sex, levels = c("female", "male")))

# Violin plots Proteins that drive enrichement of that pathway-----------------------------------------------
for (reactom in unique(reactome_df$Description)) {
  
  proteins_reactom <- reactome_df %>%
    filter(Description == reactom) %>%
    pull(geneID)
  
  plot_df <- protein_long %>%
    filter(Genes %in% proteins_reactom) %>%
    mutate(Genes = factor(Genes, levels = proteins_reactom))
  
  p <- ggplot(plot_df, aes(x = Treatment, y = ProteinValue, fill = Treatment, shape=Sex)) +
    geom_violin(trim = FALSE, alpha = 0.4) +
    geom_boxplot(width = 0.15, outlier.shape = NA) +
    geom_jitter(width = 0.08, size = 1.5) +
    scale_fill_manual(values=Treatment_colors[c("Ctrl","TAM")])+
    scale_shape_manual(values=Sex_shape)+
    facet_wrap(~ Genes, scales = "free_y") +
    theme_classic() +
    labs( title = reactom,x = NULL, y = "Log2-normalized protein abundance" )+
    geom_text(data = distinct(plot_df, Genes, adj_pvalue_Treatment),aes(x = 1.5,y = Inf,
          label = paste0("adj. p = ", format.pval(adj_pvalue_Treatment, digits = 2))),
      vjust = 1.5,inherit.aes = FALSE  )
  
  ggsave( file.path(proteom_output_pwd,  paste0("Violin_Reactome_", make.names(reactom), ".png")),
    p,  width = 12,    height = 8,  dpi = 300 )
}
# ALL SHITTTIYYY ------
# Reactome pathway heatmaps ---------------------------------------------------
# Top 10 Reactome pathways ----------------------------------------------------

top_reactome <- reactome_all %>%
  filter(p.adjust < 0.05) %>%
  arrange(p.adjust) %>%
  slice_head(n = 10) %>%
  dplyr::select(ID, Description)


# Get ALL genes from each Reactome pathway -----------------------------------

reactome_proteins <- lapply(seq_len(nrow(top_reactome)), function(i) {
    x <- getParticipants(top_reactome$ID[i],  retrieval = "ReferenceEntities" )
    x %>% dplyr::select(refDbName, refIdentifier, displayName) %>%
    mutate(  ID = top_reactome$ID[i], Description = top_reactome$Description[i] )
}) %>%
  bind_rows()

# Keep ALL pathway proteins that are present in the complete protein matrix --

pathway_genes <- reactome_proteins %>%
  filter(Genes %in% Proteins$Genes) %>%
  dplyr::select(ID, Description, Genes) %>%
  distinct()

# Heatmap data ---------------------------------------------------------------

heatmap_data <- Proteins %>%
  filter(Genes %in% pathway_genes$Genes) %>%
  dplyr::select(Genes, starts_with("F_"), starts_with("M_")) %>%
  tibble::column_to_rownames("Genes") %>%
  as.matrix()


# Row annotation: Reactome pathways ------------------------------------------

annotation_row <- pathway_genes %>%
  group_by(Genes) %>%
  summarise(
    Reactome = paste(unique(Description), collapse = " | "),
    .groups = "drop"
  ) %>%
  tibble::column_to_rownames("Genes")

annotation_row <- annotation_row[
  rownames(heatmap_data),
  ,
  drop = FALSE
]


# Column annotation -----------------------------------------------------------

annotation_col <- data.frame(
  Treatment = ifelse(
    grepl("_EtOH_", colnames(heatmap_data)),
    "Ctrl",
    "TAM"
  ),
  Sex = ifelse(
    grepl("^F_", colnames(heatmap_data)),
    "Female",
    "Male"
  ),
  row.names = colnames(heatmap_data)
)

annotation_col$Sex <- factor(
  annotation_col$Sex,
  levels = c("Female", "Male")
)

annotation_col$Treatment <- factor(
  annotation_col$Treatment,
  levels = c("Ctrl", "TAM")
)


# Heatmap ---------------------------------------------------------------------

pheatmap::pheatmap(
  heatmap_data,
  scale = "row",
  annotation_col = annotation_col,
  annotation_row = annotation_row,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  border_color = NA,
  fontsize_row = 8,
  fontsize_col = 8
)


# Save pathway-protein mapping -----------------------------------------------

write.csv2(
  reactome_proteins,
  file.path(
    stats_pwd,
    "Reactome_top10_all_pathway_proteins.csv"
  ),
  row.names = FALSE
)

write.csv2(
  pathway_genes,
  file.path(
    stats_pwd,
    "Reactome_top10_proteins_in_matrix.csv"
  ),
  row.names = FALSE
)