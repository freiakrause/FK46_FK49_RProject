rm(list = ls())
gc()

library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(circlize)
library(pheatmap)
library(gridExtra)
source("FK49_Definitions.R")
# Read Input -----
proteom_input_pwd  <- PATHS$proteomics$input
proteom_output_pwd <- PATHS$proteomics$output
stats_pwd <- file.path(proteom_output_pwd, "Statistics")

# Load protein statistics
Proteins <- read.csv2(file.path(stats_pwd, "02_LIMMA_combined_stats.csv"))
reactome_proteins <- read.csv2(file.path(proteom_output_pwd, "Data/02_Reactome_allproteins.csv"))
reactome_compounds <- read.csv2(file.path(proteom_output_pwd, "Data/02_Reactome_allcompounds.csv"))

# Significant proteins
sig_proteins <- Proteins %>%
  filter(adj_pvalue_Treatment < 0.05, abs(logFC_Treatment) > 1)

up_proteins <- sig_proteins %>% filter(logFC_Treatment > 1) %>% pull(Genes) %>% unique()

down_proteins <- sig_proteins %>% filter(logFC_Treatment < -1)%>%pull(Genes) %>%unique()

all_proteins <- sig_proteins %>%pull(Genes) %>% unique()

background <- Proteins %>%pull(Genes) %>% unique()

ORA_GO <- read.csv2(file.path(stats_pwd, "04_ORA_GO_ALL_significant.csv"))
ORA_GO_simplified <- read.csv2(file.path(stats_pwd, "04_ORA_GO_ALL_significant_simplified.csv"))
ORA_KEGG <- read.csv2(file.path(stats_pwd, "04_ORA_KEGG_ALL_significant.csv"))
ORA_reactome <- read.csv2(file.path(stats_pwd, "04_ORA_Reactome_ALL_significant.csv"))


str(ORA_reactome)
head(ORA_reactome)
# Chord plot function ---------------------------------------------------------
save_chord <- function(enrich, filename, n = 10, pathway_parents = NULL, children_of_parents = NULL,
                       pathway_lowest = NULL,  mode = "all", folder) {
  
  if (is.null(enrich) || nrow(as.data.frame(enrich)) == 0) return(NULL)
  df <- as.data.frame(enrich) %>%  arrange(p.adjust)
  
  # Reactome hierarchy filtering
  if (mode == "parents" && !is.null(pathway_parents)) {
    df <- df %>%filter(Description %in% pathway_parents)
    
  } else if (mode == "children" && !is.null(children_of_parents)) {
    df <- df %>%filter(Description %in% children_of_parents)
    
  } else if (mode == "reduced" && !is.null(pathway_lowest)) {
    df <- df %>%filter(Description %in% pathway_lowest)
  }
  
  df <- df %>%slice_head(n = n)
  
  if (nrow(df) == 0) return(NULL)
  
  # Pathway -> adjusted p-value für Beschriftung
  pathway_p <- df %>% dplyr::select(Description, p.adjust)
  
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
  
  png(paste0(file.path(proteom_output_pwd, "Pathways",folder, filename)),
    width = 14, height = 14, units = "in", res = 300)
  circos.clear()
 # circos.par(gap.degree = 1)()
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
      sector.name <- get.cell.meta.data("sector.index")
      
      if (sector.name %in% pathway_p$Description) {
          pval <- pathway_p$p.adjust[pathway_p$Description == sector.name ]
        
        print(paste0( stringr::str_wrap(sector.name, width = 25),
          "\nFDR = ", format.pval(pval, digits = 2, eps = 0.001)))
      }
      
      circos.text(
        CELL_META$xcenter,
        mean(CELL_META$ylim) + mm_y(1),
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
  
  proteins <- df %>%
    dplyr::select(Description, geneID) %>%
    tidyr::separate_rows(geneID, sep = "/") %>%
    filter(!is.na(geneID), geneID != "") %>%
    distinct(Description, geneID)
  proteins
}

# Chord plots --------------------------------------------------------------
save_chord(ORA_GO,"06_Chord_GO_ALL_significant_top10.png",n = 10,folder= "Parents")
save_chord( ORA_KEGG,  "06_Chord_KEGG_ALL_significant_top10.png",  n = 10,folder= "Parents")
chord_all <- save_chord(ORA_reactome,  "06_Chord_Reactome_ALL_significant_topall.png",n = Inf,folder="Parents")

chord_parents <- save_chord(ORA_reactome, "06_Reactome_All_parents_chord.png",  pathway_parents = PARAMETERS$Proteom$Pathway_parents,
  mode = "parents", n = Inf,folder= "Parents")

chord_children <- save_chord( ORA_reactome, "06_Reactome_All_children_chord.png", children_of_parents = PARAMETERS$Proteom$Children_of_parents,
  mode = "children",n = Inf,folder= "Children")

chord_children_top10 <- save_chord(ORA_reactome,"06_Reactome_All_children_chord_top10.png", children_of_parents = PARAMETERS$Proteom$Children_of_parents,
  mode = "children",n = 10,folder= "Children")

chord_reduced <- save_chord( ORA_reactome, "06_Reactome_All_reduced_chord.png",
  pathway_lowest = PARAMETERS$Proteom$Pathway_lowest, mode = "reduced", n = Inf,folder= "Reduced")

# Summary: number of unique proteins per chord plot --------------------------
chord_summary <- data.frame(
  Chordplot = c("All significant",  "Parents",
    "Children",   "Children top10",   "Reduced" ),
  Pathways = c(
    length(unique(chord_all$Description)),
    length(unique(chord_parents$Description)),
    length(unique(chord_children$Description)),
    length(unique(chord_children_top10$Description)),
    length(unique(chord_reduced$Description)) ),
  
  Proteins = c(  
    length(unique(chord_all$geneID)),
    length(unique(chord_parents$geneID)),
    length(unique(chord_children$geneID)),
    length(unique(chord_children_top10$geneID)),
    length(unique(chord_reduced$geneID))  )
)

print(chord_summary)


# Proteins per pathway --------------------------------------------------------
proteins_per_pathway <- chord_parents %>%
  group_by(Description) %>%
  summarise(  n_proteins = n_distinct(geneID),
    Proteins = paste(sort(unique(geneID)), collapse = ", "),  .groups = "drop") %>%
  arrange(desc(n_proteins))

print(proteins_per_pathway)

# Check: are all proteins from the significant parent pathways included? -----
significant_parent_proteins <- ORA_reactome %>%
  as.data.frame() %>%
  filter(Description %in% PARAMETERS$Proteom$Pathway_parents) %>%
  dplyr::select(Description, geneID) %>%
  tidyr::separate_rows(geneID, sep = "/") %>%
  filter(!is.na(geneID), geneID != "") %>%
  distinct(geneID)

parent_chord_proteins <- chord_parents %>%distinct(geneID)
missing_from_parent_chord <- setdiff( significant_parent_proteins$geneID, parent_chord_proteins$geneID)
extra_in_parent_chord <- setdiff( parent_chord_proteins$geneID, significant_parent_proteins$geneID)
# Save detailed tables --------------------------------------------------------
write.csv2(chord_summary,file.path(  proteom_output_pwd,  "Data/03_Reactome_Chord_protein_summary.csv"),row.names = FALSE)
write.csv2(proteins_per_pathway, file.path(proteom_output_pwd, "Data/03_Reactome_Chord_Parent_proteins_per_pathway.csv" ),
  row.names = FALSE)

# Reactome proteins -----------------------------------------------------------
reactome_df <- as.data.frame(ORA_reactome) %>%
  arrange(p.adjust) %>%
  dplyr::select(ID,Description, geneID) %>%
  separate_rows(geneID, sep = "/")


# Compare proteins in each chordplot with ALL significant ---------------------

all_chord_proteins <- chord_all %>% distinct(geneID)
compare_chord_proteins <- function(chord_proteins, name) {
  current <- chord_proteins %>%  distinct(geneID)
  missing <- setdiff( all_chord_proteins$geneID,current$geneID)
  additional <- setdiff(current$geneID, all_chord_proteins$geneID )
  common <- intersect(all_chord_proteins$geneID,current$geneID)
  data.frame( Chordplot = name,  Proteins_in_all = length(all_chord_proteins$geneID),
    Proteins_in_current = length(current$geneID),  Common = length(common),
    Missing_vs_all = length(missing),   Additional_vs_all = length(additional),
    Missing_proteins = paste(sort(missing), collapse = ", "),
    Additional_proteins = paste(sort(additional), collapse = ", "))
}


chord_comparison <- bind_rows(
  compare_chord_proteins( chord_parents,  "Parents" ),
  compare_chord_proteins( chord_children, "Children"),
  compare_chord_proteins( chord_children_top10, "Children top10"),
  compare_chord_proteins( chord_reduced,"Reduced" ))

print(chord_comparison)

write.csv2( chord_comparison, file.path(proteom_output_pwd, "Data/03_Reactome_Chord_protein_comparison_vs_all.csv" ), 
            row.names = FALSE)
write.csv2(reactome_df,file.path(proteom_output_pwd,"Data/03_Reactome_top10_df.csv"),  row.names = FALSE)


# Protein abundances in long format for violin plots ------------------------------------------
protein_long <- Proteins %>%
  dplyr::select(Name,Genes,adj_pvalue_Treatment, logFC_Treatment,starts_with("F_"), starts_with("M_")) %>%
  pivot_longer(c(-Genes, -Name,-adj_pvalue_Treatment,-logFC_Treatment),  
                 names_to = "Sample",values_to = "ProteinValue") %>%
  separate(  Sample,into = c("Sex", "Treatment", "Replicate"),sep = "_" ) %>%
  dplyr::select(-Replicate) %>%
  mutate( Treatment = case_when(Treatment == "EtOH" ~ "Ctrl", TRUE ~Treatment),
          Treatment = factor(Treatment, levels=c("Ctrl","TAM")),
          Sex = case_when(Sex == "F" ~ "female", Sex == "M" ~ "male"),
          Sex= factor(Sex, levels = c("female", "male")))


# Reactome proteins that are present in my protein matrix
reactome_prots_in_my_prots <- reactome_proteins %>%
  filter(identifier %in% Proteins$Protein.Group) %>%
  dplyr::select(ID, identifier) %>%
  distinct()%>%
  left_join(ORA_reactome%>%dplyr::select(ID, Description,Parent,Is_parent_of_significant)%>%distinct(), by = "ID") %>%
  filter(!is.na(Description))%>%
  dplyr::select(ID, Description, identifier,Parent,Is_parent_of_significant) %>%
  left_join(Proteins %>%dplyr::select(Protein.Group, Genes)%>%distinct(),   by = c("identifier" = "Protein.Group"))
head(reactome_prots_in_my_prots)
# Violin plots Proteins that drive enrichement of that pathway-----------------------------------------------
# for (reactom in unique(reactome_df$Description)) 
print_violins <- function(pathway, folder, mode= "ALL_proteins") {
  
    proteins_to_plot <- reactome_prots_in_my_prots %>%
      filter(Description == pathway | Parent == pathway ) %>%
      dplyr::select(Genes) %>%
      filter(!is.na(Genes), Genes != "") %>%
      pull(Genes) %>%
      unique()
    
   plot_df <- protein_long %>%
    filter(Genes %in% proteins_to_plot) %>%
    mutate(Genes = factor(Genes, levels = proteins_to_plot),
      Name = factor(Name))
  if(mode == "ALL_proteins"){plot_df <- plot_df%>%filter( !is.na(adj_pvalue_Treatment))}
  else if (mode == "SIG_proteins"){plot_df <- plot_df%>%filter( adj_pvalue_Treatment<0.05 & abs(logFC_Treatment) > 1)}
  else(print(" You did not specifiy how to filter proteins for plotting"))
  if (nrow(plot_df) == 0) return(NULL)
  
  n_facets <- dplyr::n_distinct(plot_df$Name) # Anzahl Spalten ungefähr quadratisch wählen 
  ncol_plot <- ceiling(sqrt(n_facets)) # Anzahl Reihen daraus berechnen 
  nrow_plot <- ceiling(n_facets / ncol_plot) # Feste Größe eines einzelnen Facets
  facet_width <- 2.2 
  facet_height <- 2.0 # Gesamtgröße des PNG 
  plot_width <- ncol_plot * facet_width
  plot_height <- nrow_plot * facet_height
  
  p <- ggplot(plot_df,aes(x = Treatment, y = ProteinValue, fill = Treatment)) +
    geom_violin(trim = FALSE, alpha = 0.4) +
    geom_point(position = position_dodge2(width = 0.5),
      size = 1.5, aes(shape = Sex), alpha = 0.5) +
    geom_boxplot( width = 0.15,  outlier.shape = NA,alpha = 0.5 ) +
    scale_fill_manual(values = Treatment_colors[c("Ctrl", "TAM")] ) +
    scale_shape_manual(values = Sex_shape) +
    facet_wrap(~ Name, scales = "free_y",ncol = ncol_plot) +
    theme_classic() +
    theme( text = element_text(size = 11), 
           strip.text = element_text(size = 10), 
           axis.text = element_text(size = 8), 
           axis.title = element_text(size = 10) )+
    labs(title = pathway,  x = NULL,   y = "Log2-normalized protein abundance"  ) +
    geom_text(  data = distinct(plot_df, Name, adj_pvalue_Treatment),
      aes(x = 1.5,  y = Inf,
        label = case_when(
          adj_pvalue_Treatment < 0.001 ~ "***",
          adj_pvalue_Treatment < 0.01  ~ "**",
          adj_pvalue_Treatment < 0.05  ~ "*",
          TRUE ~ "ns" )),      
      vjust = 1.5,inherit.aes = FALSE )
  
  ggsave( file.path(proteom_output_pwd, paste0("Pathways/", folder,"/07_Violin_",substr(make.names(pathway), 1, 15), ".png")),
    p,width = plot_width,  height = plot_height,dpi = 300 ,limitsize=FALSE)
}
for (pathway in PARAMETERS$Proteom$Pathway_parents) {
  print_violins(pathway = pathway,  folder = "Parents",mode= "SIG_proteins")
}

for (pathway in PARAMETERS$Proteom$Pathway_parents) {
  print_violins(pathway = pathway,  folder = "Parents_all",mode= "ALL_proteins")
}

for (pathway in PARAMETERS$Proteom$Children_of_parents) {
  print_violins(pathway = pathway,  folder = "Children",mode= "SIG_proteins")
}

for (pathway in PARAMETERS$Proteom$Children_of_parents) {
  print_violins(pathway = pathway,  folder = "Children_all",mode= "ALL_proteins")
}

for (pathway in PARAMETERS$Proteom$Pathway_lowest) {
  print_violins(pathway = pathway,  folder = "Reduced",mode= "SIG_proteins")
}


for (pathway in PARAMETERS$Proteom$Pathway_lowest) {
  print_violins(pathway = pathway,  folder = "Reduced_all",mode= "ALL_proteins")
}

# Heatmaps Proteins in encirhed Reactom Pathways -----

make_reactome_heatmap <- function(pathway, filename, annotation_name = NULL,mode= "SIG_proteins",folder="") {
  
     proteins_to_plot <- reactome_prots_in_my_prots %>%
    filter(Description == pathway | Parent == pathway ) %>%
    dplyr::select(Genes) %>%
    filter(!is.na(Genes), Genes != "") %>%
    pull(Genes) %>%
    unique()
    
  
  if (length(proteins_to_plot) == 0) {
    warning("Keine Proteine für ", paste(pathway, collapse = ", "))
    return(NULL)
  }
  # Matrix -----
 
  if(mode == "ALL_proteins"){Proteins <- Proteins%>%filter( !is.na(adj_pvalue_Treatment))}
  else if (mode == "SIG_proteins"){Proteins <- Proteins%>%filter( adj_pvalue_Treatment<0.05 & abs(logFC_Treatment) > 1)}
  else(print(" You did not specifiy how to filter proteins for plotting"))
  
  mat <- Proteins %>%
    filter(Genes %in% proteins_to_plot) %>%
    dplyr::select(Name, starts_with("F_"), starts_with("M_")) %>%
    distinct(Name, .keep_all = TRUE) %>%
    tibble::column_to_rownames("Name") %>%
    as.matrix()
  
  if (nrow(mat) < 2) {
    warning(  paste(pathway, collapse = ", ")," hat weniger als 2 Proteine - übersprungen"  )
    return(NULL)
  }
  
  # Row annotation -----
  annotation_row <- Proteins %>%
    filter(Name %in% rownames(mat)) %>%
    dplyr::select(Name, adj_pvalue_Treatment) %>% 
    distinct(Name, .keep_all = TRUE) %>%
    mutate( adj.p = case_when( is.na(adj_pvalue_Treatment) ~ "NA", 
                               adj_pvalue_Treatment < 0.0001 ~ "<0.0001", 
                               adj_pvalue_Treatment < 0.001 ~ "0.0001-0.001", 
                               adj_pvalue_Treatment < 0.01 ~ "0.001-0.01", 
                               adj_pvalue_Treatment < 0.05 ~ "0.01-0.05", 
                               TRUE ~ ">0.05" ) ) %>% 
    dplyr::select(Name, adj.p) %>% 
    tibble::column_to_rownames("Name") 
  annotation_row$adj.p <- factor( annotation_row$adj.p, levels = c( "<0.0001", "0.0001-0.001", "0.001-0.01", "0.01-0.05", ">0.05", "NA" ) )
  annotation_row <- annotation_row[ rownames(mat), , drop = FALSE ] 
  
  # Farben der adj. p-value Annotation ----
  p_color_list <- c( "<0.0001" = "darkgreen", 
                     "0.0001-0.001" = "forestgreen", 
                     "0.001-0.01" = "green3", 
                     "0.01-0.05" = "green1", 
                     ">0.05" = "white", 
                     "NA" = "grey80" )
    
  
  
  # Column annotation -----
  annotation_col <- data.frame(
    Treatment = ifelse( grepl("_EtOH_", colnames(mat)), "Ctrl", "TAM" ),
    Sex = ifelse( grepl("^F_", colnames(mat)), "female", "male"),    
    row.names = colnames(mat) )
  
  annotation_col$Sex <- factor( annotation_col$Sex, levels = c("female", "male"))
  annotation_col$Treatment <- factor( annotation_col$Treatment,levels = c("Ctrl", "TAM"))
  
  # Heatmap -----
  p <- pheatmap(
    mat,
    scale = "row",
    annotation_col = annotation_col,
    annotation_row = annotation_row,
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    border_color = NA,
    fontsize_row = 7.5,
    fontsize_col = 7.5,
    show_colnames = FALSE,
   # main = stringr::str_wrap(pathway, width = 50),
    annotation_colors = list( Sex = Sex_colors, Treatment = Treatment_colors[c("Ctrl", "TAM")], adj.p = p_color_list),
    annotation_legend = FALSE,
    silent = TRUE )
  
  n_rows <- nrow(mat)
  n_cols <- ncol(mat)
  
  # Heatmap speichern
  ggsave(  file.path(proteom_output_pwd,  paste0("Pathways/", folder,"/",filename)  ),
    p,  width = 4 + 0.045 * n_cols,   height = 3 + 0.08 * n_rows,   dpi = 300,  bg = "white",   limitsize = FALSE )
  print(p)
  
  invisible(p)
  
  # Set up legend sperately -----
  legend_df <- data.frame(
    Sex = factor( c("female", "male"),levels = c("female", "male")),
    Treatment = factor( c("Ctrl", "TAM"), levels = c("Ctrl", "TAM") ),
    adj.p = factor(c("<0.0001","0.0001-0.001",   "0.001-0.01", "0.01-0.05",  ">0.05",  "NA"),
                   levels = c( "<0.0001","0.0001-0.001",  "0.001-0.01",  "0.01-0.05",  ">0.05","NA"  )  ))
  
  
  
  ## Sex legend -----
  p_sex_legend <- ggplot(legend_df,aes(x = Sex,y = 1,fill = Sex)) +
    geom_tile() +
    scale_fill_manual(values = Sex_colors,name = "Sex") +
    theme_void() +
    theme(legend.position = "bottom")
  
  ## Treatment legend -----
  p_treatment_legend <- ggplot(legend_df,aes(x = Treatment,y = 1,fill = Treatment) ) +
    geom_tile() +
    scale_fill_manual(values = Treatment_colors[c("Ctrl", "TAM")],name = "Treatment") +   
    theme_void() +
    theme(legend.position = "bottom")
  
  ## adj.p legend -----
    p_pvalue_legend <- ggplot(legend_df,aes(x = adj.p,y = 1,fill = adj.p)) +
    geom_tile() +
    scale_fill_manual(values = p_color_list,name = "adj. p-value"  ) +
    theme_void() +
    theme(legend.position = "bottom")
  
 
  ## Extract ggplot legends -----
  get_legend <- function(plot) {ggplotGrob(plot)$grobs[ sapply( ggplotGrob(plot)$grobs, function(x) x$name) == "guide-box"][[1]]}
  
  
  legend_sex <- get_legend(p_sex_legend)
  legend_treatment <- get_legend(p_treatment_legend)
  legend_pvalue <- get_legend(p_pvalue_legend)
  legend_combined <- gridExtra::arrangeGrob(legend_sex,   legend_treatment,   legend_pvalue,  ncol = 1)
  
  ## Save legend -----
  ggsave(  file.path(proteom_output_pwd,paste0("Pathways/",folder, "/08_1_Heatmap_Legend.png") ),
    legend_combined,width = 6,height = 4,  dpi = 300,  bg = "white" )

}

# Call Heatmap Function -----
for (parent in PARAMETERS$Proteom$Pathway_parents) {
  filename_safe <- gsub("[^A-Za-z0-9]+", "_", parent)
  filename_safe <-substr(filename_safe,1, 15)
  # Significant proteins: Parent + alle Children
  make_reactome_heatmap(parent, filename=paste0("08_Heatmap_Reactome_",filename_safe,".png") ,mode = "SIG_proteins",folder="Parents")
  make_reactome_heatmap(parent, filename=paste0("08_Heatmap_Reactome_",filename_safe,".png") ,mode = "ALL_proteins",folder="Parents_all")
  
  }

for (parent in PARAMETERS$Proteom$Pathway_lowest) {
  filename_safe <- gsub("[^A-Za-z0-9]+", "_", parent)
  filename_safe <-substr(filename_safe,1, 15)
  # Significant proteins
  make_reactome_heatmap(parent, filename=paste0("08_Heatmap_Reactome_",filename_safe,".png") ,mode = "SIG_proteins",folder="Reduced")
  # ALL proteins
  make_reactome_heatmap(parent, filename=paste0("08_Heatmap_Reactome_",filename_safe,".png") ,mode = "ALL_proteins",folder="Reduced_all")
}

for (parent in PARAMETERS$Proteom$Children_of_parents) {
  filename_safe <- gsub("[^A-Za-z0-9]+", "_", parent)
  filename_safe <-substr(filename_safe,1, 15)
  # Significant proteins
  make_reactome_heatmap(parent, filename=paste0("08_Heatmap_Reactome_",filename_safe,".png") ,mode = "SIG_proteins",folder="Children")
  # All measured proteins
  make_reactome_heatmap(parent, filename=paste0("08_Heatmap_Reactome_",filename_safe,".png") ,mode = "ALL_proteins",folder="Children_all")
  
}
write.csv2(reactome_proteins,file.path(proteom_output_pwd,"Data/04_Reactome_all_pathway_proteins.csv"),  row.names = FALSE)
write.csv2(reactome_prots_in_my_prots, file.path(proteom_output_pwd,"Data/04_React_PWD_prots_in_myprots.csv"),row.names = FALSE)


