#####
# maybe i can to a combined heatmpt or also heatmpa of prots and metabs of the same pathwas enxt to each others with FC and p values or sth
# ich habe jetzte auch meine kegg pathways aus den metaboliten mit den reactome pathways verglichen, um später sagen zu können: die pathways die wir mit metbaoliten
# abdecken überschneiden sich nicht mit den pathwas die durch sig diff proteins enriched sind. oder sowas.
# und dann muss das ganze skript noch schick werden oder sogar auf meherer aufgetiet werden
rm(list=ls())
gc()
library(tidyverse)
library(KEGGREST)
library(org.Mm.eg.db)
library(ReactomePA)
library(VennDiagram)
source("FK49_Definitions.R")
#Set Pathways -----
proteom_output_pwd <- PATHS$proteomics$output
protein_stats_pwd <- file.path(proteom_output_pwd, "Statistics")
targeted_rawdata_pwd <- dirname(PATHS$metabolomics$rawdata)

#Load INput ----
Pathways_covered <- read.csv(file.path (PATHS$metabolomics$output,"CDHFD/targete_pathway_coverage/pathway_results.csv"))
Proteins <- read.csv2(file.path(protein_stats_pwd, "02_LIMMA_combined_stats.csv"))
tar <- readRDS(file.path(targeted_rawdata_pwd, "FK49_Analysis/01_RawData/FK49_metabolome_targeted_processed.rds"))
kegg_id_of_tar <- read.csv2(file.path(PATHS$metabolomics$rawdata,"KEGG_MetaboAnalyst.csv"))

# Find my kegg pathways from targeted metabolomics in mouse kegg -----
mouse_kegg_pathways <- as.data.frame(keggList("pathway", "mmu")) %>%
  tibble::rownames_to_column("ID") %>%
  dplyr::rename(pathway = 2) %>%
  dplyr::mutate(pathway = gsub(" - Mus musculus \\(house mouse\\)", "", pathway))

Pathways_covered$X[Pathways_covered$X == "Glycolysis or Gluconeogenesis"] <-  "Glycolysis / Gluconeogenesis"
my_kegg_pathways <- mouse_kegg_pathways%>%filter(pathway %in% Pathways_covered$X)

# Get all genes (Proteins) in my pathways -----
pathway_ids <- unique(my_kegg_pathways$ID)

pw <- unlist( lapply(seq(1, length(pathway_ids), by = 10), function(i) {keggGet(pathway_ids[i:min(i + 9, length(pathway_ids))])}),recursive = FALSE)

names(pw) <- pathway_ids

kegg_genes <- do.call(rbind, lapply(seq_along(pw), function(i) {
  genes <- pw[[i]]$GENE
  data.frame(
    pathway_id = names(pw)[i],
    kegg_gene = genes[seq(1, length(genes), 2)],
    gene_name = genes[seq(2, length(genes), 2)],
    stringsAsFactors = FALSE)
}))

kegg_genes <- kegg_genes %>%
  tidyr::extract(
    gene_name,
    into = c("Gene", "Description", "KO", "EC"),
    regex = "^([^;]+);\\s*(.*?)\\s*\\[KO:([^]]+)\\]\\s*\\[EC:([^]]+)\\]$",
    remove = TRUE ) %>%
  filter(!is.na(Gene), Gene != "")

kegg_genes <- left_join(kegg_genes,mouse_kegg_pathways,by = c("pathway_id" = "ID"))
my_prots_in_KEGG <- Proteins %>%
  filter(Genes %in% kegg_genes$Gene) %>%
  left_join(kegg_genes %>% dplyr::select(pathway_id, pathway, Gene), by = c("Genes" = "Gene"),relationship = "many-to-many")
write.csv2(my_prots_in_KEGG,file.path(proteom_output_pwd,"Data/06_my_prots_in_KEGG.csv"))
saveRDS(my_prots_in_KEGG,file.path(proteom_output_pwd,"Data/06_my_prots_in_KEGG.rds"))
# Prepare Protein Expression Data for violin plotting -----
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

# Violin Plots of Proteins in targeted KEGG Pathways
print_violins <- function(pathway, mode= "ALL_proteins") {
 
   print(paste("PATHWAY:", pathway))
  proteins_to_plot <- my_prots_in_KEGG %>%
    dplyr::filter(.data$pathway == .env$pathway) %>%
    dplyr::pull(Genes) %>%
    unique()
    
  
  print(paste("N proteins:", length(proteins_to_plot)))
  print(head(proteins_to_plot))
  
  plot_df <- protein_long %>%
    dplyr::filter(Genes %in% proteins_to_plot) %>%
    mutate(Genes = factor(Genes),
           Name = factor(Name))
  
  print(paste("N plot rows:", nrow(plot_df)))
  if(mode == "ALL_proteins"){plot_df <- plot_df%>%filter( !is.na(adj_pvalue_Treatment))}
  else if (mode == "SIG_proteins"){plot_df <- plot_df%>%filter( adj_pvalue_Treatment<0.05 & abs(logFC_Treatment) > 1)}
  else(print(" You did not specifiy how to filter proteins for plotting"))
  if (nrow(plot_df) == 0) return(NULL)
  
  n_facets <- dplyr::n_distinct(plot_df$Name) # Anzahl Spalten ungefähr quadratisch 
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
    labs(title = my_kegg_pathways$pathway[my_kegg_pathways$ID == pathway],
         x = NULL,   y = "Log2-normalized protein abundance"  ) +
    geom_text(  data = distinct(plot_df, Name, adj_pvalue_Treatment),
                aes(x = 1.5,  y = Inf,
                    label = case_when(
                      adj_pvalue_Treatment < 0.001 ~ "***",
                      adj_pvalue_Treatment < 0.01  ~ "**",
                      adj_pvalue_Treatment < 0.05  ~ "*",
                      TRUE ~ "ns" )),      
                vjust = 1.5,inherit.aes = FALSE )
  
  ggsave( file.path(proteom_output_pwd, paste0("Pathways_Metabolomics/01_Violin_",substr(make.names(pathway), 1, 15), ".png")),
          p,width = plot_width,  height = plot_height,dpi = 300 ,limitsize=FALSE)
}

# for (pathway in my_kegg_pathways$pathway) {
#   print_violins(pathway = pathway,mode= "SIG_proteins")
# }
for (pathway in my_kegg_pathways$pathway) {
  print_violins(pathway = pathway,mode= "ALL_proteins")
}

#Prepare Metabolite abundance data for violin plotting -----
metabolite_long <- as.data.frame(tar$log_values) %>%
  tibble::rownames_to_column("Animal") %>%
  tidyr::pivot_longer(
    cols = -Animal,
    names_to = "Metabolite",
    values_to = "MetaboliteValue" ) %>%
  left_join(tar$metadata %>%  dplyr::mutate(Animal = as.character(Animal)) %>%
                              dplyr::select(Animal, Sex, Treatment),  by = "Animal" )

# Get my metbolites in the KEGG Pathways targeted metabolomics covers -----
my_metabolites_in_KEGG <- lapply(seq_len(nrow(my_kegg_pathways)), function(i) {
  
  pathway_id <- my_kegg_pathways$ID[i]
  pathway_name <- my_kegg_pathways$pathway[i]
  
  pathway_data <- KEGGREST::keggGet(pathway_id)[[1]]
  
  if (is.null(pathway_data$COMPOUND)) return(NULL)
  
  compound_ids <- names(pathway_data$COMPOUND)
  
  data.frame(
    pathway = pathway_name,
    pathway_ID = pathway_id,
    KEGG = compound_ids,
    stringsAsFactors = FALSE
  )
}) %>%dplyr::bind_rows()

my_metabolites_in_KEGG <- my_metabolites_in_KEGG %>%
  dplyr::inner_join(kegg_id_of_tar %>%dplyr::select(Query, KEGG) %>%dplyr::rename(Metabolite = Query),by = "KEGG" ) %>%
  dplyr::distinct()
write.csv2(my_metabolites_in_KEGG,file.path(proteom_output_pwd,"Data/07_my_metabolites_in_KEGG.csv"))
saveRDS(my_metabolites_in_KEGG,file.path(proteom_output_pwd,"Data/07_my_metabolites_in_KEGG.rds"))
# Violin Plots of all my metaboites in the targeted KEGG Pathways -----
print_metabolite_violins <- function(pathway, mode = "ALL_metabolites") {
  
  print(paste("PATHWAY:", pathway))
  
  metabolites_to_plot <- my_metabolites_in_KEGG %>%
    dplyr::filter(.data$pathway == .env$pathway) %>%
    dplyr::pull(Metabolite) %>%
    unique()
  
  print(paste("N metabolites:", length(metabolites_to_plot)))
  print(metabolites_to_plot)
  
  plot_df <- metabolite_long %>%
    dplyr::filter(Metabolite %in% metabolites_to_plot) %>%
    dplyr::mutate(
      Metabolite = factor(Metabolite),
      Animal = factor(Animal)
    )
  
  print(paste("N plot rows:", nrow(plot_df)))
  
  if (nrow(plot_df) == 0) return(NULL)
  
  n_facets <- dplyr::n_distinct(plot_df$Metabolite)
  ncol_plot <- ceiling(sqrt(n_facets))
  nrow_plot <- ceiling(n_facets / ncol_plot)
  
  facet_width <- 2.2
  facet_height <- 2.0
  
  plot_width <- ncol_plot * facet_width
  plot_height <- nrow_plot * facet_height
  
  p <- ggplot(
    plot_df,
    aes(
      x = Treatment,
      y = MetaboliteValue,
      fill = Treatment
    )
  ) +
    geom_violin(
      trim = FALSE,
      alpha = 0.4
    ) +
    geom_point(
      position = position_dodge2(width = 0.5),
      size = 1.5,
      aes(shape = Sex),
      alpha = 0.5
    ) +
    geom_boxplot(
      width = 0.15,
      outlier.shape = NA,
      alpha = 0.5
    ) +
    scale_fill_manual(
      values = Treatment_colors[c("Ctrl", "TAM")]
    ) +
    scale_shape_manual(
      values = Sex_shape
    ) +
    facet_wrap(
      ~ Metabolite,
      scales = "free_y",
      ncol = ncol_plot
    ) +
    theme_classic() +
    theme(
      text = element_text(size = 11),
      strip.text = element_text(size = 10),
      axis.text = element_text(size = 8),
      axis.title = element_text(size = 10)
    ) +
    labs(
      title = my_kegg_pathways$pathway[
        my_kegg_pathways$pathway == pathway
      ],
      x = NULL,
      y = "Log2-normalized metabolite abundance"
    )
  
  ggsave(
    file.path(
      proteom_output_pwd,
      paste0(
        "Pathways_Metabolomics/02_Metabolites_",
        substr(make.names(pathway), 1, 15),
        ".png"
      )
    ),
    p,
    width = plot_width,
    height = plot_height,
    dpi = 300,
    limitsize = FALSE
  )
}
for (pathway in my_kegg_pathways$pathway) {
  print_metabolite_violins(pathway = pathway,mode= "ALL_metabolites")
}
# Compare KEGG and Reactome pathways -----
# Compare the KEGG pathways that where targeted by targeted Metbolomoics with reactome pathways (by enrichment)
# Ran already, takes a lot of time
# run_Reactome <- function(genes) {
#   mapped <- bitr(genes, fromType = "SYMBOL", toType = "ENTREZID",
#                  OrgDb = org.Mm.eg.db)
#   enrichPathway(
#     gene = mapped$ENTREZID,
#     universe = mapped_background$ENTREZID,
#     organism = "mouse",
#     pAdjustMethod = "BH",
#     pvalueCutoff = 0.05,
#     qvalueCutoff = 0.05,
#     readable = TRUE
#   )
# }
# reactome_all  <- run_Reactome(my_prots_in_KEGG)
# reactome_all_df <- as.data.frame(reactome_all) %>%arrange(p.adjust)



# kegg_to_reactome <- lapply(seq_len(nrow(my_kegg_pathways)), function(i) {
#   
#   pathway_id <- my_kegg_pathways$ID[i]
#   pathway_name <- my_kegg_pathways$pathway[i]
#   
#   # KEGG pathway abrufen
#   pathway_data <- KEGGREST::keggGet(pathway_id)[[1]]
#   
#   if (is.null(pathway_data$GENE)) return(NULL)
#   
#   # GENE enthält abwechselnd Entrez-ID und Beschreibung
#   entrez_genes <- pathway_data$GENE[seq(1, length(pathway_data$GENE), by = 2)]
#   entrez_genes <- unique(entrez_genes)
#   
#   if (length(entrez_genes) < 2) return(NULL)
#   
#   # Reactome enrichment
#   reactome_result <- ReactomePA::enrichPathway(
#     gene = entrez_genes,
#     organism = "mouse",
#     pvalueCutoff = 1,
#     qvalueCutoff = 1,
#     minGSSize = 5,
#     maxGSSize = 5000,
#     readable = FALSE
#   )
#   
#   if (nrow(as.data.frame(reactome_result)) == 0) return(NULL)
#   
#   result <- as.data.frame(reactome_result)
#   
#   result %>%
#     dplyr::mutate(KEGG_ID = pathway_id,
#                   KEGG_Pathway = pathway_name ) %>%
#     dplyr::select(
#       KEGG_ID,
#       KEGG_Pathway,
#       Reactome_ID = ID,
#       Reactome_Pathway = Description,
#       GeneRatio,
#       BgRatio,
#       pvalue,
#       p.adjust,
#       qvalue,
#       geneID,
#       Count
#     )
#   
# }) %>%
#   dplyr::bind_rows()
# kegg_to_reactome_best <- kegg_to_reactome %>%
#   dplyr::group_by(KEGG_ID, KEGG_Pathway) %>%
#   dplyr::arrange(p.adjust, desc(Count), .by_group = TRUE) %>%
#   dplyr::slice_head(n = 3) %>%
#   dplyr::ungroup()

#kegg_to_reactome_best
#write.csv2(kegg_to_reactome,file.path(proteom_output_pwd,"KEGG_to_REACTOME.csv"))
#saveRDS(kegg_to_reactome,file.path(proteom_output_pwd,"KEGG_to_REACTOME.rds"))
#write.csv2(kegg_to_reactome_best,file.path(proteom_output_pwd,"KEGG_to_REACTOME_best.csv"))
#saveRDS(kegg_to_reactome_best,file.path(proteom_output_pwd,"KEGG_to_REACTOME_best.rds"))
#These are the Reactome Pathways that are enriched by the KEGG IDs -----
KEGG_to_REACTOME <- readRDS(file.path(proteom_output_pwd,"KEGG_to_REACTOME.rds"))
KEGG_to_REACTOME_best <- readRDS(file.path(proteom_output_pwd,"KEGG_to_REACTOME_best.rds"))
KEGG_to_REACTOME_best2 <- KEGG_to_REACTOME %>%
    dplyr::group_by(KEGG_ID, KEGG_Pathway) %>%
    dplyr::arrange(p.adjust, desc(Count), .by_group = TRUE) %>%
    dplyr::slice_head(n = 1) %>%
    dplyr::ungroup()
enriched_reactome_byprots <- read.csv2(file.path(protein_stats_pwd, "04_ORA_Reactome_ALL_significant.csv"))
enriched_KEGG_byprots <- read.csv2(file.path(protein_stats_pwd, "04_ORA_KEGG_ALL_significant.csv"))

# Venn Diagram -----
# Pathways that are enriched wenn taking KEGG pathway participants and usind enrich with their ids on reactome
venn.diagram(
  x = list(KEGG_to_REACTOME_best$Reactome_Pathway, enriched_reactome_byprots$Description),
  category.names = c("Pathways enriched  \n  by Metabolites" , "Pathways enriched  \n by Proteins " ),
  filename = file.path(proteom_output_pwd, "Pathways_metabolomics/03_Reactomepwds.png"),
  output=FALSE,
  
  # Output features
  imagetype="png" ,
  height = 550 , 
  width = 550, 
  resolution = 300,
  compression = "lzw",
  
  # Circles
  lwd = 2,
  lty = 'blank',
  fill = c("#9ECAE1","#80CBC4"),
  
  # Numbers
  cex = 0.25,
  fontface = "bold",
  fontfamily = "sans",
  
  # Set names
  cat.cex = 0.25,
  cat.fontface = "bold",
  cat.default.pos = "outer",
  cat.pos = c(-27, 27),
  cat.dist = c(0.045, 0.045) ,
  cat.fontfamily = "sans"
)


venn.diagram(
  x = list(my_kegg_pathways$ID, enriched_KEGG_byprots$ID),
  category.names = c("Pathways enriched  \n  by Metabolites" , "Pathways enriched  \n by Proteins " ),
  filename = file.path(proteom_output_pwd, "Pathways_metabolomics/03_KEGG.png"),
  output=FALSE,
  
  # Output features
  imagetype="png" ,
  height = 550 , 
  width = 550, 
  resolution = 300,
  compression = "lzw",
  
  # Circles
  lwd = 2,
  lty = 'blank',
  fill = c("#9ECAE1","#80CBC4"),
  
  # Numbers
  cex = 0.25,
  fontface = "bold",
  fontfamily = "sans",
  ext.pos = 0,
  
  # Set names
  cat.cex = 0.25,
  cat.fontface = "bold",
  cat.default.pos = "outer",
  cat.pos = c(-27, 27),
  cat.dist = c(0.045, 0.045) ,
  cat.fontfamily = "sans"
)






