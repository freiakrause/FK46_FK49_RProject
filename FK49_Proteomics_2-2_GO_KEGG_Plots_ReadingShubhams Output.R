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
library(pathview)
library(AnnotationDbi)
source("FK49_Definitions.R")

proteom_input_pwd  <- PATHS$proteomics$input
proteom_output_pwd <- PATHS$proteomics$output
stats_pwd <- file.path(proteom_output_pwd, "Statistics")

# Load protein statistics
Proteins <- read.csv2(file.path(stats_pwd, "FK49_Proteomics_statistics.csv"))
reactome_proteins <- read.csv2(file.path(proteom_output_pwd, "Data/Reactome_Top10pwd_proteins.csv"))
reactome_compounds <- read.csv2(file.path(proteom_output_pwd, "Data/Reactome_Top10pwd_compounds.csv"))

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

# ego_up <- read.csv2(file.path(stats_pwd, "GO_ALL_up.csv"))
# ego_up_simplified <- read.csv2(file.path(stats_pwd, "GO_ALL_up_simplified.csv"))
# 
# ego_down <- read.csv2(file.path(stats_pwd, "GO_ALL_down.csv"))
# ego_down_simplified <- read.csv2(file.path(stats_pwd, "GO_ALL_down_simplified.csv"))

ego_all <- read.csv2(file.path(stats_pwd, "GO_ALL_significant.csv"))
ego_all_simplified <- read.csv2(file.path(stats_pwd, "GO_ALL_significant_simplified.csv"))

# kegg_up <- read.csv2(file.path(stats_pwd, "KEGG_up.csv"))
# kegg_down <- read.csv2(file.path(stats_pwd, "KEGG_down.csv"))
kegg_all <- read.csv2(file.path(stats_pwd, "KEGG_ALL_significant.csv"))

# reactome_up <- read.csv2(file.path(stats_pwd, "Reactome_up.csv"))
# reactome_down <- read.csv2(file.path(stats_pwd, "Reactome_down.csv"))
reactome_all <- read.csv2(file.path(stats_pwd, "Reactome_ALL_significant.csv"))


str(reactome_all)
head(reactome_all)
# Chord plot function ---------------------------------------------------------
save_chord <- function(enrich, filename, n = 10,
                       pathway_parents = NULL,
                       children_of_parents = NULL,
                       pathway_lowest = NULL,
                       mode = "all") {
  
  if (is.null(enrich) || nrow(as.data.frame(enrich)) == 0) return(NULL)
  
  df <- as.data.frame(enrich) %>%
    filter(p.adjust < 0.05) %>%
    arrange(p.adjust)
  
  # Reactome hierarchy filtering
  if (mode == "parents" && !is.null(pathway_parents)) {
    df <- df %>%
      filter(Description %in% pathway_parents)
    
  } else if (mode == "children" && !is.null(children_of_parents)) {
    df <- df %>%
      filter(Description %in% children_of_parents)
    
  } else if (mode == "reduced" && !is.null(pathway_lowest)) {
    df <- df %>%
      filter(Description %in% pathway_lowest)
  }
  
  df <- df %>%
    slice_head(n = n)
  
  if (nrow(df) == 0) return(NULL)
  
  # Pathway -> adjusted p-value für Beschriftung
  pathway_p <- df %>%
    dplyr::select(Description, p.adjust)
  
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
    paste0(file.path(proteom_output_pwd, "Plots/", filename)),
    width = 13, height = 13, units = "in", res = 300
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
      
      sector.name <- get.cell.meta.data("sector.index")
      
      if (sector.name %in% pathway_p$Description) {
        
        pval <- pathway_p$p.adjust[
          pathway_p$Description == sector.name
        ]
        
        print(paste0(
          stringr::str_wrap(sector.name, width = 25),
          "\nFDR = ",
          format.pval(pval, digits = 2, eps = 0.001)
        ))
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
}

# GO chord plots --------------------------------------------------------------
# save_chord(ego_up,   "06_Chord_GO_up.png")
# save_chord(ego_down, "06_Chord_GO_down.png")
save_chord(
  ego_all,
  "06_Chord_GO_ALL_significant_top10.png",
  n = 10
)


# KEGG chord plots ------------------------------------------------------------
# save_chord(kegg_up,   "06_Chord_KEGG_up.png")
# save_chord(kegg_down, "06_Chord_KEGG_down.png")
save_chord(
  kegg_all,
  "06_Chord_KEGG_ALL_significant_top10.png",
  n = 10
)


# Reactome chord plots --------------------------------------------------------
# save_chord(reactome_up,   "06_Chord_Reactome_up.png")
# save_chord(reactome_down, "06_Chord_Reactome_down.png")

save_chord(
  reactome_all,
  "06_Chord_Reactome_ALL_significant_topall.png",
  n = 28
)

save_chord(
  reactome_all,
  "06_Reactome_All_parents_chord.png",
  pathway_parents = PARAMETERS$Proteom$Pathway_parents,
  mode = "parents",
  n = 28
)

save_chord(
  reactome_all,
  "06_Reactome_All_children_chord.png",
  children_of_parents = PARAMETERS$Proteom$Children_of_parents,
  mode = "children",
  n = 28
)

save_chord(
  reactome_all,
  "06_Reactome_All_children_chord_top10.png",
  children_of_parents = PARAMETERS$Proteom$Children_of_parents,
  mode = "children",
  n = 10
)

save_chord(
  reactome_all,
  "06_Reactome_All_reduced_chord.png",
  pathway_lowest = PARAMETERS$Proteom$Pathway_lowest,
  mode = "reduced",
  n = 28
)
# Reactome proteins -----------------------------------------------------------
reactome_df <- as.data.frame(reactome_all) %>%
  filter(p.adjust < 0.05) %>%
  arrange(p.adjust) %>%
  slice_head(n = 10) %>%
  dplyr::select(ID,Description, geneID) %>%
  separate_rows(geneID, sep = "/")
write.csv2(reactome_df,file.path(proteom_output_pwd,"Data/Reactome_top10_df.csv"),  row.names = FALSE)


# Protein abundances in long format ------------------------------------------
protein_long <- Proteins %>%
  dplyr::select(Genes,adj_pvalue_Treatment, starts_with("F_"), starts_with("M_")) %>%
  pivot_longer(c(-Genes, -adj_pvalue_Treatment),  
                 names_to = "Sample",values_to = "ProteinValue") %>%
  separate(  Sample,into = c("Sex", "Treatment", "Replicate"),sep = "_" ) %>%
  dplyr::select(-Replicate) %>%
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
  
  ggsave( file.path(proteom_output_pwd,  paste0("Plots/Violin_Reactome_", make.names(reactom), ".png")),
    p,  width = 12,    height = 8,  dpi = 300 )
}
# ALLlittle less shittyyyy  ------
# #maybe reanem stuff so that it is not that ocnfusing
# mabe in clude p value of tested proteins as row annotation 
# maybe do violins of proteins
# how can i vizualise the pathyw itself?
# Reactome pathway heatmaps ---------------------------------------------------


# Keep ALL pathway proteins that are present in the complete protein matrix --
# Notes:
# reactome_df<- (signifinaclt diff )proteins of my matric that drive the enrichemnt of the rectome pwd (top10 pwd)
# reactome_proteins = all proteins in my top10 enriched reactome oathways downlodad from reactome
# reactome_prots_in_my_prots = overlap of reactom_proteins and my protein matrix
# pathway_id: pathway id of my significantly enrichet top 10 reactome pathways
# # pathway_names: pathway id of my significantly enrichet top 10 reactome pathways
# 
# reactome_prots_in_my_prots <- reactome_proteins %>%
#   filter(identifier %in% Proteins$Protein.Group) %>%
#   dplyr::select(ID, identifier) %>%
#   distinct()
# 
# pathway_ids <- unique(reactome_prots_in_my_prots$ID)
# pathway_names <- reactome_df %>%
#   dplyr::filter(ID %in% pathway_ids) %>%
#   dplyr::select(ID, Description)
# # Heatmap data ---------------------------------------------------------------
# heatmap_list <- list()
# 
# for (pw in pathway_ids) {
#   
#   ##Proteins of this pathway -----
#   ids_this_pathway <- reactome_prots_in_my_prots %>%
#     dplyr::filter(ID == pw) %>%
#     dplyr::pull(identifier)
#   
#   ## matrix for this proteins ------
#   mat <- Proteins %>%
#     dplyr::filter(Protein.Group %in% ids_this_pathway) %>%
#     dplyr::select(Genes, starts_with("F_"), starts_with("M_")) %>%
#     tibble::column_to_rownames("Genes") %>%
#     as.matrix()
#   ## Colannotation -----------------------------------------------------------
#   annotation_col <- data.frame( Treatment = ifelse(grepl("_EtOH_", colnames(mat)), "Ctrl","TAM" ),
#     Sex = ifelse( grepl("^F_", colnames(mat)),"female","male"),
#     row.names = colnames(mat))
#   annotation_col$Sex <- factor( annotation_col$Sex, levels = c("female", "male"))
#   annotation_col$Treatment <- factor( annotation_col$Treatment,  levels = c("Ctrl", "TAM"))
#   ## if only one protein dont do heatmap -----
#   if (nrow(mat) < 2) {
#     warning(paste("Pathway", pw, "hat weniger als 2 Proteine - übersprungen"))
#     next
#   }
#     title_txt <- pathway_names$Description[pathway_names$ID == pw]
#   
#   ##generated the heatmap -----
#   heatmap_list[[pw]] <- pheatmap(
#     mat,
#     scale = "row",
#     annotation_col = annotation_col,
#     # annotation_row = annotation_row,
#     cluster_rows = TRUE,
#     cluster_cols = TRUE,
#     border_color = NA,
#     fontsize_row = 8,
#     fontsize_col = 8,
#     show_colnames = FALSE,
#     annotation_colors = list( Sex = Sex_colors,Treatment = Treatment_colors[c("Ctrl", "TAM")])
#     )
# }
# for (pw in pathway_ids) {
#   if (!is.null(heatmap_list[[pw]])) {
#     title_txt <- pathway_names$Description[pathway_names$ID == pw]
#     filename_safe <- gsub("[^A-Za-z0-9]", "_", title_txt)  # Sonderzeichen raus für Dateinamen
#     
#     png(file.path(proteom_output_pwd, paste0("Plots/08_Heatmap_",filename_safe, ".png")),
#         width = 1200, height = 1000, res = 300)
#     print(heatmap_list[[pw]])
#     dev.off()
#   }
# }

# Reactome pathway heatmaps ---------------------------------------------------
str(reactome_all)
head(reactome_all)
# Significant Reactome pathways
reactome_sig <- as.data.frame(reactome_all) %>%
  dplyr::select(ID, Description) %>%
  distinct()
str(reactome_sig)
head(reactome_sig)
# Reactome proteins that are present in my protein matrix
reactome_prots_in_my_prots <- reactome_proteins %>%
  filter(identifier %in% Proteins$Protein.Group) %>%
  dplyr::select(ID, identifier) %>%
  distinct()
str(reactome_prots_in_my_prots)
head(reactome_prots_in_my_prots)
# Add pathway names
reactome_prots_in_my_prots <- reactome_prots_in_my_prots %>%
  left_join(reactome_sig, by = "ID") %>%
  filter(!is.na(Description))
str(reactome_prots_in_my_prots)
head(reactome_prots_in_my_prots)

# Function to generate heatmap -----------------------------------------------
make_reactome_heatmap <- function(pathways, filename, annotation_name = NULL) {
  
  pathway_proteins <- reactome_prots_in_my_prots %>%
    filter(Description %in% pathways)
  
  if (nrow(pathway_proteins) == 0) {
    warning("Keine Proteine für ", pathways)
    return(NULL)
  }
  
  # Proteins = union of all selected pathways
  protein_ids <- unique(pathway_proteins$identifier)
  
  mat <- Proteins %>%
    filter(Protein.Group %in% protein_ids) %>%
    dplyr::select(Genes, starts_with("F_"), starts_with("M_")) %>%
    tibble::column_to_rownames("Genes") %>%
    as.matrix()
  
  if (nrow(mat) < 2) {
    warning(annotation_name, " hat weniger als 2 Proteine - übersprungen")
    return(NULL)
  }
  
  # Row annotation --------------------------------------------------------------
  if (!is.null(annotation_name)) {
    
    annotation_row <- pathway_proteins %>%
      filter(identifier %in% Proteins$Protein.Group) %>%
      left_join(
        Proteins %>% dplyr::select(Protein.Group, Genes),
        by = c("identifier" = "Protein.Group")
      ) %>%
      filter(Genes %in% rownames(mat)) %>%
      group_by(Genes) %>%
      summarise(
        !!annotation_name := paste(unique(Description), collapse = "; "),
        .groups = "drop"
      ) %>%
      tibble::column_to_rownames("Genes")
    
    annotation_row <- annotation_row[rownames(mat), , drop = FALSE]
    
  } else {
    annotation_row <- NULL
  }
  
  # Column annotation
  annotation_col <- data.frame(
    Treatment = ifelse(
      grepl("_EtOH_", colnames(mat)), "Ctrl", "TAM"
    ),
    Sex = ifelse(
      grepl("^F_", colnames(mat)), "female", "male"
    ),
    row.names = colnames(mat)
  )
  
  annotation_col$Sex <- factor(
    annotation_col$Sex,
    levels = c("female", "male")
  )
  
  annotation_col$Treatment <- factor(
    annotation_col$Treatment,
    levels = c("Ctrl", "TAM")
  )
  
  # Keep annotation in same order as matrix
  annotation_row <- annotation_row[rownames(mat), , drop = FALSE]
  
  # Heatmap
  p <- pheatmap(
    mat,
    scale = "row",
    annotation_col = annotation_col,
    annotation_row = annotation_row,
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    border_color = NA,
    fontsize_row = 7.9,
    fontsize_col = 7.9,
    show_colnames = FALSE,
    annotation_colors = list(
      Sex = Sex_colors,
      Treatment = Treatment_colors[c("Ctrl", "TAM")]
    )
  )
  n_rows <- nrow(mat)
  n_cols <- ncol(mat)
  ggsave( file.path(proteom_output_pwd,  paste0("Plots/", filename)),
          p,  width = 4 + 0.045 * n_rows,   height =1.5 + 0.08 * n_rows,   dpi= 300, bg = "white", limitsize = FALSE)
  
  print(p)

  
  invisible(p)
}


# Parent heatmap --------------------------------------------------------------
for (parent in PARAMETERS$Proteom$Pathway_parents) {
  
  filename_safe <- gsub(  "[^A-Za-z0-9]+",  "_",  parent )
  make_reactome_heatmap(  parent,  paste0("08_Heatmap_Reactome_Parent_",  filename_safe,  ".png") )
}

# Reduced pathway heatmap -----------------------------------------------------

make_reactome_heatmap(PARAMETERS$Proteom$Pathway_lowest,"08_Heatmap_Reactome_Reduced.png","Pathway")


# Individual child heatmaps ---------------------------------------------------

for (child in PARAMETERS$Proteom$Children_of_parents) {
  
  filename_safe <- gsub(  "[^A-Za-z0-9]+",  "_",  child )
  
  make_reactome_heatmap(  child,  paste0("08_Heatmap_Reactome_Child_",  filename_safe,  ".png") )
}
write.csv2( reactome_proteins,file.path(proteom_output_pwd,"Data/Reactome_top10_all_pathway_proteins.csv"),  row.names = FALSE)
write.csv2(reactome_prots_in_my_prots, file.path(proteom_output_pwd,"Data/Reactome_top10_proteins_in_myproteins.csv"),row.names = FALSE)
# Reactome Pathway Diagramms ind hässlich XX

# Connect Pathway Diagramm with expression data -----
##Prepare Exression data for upload -----

abundance_data_metabos <- read.csv2(file.path(PATHS$metabolomics$output,"CDHFD/FK49_metabolome_statistics.csv"))
KEGG_ID_targeted <- read.csv2(file.path(PATHS$metabolomics$rawdata,"KEGG_MetaboAnalyst.csv"))%>%mutate(ChEBI = NA,METLIN = NA)

untargeted_pwd<-paste0(parent,"/OneDrive - Universität Salzburg/AG_Tumorimmunologie - Dokumente/Data/Freia Krause/01_Experiments/FK49_CD-HFD_13wks/FK49_Analysis/02_GeneratedData")
ID_untargeted_pos <- read.csv2(file.path(untargeted_pwd,"Archive/untargetedLivMetabolome/CDHFD_a/CDHFD_a_p_trend_METABOANALYST.csv"))%>%  dplyr::rename(METLIN = MetLin)
ID_untargeted_neg <- read.csv2(file.path(untargeted_pwd,"Archive/untargetedLivMetabolome/CDHFD_a/CDHFD_a_n_trend_METABOANALYST.csv"))
ID_ALL <- dplyr::bind_rows(
  KEGG_ID_targeted %>% dplyr::select(Query, HMDB, PubChem, ChEBI, KEGG, METLIN),
  ID_untargeted_pos %>% dplyr::select(Query, HMDB, PubChem, ChEBI, KEGG, METLIN),
  ID_untargeted_neg %>% dplyr::select(Query, HMDB, PubChem, ChEBI, KEGG, METLIN)) %>%
  dplyr::distinct(Query, .keep_all = TRUE)

abundance_data_metabos <- abundance_data_metabos %>%
  dplyr::left_join(   ID_ALL,    by = c("Metabolite" = "Query")  )

expr_data_proteins <- Proteins %>%
  dplyr::filter(Protein.Group %in% reactome_prots_in_my_prots$identifier) %>%
  dplyr::select(ID = Protein.Group, logFC)

expr_data_compounds <- abundance_data_metabos %>%
  dplyr::filter(!is.na(KEGG), KEGG != "") %>%
  dplyr::select(ID = KEGG, logFC = log2FC)

expr_data <- dplyr::bind_rows(
  expr_data_proteins,
  expr_data_compounds
)
saveRDS(expr_data_proteins, file.path(proteom_output_pwd, "Data/Expr_Data_Proteins.rds"))
saveRDS(expr_data_compounds, file.path(proteom_output_pwd, "Data/Expr_Data_Compounds.rds"))
saveRDS(expr_data,           file.path(proteom_output_pwd, "Data/Expr_Data_ALL.rds"))# Header mit # davor, wie von Reactome gefordert
header_line <- paste0("#", paste(colnames(expr_data), collapse = "\t"))
data_lines <- apply(expr_data, 1, function(row) {paste(trimws(row), collapse = "\t")})
tsv_payload <- paste(c(header_line, data_lines), collapse = "\n")

## Do upload and load toke ---
analysis_res <- httr::POST(
  url = "https://reactome.org/AnalysisService/identifiers/?species=10090",
  httr::content_type("text/plain"), body = tsv_payload)

if (httr::status_code(analysis_res) != 200) {
  stop(paste("Analyse fehlgeschlagen - Status:", httr::status_code(analysis_res)))
}

analysis_json <- httr::content(analysis_res, as = "text", encoding = "UTF-8") %>%
  jsonlite::fromJSON(flatten = TRUE)

token <- analysis_json$summary$token
print(token)
# Wichtig: prüfen, wie viele deiner IDs NICHT gemapped werden konnten
not_found <- httr::GET(paste0("https://reactome.org/AnalysisService/token/", token, "/notFound")) %>%
  httr::content(as = "text", encoding = "UTF-8") %>%
  jsonlite::fromJSON()
analysis_json$summary
not_found
length(expr_data$ID)
print(paste(length(not_found), "IDs konnten nicht gemappt werden"))
##### Download diagram with overaly
download_reactome_diagram_with_overlay <- function(pathway_id, description, token, out_dir, 
                                                   ext = "png", quality = 10,
                                                   diagram_profile = "Modern",
                                                   analysis_profile = "Strosobar",
                                                   exp_column = NULL, ehld = FALSE,
                                                   flg_ids = NULL, flg_interactors = FALSE,coverage = FALSE) {
  
  url <- paste0("https://reactome.org/ContentService/exporter/diagram/",
                pathway_id, ".", ext,
                "?quality=", quality,
                "&token=", token,
                "&resource=TOTAL",
                "&diagramProfile=", diagram_profile,
                "&analysisProfile=", utils::URLencode(analysis_profile),
                "&ehld=", tolower(as.character(ehld)))
  
  if (!is.null(exp_column)) url <- paste0(url, "&expColumn=", exp_column)
  if (!is.null(coverage))   url <- paste0(url, "&coverage=", tolower(as.character(coverage)))
  
  if (!is.null(flg_ids)) {
    url <- paste0(url, "&flg=", paste(flg_ids, collapse = ","),
                  "&flgInteractors=", tolower(as.character(flg_interactors)))
  }

  res <- httr::GET(url)
  
  if (httr::status_code(res) != 200) {
    warning(paste("Fehler bei", pathway_id, "- Status:", httr::status_code(res)))
    return(NULL)
  }
  
  desc_safe <- gsub("[^A-Za-z0-9]+", "_", description)
  desc_safe <- gsub("_+$", "", desc_safe)
  
  out_file <- file.path(out_dir, paste0(desc_safe, "_", analysis_profile, ".", ext))
  writeBin(httr::content(res, "raw"), out_file)
  message(paste("Gespeichert:", out_file))
}
# Pathway-ID -> Description Mapping (aus top_reactome)
pathway_lookup <- reactome_df %>%
  dplyr::select(ID, Description) %>%
  dplyr::distinct()

# Schleife über alle Top-Pathways
for (pw in pathway_ids) {
  desc <- pathway_lookup$Description[pathway_lookup$ID == pw]
  download_reactome_diagram_with_overlay(
    pathway_id = pw, description = desc,
    token = token, out_dir = file.path(proteom_output_pwd, "Pathways"),
    flg_ids = expr_data$ID  # deine kompletten Proteine+Metabolite IDs
  )
}
# 
# #######
library(ReactomeContentService4R)
library(dplyr)
library(purrr)
library(tidyr)
library(igraph)
library(ggraph)
library(ggplot2)
node_selection <- "Protein"   # "Protein", "Metabolite" oder "Both"

for (i in seq_len(min(10, length(pathway_ids)))) {
  i <- 2
  pw <- pathway_ids[i]
  pw_name <- pathway_names$Description[pathway_names$ID==pw]
  pw_label <- gsub(" \\[.*", "", pw_name)
  pw_label <- gsub("[^[:alnum:] _/-]", "", pw_label)
  
  events <- getParticipants(pw, retrieval = "EventsInPathways")
  reactions <- events[purrr::map_lgl(events, ~ is.list(.x) && identical(.x$className, "Reaction"))]
  
  reaction_data <- purrr::map_dfr(reactions, function(r) {
    x <- tryCatch(getParticipants(r$stId, retrieval = "AllInstances"), error = function(e) NULL)
    if (is.null(x) || nrow(x) == 0) return(NULL)
    x %>% mutate(reaction_id = r$stId, reaction_name = r$displayName)
  })
  
  # NODES
  reaction_nodes <- reactions %>%
    purrr::map_dfr(~ tibble(
      id = paste0("R_", .x$stId),
      label = .x$displayName,
      type = "Reaction"
    ))
  
  physical_nodes <- reaction_data %>%
    transmute(
      id = paste0("PE_", peDbId),
      label = displayName,
      type = case_when(
        schemaClass == "SimpleEntity" ~ "Metabolite",
        schemaClass %in% c("DefinedSet", "Complex") ~ "Protein",
        TRUE ~ schemaClass
      ),
      peDbId = peDbId,
      refEntities = refEntities
    ) %>%
    distinct(id, .keep_all = TRUE)
  
  # NUR GEWÄHLTE MOLEKÜLKLASSE
  if (node_selection == "Protein") {
    physical_nodes <- physical_nodes %>%
      filter(type == "Protein")
  }
  
  if (node_selection == "Metabolite") {
    physical_nodes <- physical_nodes %>%
      filter(type == "Metabolite")
  }
  
  nodes <- bind_rows(reaction_nodes, physical_nodes)
  
  # EDGES
  edges <- reaction_data %>%
    filter(type %in% c("input", "output", "catalyst")) %>%
    transmute(
      from = case_when(
        type == "input" ~ paste0("PE_", peDbId),
        type == "output" ~ paste0("R_", reaction_id),
        type == "catalyst" ~ paste0("PE_", peDbId)
      ),
      to = case_when(
        type == "input" ~ paste0("R_", reaction_id),
        type == "output" ~ paste0("PE_", peDbId),
        type == "catalyst" ~ paste0("R_", reaction_id)
      ),
      relation = type
    )
  
  # Nur Kanten behalten, deren Molekülknoten ausgewählt wurden
  edges <- edges %>%
    filter(from %in% nodes$id, to %in% nodes$id)
  
  # REFERENCE TABLE
  ref_table <- purrr::map2_dfr(
    physical_nodes$id,
    physical_nodes$refEntities,
    function(x, y) {
      if (is.null(y) || nrow(y) == 0) {
        NULL
      } else {
        y %>% mutate(node_id = x)
      }
    }
  )
  
  # PROTEIN FC
  protein_fc <- expr_data_proteins %>%
    mutate(identifier = toupper(ID)) %>%
    dplyr::select(identifier, logFC) %>%
    dplyr::rename(log2FC_protein = logFC)
  
  nodes <- nodes %>%
    left_join(
      ref_table %>%
        filter(schemaClass == "ReferenceGeneProduct") %>%
        dplyr::select(node_id, identifier) %>%
        left_join(protein_fc, by = "identifier") %>%
        group_by(node_id) %>%
        summarise(
          log2FC_protein = dplyr::first(log2FC_protein[!is.na(log2FC_protein)]),
          .groups = "drop"
        ),
      by = c("id" = "node_id")
    )
  
  # METABOLITE FC
  metab_fc <- abundance_data_metabos %>%
    mutate(ChEBI = as.character(ChEBI)) %>%
    filter(!is.na(ChEBI), ChEBI != "") %>%
    dplyr::select(ChEBI, log2FC)
  
  nodes <- nodes %>%
    left_join(
      ref_table %>%
        filter(schemaClass == "ReferenceMolecule") %>%
        mutate(ChEBI = as.character(identifier)) %>%
        dplyr::select(node_id, ChEBI) %>%
        left_join(metab_fc, by = "ChEBI") %>%
        group_by(node_id) %>%
        summarise(
          log2FC_metab = dplyr::first(log2FC[!is.na(log2FC)]),
          .groups = "drop"
        ),
      by = c("id" = "node_id")
    ) %>%
    mutate(
      log2FC = coalesce(log2FC_protein, log2FC_metab),
      plot_label = gsub(" \\[.*", "", label),
      plot_label = case_when(
        type == "Reaction" ~ paste0("R: ", plot_label),
        type %in% c("Protein", "CandidateSet", "EntityWithAccessionedSequence") ~
          sub(" .*", "", plot_label),
        TRUE ~ plot_label
      )
    ) %>%
    dplyr::select(-log2FC_protein, -log2FC_metab)
  
  nodes <- nodes %>%
    mutate(
      in_data = !is.na(log2FC),
      
      node_shape = case_when(
        type == "Reaction"   ~ 15L,
        type == "Protein"    ~ 22L,
        type == "Metabolite" ~ 21L,
        TRUE                 ~ 20L
      ),
      
      border_col = if_else(in_data, "black", "grey70"),
      fill_val = log2FC
    )
  
  # GRAPH
  connected_nodes <- nodes %>%
    filter(id %in% unique(c(edges$from, edges$to)))
  
  g <- igraph::graph_from_data_frame(
    edges,
    vertices = connected_nodes,
    directed = TRUE
  )
  
  # HIER explizit igraph verwenden
  comp <- igraph::components(g)
  
  # PLOT
  plot_component <- function(component_id) {
    
    g_sub <- igraph::induced_subgraph(
      g,
      vids = igraph::V(g)[comp$membership == component_id]
    )
    
    ggraph(g_sub, layout = "stress") +
      
      # Kanten
      geom_edge_link(
        aes(linetype = relation),
        linewidth = 0.2,
        alpha = 0.5,
        arrow = grid::arrow(
          length = grid::unit(1, "mm"),
          type = "closed"
        ),
        end_cap = circle(2, "mm")
      ) +
      
      # Reactions
      geom_node_point(
        data = function(x) x %>% filter(type == "Reaction"),
        shape = 15,
        size = 3,
        colour = "grey40"
      ) +
      
      # Proteine / Metabolite
      geom_node_point(
        data = function(x) x %>% filter(type != "Reaction"),
        aes(
          shape = node_shape,
          fill = fill_val,
          colour = border_col
        ),
        size = 6,
        stroke = 0.8
      ) +
      
      scale_shape_identity() +
      scale_colour_identity() +
      
      scale_fill_gradient2(
        low = "blue",
        mid = "white",
        high = "red",
        midpoint = 0,
        na.value = "grey70",
        name = "log2FC"
      ) +
      
      # Labels
      geom_node_text(
        data = function(x) x %>% filter(type != "Reaction"),
        aes(label = plot_label),
        repel = TRUE,
        size = 2.5,
        colour = "grey20"
      ) +
      
      geom_node_text(
        data = function(x) x %>% filter(type == "Reaction"),
        aes(label = plot_label),
        repel = TRUE,
        size = 2.2,
        colour = "grey40"
      ) +
      
      labs(title = pw_label) +
      
      theme_void() +
      
      theme(
        legend.position = "right",
        plot.title = element_text(
          size = 16,
          face = "bold",
          hjust = 0.5
        )
      )
  }
  
  # Alle Komponenten speichern
  for (j in sort(unique(comp$membership))) {
    
    p <- plot_component(j)
    
    ggsave(
      file.path(
        proteom_output_pwd,
        paste0(
          "/Pathways/Reactome_",
          i,
          "_",
          pw_label,
          "_",
          node_selection,
          "_Component_",
          j,
          ".png"
        )
      ),
      p,
      width = 25,
      height = 25,
      units = "in",
      dpi = 100,
      bg = "white",
      limitsize = FALSE
    )
  }
}
### Notiz zum weitermachen: Sind alle paorents pathway significant? sinfd alle child pathways
###  und alle reduced pathways significant?
### wie filtere ich im stats skript um alle proteinen für alle möglichen geplotteten pathways zu finden? 
### sollte ich dort evtl nicht das slicen machen? Ich sollte evtl auch verschieden csv speichern in denen
###  genua aufgelistet ist, welcher parent pathways welchen p value hat und wer siene chidlren sind 
###  genauso wie jedes childeren seinen parents und weteren childere haben sollte als übersicht
###  
###  Dann uach chekcne, dass die heatmaps für alles pathways proteien plotten.
###  