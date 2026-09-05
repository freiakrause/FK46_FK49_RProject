
library(dplyr)
library(httr)
library(jsonlite)
library(org.Mm.eg.db)
library(pheatmap)
library(pathview)
library(AnnotationDbi)
# Connect Pathway Diagramm with expression data -----
# #### Notizen zum weitermachen
# Kann ich Reactome files für alle 28 enriched pathway herunterladen? im moment habe ich nur einen anteil
# generall funktion checkn un säubern
##Prepare Exression data for upload -----
reactome_all <- read.csv2(file.path(stats_pwd, "Reactome_ALL_significant.csv"))

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

expr_data <- dplyr::bind_rows(expr_data_proteins, expr_data_compounds)

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
# # #######
# library(ReactomeContentService4R)
# library(dplyr)
# library(purrr)
# library(tidyr)
# library(igraph)
# library(ggraph)
# library(ggplot2)
# node_selection <- "Protein"   # "Protein", "Metabolite" oder "Both"
# 
# for (i in seq_len(min(10, length(pathway_ids)))) {
#   i <- 2
#   pw <- pathway_ids[i]
#   pw_name <- pathway_names$Description[pathway_names$ID==pw]
#   pw_label <- gsub(" \\[.*", "", pw_name)
#   pw_label <- gsub("[^[:alnum:] _/-]", "", pw_label)
#   
#   events <- getParticipants(pw, retrieval = "EventsInPathways")
#   reactions <- events[purrr::map_lgl(events, ~ is.list(.x) && identical(.x$className, "Reaction"))]
#   
#   reaction_data <- purrr::map_dfr(reactions, function(r) {
#     x <- tryCatch(getParticipants(r$stId, retrieval = "AllInstances"), error = function(e) NULL)
#     if (is.null(x) || nrow(x) == 0) return(NULL)
#     x %>% mutate(reaction_id = r$stId, reaction_name = r$displayName)
#   })
#   
#   # NODES
#   reaction_nodes <- reactions %>%
#     purrr::map_dfr(~ tibble(
#       id = paste0("R_", .x$stId),
#       label = .x$displayName,
#       type = "Reaction"
#     ))
#   
#   physical_nodes <- reaction_data %>%
#     transmute(
#       id = paste0("PE_", peDbId),
#       label = displayName,
#       type = case_when(
#         schemaClass == "SimpleEntity" ~ "Metabolite",
#         schemaClass %in% c("DefinedSet", "Complex") ~ "Protein",
#         TRUE ~ schemaClass
#       ),
#       peDbId = peDbId,
#       refEntities = refEntities
#     ) %>%
#     distinct(id, .keep_all = TRUE)
#   
#   # NUR GEWÄHLTE MOLEKÜLKLASSE
#   if (node_selection == "Protein") {
#     physical_nodes <- physical_nodes %>%
#       filter(type == "Protein")
#   }
#   
#   if (node_selection == "Metabolite") {
#     physical_nodes <- physical_nodes %>%
#       filter(type == "Metabolite")
#   }
#   
#   nodes <- bind_rows(reaction_nodes, physical_nodes)
#   
#   # EDGES
#   edges <- reaction_data %>%
#     filter(type %in% c("input", "output", "catalyst")) %>%
#     transmute(
#       from = case_when(
#         type == "input" ~ paste0("PE_", peDbId),
#         type == "output" ~ paste0("R_", reaction_id),
#         type == "catalyst" ~ paste0("PE_", peDbId)
#       ),
#       to = case_when(
#         type == "input" ~ paste0("R_", reaction_id),
#         type == "output" ~ paste0("PE_", peDbId),
#         type == "catalyst" ~ paste0("R_", reaction_id)
#       ),
#       relation = type
#     )
#   
#   # Nur Kanten behalten, deren Molekülknoten ausgewählt wurden
#   edges <- edges %>%
#     filter(from %in% nodes$id, to %in% nodes$id)
#   
#   # REFERENCE TABLE
#   ref_table <- purrr::map2_dfr(
#     physical_nodes$id,
#     physical_nodes$refEntities,
#     function(x, y) {
#       if (is.null(y) || nrow(y) == 0) {
#         NULL
#       } else {
#         y %>% mutate(node_id = x)
#       }
#     }
#   )
#   
#   # PROTEIN FC
#   protein_fc <- expr_data_proteins %>%
#     mutate(identifier = toupper(ID)) %>%
#     dplyr::select(identifier, logFC) %>%
#     dplyr::rename(log2FC_protein = logFC)
#   
#   nodes <- nodes %>%
#     left_join(
#       ref_table %>%
#         filter(schemaClass == "ReferenceGeneProduct") %>%
#         dplyr::select(node_id, identifier) %>%
#         left_join(protein_fc, by = "identifier") %>%
#         group_by(node_id) %>%
#         summarise(
#           log2FC_protein = dplyr::first(log2FC_protein[!is.na(log2FC_protein)]),
#           .groups = "drop"
#         ),
#       by = c("id" = "node_id")
#     )
#   
#   # METABOLITE FC
#   metab_fc <- abundance_data_metabos %>%
#     mutate(ChEBI = as.character(ChEBI)) %>%
#     filter(!is.na(ChEBI), ChEBI != "") %>%
#     dplyr::select(ChEBI, log2FC)
#   
#   nodes <- nodes %>%
#     left_join(
#       ref_table %>%
#         filter(schemaClass == "ReferenceMolecule") %>%
#         mutate(ChEBI = as.character(identifier)) %>%
#         dplyr::select(node_id, ChEBI) %>%
#         left_join(metab_fc, by = "ChEBI") %>%
#         group_by(node_id) %>%
#         summarise(
#           log2FC_metab = dplyr::first(log2FC[!is.na(log2FC)]),
#           .groups = "drop"
#         ),
#       by = c("id" = "node_id")
#     ) %>%
#     mutate(
#       log2FC = coalesce(log2FC_protein, log2FC_metab),
#       plot_label = gsub(" \\[.*", "", label),
#       plot_label = case_when(
#         type == "Reaction" ~ paste0("R: ", plot_label),
#         type %in% c("Protein", "CandidateSet", "EntityWithAccessionedSequence") ~
#           sub(" .*", "", plot_label),
#         TRUE ~ plot_label
#       )
#     ) %>%
#     dplyr::select(-log2FC_protein, -log2FC_metab)
#   
#   nodes <- nodes %>%
#     mutate(
#       in_data = !is.na(log2FC),
#       
#       node_shape = case_when(
#         type == "Reaction"   ~ 15L,
#         type == "Protein"    ~ 22L,
#         type == "Metabolite" ~ 21L,
#         TRUE                 ~ 20L
#       ),
#       
#       border_col = if_else(in_data, "black", "grey70"),
#       fill_val = log2FC
#     )
#   
#   # GRAPH
#   connected_nodes <- nodes %>%
#     filter(id %in% unique(c(edges$from, edges$to)))
#   
#   g <- igraph::graph_from_data_frame(
#     edges,
#     vertices = connected_nodes,
#     directed = TRUE
#   )
#   
#   # HIER explizit igraph verwenden
#   comp <- igraph::components(g)
#   
#   # PLOT
#   plot_component <- function(component_id) {
#     
#     g_sub <- igraph::induced_subgraph(
#       g,
#       vids = igraph::V(g)[comp$membership == component_id]
#     )
#     
#     ggraph(g_sub, layout = "stress") +
#       
#       # Kanten
#       geom_edge_link(
#         aes(linetype = relation),
#         linewidth = 0.2,
#         alpha = 0.5,
#         arrow = grid::arrow(
#           length = grid::unit(1, "mm"),
#           type = "closed"
#         ),
#         end_cap = circle(2, "mm")
#       ) +
#       
#       # Reactions
#       geom_node_point(
#         data = function(x) x %>% filter(type == "Reaction"),
#         shape = 15,
#         size = 3,
#         colour = "grey40"
#       ) +
#       
#       # Proteine / Metabolite
#       geom_node_point(
#         data = function(x) x %>% filter(type != "Reaction"),
#         aes(
#           shape = node_shape,
#           fill = fill_val,
#           colour = border_col
#         ),
#         size = 6,
#         stroke = 0.8
#       ) +
#       
#       scale_shape_identity() +
#       scale_colour_identity() +
#       
#       scale_fill_gradient2(
#         low = "blue",
#         mid = "white",
#         high = "red",
#         midpoint = 0,
#         na.value = "grey70",
#         name = "log2FC"
#       ) +
#       
#       # Labels
#       geom_node_text(
#         data = function(x) x %>% filter(type != "Reaction"),
#         aes(label = plot_label),
#         repel = TRUE,
#         size = 2.5,
#         colour = "grey20"
#       ) +
#       
#       geom_node_text(
#         data = function(x) x %>% filter(type == "Reaction"),
#         aes(label = plot_label),repel = TRUE, size = 2.2, colour = "grey40" ) +
#       
#       labs(title = pw_label) +
#       theme_void() +
#       theme(
#         legend.position = "right",
#         plot.title = element_text(size = 16,  face = "bold",    hjust = 0.5 )
#       )
#   }
#   
#   # Alle Komponenten speichern
#   for (j in sort(unique(comp$membership))) {
#     
#     p <- plot_component(j)
#     ggsave( file.path(proteom_output_pwd,
#         paste0( "/Pathways/Reactome_",i,"_",pw_label, "_",node_selection,"_Component_",j,".png")),
#       p,width = 25,  height = 25,units = "in", dpi = 100,bg = "white", limitsize = FALSE   )
#   }
# }