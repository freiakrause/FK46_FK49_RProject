###############################################################################
# FK49_Proteomics_3_Pathway_Analysis.R
#
# GO, KEGG and Reactome enrichment analysis of significantly regulated proteins
###############################################################################

rm(list = ls())
gc()
library(dplyr)
library(clusterProfiler)
library(org.Mm.eg.db)
library(ReactomePA)
library(ggplot2)
library(httr)
library(jsonlite)
source("FK49_Definitions.R")

proteom_input_pwd  <- PATHS$proteomics$input
proteom_output_pwd <- PATHS$proteomics$output
stats_pwd <- file.path(proteom_output_pwd, "Statistics")

# Load protein statistics
Proteins <- read.csv2(file.path(stats_pwd, "FK49_Proteomics_statistics.csv"))

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


# Enrichment function ---------------------------------------------------------
run_GO <- function(genes) enrichGO(
  gene = genes, universe = background, OrgDb = org.Mm.eg.db,
  keyType = "SYMBOL", ont = "ALL", pAdjustMethod = "BH",
  pvalueCutoff = 0.05, qvalueCutoff = 0.05, readable = TRUE
)

run_KEGG <- function(genes) {
  mapped <- bitr(genes, fromType = "SYMBOL", toType = "ENTREZID",
                 OrgDb = org.Mm.eg.db)
  enrichKEGG(
    gene = mapped$ENTREZID,
    universe = mapped_background$ENTREZID,
    organism = "mmu",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05
  )
}

run_Reactome <- function(genes) {
  mapped <- bitr(genes, fromType = "SYMBOL", toType = "ENTREZID",
                 OrgDb = org.Mm.eg.db)
  enrichPathway(
    gene = mapped$ENTREZID,
    universe = mapped_background$ENTREZID,
    organism = "mouse",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.05,
    readable = TRUE
  )
}


# KEGG background -------------------------------------------------------------
mapped_background <- bitr(
  background, fromType = "SYMBOL", toType = "ENTREZID",
  OrgDb = org.Mm.eg.db
)


# GO enrichment ---------------------------------------------------------------
ego_up   <- run_GO(up_proteins)
ego_down <- run_GO(down_proteins)
ego_all  <- run_GO(all_proteins)

ego_up_simplified   <- simplify(ego_up,   cutoff = 0.7, by = "p.adjust", select_fun = min)
ego_down_simplified <- simplify(ego_down, cutoff = 0.7, by = "p.adjust", select_fun = min)
ego_all_simplified  <- simplify(ego_all,  cutoff = 0.7, by = "p.adjust", select_fun = min)


# KEGG enrichment -------------------------------------------------------------
kegg_up   <- run_KEGG(up_proteins)
kegg_down <- run_KEGG(down_proteins)
kegg_all  <- run_KEGG(all_proteins)

kegg_up <- setReadable(kegg_up, OrgDb = org.Mm.eg.db, keyType = "ENTREZID")
kegg_down <- setReadable(kegg_down, OrgDb = org.Mm.eg.db, keyType = "ENTREZID")
kegg_all <- setReadable(kegg_all, OrgDb = org.Mm.eg.db, keyType = "ENTREZID")


# Reactome enrichment ---------------------------------------------------------
reactome_up   <- run_Reactome(up_proteins)
reactome_down <- run_Reactome(down_proteins)
reactome_all  <- run_Reactome(all_proteins)


reactome_all_df <- as.data.frame(reactome_all) %>%
  filter(p.adjust < 0.05) %>%
  arrange(p.adjust)
# 
# get_reactome_parents <- function(pathway_id) {
#   url <- paste0(  "https://reactome.org/ContentService/data/pathway/",
#                   pathway_id,   "/containedEvents")
#   
#   res <- GET(url)
#   
#   if (status_code(res) != 200) return(NULL)
#   
#   fromJSON(content(res, "text", encoding = "UTF-8"))
# }
# reactome_parents <- lapply(seq_len(nrow(reactome_all_df)), function(i) {
#   get_reactome_parents(reactome_all_df$ID[i])}) 
# 
# 
# hierarchy <- lapply(reactome_parents, function(x) {
#   
#   do.call(rbind, lapply(x, function(event) {
#     
#     if (!is.list(event) ||
#         is.null(event$eventOf) ||
#         event$schemaClass != "Pathway") {
#       return(NULL)
#     }
#     
#     parent <- event$eventOf
#     
#     if (nrow(parent) == 0 || parent$schemaClass != "Pathway") {
#       return(NULL)
#     }
#     
#     data.frame(
#       Child = event$displayName,
#       Child_ID = event$stId,
#       Parent = parent$displayName,
#       Parent_ID = parent$stId,
#       stringsAsFactors = FALSE
#     )
#   }))
# })
# 
# hierarchy <- dplyr::bind_rows(hierarchy) %>%
#   distinct()
# reactome_all_df_hierachy <- reactome_all_df %>%
#   left_join(  hierarchy %>%   dplyr::select(Child, Parent),
#               by = c("Description" = "Child") )%>%
#   mutate( Is_parent_of_significant = Description %in% Parent)
get_reactome_parents <- function(pathway_id) {
  
  url <- paste0(
    "https://reactome.org/ContentService/data/event/",
    pathway_id,
    "/ancestors"
  )
  
  res <- GET(url)
  
  if (status_code(res) != 200) return(NULL)
  
  fromJSON(
    content(res, "text", encoding = "UTF-8"),
    simplifyVector = FALSE
  )
}


reactome_parents <- lapply(
  reactome_all_df$ID,
  get_reactome_parents
)


hierarchy <- lapply(seq_along(reactome_parents), function(i) {
  
  paths <- reactome_parents[[i]]
  
  if (is.null(paths) || length(paths) == 0) {
    return(NULL)
  }
  
  child_id <- reactome_all_df$ID[i]
  child_name <- reactome_all_df$Description[i]
  
  do.call(rbind, lapply(paths, function(path) {
    
    if (!is.list(path) || length(path) == 0) {
      return(NULL)
    }
    
    # Nur Pathways behalten
    path <- path[
      vapply(
        path,
        function(x) {
          is.list(x) &&
            !is.null(x$schemaClass) &&
            x$schemaClass == "Pathway"
        },
        logical(1)
      )
    ]
    
    if (length(path) == 0) {
      return(NULL)
    }
    
    path_df <- do.call(rbind, lapply(path, function(x) {
      
      data.frame(
        Child = child_name,
        Child_ID = child_id,
        Parent = x$displayName,
        Parent_ID = x$stId,
        stringsAsFactors = FALSE
      )
    }))
    
    path_df
  }))
  
})


hierarchy <- dplyr::bind_rows(hierarchy) %>%
  dplyr::filter(Child_ID != Parent_ID) %>%
  dplyr::distinct()
str(hierarchy)
head(hierarchy)
hierarchy %>%
  arrange(Child, Parent) %>%
  print()
hierarchy %>%
  filter(
    Child %in% c(
      "Phase I - Functionalization of compounds",
      "Cytochrome P450 - arranged by substrate type",
      "Phase II - Conjugation of compounds",
      "Glutathione conjugation",
      "Metabolism of steroids",
      "Bile acid and bile salt metabolism",
      "Recycling of bile acids and salts",
      "Fibrin formation",
      "Biosynthesis of DHA-derived SPMs"
    )
  ) %>%
  arrange(Child, Parent)

reactome_all_df_hierachy <- reactome_all_df %>%
    left_join(  hierarchy %>%   dplyr::select(Child, Parent),
                by = c("Description" = "Child") )%>%
    mutate( Is_parent_of_significant = Description %in% Parent)
# Save enrichment results -----------------------------------------------------
write.csv2(as.data.frame(ego_up),
           file.path(stats_pwd, "GO_ALL_up.csv"), row.names = FALSE)
write.csv2(as.data.frame(ego_up_simplified),
           file.path(stats_pwd, "GO_ALL_up_simplified.csv"), row.names = FALSE)

write.csv2(as.data.frame(ego_down),
           file.path(stats_pwd, "GO_ALL_down.csv"), row.names = FALSE)
write.csv2(as.data.frame(ego_down_simplified),
           file.path(stats_pwd, "GO_ALL_down_simplified.csv"), row.names = FALSE)

write.csv2(as.data.frame(ego_all),
           file.path(stats_pwd, "GO_ALL_significant.csv"), row.names = FALSE)
write.csv2(as.data.frame(ego_all_simplified),
           file.path(stats_pwd, "GO_ALL_significant_simplified.csv"), row.names = FALSE)

write.csv2(as.data.frame(kegg_up),
           file.path(stats_pwd, "KEGG_up.csv"), row.names = FALSE)
write.csv2(as.data.frame(kegg_down),
           file.path(stats_pwd, "KEGG_down.csv"), row.names = FALSE)
write.csv2(as.data.frame(kegg_all),
           file.path(stats_pwd, "KEGG_ALL_significant.csv"), row.names = FALSE)

write.csv2(as.data.frame(reactome_up),
           file.path(stats_pwd, "Reactome_up.csv"), row.names = FALSE)
write.csv2(as.data.frame(reactome_down),
           file.path(stats_pwd, "Reactome_down.csv"), row.names = FALSE)
write.csv2(as.data.frame(reactome_all_df_hierachy),
           file.path(stats_pwd, "Reactome_ALL_significant.csv"), row.names = FALSE)
# Dotplots --------------------------------------------------------------------
# ggsave(file.path(proteom_output_pwd, "Plots/05_Dotplot_KEGG_up.png"),
#        dotplot(kegg_up, showCategory = 15), width = 8, height = 6, dpi = 300)
# 
# ggsave(file.path(proteom_output_pwd, "Plots/05_Dotplot_KEGG_down.png"),
#        dotplot(kegg_down, showCategory = 15), width = 8, height = 6, dpi = 300)

ggsave(file.path(proteom_output_pwd, "Plots/05_Dotplot_KEGG_ALL_significant.png"),
       dotplot(kegg_all, showCategory = 15), width = 8, height = 6, dpi = 300)

# ggsave(file.path(proteom_output_pwd, "Plots/05_Dotplot_GO_up.png"),
#        dotplot(ego_up, showCategory = 15), width = 8, height = 6, dpi = 300)
# 
# ggsave(file.path(proteom_output_pwd, "Plots/05_Dotplot_GO_up_simplified.png"),
#        dotplot(ego_up_simplified, showCategory = 15), width = 8, height = 6, dpi = 300)
# 
# ggsave(file.path(proteom_output_pwd, "Plots/05_Dotplot_GO_down.png"),
#        dotplot(ego_down, showCategory = 15), width = 8, height = 6, dpi = 300)
# 
# ggsave(file.path(proteom_output_pwd, "Plots/05_Dotplot_GO_down_simplified.png"),
#        dotplot(ego_down_simplified, showCategory = 15), width = 8, height = 6, dpi = 300)

ggsave(file.path(proteom_output_pwd, "Plots/05_Dotplot_GO_ALL_significant.png"),
       dotplot(ego_all, showCategory = 15), width = 8, height = 6, dpi = 300)

ggsave(file.path(proteom_output_pwd, "Plots/05_Dotplot_GO_ALL_significant_simplified.png"),
       dotplot(ego_all_simplified, showCategory = 15), width = 8, height = 6, dpi = 300)

# ggsave(file.path(proteom_output_pwd, "Plots/05_Dotplot_Reactome_up.png"),
#        dotplot(reactome_up, showCategory = 15), width = 8, height = 6, dpi = 300)
# 
# ggsave(file.path(proteom_output_pwd, "Plots/05_Dotplot_Reactome_down.png"),
#        dotplot(reactome_down, showCategory = 15), width = 8, height = 6, dpi = 300)

ggsave(file.path(proteom_output_pwd, "Plots/05_Dotplot_Reactome_ALL_significant.png"),
       dotplot(reactome_all, showCategory = 15), width = 8, height = 6, dpi = 300)

#Get ALL genes componentd of Reactome pathway -----------------------------------
# Top 10 Reactome pathways ----------------------------------------------------
# only use top 10


top_reactome <- as.data.frame(reactome_all) %>%
  filter(p.adjust < 0.05) %>%
  arrange(p.adjust) %>%
  dplyr::slice_head(n = 10) %>%
  dplyr::select(ID, Description)

get_reactome_participants <- function(pathway_id, description) {
  url <- paste0("https://reactome.org/ContentService/data/participants/",pathway_id, "/referenceEntities")
  res <- GET(url)
  if (status_code(res) != 200) {
    warning(paste("Fehler bei", pathway_id, "- Status:", status_code(res)))
    return(NULL)
  }
  content_raw <- httr::content(res, as = "text", encoding = "UTF-8")
  df <- fromJSON(content_raw, flatten = TRUE)
  if (length(df) == 0 || nrow(df) == 0) return(NULL)
  df <- df %>%
    dplyr::mutate(dplyr::across(where(is.list), ~ sapply(., function(x) paste(unlist(x), collapse = "; "))))
  df %>%as_tibble() %>% mutate(ID = pathway_id, Description = description)
}

reactome_participants <- lapply(seq_len(nrow(top_reactome)), function(i) {
  get_reactome_participants(top_reactome$ID[i], top_reactome$Description[i])}) %>%bind_rows()

reactome_participants <- reactome_participants%>%dplyr::select(dbId,displayName,stId,databaseName,identifier,moleculeType,name,geneName,ID)
reactome_proteins <- reactome_participants%>%filter(moleculeType == "Protein")
reactome_compounds <- reactome_participants%>%filter(moleculeType == "Chemical" | moleculeType =="ChemicalDrug")
saveRDS(reactome_participants, file.path(proteom_output_pwd,"Data/Reactome_Top10pwd_allparticipants.rds"))
saveRDS(reactome_proteins, file.path(proteom_output_pwd,"Data/Reactome_Top10pwd_proteins.rds"))
saveRDS(reactome_compounds, file.path(proteom_output_pwd,"Data/Reactome_Top10pwd_compounds.rds"))

write.csv2(reactome_participants, file = file.path(proteom_output_pwd,"Data/Reactome_Top10pwd_allparticipants.csv"))
write.csv2(reactome_proteins, file = file.path(proteom_output_pwd,"Data/Reactome_Top10pwd_proteins.csv"))
write.csv2(reactome_compounds, file = file.path(proteom_output_pwd,"Data/Reactome_Top10pwd_compounds.csv"))




