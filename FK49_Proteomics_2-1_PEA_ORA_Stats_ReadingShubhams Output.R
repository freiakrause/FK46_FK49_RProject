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
Proteins <- read.csv2(file.path(stats_pwd, "02_LIMMA_combined_stats.csv"))

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
mapped_background <- bitr(background, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db)


# GO enrichment ---------------------------------------------------------------
ego_all  <- run_GO(all_proteins)
ego_all_simplified  <- simplify(ego_all,  cutoff = 0.7, by = "p.adjust", select_fun = min)

# KEGG enrichment -------------------------------------------------------------
kegg_all  <- run_KEGG(all_proteins)
kegg_all <- setReadable(kegg_all, OrgDb = org.Mm.eg.db, keyType = "ENTREZID")

# Reactome enrichment ---------------------------------------------------------
reactome_all  <- run_Reactome(all_proteins)
reactome_all_df <- as.data.frame(reactome_all) %>%arrange(p.adjust)

hierarchy <- PARAMETERS$Proteom$hierarchy # hierachy defined in FK49_Definitions, from reactome images, former parent searches and AI 
                                          # it only includes root, and child parent realtions from pathways that are sig enriched
reactome_all_df_hierachy <- reactome_all_df %>% left_join(hierarchy %>% dplyr::select(Child, Parent), by = c("Description" = "Child") )%>%
    mutate(Is_parent_of_significant = Description %in% Parent)

write.csv2(hierarchy,file.path(proteom_output_pwd, "Data/02_Pathway_Hierarchy_used.csv"), row.names = FALSE)
# Save enrichment results -----------------------------------------------------
write.csv2(as.data.frame(ego_all),             file.path(stats_pwd, "04_ORA_GO_ALL_significant.csv"), row.names = FALSE)
write.csv2(as.data.frame(ego_all_simplified),  file.path(stats_pwd, "04_ORA_GO_ALL_significant_simplified.csv"), row.names = FALSE)
write.csv2(as.data.frame(kegg_all),    file.path(stats_pwd, "04_ORA_KEGG_ALL_significant.csv"), row.names = FALSE)
write.csv2(as.data.frame(reactome_all_df_hierachy),  file.path(stats_pwd, "04_ORA_Reactome_ALL_significant.csv"), row.names = FALSE)
# Dotplots --------------------------------------------------------------------

ggsave(file.path(proteom_output_pwd, "Plots/05_Dotplot_KEGG_ALL_significant.png"),
       dotplot(kegg_all, showCategory = 15), width = 8, height = 6, dpi = 300)

ggsave(file.path(proteom_output_pwd, "Plots/05_Dotplot_GO_ALL_significant.png"),
       dotplot(ego_all, showCategory = 15), width = 8, height = 6, dpi = 300)

ggsave(file.path(proteom_output_pwd, "Plots/05_Dotplot_GO_ALL_significant_simplified.png"),
       dotplot(ego_all_simplified, showCategory = 15), width = 8, height = 6, dpi = 300)


ggsave(file.path(proteom_output_pwd, "Plots/05_Dotplot_Reactome_ALL_significant.png"),
       dotplot(reactome_all, showCategory = 28), width = 8, height = 12, dpi = 300)

#Get ALL genes componentd of Reactome pathway -----------------------------------
# Enriched Reactome pathways ----------------------------------------------------

 top_reactome <- as.data.frame(reactome_all) %>%
   filter(p.adjust < 0.05) %>%
   arrange(p.adjust) %>%
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
   get_reactome_participants(as.data.frame(top_reactome)$ID[i], as.data.frame(top_reactome)$Description[i])}) %>%bind_rows()

reactome_participants <- reactome_participants%>%dplyr::select(dbId,displayName,stId,databaseName,identifier,moleculeType,name,geneName,ID)
reactome_proteins <- reactome_participants%>%filter(moleculeType == "Protein")
reactome_compounds <- reactome_participants%>%filter(moleculeType == "Chemical" | moleculeType =="ChemicalDrug")

saveRDS(reactome_participants, file.path(proteom_output_pwd,"Data/02_Reactome_allparticipants.rds"))
saveRDS(reactome_proteins, file.path(proteom_output_pwd,"Data/02_Reactome__proteins.rds"))
saveRDS(reactome_compounds, file.path(proteom_output_pwd,"Data/02_Reactome_compounds.rds"))

write.csv2(reactome_participants, file = file.path(proteom_output_pwd,"Data/02_Reactome_allparticipants.csv"))
write.csv2(reactome_proteins, file = file.path(proteom_output_pwd,"Data/02_Reactome_allproteins.csv"))
write.csv2(reactome_compounds, file = file.path(proteom_output_pwd,"Data/02_Reactome_allcompounds.csv"))




