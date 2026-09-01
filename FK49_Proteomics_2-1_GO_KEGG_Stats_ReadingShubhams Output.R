###############################################################################
# FK49_Proteomics_3_Pathway_Analysis.R
#
# GO, KEGG and Reactome enrichment analysis of significantly regulated proteins
###############################################################################

rm(list = ls())
gc()
ibrary(dplyr)
library(clusterProfiler)
library(org.Mm.eg.db)
library(ReactomePA)
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
write.csv2(as.data.frame(reactome_all),
           file.path(stats_pwd, "Reactome_ALL_significant.csv"), row.names = FALSE)
# Dotplots --------------------------------------------------------------------
ggsave(file.path(proteom_output_pwd, "Dotplot_KEGG_up.png"),
       dotplot(kegg_up, showCategory = 15), width = 8, height = 6, dpi = 300)

ggsave(file.path(proteom_output_pwd, "Dotplot_KEGG_down.png"),
       dotplot(kegg_down, showCategory = 15), width = 8, height = 6, dpi = 300)

ggsave(file.path(proteom_output_pwd, "Dotplot_KEGG_ALL_significant.png"),
       dotplot(kegg_all, showCategory = 15), width = 8, height = 6, dpi = 300)

ggsave(file.path(proteom_output_pwd, "Dotplot_GO_up.png"),
       dotplot(ego_up, showCategory = 15), width = 8, height = 6, dpi = 300)

ggsave(file.path(proteom_output_pwd, "Dotplot_GO_up_simplified.png"),
       dotplot(ego_up_simplified, showCategory = 15), width = 8, height = 6, dpi = 300)

ggsave(file.path(proteom_output_pwd, "Dotplot_GO_down.png"),
       dotplot(ego_down, showCategory = 15), width = 8, height = 6, dpi = 300)

ggsave(file.path(proteom_output_pwd, "Dotplot_GO_down_simplified.png"),
       dotplot(ego_down_simplified, showCategory = 15), width = 8, height = 6, dpi = 300)

ggsave(file.path(proteom_output_pwd, "Dotplot_GO_ALL_significant.png"),
       dotplot(ego_all, showCategory = 15), width = 8, height = 6, dpi = 300)

ggsave(file.path(proteom_output_pwd, "Dotplot_GO_ALL_significant_simplified.png"),
       dotplot(ego_all_simplified, showCategory = 15), width = 8, height = 6, dpi = 300)

ggsave(file.path(proteom_output_pwd, "Dotplot_Reactome_up.png"),
       dotplot(reactome_up, showCategory = 15), width = 8, height = 6, dpi = 300)

ggsave(file.path(proteom_output_pwd, "Dotplot_Reactome_down.png"),
       dotplot(reactome_down, showCategory = 15), width = 8, height = 6, dpi = 300)

ggsave(file.path(proteom_output_pwd, "Dotplot_Reactome_ALL_significant.png"),
       dotplot(reactome_all, showCategory = 15), width = 8, height = 6, dpi = 300)

