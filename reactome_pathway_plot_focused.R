# =============================================================================
# Reactome Pathway Plotting - Focused + Zoomable
# - Filtert Komponenten auf >= MIN_MEASURED_PROTEINS gemessene (eindeutige!) Proteine
# - Metabolit-Matching: Prioritaet KEGG -> ChEBI -> PubChem
# - Speichert Pathway-Ancestry (Hierarchie) als Tabelle
# - Zoom: Overview-Plot + Detail-Plots pro Community + GraphML fuer Cytoscape
# =============================================================================
rm(list=ls())
gc()
library(ReactomeContentService4R)
library(dplyr)
library(purrr)
library(tidyr)
library(igraph)
library(ggraph)
library(ggplot2)
source("FK49_Definitions.R")

# --- Konfiguration -----------------------------------------------------------
MIN_MEASURED_PROTEINS <- 5   # mindestanzahl gemessener (eindeutiger) Proteine pro Komponente
DETAIL_THRESHOLD      <- 25    # Komponenten > N Knoten werden in Communities gezoomt
DETAIL_MAX_NODES      <- 40    # Communities mit > N Knoten werden nicht einzeln geplottet
LAYOUT                <- "stress"
proteom_output_pwd <- PATHS$proteomics$output
out_dir <- file.path(proteom_output_pwd, "Pathways")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

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

reactome_prots_in_my_prots <- read.csv2( file.path(proteom_output_pwd,"Data/Reactome_top10_proteins_in_myproteins.csv"))
pathway_ids <- unique(reactome_prots_in_my_prots$ID)
reactome_df <- read.csv2(file.path(proteom_output_pwd,"Data/Reactome_top10_df.csv"))
pathway_names <- reactome_df %>%
  dplyr::filter(ID %in% pathway_ids) %>%
  dplyr::select(ID, Description)
expr_data_proteins <- readRDS( file.path(proteom_output_pwd, "Data/Expr_Data_Proteins.rds"))

# --- KEGG -> ChEBI Uebersetzungstabelle (einmalig laden, dann cachen) --------
kegg_chebi_map_path <- file.path(proteom_output_pwd, "Data/kegg_chebi_map.rds")

if (!file.exists(kegg_chebi_map_path)) {
  kegg_chebi_raw <- httr::GET("https://rest.kegg.jp/conv/chebi/compound") %>%
    httr::content(as = "text", encoding = "UTF-8")
  
  kegg_chebi_map <- read.delim(text = kegg_chebi_raw, header = FALSE,
                               col.names = c("kegg_raw", "chebi_raw")) %>%
    dplyr::mutate(
      KEGG  = sub("^cpd:", "", kegg_raw),
      ChEBI = sub("^chebi:", "", chebi_raw)
    ) %>%
    dplyr::select(KEGG, ChEBI) %>%
    dplyr::distinct()
  
  saveRDS(kegg_chebi_map, kegg_chebi_map_path)
} else {
  kegg_chebi_map <- readRDS(kegg_chebi_map_path)
}
# --- Metabolomik-Lookup-Tabellen (Prioritaet: KEGG -> ChEBI -> PubChem) ------
# WICHTIG: KEGG ist bei euch fast vollstaendig vorhanden (0 NAs), ChEBI fast nie
# (640/643 NA), PubChem teilweise (572/643 NA) - KEGG ist daher die wichtigste
# Quelle und wird als erstes geprueft.
# Native ChEBI-Werte (selten, aber falls vorhanden direkt nutzen)
metab_native_chebi <- abundance_data_metabos %>%
  dplyr::filter(!is.na(ChEBI)) %>%
  dplyr::mutate(ChEBI_final = as.character(as.integer(ChEBI))) %>%
  dplyr::select(Metabolite, ChEBI_final, log2FC)

# Ueber KEGG->ChEBI uebersetzte Werte (deckt fast alles ab bei euch)
metab_via_kegg <- abundance_data_metabos %>%
  dplyr::filter(!is.na(KEGG), KEGG != "") %>%
  dplyr::mutate(KEGG = as.character(KEGG)) %>%select(-ChEBI)%>%
  dplyr::left_join(kegg_chebi_map, by = "KEGG") %>%
  dplyr::filter(!is.na(ChEBI)) %>%
  dplyr::select(Metabolite, ChEBI_final = ChEBI, log2FC)

# Kombinieren: native ChEBI hat Vorrang, sonst die KEGG-Uebersetzung
metab_fc_chebi_unified <- dplyr::bind_rows(metab_native_chebi, metab_via_kegg) %>%
  dplyr::distinct(Metabolite, .keep_all = TRUE) %>%
  dplyr::select(ChEBI = ChEBI_final, log2FC) %>%
  dplyr::distinct(ChEBI, .keep_all = TRUE)

# --- Hilfsfunktion: Pathway-Ancestry -----------------------------------------
get_pathway_ancestry <- function(stId) {
  anc <- tryCatch(getEventAncestors(stId), error = function(e) NULL)
  if (is.null(anc) || length(anc) == 0) {
    return(tibble(pathway_stId = stId, level = NA_integer_,
                  ancestor_stId = NA_character_, ancestor_name = NA_character_,
                  ancestor_type = NA_character_))
  }
  purrr::imap_dfr(anc, function(a, lvl) {
    if (is.list(a) && !is.data.frame(a)) {
      tibble(pathway_stId = stId, level = as.integer(lvl),
             ancestor_stId  = a$stId %||% NA_character_,
             ancestor_name  = a$displayName %||% NA_character_,
             ancestor_type  = a$type %||% NA_character_)
    } else if (is.data.frame(a)) {
      a %>% mutate(pathway_stId = stId, level = as.integer(lvl)) %>%
        dplyr::rename(ancestor_stId = stId, ancestor_name = displayName,
                      ancestor_type = type) %>%
        dplyr::select(pathway_stId, level, ancestor_stId, ancestor_name, ancestor_type)
    } else {
      tibble(pathway_stId = stId, level = as.integer(lvl),
             ancestor_stId = NA_character_, ancestor_name = as.character(a),
             ancestor_type = NA_character_)
    }
  })
}

`%||%` <- function(a, b) if (is.null(a) || length(a) == 0) b else a

# --- Hilfsfunktion: Fokus-Filter ---------------------------------------------
build_focused_graph <- function(nodes, edges,
                                measured_protein_ids, measured_metab_ids,
                                min_proteins = MIN_MEASURED_PROTEINS) {
  
  measured_all <- union(measured_protein_ids, measured_metab_ids)
  edges_meas   <- edges %>% dplyr::filter(from %in% measured_all | to %in% measured_all)
  
  kept_reactions <- unique(c(
    edges_meas$to[grepl("^R_", edges_meas$to)],
    edges_meas$from[grepl("^R_", edges_meas$from)]
  ))
  
  kept_entities <- unique(c(
    edges %>% dplyr::filter(to %in% kept_reactions) %>% pull(from),
    edges %>% dplyr::filter(from %in% kept_reactions) %>% pull(to)
  ))
  
  keep_node_ids <- unique(c(kept_reactions, kept_entities, measured_all))
  keep_node_ids <- intersect(keep_node_ids, nodes$id)
  
  nodes_f <- nodes %>% dplyr::filter(id %in% keep_node_ids)
  edges_f <- edges  %>% dplyr::filter(from %in% keep_node_ids & to %in% keep_node_ids)
  
  if (nrow(edges_f) == 0) return(NULL)
  
  # nur Knoten behalten, die tatsaechlich mind. eine Kante in edges_f haben
  # (verhindert isolierte Ein-Knoten-"Komponenten" ohne edge-Attribute)
  nodes_with_edges <- union(edges_f$from, edges_f$to)
  nodes_f <- nodes_f %>% dplyr::filter(id %in% nodes_with_edges)
  
  g <- graph_from_data_frame(
    edges_f,
    vertices = nodes_f %>% dplyr::select(-any_of("refEntities")),
    directed = TRUE)
  comp <- components(g)
  
  # eindeutige gemessene Proteine pro Komponente zaehlen (statt Knoten!)
  id_to_identifier <- setNames(nodes_f$identifier, nodes_f$id)
  
  comp_protein_count <- sapply(seq_len(comp$no), function(k) {
    vids <- V(g)[comp$membership == k]$name
    prot_vids <- vids[vids %in% measured_protein_ids]
    length(unique(na.omit(id_to_identifier[prot_vids])))
  })
  
  # zusaetzliche Sicherheit - Komponenten ohne Kanten explizit raus
  comp_ecount <- sapply(seq_len(comp$no), function(k) {
    vids <- V(g)[comp$membership == k]$name
    sum(edges_f$from %in% vids & edges_f$to %in% vids)
  })
  
  keep_comps <- which(comp_protein_count >= min_proteins & comp_ecount > 0)
  
  if (length(keep_comps) == 0) return(NULL)
  
  list(graph = g, comp = comp, keep_comps = keep_comps,
       comp_protein_count = comp_protein_count, nodes = nodes_f, edges = edges_f)
}

# Hilfsfunktion: eindeutige gemessene Proteine in einem (Teil-)Graphen zaehlen
count_unique_measured_proteins <- function(g_x, measured_protein_ids, id_to_identifier) {
  vids <- V(g_x)$name
  prot_vids <- vids[vids %in% measured_protein_ids]
  length(unique(na.omit(id_to_identifier[prot_vids])))
}

# --- Plot-Funktion: gerichtete Pathway-Map -----------------------------------
make_plot <- function(g_sub, title) {
  
  if (igraph::ecount(g_sub) == 0) {
    warning("make_plot: g_sub hat keine Kanten, ueberspringe: ", title)
    return(NULL)
  }
  
  # Kantentyp fuer die Legende vereinheitlichen (VOR create_layout setzen!)
  E(g_sub)$edge_kind <- dplyr::if_else(
    E(g_sub)$relation %in% c("input", "output"),
    "Stofffluss (input/output)",
    "Katalysator"
  )
  
  lay <- create_layout(
    g_sub,
    layout  = "sugiyama",
    hgap    = 1.5,
    vgap    = 1,
    maxiter = 100,
    weights = NA
  )
  
  p <- ggraph(lay) +
    
    # Stofffluss-Kanten (durchgezogen, mit Pfeil)
    geom_edge_link(
      aes(filter = relation %in% c("input", "output"),
          edge_linetype = edge_kind),
      linewidth = 0.35,
      alpha = 0.55,
      arrow = grid::arrow(length = grid::unit(2, "mm"), type = "closed"),
      end_cap = circle(2, "mm")
    ) +
    
    # Katalysator-Kanten (gestrichelt)
    geom_edge_link(
      aes(filter = relation == "catalyst",
          edge_linetype = edge_kind),
      linewidth = 0.3,
      alpha = 0.45,
      end_cap = circle(2, "mm")
    ) +
    
    scale_edge_linetype_manual(
      name   = "Kantentyp",
      values = c("Stofffluss (input/output)" = "solid",
                 "Katalysator" = "dashed")
    ) +
    
    geom_node_point(
      data = function(x) x %>% dplyr::filter(type == "Reaction"),
      shape = 15, size = 3, colour = "grey30"
    ) +
    
    geom_node_point(
      data = function(x) x %>% dplyr::filter(type != "Reaction", type != "Metabolite"),
      aes(fill = log2FC),
      shape = 22, size = 9, stroke = 1, colour = "black"
    ) +
    
    geom_node_point(
      data = function(x) x %>% dplyr::filter(type == "Metabolite"),
      aes(fill = log2FC),
      shape = 21, size = 5, stroke = 0.5, colour = "black"
    ) +
    
    scale_fill_gradient2(
      low = "blue", mid = "white", high = "red", midpoint = 0,
      na.value = "grey85", name = "log2FC"
    ) +
    
    geom_node_text(
      data = function(x) x %>% dplyr::filter(type != "Reaction", !is.na(plot_label)),
      aes(label = plot_label),
      size = 2.8, repel = TRUE
    ) +
    
    labs(title = title) +
    theme_void() +
    theme(
      plot.title = element_text(size = 13, face = "bold", hjust = 0.5),
      legend.position = "right"
    )
  
  p
}

# --- Ancestry-Tabelle ueber alle Top-10 Pathways -----------------------------
ancestry_all <- purrr::map_dfr(pathway_ids, get_pathway_ancestry)
write.csv(ancestry_all,  file.path(out_dir, "pathway_ancestry.csv"),row.names = FALSE)
pathway_ids <- pathway_ids[3]

# --- Hauptschleife: Top 10 Pathways ------------------------------------------
for (i in seq_len(min(10, length(pathway_ids)))) {
  i<-1
  pw      <- pathway_ids[i]
  pw_name <- pathway_names$Description[match(pw, pathway_names$ID)]
  pw_label <- gsub(" \\[.*", "", pw_name[1])
  pw_label <- gsub("[^[:alnum:] _/-]", "", pw_label)
  
  cat(sprintf("[%d] %s\n", i, pw_label))
  
  events <- getParticipants(pw, retrieval = "EventsInPathways")
  
  reactions <- events[
    purrr::map_lgl(events, function(x) {
      is.list(x) &&
        identical(x$className, "Reaction") &&
        !is.null(x$stId) &&
        length(x$stId) == 1 &&
        !is.na(x$stId)
    })
  ]
  
  reaction_data <- purrr::map_dfr(reactions, function(r) {
    x <- tryCatch(
      getParticipants(r$stId, retrieval = "AllInstances"),
      error = function(e) NULL
    )
    if (is.null(x) || !is.data.frame(x) || nrow(x) == 0) return(NULL)
    x %>% mutate(reaction_id = r$stId, reaction_name = r$displayName)
  })
  
  if (nrow(reaction_data) == 0 || !"peDbId" %in% names(reaction_data)) {
    cat("  -> keine gültigen Participants gefunden, überspringe.\n")
    next
  }
  
  # NODES
  reaction_nodes <- reactions %>%
    purrr::map_dfr(~ tibble(id = paste0("R_", .x$stId),
                            label = .x$displayName, type = "Reaction"))
  
  physical_nodes <- reaction_data %>%
    transmute(
      id = paste0("PE_", peDbId), label = displayName,
      type = case_when(
        schemaClass == "SimpleEntity"                  ~ "Metabolite",
        schemaClass %in% c("DefinedSet", "Complex")    ~ "Protein",
        TRUE                                           ~ schemaClass
      ),
      peDbId = peDbId, refEntities = refEntities
    ) %>%
    distinct(id, .keep_all = TRUE)
  
  nodes <- bind_rows(reaction_nodes, physical_nodes)
  
  # EDGES
  edges <- reaction_data %>%
    dplyr::filter(type %in% c("input", "output", "catalyst")) %>%
    transmute(
      from = case_when(
        type == "input"    ~ paste0("PE_", peDbId),
        type == "output"   ~ paste0("R_",  reaction_id),
        type == "catalyst" ~ paste0("PE_", peDbId)
      ),
      to = case_when(
        type == "input"    ~ paste0("R_",  reaction_id),
        type == "output"   ~ paste0("PE_", peDbId),
        type == "catalyst" ~ paste0("R_",  reaction_id)
      ),
      relation = type
    )
  
  # REFERENCE TABLE
  ref_table <- purrr::map2_dfr(
    physical_nodes$id, physical_nodes$refEntities,
    function(x, y) {
      if (is.null(y) || nrow(y) == 0) NULL else y %>% mutate(node_id = x)
    })
  
  # ---------------------------------------------------------------
  # OPTIONAL, aber empfohlen VOR dem ersten Lauf einmal pruefen:
  # ref_table %>% dplyr::filter(schemaClass == "ReferenceMolecule") %>%
  #   dplyr::mutate(db = tolower(sub(":.*", "", stId))) %>% dplyr::count(db)
  # -> zeigt, welche DB-Praefixe (chebi/kegg/...) bei euch tatsaechlich
  #    in Reactomes ReferenceMolecule-Eintraegen vorkommen.
  # ---------------------------------------------------------------
  
  # PROTEIN FC
  protein_fc <- expr_data_proteins %>%
    mutate(identifier = toupper(ID)) %>%
    dplyr::select(identifier, logFC) %>%
    dplyr::rename(log2FC_protein = logFC)
  
  protein_identifier_map <- ref_table %>%
    dplyr::filter(schemaClass == "ReferenceGeneProduct") %>%
    dplyr::select(node_id, identifier) %>%
    distinct(node_id, .keep_all = TRUE)
  
  nodes <- nodes %>%
    left_join(protein_identifier_map, by = c("id" = "node_id"))
  
  nodes <- nodes %>%
    left_join(
      ref_table %>%
        dplyr::filter(schemaClass == "ReferenceGeneProduct") %>%
        dplyr::select(node_id, identifier) %>%
        left_join(protein_fc, by = "identifier") %>%
        group_by(node_id) %>%
        summarise(log2FC_protein =
                    dplyr::first(log2FC_protein[!is.na(log2FC_protein)]),
                  .groups = "drop"),
      by = c("id" = "node_id")
    )
  
  # --- METABOLITE FC: Prioritaet KEGG -> ChEBI -> PubChem --------------------
  ref_ids_long <- ref_table %>%
    dplyr::filter(schemaClass == "ReferenceMolecule") %>%
    dplyr::mutate(db = tolower(sub(":.*", "", stId))) %>%
    dplyr::select(node_id, db, identifier) %>%
    dplyr::distinct(node_id, db, .keep_all = TRUE)
  
  ref_ids_wide <- ref_ids_long %>%
    tidyr::pivot_wider(names_from = db, values_from = identifier)
  
  # Fehlende Spalten sicher ergaenzen, falls eine DB gar nicht vorkommt
  for (col in c("kegg", "chebi", "pubchem.compound")) {
    if (!col %in% names(ref_ids_wide)) ref_ids_wide[[col]] <- NA_character_
  }
  
  ref_ids_wide <- ref_ids_wide %>%
    dplyr::rename(KEGG = kegg, ChEBI = chebi, PubChem = `pubchem.compound`) %>%
    dplyr::mutate(KEGG = as.character(KEGG),
                  ChEBI = as.character(ChEBI),
                  PubChem = as.character(PubChem))
  
  # metab_match <- ref_ids_wide %>%
  #   dplyr::left_join(metab_fc_kegg,    by = "KEGG")    %>% dplyr::rename(log2FC_kegg = log2FC) %>%
  #   dplyr::left_join(metab_fc_chebi,   by = "ChEBI")   %>% dplyr::rename(log2FC_chebi = log2FC) %>%
  #   dplyr::left_join(metab_fc_pubchem, by = "PubChem") %>% dplyr::rename(log2FC_pubchem = log2FC) %>%
  #   dplyr::mutate(
  #     log2FC_metab = dplyr::coalesce(log2FC_kegg, log2FC_chebi, log2FC_pubchem)
  #   ) %>%
  #   dplyr::select(node_id, log2FC_metab)
  nodes2 <- nodes %>%
    left_join(
      x<-ref_table %>%
        dplyr::filter(schemaClass == "ReferenceMolecule") %>%
        mutate(ChEBI = as.character(identifier)) %>%
        dplyr::select(node_id, ChEBI) %>%
        left_join(metab_fc_chebi_unified, by = "ChEBI") %>%
        group_by(node_id) %>%
        summarise(log2FC_metab = dplyr::first(log2FC[!is.na(log2FC)]), .groups = "drop"),
      by = c("id" = "node_id")
    ) %>%
    mutate(
      log2FC = coalesce(log2FC_protein, log2FC_metab),
  # nodes <- nodes %>%
  #   dplyr::left_join(metab_match, by = c("id" = "node_id")) %>%
  #   dplyr::mutate(
  #     log2FC = coalesce(log2FC_protein, log2FC_metab),
      plot_label = gsub(" \\[.*", "", label),
      plot_label = case_when(
        type == "Reaction"                                              ~ paste0("R: ", plot_label),
        type %in% c("Protein", "CandidateSet", "EntityWithAccessionedSequence") ~ sub(" .*", "", plot_label),
        TRUE                                                            ~ plot_label
      )
    ) %>%
    dplyr::select(-log2FC_protein, -log2FC_metab) %>%
    mutate(
      in_data     = !is.na(log2FC),
      node_shape  = case_when(
        type == "Reaction"   ~ 15L,
        type == "Protein"    ~ 22L,
        type == "Metabolite" ~ 21L,
        TRUE                 ~ 20L
      ),
      border_col  = if_else(in_data, "black", "grey70"),
      fill_val    = log2FC
    )
  
  # --- gemessene Protein-/Metabolit-IDs -------------------------------------
  measured_protein_ids <- nodes %>%
    dplyr::filter(
      type != "Reaction",
      type != "Metabolite",
      in_data
    ) %>%
    pull(id)
  measured_metab_ids   <- nodes %>%
    dplyr::filter(type == "Metabolite", in_data) %>% pull(id)
  
  # --- Fokus-Graph ----------------------------------------------------------
  fg <- build_focused_graph(nodes, edges,
                            measured_protein_ids, measured_metab_ids,
                            min_proteins = MIN_MEASURED_PROTEINS)
  if (is.null(fg)) {
    cat("  -> keine Komponente mit >= ", MIN_MEASURED_PROTEINS,
        " gemessenen Proteinen, ueberspringe.\n", sep = "")
    next
  }
  
  g        <- fg$graph
  comp     <- fg$comp
  keep_comps <- fg$keep_comps
  
  # --- GraphML fuer Cytoscape (interaktiver Zoom) ---------------------------
  graphml_path <- file.path(out_dir,
                            sprintf("Reactome_%d_%s.graphml", i, pw_label))
  write_graph(g, as.character(graphml_path[1]), format = "graphml")
  
  # --- Overview + Detail-Plots pro Komponente -------------------------------
  id_to_identifier <- setNames(fg$nodes$identifier, fg$nodes$id)
  
  for (k in keep_comps) {
    g_sub <- induced_subgraph(g, vids = V(g)[comp$membership == k])
    n_nodes <- vcount(g_sub)
    n_prot  <- fg$comp_protein_count[k]
    
    ov_title <- sprintf("%s | Komp. %d | %d Knoten | %d Proteine",
                        pw_label, k, n_nodes, n_prot)
    p <- make_plot(g_sub, ov_title)
    if (is.null(p)) next
    ggsave(file.path(out_dir,
                     sprintf("Reactome_%d_%s_Overview_K%d.png", i, pw_label, k)),
           p, width = 16, height = 12, units = "in", dpi = 200,
           bg = "white", limitsize = FALSE)
    
    # --- ZOOM: grosse Komponenten in Communities zerlegen -------------------
    if (n_nodes > DETAIL_THRESHOLD) {
      cl <- cluster_walktrap(as.undirected(g_sub, mode = "collapse"))
      mem <- membership(cl)
      comm_ids <- unique(mem)
      
      for (cm in comm_ids) {
        g_comm <- induced_subgraph(g_sub, vids = V(g_sub)[mem == cm])
        if (vcount(g_comm) < 3 || vcount(g_comm) > DETAIL_MAX_NODES) next
        
        comm_prot_count <- count_unique_measured_proteins(g_comm, measured_protein_ids, id_to_identifier)
        if (comm_prot_count < MIN_MEASURED_PROTEINS) next
        
        cm_title <- sprintf("%s | Komp.%d -> Cluster %d | %d Knoten | %d Proteine",
                            pw_label, k, cm, vcount(g_comm), comm_prot_count)
        p2 <- make_plot(g_comm, cm_title)
        if (is.null(p2)) next
        ggsave(file.path(out_dir,
                         sprintf("Reactome_%d_%s_Detail_K%d_C%d.png", i, pw_label, k, cm)),
               p2, width = 12, height = 10, units = "in", dpi = 200,
               bg = "white", limitsize = FALSE)
      }
    }
  }
}
