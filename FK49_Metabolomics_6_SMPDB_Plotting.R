gc()
rm(list = ls())
home <- normalizePath("~") # bc Windows does not start inuser dir but in user/documents dir
parent <- dirname(home)
setwd(paste0(parent,"/OneDrive - Universität Salzburg/AG_Tumorimmunologie - Freia Krause/01_Experiments/FK46_FK49_RProject"))

library(ggplot2)
library(dplyr)
library(RColorBrewer)
library(tidyr)
library(readr)
library(stringr)
source("FK49_Definitions.R")

ExpId = "FK49"
path_smpdb <- paste0( parent, "/OneDrive - Universität Salzburg/AG_Tumorimmunologie - Freia Krause/01_Experiments/FK49_CD-HFD_13wks/FK49_HILIC09/smpdb_metabolites")
# Define recoding vector once
metabolite_recoding <- c(
  # Nucleotides
  "Adenosine triphosphate" = "ATP",
  "Adenosine monophosphate" = "AMP",
  "Adenosine diphosphate" = "ADP",   # NOTE: not present in your data
  "Cytidine diphosphate" = "CDP",
  "Cytidine monophosphate" = "CMP",
  "Guanosine diphosphate" = "GDP",
  "Guanosine triphosphate" = "GTP",
  "D-Glucose" = "Glucose",
  "Oxoglutaric acid" = "alpha-Ketoglutaric acid",
  "L-Glutamic acid" = "L-Glutamic Acid",
  "L-Aspartic acid" = "L-Aspartic Acid",
  "Aminoadipic acid" = "Aminoadipic Acid",
  "Sarcosine" = "L-Sarcosine"
)


# Read files 
# Pathway search means: enriched pathways were taken and our data file od targete metabolomics was search for these pathways and all metabolites of these pathway we have also data on are taken 
#  Metabolite search means SMPDB was search for trend (pvalus, FC) metabolites, pathways were these metabos were included where searched for additionaö metabolites that are in our targeted metabo data
PWD_NDvsCDHFD <- readRDS(paste0(targeted_pwd, "/NDvsCDHFD/SMPDB_Pathway_search_NDvsCDHFD.rds")) %>%
  mutate(`Metabolite Name` = recode(`Metabolite Name`, !!!metabolite_recoding))

CDHFD_f <- read.csv(paste0(targeted_pwd, "/CDHFD_f/SMPDBsearch_metabos_manual_cleaned.csv"))%>%
  mutate(`Metabolite` = recode(`Metabolite`, !!!metabolite_recoding))
ND_f <- read.csv(paste0(targeted_pwd, "/ND_f/SMPDBsearch_metabos_ND_f_manual.csv"),sep=";")%>%
  mutate(`Metabolite` = recode(`Metabolite`, !!!metabolite_recoding))

Met_NDvsCDHFD <- readRDS(paste0(targeted_pwd, "/NDvsCDHFD/SMPDB_Metabolite_search_NDvsCDHFD.rds")) %>%
  mutate(`Metabolite Name` = recode(`Metabolite`, !!!metabolite_recoding))

metabolome_targeted <-readRDS(
                      paste0(parent,
                      "/OneDrive - Universität Salzburg/AG_Tumorimmunologie - Freia Krause/01_Experiments/FK49_CD-HFD_13wks/FK49_Analysis/01_RawData/FK49_metabolome_targeted.rds"
                      ))%>%
                      rename_with(~ recode(.x, !!!metabolite_recoding))
str(metabolome_targeted)

plot_violins <- function(data,
                         metabolites,
                         comparison = "Treatment",
                         fill = "Treatment",
                         Sex_filter = NULL,
                         Diet_filter = NULL,
                         ExpID_filter = NULL,
                         method = NULL,
                         pwd_title = "u gave no title idiot",
                         PATH=NULL,
                         FOLDER=NULL) {
  
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(stringr)
  
  # -------------------------
  # 1️⃣ Subset
  # -------------------------
  df <- subset_data(data, 
                    Sex_filter = Sex_filter, 
                    Diet_filter = Diet_filter, 
                    ExpID_filter = ExpID_filter,
                    method = method
  )
  
  # -------------------------
  # 2️⃣ Choose data type
  # -------------------------
  if (method == "untargeted") {
    df_values <- df$norm_values
    y_label_method <- "relative abundance [%]"
    
  } else if (method == "targeted") {
    df_values <- df$raw_values
    y_label_method <- "Abundance [pM/mg]"
    
  } else {
    stop("I need method 'untargeted' or 'targeted'.")
  }
  
  # -------------------------
  # 3️⃣ Choose color map
  # -------------------------
  if (comparison == "Sex") {
    color_map <- Sex_colors
  } else if (comparison == "Treatment") {
    color_map <- Treatment_colors
  } else if (comparison == "T_D_S") {
    color_map <- T_D_S_colors
  } else if (comparison == "T_S") {
    color_map <- T_S_colors
  } else if (comparison == "Diet") {
    color_map <- Diet_colors
  } else {
    stop("Unknown plot variable.")
  }
  
  # -------------------------
  # 4️⃣ Wide → Long
  # -------------------------
  valid_metabs <- intersect(metabolites, colnames(df_values))
  
  if (length(valid_metabs) == 0) {
    # No metabolites found → make empty plot with text
    message("No metabolite of this pathway found in dataset: ", pwd_title)
    
    empty_plot <- ggplot() +
      annotate("text", x = 1, y = 1, label = "No metabolite of this pathway found in dataset", size = 6) +
      theme_void() +
      labs(title = pwd_title)
    
    # make folder if missing
    dir.create(file.path(PATH, FOLDER), showWarnings = FALSE, recursive = TRUE)
    safe_title <- stringr::str_replace_all(pwd_title, "[/\\?<>\\:*|\"^]", "_")
      safe_title <- stringr::str_trunc(safe_title, width = 60, side = "right", ellipsis = "")
    print(safe_title)
    
    # save empty plot
    ggsave(filename = file.path(PATH, FOLDER, paste0("Vio_Out/Violin_", safe_title, ".png")),
           plot = empty_plot, width = 10, height = 6)
    
    return(list(violin = empty_plot))
  }
  
  missing_metabs <- setdiff(metabolites, colnames(df_values))
  if (length(missing_metabs) > 0) {
    warning("Missing metabolites: ", paste(missing_metabs, collapse = ", "))
  }
  
  df_plot <- cbind(df$metadata,
                   df_values[, valid_metabs, drop = FALSE])
  
  df_long <- df_plot %>%
    pivot_longer(cols = all_of(valid_metabs),
                 names_to = "Metabolite",
                 values_to = "value")
  
  # -------------------------
  # 5️⃣ Plot
  # -------------------------
  violin <- ggplot(df_long,
                   aes(x = .data[[comparison]],
                       y = value,
                       fill = .data[[fill]])) +
    geom_violin(trim = FALSE, alpha = 0.5) +
    geom_boxplot(width = 0.1, outlier.shape = NA) +
    geom_jitter(width = 0.1, size = 1.5, alpha = 0.7) +
    facet_wrap(~Metabolite, scales = "free_y", ncol = 4) +
    theme_classic() +
    scale_fill_manual(values = color_map) +
    labs(title = pwd_title,
         x = comparison,
         y = y_label_method)
  
  # -------------------------
  # 6️⃣ Save plot
  # -------------------------
  dir.create(file.path(PATH, FOLDER), showWarnings = FALSE, recursive = TRUE)
  safe_title <- stringr::str_replace_all(pwd_title, "[/\\?<>\\:*|\"^]", "_")
  safe_title <- stringr::str_trunc(safe_title, width = 60, side = "right", ellipsis = "")
  # truncate to max 50 characters to avoid Windows/OneDrive issues
  
  # save plot
  ggsave(
    filename = file.path(PATH, FOLDER, paste0("Vio_Out/Violin_", safe_title, ".png")),
    plot = violin,
    width = 10,
    height = 6
  )
  return(list(violin = violin))
  
}

pathways <- unique(CDHFD_f$pathway_name)
results <- lapply(pathways, function(pw) {
  
  # get metabolites for this pathway
  mets <- CDHFD_f %>%
    dplyr::ungroup() %>%
    dplyr::filter(pathway_name == pw) %>%
    dplyr::pull(Metabolite)
  
  # run your function
  plot_violins(
    data = metabolome_targeted,
    metabolites = mets,
    comparison = "Treatment",
    fill = "Treatment",
    Sex_filter = c("female"),
    Diet_filter = c("CDHFD13"),
    ExpID_filter = c("FK49"),
    method = "targeted",
    PATH=targeted_pwd,
    FOLDER="CDHFD_f",
    pwd_title = pw
  )
})

pathways <- unique(ND_f$pathway_name)
results <- lapply(pathways, function(pw) {
  
  # get metabolites for this pathway
  mets <- ND_f %>%
    dplyr::ungroup() %>%
    dplyr::filter(pathway_name == pw) %>%
    dplyr::pull(Metabolite)
  
  # run your function
  plot_violins(
    data = metabolome_targeted,
    metabolites = mets,
    comparison = "Treatment",
    fill = "Treatment",
    Sex_filter = c("female"),
    Diet_filter = c("ND"),
    ExpID_filter = c("BH"),
    method = "targeted",
    PATH=targeted_pwd,
    FOLDER="ND_f",
    pwd_title = pw
  )
})


#If I only want to plot specific single pathways
pwd_CDHFD_f<-CDHFD_f%>%
  dplyr::filter(pathway_name == "Isovaleric Aciduria") %>%
  dplyr::pull(Metabolite)

plot_violins( data = metabolome_targeted,
              metabolites = pwd_CDHFD_f,  comparison = "Treatment",
              fill = "Treatment",  Sex_filter  = c("female"),  Diet_filter = c("CDHFD13"),
              ExpID_filter = c("FK49"),  method = "targeted",  pwd_title="Isovaleric Aciduria")
