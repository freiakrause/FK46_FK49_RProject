#### log transformation, normalization sclaing,
#PCA
# test linear model, FC , p adjust, volcano and violins of trend or significant metabolites perforem in FK49_untartgeted.R
# Also Metbaoanalyst was used to perfomr over representation analysis (ora) on the trend metabolites. USING smpdb pathways (KEGG did not work with metaboanalyst)
#From the recevied plots pathway names where manually read and included here
# Here I want to search the downloaded SMPDB csv files for Pathways that are matching my ora pathways
# In the next step, metabolites from the matched ORA pathways are compared to the list of metabolites I receveid from MetaboNet (comparison on KEGG ID level, not name, ID is more precise)
# I save this df of pathways and metabolites and read it in enxt script
# Trend Metabolites come from analyssed sets: ND f vs CDHFD f (ttest+FDR); CDHFD a TAM vs EtOH (linear model+FDR); CDHFD female TAM vs EtOH (ttest+FDR); CDHFD male TAM vs EtOH (ttest+FDR)
# ND vs CDHFD is only to be bale to generally explain " in ND this an that pathway is used in CDHFD this and that pathway is used" 
# ND and CDHFD were not performed together an ND is not a planned controlled ctrl group for CDHFD. ND is young female mice that were used in Birgits Experimetns. CDHFD is my FK49 experiment
# in next the script all measured metabolites that are in an ora pathways should be plotted as violins (abundance) to check potential regulation of pathway
gc()
rm(list = ls())
home <- normalizePath("~") # bc Windows does not start inuser dir but in user/documents dir
parent <- dirname(home)
setwd(paste0(parent,"/OneDrive - Universität Salzburg/Freigegebene Dokumente - AG_Tumorimmunologie/Data/Freia Krause/01_Experiments/FK46_FK49_RProject"))

library(data.table)
library(future.apply)
library(tools)
library(future)
library(dplyr)
source("FK49_Definitions.R")
home <- normalizePath("~") # bc Windows does not start inuser dir but in user/documents dir
parent <- dirname(home) # goves users dir works for me since i aways have the data from work in the onedrive at the same location and only local home dir changes
pwd_to_experiments <-paste0("/OneDrive - Universität Salzburg/Freigegebene Dokumente - AG_Tumorimmunologie/Data/Freia Krause/01_Experiments/")
ExpId = "FK49"
# ---- Define workers ----
# IDs for metabolites analysed in targeted HILIC09 from MetaboNet
targeted_metabolites<-read.csv(file= paste0(parent,pwd_to_experiments,"FK49_CD-HFD_13wks/FK49_HILIC09/KEGG_MetaboAnalyst.csv"),sep= ";")
untargeted_metabolites<-
output_folders <- c("CDHFD_a",   "CDHFD_m", "CDHFD_f",  "ND_f",          "NDvsCDHFD")
names <- c("CDHFD_a", 
           "CDHFD_a_p", 
           "CDHFD_a_n", 
           "CDHFD_m",  
           "CDHFD_m_p",
           "CDHFD_M_n",
           "CDHFD_f", 
           "CDHFD_f_p",
           "CDHFD_f_n",
           "ND_f",
           "ND_f_p",
           "ND_f_n",
           "NDvsCDHFD",
           "NDvsCDHFD_p",
           "NDvsCDHFD_n")

read_metabos<-function(PATH= NULL, FOLDER = NULL, NAME = NULL){
  setwd(paste0(PATH,"/",FOLDER,"/"))
  trend_metabolites <-read.csv(paste0(NAME,"_trend_METABOANALYST.csv"), sep= ";")
  setwd(paste0(parent,"/OneDrive - Universität Salzburg/Freigegebene Dokumente - AG_Tumorimmunologie/Data/Freia Krause/01_Experiments/FK46_FK49_RProject"))
  message("I read ",paste0(NAME,"_trend_METABOANALYST.csv"))
  HMDB_KEGG <- na.omit(c(trend_metabolites$Query, trend_metabolites$KEGG, trend_metabolites$HMDB))
  return(HMDB_KEGG)
  }
  # metabolite_vector<-readRDS(paste0(PATH,"/",FOLDER,"/",NAME,".rds")) # this was reading the rds but now i manually colected all HMDB or KEGG ID for the metabolutes. Thae are inthe files with addiotn METABOANALYST
  # message("I read ",paste0(PATH,"/",FOLDER,"/",NAME,".rds"))
  # trend_metabolites<-metabolite_vector%>%filter(trend == TRUE)
  # trend_metabolites<-trend_metabolites[["Metabolite"]]
  # print(head(trend_metabolites))


# met_ND_CDHFD<-read_metabos(PATH=targeted_pwd, FOLDER="NDvsCDHFD" , NAME="NDvsCDHFD") #vector of metabolites significantöy (or trend) diff in ND vs CDHFD
# matched_rows <- target_metabolites$Query %in% met_ND_CDHFD
# met_diff_ND_CDHFD <- c( target_metabolites$HMDB[matched_rows],  target_metabolites$KEGG[matched_rows]) # get KEGG ID and HMDB ID s for these metabolites
# rm(met_ND_CDHFD)
CDHFD_a_p<-read_metabos(PATH=untargeted_pwd , FOLDER="CDHFD_a" , NAME="CDHFD_a_p")
CDHFD_a_n<-read_metabos(PATH=untargeted_pwd , FOLDER="CDHFD_a" , NAME="CDHFD_a_n")
CDHFD_f  <-read_metabos(PATH=targeted_pwd , FOLDER="CDHFD_f" , NAME="CDHFD_f") # need do do files with HMDB IDs
#ND_f_p  <-read_metabos(PATH=untargeted_pwd , FOLDER="ND_f" , NAME="ND_f_p") # need do do files with HMDB IDs
#ND_f_n  <-read_metabos(PATH=untargeted_pwd , FOLDER="ND_f" , NAME="ND_f_n") # need do do files with HMDB IDs
ND_f     <-read_metabos(PATH=targeted_pwd, FOLDER="ND_f" , NAME="ND_f") 
CDHFD_m_p<-read_metabos(PATH=untargeted_pwd , FOLDER="CDHFD_m" , NAME="CDHFD_m_p")
CDHFD_m_n<-read_metabos(PATH=untargeted_pwd , FOLDER="CDHFD_m" , NAME="CDHFD_m_n")
CDHFD_f_p<-read_metabos(PATH=untargeted_pwd , FOLDER="CDHFD_f" , NAME="CDHFD_f_p")
CDHFD_f_n<-read_metabos(PATH=untargeted_pwd , FOLDER="CDHFD_f" , NAME="CDHFD_f_n")


plan(multisession, workers = 5)
#### from targeted Metabolomics
# ---- Your pathway list ----
pwd_ora_diff_ND_CDHFD <- c(
  "Warburg Effect",
  "Gluconeogenesis",
  "Glucose-Alanine-Cycle", 
  "Lactose Synthesis", 
  "Transfer of Acetyl Groups into Mitochondria", 
  "Glycolysis", 
  "Glycine and Serine Metabolism")


path_smpdb <- paste0(
  parent,
  "/OneDrive - Universität Salzburg/Freigegebene Dokumente - AG_Tumorimmunologie/Data/Freia Krause/01_Experiments/FK49_CD-HFD_13wks/FK49_HILIC09/smpdb_metabolites" #
)

all_files <- list.files(path_smpdb, pattern = "\\.csv$", full.names = TRUE)
### Find ora pathways in my smpdb "database" 
# if path is in ora list and  database, read csv and compare with my total list of analysed metabolits (from MetaboNet, all metabolites where i have values)
#save pathways and all metabolites that ihave values for that elngto these pathways in an object/df
# Let this part run on 12.02. 16:35. saved RDS, for useage rds can now be read without the calculations
# NDvsCDHFD <- future_lapply(all_files, 
#                            function(file) {
#   target_pathways<-pwd_ora_diff_ND_CDHFD
#     # Read first rows to check pathway name
#    header_check <- tryCatch(
#      fread(file, nrows = 2, showProgress = FALSE),
#      error = function(e) return(NULL)
#     )
#   
#    if (is.null(header_check) || !"Pathway Name" %in% colnames(header_check)) {
#       return(NULL)
#     }
#   
#     pathway_name <- header_check[["Pathway Name"]][1]
#   
#     if (!pathway_name %in% target_pathways) {
#       return(NULL)
#     }
#     
#     # Read full file
#     df <- fread(file, showProgress = FALSE)
#     if (!'KEGG ID' %in% colnames(df)) {
#       # Return pathway-only row
#       return(data.table(
#         pathway_name = pathway_name,
#         Metabolite = NA_character_,
#         pathway_file = tools::file_path_sans_ext(basename(file))
#       ))
#     }
#     print(df$'KEGG ID')
#     print(target_metabolites)
#     df_filtered <- df[`KEGG ID` %in% target_metabolites]
#   
#     if (nrow(df_filtered) == 0) {
#       # Return pathway-only row when no metabolites matched
#       return(data.table(
#         pathway_name = pathway_name,
#         Metabolite = NA_character_,
#        pathway_file = tools::file_path_sans_ext(basename(file))
#      ))
#     }
#   
#     # Add pathway info to matched metabolites
#     df_filtered[, pathway_name := pathway_name]
#    df_filtered[, pathway_file := tools::file_path_sans_ext(basename(file))]
#   
#     return(df_filtered)
#   }) 
# 
# # Combine
# final_NDvsCDHFD <- rbindlist(NDvsCDHFD, fill = TRUE)
# # Optional cleanup
# rm(NDvsCDHFD)
# gc()
#saveRDS(final_NDvsCDHFD,  file = paste0(targeted_pwd,"/SMPDB_Pathway_search_NDvsCDHFD.rds"))
#write.csv(final_NDvsCDHFD,  file = paste0(targeted_pwd,"/NDvsCDHFD/SMPDB_Pathway_search_NDvsCDHFD.csv"), row.names = FALSE)
#final_NDvsCDHFD <- readRDS(paste0(targeted_pwd, "/NDvsCDHFD/SMPDB_Pathway_search_NDvsCDHFD.rds"))

### Search single different meatbolite in "my" smpdb database and give back df/rds with pathwa that contains that metabolite ans also all other metabolite of that pathways whcih match with my total metabolite list
# Function to search for a specific metabolite in the SMPDB database
search_metabolite_in_smpdb <- function(PATH, FOLDER,NAME, search_metabolites, 
                                       target_metabolites, 
                                       path_library) {
  
  all_files <- list.files(path_library, pattern = "\\.csv$", full.names = TRUE)
  
  results <- future_lapply(all_files, function(file) {
    
    header_check <- tryCatch(
      fread(file, nrows = 2, showProgress = FALSE),
      error = function(e) return(NULL)
    )
    
    if (is.null(header_check) || !"Pathway Name" %in% colnames(header_check)) {
      return(NULL)
    }
    
    pathway_name <- header_check[["Pathway Name"]][1]
    
    df <- fread(file, showProgress = FALSE)
    
    # Identify available ID columns
    id_cols <- intersect(c("Metabolite Name","KEGG ID", "HMDB ID"), colnames(df))
    
    if (length(id_cols) == 0) {
      return(NULL)
    }
    
    # Check if ANY search metabolite matches KEGG or HMDB
    matched <- Reduce(`|`, lapply(id_cols, function(col) {
      df[[col]] %in% search_metabolites
    }))
    
    if (any(matched)) {
      
      # Filter using target metabolites across both ID types
      df_filtered <- df[
        ( "KEGG ID" %in% id_cols  & df[["KEGG ID"]] %in% target_metabolites$KEGG ) |
          ( "HMDB ID" %in% id_cols  & df[["HMDB ID"]] %in% target_metabolites$HMDB |
            ("Metabolite Name" %in% id_cols  & df[["Metabolite Name"]] %in% target_metabolites$Query))
      ]
      
      if (nrow(df_filtered) == 0) return(NULL)
      
      df_filtered[, pathway_name := pathway_name]
      df_filtered[, pathway_file := tools::file_path_sans_ext(basename(file))]
      
      # Create unified ID column (KEGG preferred if available)
      df_filtered[, is_significant := 
                    (`KEGG ID` %in% search_metabolites) |
                    (`HMDB ID` %in% search_metabolites) |
                    (`Metabolite Name` %in% search_metabolites)
      ]
      df_filtered <- df_filtered[, .(
        pathway_name,
        is_significant,
        KEGG_ID = `KEGG ID`,
        HMDB_ID = `HMDB ID`,
        Metabolite = `Metabolite Name`,
        pathway_file
      )]
      
      return(df_filtered)
      
    } else {
      return(NULL)
    }
  })
  
  final_results <- rbindlist(results, fill = TRUE)
  saveRDS(final_results, file = paste0(PATH,"/",FOLDER,"/SMPDBsearch_metabos_",NAME, ".rds"))
  write.csv(final_results, file = paste0(PATH,"/",FOLDER,"/SMPDBsearch_metabos_",NAME, ".csv"))
  return(final_results)
  
  
  
}



# result_df <- search_metabolite_in_smpdb(PATH= targeted_pwd, FOLDER="CDHFD_f",NAME="CDHFD_f",      search_metabolites =  CDHFD_f, 
#                                         target_metabolites= targeted_metabolites,   path_library= path_smpdb) results was 0
#atuomatic search gave 0, manual vor isovaleryl carnitine gave SMP0126748
#SMP0126748 pathwy does not directly involve isovalerylcarninitnt
#manual for buturylcarnitine gave SMPDB
#SMP0000235  	pathwya does not diretly involve buturylcarinitine but is associated with its elveation in disease
#manually put them in SMPD metabolite search file.
# It does not contain Buturylcarinitne but elevation of buturylcarinitine is linke to deficiency in this pathways.Short chain Acyl Coa dehydrogenase deficiency
result_df <- search_metabolite_in_smpdb(PATH= targeted_pwd, FOLDER="ND_f",NAME="ND_f",      search_metabolites = ND_f, 
                                         target_metabolites= targeted_metabolites,   path_library= path_smpdb) 
#results was 0 but manual smpdb search gave SMP012328 ( in my downloaded files are only until SMP0119305)
#SMP012328 Creatine Creatinine Glycine Guanidoacetic acid L-Arginine Ornithine S-Adenosylhomocysteine, S-Adenosylmethionine Water
### Result was 0
### Set up a list like i did for targteed with the metabolites in it and get ids
#result_df <- search_metabolite_in_smpdb(met_diff_ND_CDHFD, target_metabolites, path_smpdb)### 


# result_df <- readRDS(paste0(dirname(path_smpdb), "/SMDPB_Metabolite_search_NDvsCDHFD.rds"))