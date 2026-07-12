
gc()
rm(list = ls())
home <- normalizePath("~") # bc Windows does not start inuser dir but in user/documents dir
parent <- dirname(home)
setwd(paste0(parent,"/OneDrive - Universität Salzburg/Freigegebene Dokumente - AG_Tumorimmunologie/Data/Freia Krause/01_Experiments/FK46_FK49_RProject"))

library(MetaboAnalystR)
library(dplyr)
library(tidyr)
library(ggplot2)
library(readr)
library(data.table)
source("FK49_Definitions.R")

#### Find out which datatabase to choose to analyse metabolite set enrichement and analyse pathways etc.
#### Plot Metabolite Sets from interesting pathways
#### Rank based enrichment analysis
#### Code from MetaboAnalyst online 
### Metabolite liste targeted, -LOG FC ND vs CDHFD

#### Not done YET #### ----

### Enrichemt of Metabolite Sets based on the total Metabolite-Conc Table
### I was confused why this gives enriched set or what is computed behind the code so 
# So i tried enrichement analysis by only suplying oall the analysed metabolites. 
# The enriched sets differ (at least in sequence of appearence and strength in enrichmen plot)
# So it takes into account the concentrations of metabolites. Still, I am confused if this is correct.
#I submitted csv file with sample names, Treatment and Metabolites with concentrations.
# Analysis is down with all diff emtabolites, since metboalite up or down does not necessary correalte to pathway up or down. 
#So all diff metabolites could inform about diff regulated pathways
# goves users dir works for me since i aways have the data from work in the onedrive at the same location and only local home dir changes
ExpId = "FK49"
targeted_pwd        <-paste0(parent,"/OneDrive - Universität Salzburg/Freigegebene Dokumente - AG_Tumorimmunologie/Data/Freia Krause/01_Experiments/FK49_CD-HFD_13wks/FK49_Analysis/02_GeneratedData/targetedLivMetabolome")
untargeted_pwd      <-paste0(parent,"/OneDrive - Universität Salzburg/Freigegebene Dokumente - AG_Tumorimmunologie/Data/Freia Krause/01_Experiments/FK49_CD-HFD_13wks/FK49_Analysis/02_GeneratedData/untargetedLivMetabolome")
output_folders <- c("CDHFD_a", 
                    "CDHFD_m", 
                    "CDHFD_f", 
                    "ND_f",
                    "NDvsCDHFD")
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

#### 
do_MetaboMapping<-function(PATH= NULL, FOLDER = NULL, NAME = NULL){
  setwd(paste0(PATH,"/",FOLDER,"/"))
 metabolite_vector<-readRDS(paste0(PATH,"/",FOLDER,"/",NAME,".rds"))
 message("I read ",paste0(PATH,"/",FOLDER,"/",NAME,".rds"))
 trend_metabolites<-metabolite_vector%>%filter(trend == TRUE)
 trend_metabolites<-trend_metabolites[["Metabolite"]]
 print(head(trend_metabolites))
 mSet<-InitDataObjects("conc", "msetora", FALSE, 150)
 mSet<-Setup.MapData(mSet, trend_metabolites);
 mSet<-CrossReferencing(mSet, "name");
 return(mSet)
}
do_MetaboPrinting<-function(Set=NULL,PATH= NULL, FOLDER = NULL, NAME = NULL){
  setwd(paste0(PATH,"/",FOLDER,"/"))
  Set<-CreateMappingResultTable(Set)
  message("I created Mapping")
  Set<-SetMetabolomeFilter(Set, F);
  message("I set Filter")
  Set<-SetCurrentMsetLib(Set, "smpdb_pathway", 2);
  message("I defined smpdb as library")
  Set<-CalculateHyperScore(Set)
  message("I calculates HyperScore")
  Set<-PlotORA(Set, paste0(NAME,"_ora_SMPDB_1"), "net", "png", 300, width=NA)
  message("I saved image 1")
  Set<-PlotEnrichDotPlot(Set, "ora", paste0(NAME,"_ora_SMPDB_2"), "png", 300, width=NA)
  message("I saved image 2")
}
# Targeted CDHFD a ---- 
#targeted_mSET_CDHFD_a<-do_MetaboMapping(PATH = targeted_pwd, FOLDER= "CDHFD_a", NAME = "CDHFD_a") targeted_CDHFD_a has no trend or sig metabolites

# Targeted CDHFD m ---- 
#targeted_mSET_CDHFD_m<-do_MetaboMapping(PATH = targeted_pwd, FOLDER= "CDHFD_m", NAME = "CDHFD_m") #targeted_CDHFD_m has no trend or sig metabolites

# Targeted CDHFD f ---- 
targeted_mSET_CDHFD_f<-do_MetaboMapping(PATH = targeted_pwd, FOLDER= "CDHFD_f", NAME = "CDHFD_f") #2 Metbolites, names matched
targeted_mSET_CDHFD_f$name.map   #OK
do_MetaboPrinting(Set = targeted_mSET_CDHFD_f, PATH = targeted_pwd, FOLDER= "CDHFD_f", NAME = "CDHFD_f") #"No match was found to the selected metabolite set library! "SMPDB
rm(targeted_mSET_CDHFD_f)

# Targeted ND ---- 
targeted_mSET_ND_f<-do_MetaboMapping(PATH = targeted_pwd, FOLDER= "ND_f", NAME = "ND_f") #Status 1 ok check manually
targeted_mSET_ND_f$name.map      #OK
do_MetaboPrinting(Set = targeted_mSET_ND_f,    PATH = targeted_pwd, FOLDER= "ND_f", NAME = "ND_f") #No match was found to the selected metabolite set library!"SMPDB
rm(targeted_mSET_ND_f)


# Targeted ND vs CDHFD---- 
targeted_mSET_NDvsCDHFD<-do_MetaboMapping(PATH = targeted_pwd, FOLDER= "NDvsCDHFD", NAME = "NDvsCDHFD") #Status 1 ok check manually
targeted_mSET_NDvsCDHFD$name.map #OK
do_MetaboPrinting(Set = targeted_mSET_NDvsCDHFD, PATH = targeted_pwd, FOLDER= "NDvsCDHFD", NAME = "NDvsCDHFD")
rm(targeted_mSET_NDvsCDHFD)

# Untargeted CDHFD a positive -----
untargeted_mSET_CDHFD_a_p<-do_MetaboMapping(PATH = untargeted_pwd, FOLDER= "CDHFD_a", NAME = "CDHFD_a_p") 
#Over Half of List could not be matched also happens if i tryy to match cmpd id 
namesCDHFD_a_p<-untargeted_mSET_CDHFD_a_p$name.map 
#mSet<-InitDataObjects("conc", "msetora", FALSE, 150)
#cmpd.vec<-c("HMDB0001488","HMDB0000378","HMDB0248109","HMDB0257555","HMDB0254596","HMDB0258477")
#mSet<-Setup.MapData(mSet, cmpd.vec); only nicotinic acid and 2-Methylbutyroylcarnitine are matched
#mSet<-CrossReferencing(mSet, "hmdb");
#mSet<-CreateMappingResultTable(mSet)
do_MetaboPrinting(Set = untargeted_mSET_CDHFD_a_p, PATH = untargeted_pwd, FOLDER= "CDHFD_a", NAME = "CDHFD_a_p")# only Nictinate metabolism is overreprestend...
rm(untargeted_mSET_CDHFD_a_p)

# Untargeted CDHFD a negative -----
untargeted_mSET_CDHFD_a_n<-do_MetaboMapping(PATH = untargeted_pwd, FOLDER= "CDHFD_a", NAME = "CDHFD_a_n") #Status 1 ok check manually
namesCDHFD_a_n<-untargeted_mSET_CDHFD_a_n$name.map #1 matched 1 not same when i do it with HMDB IDs HMDB0000755 HMDB0255319
do_MetaboPrinting(Set = untargeted_mSET_CDHFD_a_n, PATH = untargeted_pwd, FOLDER= "CDHFD_a", NAME = "CDHFD_a_n")# #No match was found to the selected metabolite set library!"SMPDB
rm(untargeted_mSET_CDHFD_a_n)


# Untargeted CDHFD male positive -----
untargeted_mSET_CDHFD_m_p<-do_MetaboMapping(PATH = untargeted_pwd, FOLDER= "CDHFD_m", NAME = "CDHFD_m_p") # Over Half of List could not be matched
namesCDHFD_m_p<-untargeted_mSET_CDHFD_m_p$name.map 
#HMDB0251344 HMDB0133970 HMDB0259954 HMDB0245260 HMDB0000468 HMDB0014975 HMDB0250083 HMDB0249037 HMDB0252098 HMDB0250763 HMDB0245139 
#HMDB0254717 HMDB0000244 HMDB0245134 HMDB0242353 HMDB0257557 HMDB0258477
#over half of list could not be matched even if i use MEtaboanalyst webserver and HMDB IDs  
#Matched Nitrophenol, Biopterin, Pragabide, D-Histidine, Riboflavin
do_MetaboPrinting(Set = untargeted_mSET_CDHFD_m_p, PATH = untargeted_pwd, FOLDER= "CDHFD_m", NAME = "CDHFD_m_p")# #No match was found to the selected metabolite set library!"SMPDB
#"over represetend"riboflavin metabolism

# Untargeted CDHFD male negative -----
untargeted_mSET_CDHFD_m_n<-do_MetaboMapping(PATH = untargeted_pwd, FOLDER= "CDHFD_m", NAME = "CDHFD_m_n") # Over Half of List could not be matched
namesCDHFD_m_n<-untargeted_mSET_CDHFD_m_n$name.map 
#HMDB0246484 HMDB0000714 HMDB0000739 HMDB0255153 HMDB0246459 HMDB0010207HMDB0255319 same results if i match HMDB in the webserver
do_MetaboPrinting(Set = untargeted_mSET_CDHFD_m_n, PATH = untargeted_pwd, FOLDER= "CDHFD_m", NAME = "CDHFD_m_n")# #No match was found to the selected metabolite set library!"SMPDB
rm(untargeted_mSET_CDHFD_m_n)

# Untargeted CDHFD female positive -----
untargeted_mSET_CDHFD_f_p<-do_MetaboMapping(PATH = untargeted_pwd, FOLDER= "CDHFD_f", NAME = "CDHFD_f_p") #half could not be matched
namesCDHFD_f_p<-untargeted_mSET_CDHFD_f_p$name.map
#HMDB0000378 HMDB0001488 HMDB0242723 HMDB0251500 HMDB0000742 HMDB0253857 HMDB0244983 #half could not be matched also using  HMDB IDs
#HMDB0247693 HMDB0248109 HMDB0257555 HMDB0258477 HMDB0254596 
do_MetaboPrinting(Set = untargeted_mSET_CDHFD_f_p, PATH = untargeted_pwd, FOLDER= "CDHFD_f", NAME = "CDHFD_f_p")# #No match was found to the selected metabolite set library!"SMPDB

# Untargeted CDHFD female negative -----
untargeted_mSET_CDHFD_f_n<-do_MetaboMapping(PATH = untargeted_pwd, FOLDER= "CDHFD_f", NAME = "CDHFD_f_n") #Status 1 ok check manually
namesCDHFD_f_n<-untargeted_mSET_CDHFD_f_n$name.map 
#HMDB0255319 HMDB00613 HMDB0000783 HMDB0000755
do_MetaboPrinting(Set = untargeted_mSET_CDHFD_f_n, PATH = untargeted_pwd, FOLDER= "CDHFD_f", NAME = "CDHFD_f_n")# #No match was found to the selected metabolite set library!"SMPDB
rm(untargeted_mSET_CDHFD_f_n)

# Untargeted ND positive -----

#untargeted_mSET_ND_f_p<-do_MetaboMapping(PATH = untargeted_pwd, FOLDER= "ND_f", NAME = "ND_f_p") # Over Half of List could not be matched
#namesND_f_p   <-untargeted_mSET_ND_f_p$name.map 
# Untargeted ND negative -----

#untargeted_mSET_ND_f_n<-do_MetaboMapping(PATH = untargeted_pwd, FOLDER= "ND_f", NAME = "ND_f_n") # Over Half of List could not be matched
#namesND_f_n     <-untargeted_mSET_ND_f_n$name.map 
# Untargeted ND vs CDHFD positive -----

untargeted_mSET_NDvsCDHFD_p<-do_MetaboMapping(PATH = untargeted_pwd, FOLDER= "NDvsCDHFD", NAME = "NDvsCDHFD_p") ## Over Half of List could not be matched
namesNDvsCDHFD_p<-untargeted_mSET_NDvsCDHFD_p$name.map
# Untargeted ND vs CDHFD negative -----
untargeted_mSET_NDvsCDHFD_n<-do_MetaboMapping(PATH = untargeted_pwd, FOLDER= "NDvsCDHFD", NAME = "NDvsCDHFD_n") ## Over Half of List could not be matched
namesNDvsCDHFD_n<-untargeted_mSET_NDvsCDHFD_n$name.map 




 
### KEGG does not work gives this error. ####
# ND_CDHFD_down_Set<-CalculateHyperScore(ND_CDHFD_down_Set)
# [1] "Loaded files from MetaboAnalyst web-server."
# [1] "Failed to connect to the API Server!"
# [1] "Error! Mset ORA via xialab.ca/api unsuccessful!"
# ND_CDHFD_down_Set<-SetCurrentMsetLib(ND_CDHFD_down_Set, "kegg_pathway", 2);
# ND_CDHFD_down_Set<-CalculateHyperScore(ND_CDHFD_down_Set)
# ND_CDHFD_down_Set<-PlotORA(ND_CDHFD_down_Set, "NDvsCDHFD_down_ora_kegg1", "net", "png", 150, width=NA)
# ND_CDHFD_down_Set<-PlotEnrichDotPlot(ND_CDHFD_down_Set, "ora", "NDvsCDHFDdown_ora_kegg2", "png", 150, width=NA)

#### Plot Metabolites from ORA Pathways ####

# I now have DotPlot of overepresented pathways for example for CDHFD vs ND
# I downloaded https://smpdb.ca/downloads "Metabolite names linked to SMPDB pathwys CSV" 
# Now in FK49/HILIC09 I have folder with one csv per pathway with pwd realted metabolites. 
# Metabolites in Rows, infos in Cols, Pathway name as "Pathway Name" metabolite name: "Metabolite Name", "KEGG ID"
# I now want to read the pathways. On plots I have the Pathway Name but not the SMPDB ID. 
#First I want to give a lit of pathway names, that should be serach in the csv files. 
#Then per pathwy i want to access csv of the correct name, compare metabolites with 
#my metabolutes in the metabolome datafram and plot violin of abundance of these metabolites. Maybe also state how many metabolites I have of that pathwy and howw many are there 3/35

# I do these thing in FK49_Metabolomics_3 and 4 with mapping and plotting.

#### This is not working yet work on it later -----
  


#### Not done YET #### ----
setwd("C:/Users/b1084855/OneDrive - Universität Salzburg/Freigegebene Dokumente - AG_Tumorimmunologie/Data/Freia Krause/01_Experiments/FK49_CD-HFD_13wks/FK49_Analysis/02_GeneratedData/targetedLivMetabolome/Metaboanalyst")
library(MetaboAnalystR)
#Pathway Analysis with Concntration table
fname <- paste0("C:/Users/b1084855/OneDrive - Universität Salzburg/Freigegebene Dokumente - AG_Tumorimmunologie/Data/Freia Krause/01_Experiments/FK49_CD-HFD_13wks/FK49_HILIC09/Report_M086_HILIC09_20251222_forWebMetaboanalyst.csv")
mSet<-InitDataObjects("conc", "pathqea", FALSE, 150)
mSet<-Read.TextData(mSet, "Replacing_with_your_file_path", "rowu", "disc");
mSet<-SanityCheckData(mSet)
mSet<-PerformSanityClosure (mSet);
mSet<-CrossReferencing(mSet, "name");
mSet<-CreateMappingResultTable(mSet)
mSet<-PerformDetailMatch(mSet, "3-Hydroxybutyrate");
mSet<-GetCandidateList(mSet);
mSet<-SetCandidate(mSet, "3-Hydroxybutyrate", "3-Hydroxybutyric acid");
mSet<-PreparePrenormData(mSet)
mSet<-Normalization(mSet, "SumNorm", "Log2Norm", "AutoNorm", ratio=FALSE, ratioNum=20)
mSet<-PlotNormSummary(mSet, "norm_0_", "png", 150, width=NA)
mSet<-PlotSampleNormSummary(mSet, "snorm_0_", "png", 150, width=NA)
mSet<-SetSMPDB.PathLib(mSet, "hsa")
mSet<-SetOrganism(mSet, "hsa")
mSet<-SetMetabolomeFilter(mSet, F);
mSet<-CalculateQeaScore(mSet, "rbc", "gt")
mSet<-PlotPathSummary(mSet, F, "path_view_0_", "png", 150, width=NA, NA, NA )
mSet<-SaveTransformedData(mSet)
