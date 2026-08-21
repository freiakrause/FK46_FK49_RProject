#FK49_Color definitions
Sex_colors <- c("male"   = "#A2C2E5",
                "female" = "#F5A9B8")
Sex_shape <- c("male"   = 22,
                "female" = 25)

Treatment_colors <- c("none" = "black",
                      "Ctrl" = "#4D4D4DBF", # If etoh is replaced by Ctrl, which is the goal
                      "EtOH" = "#4D4D4DBF", 
                      "TAM"  = "#8B0000BF"
                      )

T_D_S_colors <- c("EtOH_CDHFD13_female"  = "#F5A9B8",
                  "EtOH_CDHFD13_male"    = "#A2C2E5",
                  "TAM_CDHFD13_female"   = "#F9D1D1",
                  "TAM_CDHFD13_male"     = "#B2A2D2",
                  "TAM_ND_female"        = "#F9B2B2",
                  "EtOH_ND_female"       = "#F5A9B8", 
                  "EtOH_ND_male"         = "#A2C2E5",  
                  "TAM_ND_male"          = "#B2A2D2"  ,
                  #if script does not spcecific CDHFD after 13wks but only CDHFD
                  "EtOH_CDHFD_female"  = "#F5A9B8",
                  "EtOH_CDHFD_male"    = "#A2C2E5",
                  "TAM_CDHFD_female"   = "#F9D1D1",
                  "TAM_CDHFD_male"     = "#B2A2D2",
                  "TAM_ND_female"        = "#F9B2B2",
                  "EtOH_ND_female"       = "#F5A9B8", 
                  "EtOH_ND_male"         = "#A2C2E5",  
                  "TAM_ND_male"          = "#B2A2D2")


T_S_colors <- c("EtOH_female" = "#F5A9B8", 
                "TAM_female"  = "#F9D1D1",  
                "EtOH_male"   = "#A2C2E5",
                "TAM_male"    = "#B2A2D2" )


Diet_colors <- c("CDHFD" = "darkviolet" ,
          "CDHFD13" = "darkviolet" , 
          "ND"      = "darkorange3")
Batch_colors<-c("1"= "blue",
                "2"= "red")

# Function to create directories
create_output_folders <- function(base_path, folders) {
  for (folder in folders) {
    # Full path to the output folder
    folder_path <- file.path(base_path, folder)
    
    # Check if the folder exists, and if not, create it
    if (!dir.exists(folder_path)) {
      dir.create(folder_path, recursive = TRUE)
      cat("Created folder:", folder_path, "\n")
    } else {
      cat("Folder already exists:", folder_path, "\n")
    }
  }
}
#Pathways
home <- normalizePath("~") # bc Windows does not start inuser dir but in user/documents dir
parent <- dirname(home) 
PATHS <- list(
  targetedMet = list(
    output = paste0(parent,"/OneDrive - Universität Salzburg/AG_Tumorimmunologie - Dokumente/Data/Freia Krause/01_Experiments/FK49_CD-HFD_13wks/FK49_Analysis/02_GeneratedData/targetedLivMetabolome")
  ),
  untargeted_Met = list(
    input = paste0(parent,"/OneDrive - Universität Salzburg/AG_Tumorimmunologie - Dokumente/Data/Freia Krause/01_Experiments/FK49_CD-HFD_13wks/FK49_Analysis/02_GeneratedData/untargetedLivMetabolome")
    ),
  BA = list(
    input =  paste0(parent,"/OneDrive - Universität Salzburg/AG_Tumorimmunologie - Dokumente/Data/Freia Krause/01_Experiments/FK49_CD-HFD_13wks/FK49_BA/Report_M087_BACID01_20251222_withMeta.csv"),
    output = paste0(parent,"/OneDrive - Universität Salzburg/AG_Tumorimmunologie - Dokumente/Data/Freia Krause/01_Experiments/FK49_CD-HFD_13wks/FK49_Analysis/02_GeneratedData/BA")
  ),
  microbiome = list(
    input  = "D:/Data/Experiment2/Input",
    output = "D:/Data/Experiment2/Output"
  ),
  proteomics = list(
    input  = "D:/Data/Experiment2/Input",
    output = "D:/Data/Experiment2/Output"
  ),
  TEM = list(
    input  = paste0(parent, "/OneDrive - Universität Salzburg/AG_Tumorimmunologie - Dokumente/Data/Freia Krause/01_Experiments/FK49_CD-HFD_13wks/FK49_Analysis/02_GeneratedData/FK49_TEM/QuPath"),
    output = paste0(parent, "/OneDrive - Universität Salzburg/AG_Tumorimmunologie - Dokumente/Data/Freia Krause/01_Experiments/FK49_CD-HFD_13wks/FK49_Analysis/02_GeneratedData/FK49_TEM")
  ),
  legendplex = list(
    FK49_input = paste0(parent,"/OneDrive - Universität Salzburg/AG_Tumorimmunologie - Dokumente/Data/Freia Krause/01_Experiments/FK49_CD-HFD_13wks/FK49_Legendplex/02_generated"),
    FK46_input = paste0(parent,"/OneDrive - Universität Salzburg/AG_Tumorimmunologie - Dokumente/Data/Freia Krause/01_Experiments/FK46_iALmice_high Fat diet 52 weeks 7d after injection/FK46_Legendplex"),
    FK46_output =paste0(parent,"/OneDrive - Universität Salzburg/AG_Tumorimmunologie - Dokumente/Data/Freia Krause/01_Experiments/FK46_iALmice_high Fat diet 52 weeks 7d after injection/Analysis/02_GeneratedData/Legendplex"),
    FK49_output =paste0(parent,"/OneDrive - Universität Salzburg/AG_Tumorimmunologie - Dokumente/Data/Freia Krause/01_Experiments/FK49_CD-HFD_13wks/FK49_Analysis/02_GeneratedData/Legendplex")
  ),
  MASH = list(
    FK49_output  =paste0(parent,"/OneDrive - Universität Salzburg/AG_Tumorimmunologie - Dokumente/Data/Freia Krause/01_Experiments/FK49_CD-HFD_13wks/FK49_Analysis/02_GeneratedData/NASH_Score"),
    FK49_input = paste0(parent,"/OneDrive - Universität Salzburg/AG_Tumorimmunologie - Dokumente/Data/Freia Krause/01_Experiments/FK49_CD-HFD_13wks/FK49_Analysis/01_RawData"),
    FK46_output = paste0(parent,"/OneDrive - Universität Salzburg/AG_Tumorimmunologie - Dokumente/Data/Freia Krause/01_Experiments/FK46_iALmice_high Fat diet 52 weeks 7d after injection/Analysis/02_GeneratedData/NASH_Score"),
    FK46_input = paste0(parent,"/OneDrive - Universität Salzburg/AG_Tumorimmunologie - Dokumente/Data/Freia Krause/01_Experiments/FK49_CD-HFD_13wks/FK49_Analysis/01_RawData")
  ),
  exigo = list(
    FK46_input = paste0(parent,"/OneDrive - Universität Salzburg/AG_Tumorimmunologie - Dokumente/Data/Freia Krause/01_Experiments/FK46_iALmice_high Fat diet 52 weeks 7d after injection/Analysis/01_RawData"),
    FK46_output = paste0(parent,"/OneDrive - Universität Salzburg/AG_Tumorimmunologie - Dokumente/Data/Freia Krause/01_Experiments/FK46_iALmice_high Fat diet 52 weeks 7d after injection/Analysis/02_GeneratedData/Exigo"),
    FK49_input = paste0(parent,"/OneDrive - Universität Salzburg/AG_Tumorimmunologie - Dokumente/Data/Freia Krause/01_Experiments/FK49_CD-HFD_13wks/FK49_Analysis/01_RawData"),
    FK49_output  =paste0(parent,"/OneDrive - Universität Salzburg/AG_Tumorimmunologie - Dokumente/Data/Freia Krause/01_Experiments/FK49_CD-HFD_13wks/FK49_Analysis/02_GeneratedData/Exigo/FK49"),
    BH15_input= paste0(parent,"/OneDrive - Universität Salzburg/AG_Tumorimmunologie - Dokumente/Data/Freia Krause/01_Experiments/FK49_CD-HFD_13wks/FK49_Analysis/01_RawData"),
    BH15_output= paste0(parent,"/OneDrive - Universität Salzburg/AG_Tumorimmunologie - Dokumente/Data/Freia Krause/01_Experiments/FK49_CD-HFD_13wks/FK49_Analysis/02_GeneratedData/Exigo/BH15")
    
    ),
  general_data = list(
    FK49_output = paste0(parent,"/OneDrive - Universität Salzburg/AG_Tumorimmunologie - Dokumente/Data/Freia Krause/01_Experiments/FK49_CD-HFD_13wks/FK49_Analysis"),
    FK46_output = paste0(parent,"/OneDrive - Universität Salzburg/AG_Tumorimmunologie - Dokumente/Data/Freia Krause/01_Experiments/FK46_iALmice_high Fat diet 52 weeks 7d after injection/Analysis"),
    FK46_input = paste0(parent,"/OneDrive - Universität Salzburg/AG_Tumorimmunologie - Dokumente/Data/Freia Krause/01_Experiments/FK46_iALmice_high Fat diet 52 weeks 7d after injection/FK46_Organ_GeneralMouseData"),
    FK49_input = paste0(parent,"/OneDrive - Universität Salzburg/AG_Tumorimmunologie - Dokumente/Data/Freia Krause/01_Experiments/FK49_CD-HFD_13wks")
    ),
 FoodIntake = list(
    output = paste0(parent,"/OneDrive - Universität Salzburg/AG_Tumorimmunologie - Dokumente/Data/Freia Krause/01_Experiments/FK49_CD-HFD_13wks/FK49_Analysis/02_GeneratedData/FoodIntake"),
    input = paste0(parent,"/OneDrive - Universität Salzburg/AG_Tumorimmunologie - Dokumente/Data/Freia Krause/01_Experiments/FK49_CD-HFD_13wks/FK49_Analysis/01_RawData")
  ),
 TAG = list(
   output = paste0(parent,"/OneDrive - Universität Salzburg/AG_Tumorimmunologie - Dokumente/Data/Freia Krause/01_Experiments/FK49_CD-HFD_13wks/FK49_Analysis/02_GeneratedData/Tumor_Ascites_Granuloma"),
   input = paste0(parent,"/OneDrive - Universität Salzburg/AG_Tumorimmunologie - Dokumente/Data/Freia Krause/01_Experiments/FK49_CD-HFD_13wks/FK49_Analysis/01_RawData")
 ),
 Organs = list(
   output = paste0(parent,"/OneDrive - Universität Salzburg/AG_Tumorimmunologie - Dokumente/Data/Freia Krause/01_Experiments/FK49_CD-HFD_13wks/FK49_Analysis/02_GeneratedData/Weight_Organs"),
   input = paste0(parent,"/OneDrive - Universität Salzburg/AG_Tumorimmunologie - Dokumente/Data/Freia Krause/01_Experiments/FK49_CD-HFD_13wks/FK49_Analysis/01_RawData")
 ),
 BH_baseline = list(
   input = paste0(parent,"/OneDrive - Universität Salzburg/AG_Tumorimmunologie - Dokumente/Data/Freia Krause/01_Experiments/FK49_CD-HFD_13wks/FK49_Analysis/01_RawData")
 )
 
)
# targeted_in_pwd 
# untargeted_pwd      
# untargeted_in_pwd
# BApwd <- paste0(parent, "/OneDrive - Universität Salzburg/AG_Tumorimmunologie - Dokumente/Data/Freia Krause/01_Experiments/FK49_CD-HFD_13wks/FK49_Analysis/02_GeneratedData/BA")
# BA_in_pwd
# exigo_in_pwd
# exigo_ou_pwd
# Legendplex_in_pwd
# Legendplex_out_pwd
# NASH_in_pwd
# NASH_out_pwd
# microbiome_in_pwd
# microbiome_out_pwd
# prot_in_pwd
# prot_out_pwd
# general_data_in
# general_data_out
# TEM_out_pwd <- paste0(parent, "/OneDrive - Universität Salzburg/AG_Tumorimmunologie - Dokumente/Data/Freia Krause/01_Experiments/FK49_CD-HFD_13wks/FK49_Analysis/02_GeneratedData/FK49_TEM")
# TEM_in_pwd <- paste0(parent, "/OneDrive - Universität Salzburg/AG_Tumorimmunologie - Dokumente/Data/Freia Krause/01_Experiments/FK49_CD-HFD_13wks/FK49_Analysis/02_GeneratedData/FK49_TEM/QuPath")

subset_data <- function(data,
                        Sex_filter  = NULL,
                        Diet_filter = NULL,
                        ExpID_filter = NULL,
                        Timepoint_filter = NULL,
                        filters = list(),
                        method = NULL
) {
  
  # Old-style filters collected into a list
  legacy_filters <- list(
    Sex = Sex_filter,
    Diet = Diet_filter,
    ExpID = ExpID_filter,
    Timepoint = Timepoint_filter
  )
  
  # Remove NULLs from legacy filters
  legacy_filters <- legacy_filters[!sapply(legacy_filters, is.null)]
  
  # Merge: new `filters` overrides old-style arguments if both provided
  filters <- modifyList(legacy_filters, filters)
  cat("Applying filters:\n")
  print(filters)
  dataset <- data
  
  for (col in names(filters)) {
    dataset <- dataset %>%
      dplyr::filter(.data[[col]] %in% filters[[col]])
  }
  
  # Metadata columns
  meta_cols <- c( "Sample", "Animal", "Sex", "Treatment", "Diet", "ExpID", "T_D_S", "T_D", "T_S","Time_Treat","Timepoint")
  
  # Extract metabolite matrix
  numerical_data <- dataset %>%
    dplyr::select(-all_of(meta_cols)) %>%
    dplyr::select(where(is.numeric))
  
  if(method == "untargeted"){
    ####### Article about Centering, Scaling, Transfor in metabolomics
    #https://link.springer.com/article/10.1186/1471-2164-7-142 van den Berg, R.A., Hoefsloot, H.C., Westerhuis, J.A. et al. Centering, scaling, and transformations: improving the biological information content of metabolomics data. BMC Genomics 7, 142 (2006). https://doi.org/10.1186/1471-2164-7-142
    ######
    nzv <- apply(numerical_data, 2, var, na.rm = TRUE) > 0
    numerical_data <- numerical_data[, nzv, drop = FALSE] # reduce dimensinality by removing features that have 0 variance
    
    rs <- rowSums(numerical_data, na.rm = TRUE)   # peak sum per sample for normalization
    rs[rs == 0] <- NA                              # when sum 0, norm fails so exchange to na so 
    
    data_norm <- numerical_data / rs
    data_log <- log2(data_norm + 1e-9) #heteroscedasticity in data, transformation removes this, log can not deal with 0 therefore small values added.
    data_scaled <- as.data.frame(scale(data_log, center= TRUE, scale = TRUE)) #method is autoscaling SD as scaling factor
    
  } else if (method == "targeted"){
    # Extract metabolite matrix
    numerical_data <- dataset %>%
      dplyr::select(-all_of(meta_cols)) %>%
      dplyr::select(where(is.numeric))
    
    # Remove zero-variance features
    nzv <- apply(numerical_data, 2, var, na.rm = TRUE) > 0
    numerical_data <- numerical_data[, nzv, drop = FALSE]
    
    # Optional scaling
    eps <- min(numerical_data[numerical_data > 0], na.rm = TRUE) / 2
    #eps hier noch addieren innerhalb des log()? Macht das einen unterschied?
    data_norm <- NULL
    data_log <- log2(numerical_data+eps)
    data_scaled <- as.data.frame(scale(data_log))
    
  } else{ print("You need to give me method 'targeted' or 'untargeted' so that I can perform correct preprocessing.")}
  
  #i am not sure if i should to 1.log transform 2 normalize 3 scale 
  # #           or if I should do 1.normalize 2. log trasnform 3 scale. 
  # # 1- log looks less skewed in the finally scaled data. But it might compress data and variances to much and i might lose signal
  raw_data <-  numerical_data%>%
    as.data.frame() %>%
    tidyr::pivot_longer(cols = everything(), names_to = "Metabolite", values_to = "Value") %>%
    dplyr::mutate(type = "raw")
  
  norm_data <- if (is.null(data_norm)) {
    message("Targeted data: no normalization applied")
    norm_data <- NULL
  } else {
    norm_data <- data_norm %>%
      as.data.frame() %>%
      tidyr::pivot_longer(cols = everything(),
                          names_to = "Metabolite",
                          values_to = "Value") %>%
      dplyr::mutate(type = "normalized")
  }
  
  log_data <- data_log %>%
    as.data.frame() %>%
    tidyr::pivot_longer(cols = everything(), names_to = "Metabolite", values_to = "Value") %>%
    dplyr::mutate(type = "log2-transformed")
  
  scaled_data <- data_scaled %>%
    as.data.frame() %>%
    tidyr::pivot_longer(cols = everything(), names_to = "Metabolite", values_to = "Value") %>%
    dplyr::mutate(type = "autoscaled")
  
  # Combine all data into one data frame
  if (is.null(data_norm)) {
    combined_data <- bind_rows(raw_data, log_data, scaled_data) %>%
      dplyr::mutate(type = factor(type,
                                  levels = c("raw", "log2-transformed", "autoscaled")))
  } else {
    combined_data <- bind_rows(raw_data, norm_data, log_data, scaled_data) %>%
      dplyr::mutate(type = factor(type,
                                  levels = c("raw", "normalized", "log2-transformed", "autoscaled")))
  }
  
  # Plot histograms
  plot_d<- ggplot(combined_data, aes(x = Value, fill = type)) +
    geom_histogram(bins = 60, alpha = 0.6, position = "identity") +
    facet_wrap(~ type, scales = "free_x", ncol = 4) +
    scale_fill_manual(values = c("gray", "blue","violet", "green")) +
    labs(title = "Distribution of Data at Different Preprocessing Stages",
         x = "Value",
         y = "Frequency") +
    theme_minimal() +
    theme(legend.position = "none")
  
  print(plot_d)
  
  return(list(Preprocessing =plot_d,
              raw_values  = numerical_data,
              log_values = data_log,
              norm_values = data_norm,
              scaled_log_values = data_scaled,
              metadata   = dataset[, meta_cols, drop = FALSE]
  ))
}


PARAMETERS<-list(
  BA = list(
    BA_sort = c("CA","GCA","TCA",                                        #primary cholic acid and conjugates
                "DCA","TDCA","LCA","GDCA","GLCA","THDCA","TLCA","TUDCA", #secodnray cholic and conjugates
                  "CDCA","GCDCA","TCDCA",                                #primary Chenodeoxycholic acid and conjugates
                  "12-ketoCDCA","HDCA","UDCA",                           #secondray Chenodeoxycholic acid and conjugates
                  "alpha-MCA","beta-MCA","omega-MCA","TMCA"),           # muri cholich acid
    BA_primary =c("CA","GCA","TCA",                                        #Cholic acid and conjugates
                  "CDCA","GCDCA","TCDCA",                                    #Chenodeoxychiolic acid and conjugates
                  "alpha-MCA","beta-MCA","omega-MCA","TMCA") ,    # Muricholic acids and conjugetes (only murine not human),
    BA_secondary =c("DCA", "HDCA","LCA","UDCA","GDCA","GLCA","TDCA","THDCA","TLCA","TUDCA","12-ketoCDCA"),
    #Status Origin con
    BA_primary_uncon =c("CA", "CDCA","alpha-MCA","beta-MCA","omega-MCA") , 
    BA_primary_con=c("GCA","TCA","GCDCA","TCDCA" ,"TMCA"),
    BA_secondary_uncon =c("DCA", "HDCA","LCA","UDCA","12-ketoCDCA"),
    BA_secondary_con =c("GDCA","GLCA","TDCA","THDCA","TLCA","TUDCA"),
    #Status con un con
    BA_uncon =c("CA","DCA","CDCA", "alpha-MCA","beta-MCA","omega-MCA","LCA","UDCA","HDCA","12-ketoCDCA"),
    BA_con  =c("GCA","TCA","GCDCA","TCDCA","GDCA","TDCA","GLCA","TLCA","TUDCA","THDCA","TMCA"),
    meta_cols =c("Animal","Sample","Sex","Treatment","Diet","ExpID","T_D_S","T_D" ,"T_S")
   ),
  EXIGO = list(
    FK49_Exigo_Comprehensive_Panel = list(
      list(value="ALB",  y_title="Alb [g/L]",   normal_range=c(20,48), lowlimit=2),
      list(value="TP",   y_title="TP [g/L]",    normal_range=c(36,66)),
      list(value="GLOB", y_title="GLOB [g/L]"),
      list(value="A.G",  y_title="A/G"),
      list(value="TB",   y_title="TB [µmol/L]", normal_range=c(1,15), lowlimit=0.1),
      list(value="GGT",  y_title="GGT [U/L]",   lowlimit=2),
      list(value="AST",  y_title="AST [U/L]",   normal_range=c(59,247), hilimit=650, lowlimit=5),
      list(value="ALT",  y_title="ALT [U/L]",   normal_range=c(28,132)),
      list(value="ALP",  y_title="ALP [U/L]",   normal_range=c(62,209), lowlimit=5),
      list(value="AMY",  y_title="AMY [U/L]",   normal_range=c(1691,3615)),
      list(value="Crea", y_title="Crea [?]",    normal_range=c(12,71)),
      list(value="UA",   y_title="UA [µmol/L]", normal_range=c(101,321), lowlimit=10),
      list(value="BUN",  y_title="BUN [mmol/L]",normal_range=c(4,11.8)),
      list(value="GLU",  y_title="GLU [mmol/L]",normal_range=c(5,10.67)),
      list(value="TC",   y_title="TC [mmol/L]", normal_range=c(0.93,4.04)),
      list(value="TG",   y_title="TG [mmol/L]", normal_range=c(0.62,1.63))),
    
    FK46_Exigo_Liver_Panel = list(
      list(value="ALB",  y_title="Alb [g/L]",     normal_range=c(20,48), lowlimit=2),
      list(value="TP",   y_title="TP [g/L]",      normal_range=c(36,66)),
      list(value="GLOB", y_title="GLOB [g/L]"),
      list(value="A.G",  y_title="A/G"),
      list(value="TB",   y_title="TB [µmol/L]",   normal_range=c(0,15),  lowlimit=0.1),
      list(value="GGT",  y_title="GGT [U/L]",     lowlimit=2),
      list(value="AST",  y_title="AST [U/L]",     normal_range=c(59,247), hilimit=650, lowlimit=5),
      list(value="ALT",  y_title="ALT [U/L]",     normal_range=c(28,132)),
      list(value="ALP",  y_title="ALP [U/L]",     normal_range=c(62,209), lowlimit=5),
      list(value="TBA",  y_title="TBA [µmol/L]",  lowlimit=1),
      list(value="TC",   y_title="TC [mmol/L]",   normal_range=c(0.93,4.04))),
  
  FK49_Exigo_cols = c("ALB", "TP", "GLOB", "A.G", "TB", "GGT","AST", "ALT", "ALP", "AMY", "Crea", "UA", "BUN", "GLU", "TC", "TG"),
  FK46_Exigo_cols = c("ALB", "TP", "GLOB", "A.G", "TB", "GGT","AST", "ALT", "ALP", "TBA", "TC")
  ),
  Legendplex = list(
    cytokine_list= list(
      list(value = "IL23",    y_title = "IL-23 [pg/mL]"),
      list(value = "IL1a",    y_title = "IL-1α [pg/mL]"),
      list(value = "IFNy",    y_title = "IFN-γ [pg/mL]"),
      list(value = "TNFa",    y_title = "TNF-α [pg/mL]"),
      list(value = "MCP1",    y_title = "MCP-1 [pg/mL]"),
      list(value = "IL12p70", y_title = "IL-12p70 [pg/mL]"),
      list(value = "IL1ß",    y_title = "IL-1β [pg/mL]"),
      list(value = "IL10",    y_title = "IL-10 [pg/mL]"),
      list(value = "IL6",     y_title = "IL-6 [pg/mL]"),
      list(value = "IL27",    y_title = "IL-27 [pg/mL]"),
      list(value = "IL17A",   y_title = "IL-17A [pg/mL]"),
      list(value = "IFNß",    y_title = "IFN-β [pg/mL]"),
      list(value = "GMCSF",   y_title = "GM-CSF [pg/mL]")
    ),
    
    cytokines = c(
      "IL23", "IL1a", "IFNy", "TNFa", "MCP1", "IL12p70",
      "IL1ß", "IL10", "IL6", "IL27", "IL17A", "IFNß", "GMCSF"
    )
  )
  )


