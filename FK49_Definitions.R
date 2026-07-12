#FK49_Color definitions
Sex_colors <- c("male"   = "#A2C2E5",
                "female" = "#F5A9B8")

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
targeted_pwd        <-paste0(parent,"/OneDrive - Universität Salzburg/AG_Tumorimmunologie - Dokumente/Data/Freia Krause/01_Experiments/FK49_CD-HFD_13wks/FK49_Analysis/02_GeneratedData/targetedLivMetabolome")
untargeted_pwd      <-paste0(parent,"/OneDrive - Universität Salzburg/AG_Tumorimmunologie - Dokumente/Data/Freia Krause/01_Experiments/FK49_CD-HFD_13wks/FK49_Analysis/02_GeneratedData/untargetedLivMetabolome")
BApwd <- paste0(parent, "/OneDrive - Universität Salzburg/AG_Tumorimmunologie - Dokumente/Data/Freia Krause/01_Experiments/FK49_CD-HFD_13wks/FK49_Analysis/02_GeneratedData/BA")
TEM_out_pwd <- paste0(parent, "/OneDrive - Universität Salzburg/AG_Tumorimmunologie - Dokumente/Data/Freia Krause/01_Experiments/FK49_CD-HFD_13wks/FK49_Analysis/02_GeneratedData/FK49_TEM")
TEM_in_pwd <- paste0(parent, "/OneDrive - Universität Salzburg/AG_Tumorimmunologie - Dokumente/Data/Freia Krause/01_Experiments/FK49_CD-HFD_13wks/FK49_Analysis/02_GeneratedData/FK49_TEM/QuPath")
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



#### Bile Acids Sorted and Sorted into Primary and Secondary ones -----
BA_sort <- c("CA","GCA","TCA",#primary cholic acid and conjugates
             "DCA","TDCA","LCA","GDCA","GLCA","THDCA","TLCA","TUDCA", #secodnray cholic and conjugates
             "CDCA","GCDCA","TCDCA",#primary Chenodeoxycholic acid and conjugates
             "12-ketoCDCA","HDCA","UDCA", #secondray Chenodeoxycholic acid and conjugates
             "alpha-MCA","beta-MCA","omega-MCA","TMCA" # muri cholich acid
)                                                       
#Status Origin
BA_primary <-c("CA","GCA","TCA",                                        #Cholic acid and conjugates
               "CDCA","GCDCA","TCDCA",                                     #Chenodeoxychiolic acid and conjugates
               "alpha-MCA","beta-MCA","omega-MCA","TMCA")  
  # Muricholic acids and conjugetes (only murine not human)
BA_secondary <-c("DCA", "HDCA","LCA","UDCA","GDCA","GLCA","TDCA","THDCA","TLCA","TUDCA","12-ketoCDCA")

#Status Origin con
BA_primary_uncon <-c("CA", "CDCA","alpha-MCA","beta-MCA","omega-MCA")  
BA_primary_con<-c("GCA","TCA","GCDCA","TCDCA" ,"TMCA")
BA_secondary_uncon <-c("DCA", "HDCA","LCA","UDCA","12-ketoCDCA")
BA_secondary_con <-c("GDCA","GLCA","TDCA","THDCA","TLCA","TUDCA")
#Status con un con

BA_uncon <-c("CA","DCA","CDCA", "alpha-MCA","beta-MCA","omega-MCA","LCA","UDCA","HDCA","12-ketoCDCA")
BA_con  <-c("GCA","TCA", 
              "GCDCA","TCDCA",
              "GDCA","TDCA","GLCA","TLCA","TUDCA","THDCA","TMCA"              )