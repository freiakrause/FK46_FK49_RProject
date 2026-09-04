#FK49_Color definitions
Sex_colors <- c("male"   = "#A2C2E5",
                "female" = "#F5A9B8")
Sex_shape <- c("male"   = 22,
                "female" = 25)

Treatment_colors <- c("none" = "black",
                      "Ctrl" = "#4D4D4DBF",
                      "TAM"  = "#8B0000BF"
                      )

T_D_S_colors <- c("Ctrl_CDHFD13_female"  = "#F5A9B8",
                  "Ctrl_CDHFD13_male"    = "#A2C2E5",
                  "TAM_CDHFD13_female"   = "#F9D1D1",
                  "TAM_CDHFD13_male"     = "#B2A2D2",
                  "TAM_ND_female"        = "#F9B2B2",
                  "Ctrl_ND_female"       = "#F5A9B8",
                  "Ctrl_ND_male"         = "#A2C2E5",
                  "TAM_ND_male"          = "#B2A2D2"  ,
                  #if script does not spcecific CDHFD after 13wks but only CDHFD
                  "Ctrl_CDHFD_female"  = "#F5A9B8",
                  "Ctrl_CDHFD_male"    = "#A2C2E5",
                  "TAM_CDHFD_female"   = "#F9D1D1",
                  "TAM_CDHFD_male"     = "#B2A2D2",
                  "TAM_ND_female"        = "#F9B2B2",
                  "Ctrl_ND_female"       = "#F5A9B8",
                  "Ctrl_ND_male"         = "#A2C2E5",
                  "TAM_ND_male"          = "#B2A2D2")


T_S_colors <- c("Ctrl_female" = "#F5A9B8",
                "TAM_female"  = "#F9D1D1",
                "Ctrl_male"   = "#A2C2E5",
                "TAM_male"    = "#B2A2D2" )


Diet_colors <- c("CDHFD" = "darkviolet" ,
          "CDHFD13" = "darkviolet" ,
          "ND"      = "darkorange3")
Batch_colors<-c("1"= "blue",
                "2"= "red")

Phylum_colors <- c(
  "Firmicutes"       = "#B58AD0",
  "Bacteroidetes"    = "#D99A55",
  "Proteobacteria"   = "#72B56B",
  "Actinobacteria"   = "#6FA8D8",
  "Tenericutes"      = "#D88989",
  "Patescibacteria"  = "#7FB8C8",
  "Cyanobacteria"    = "#82B87A"
)
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
# Function to group close timepoints (used in FK46 preprocessing)
# Since batches were sometimes measured 1-2 days apart, this groups
# timepoints within a tolerance and replaces them with the mean
group_close_timepoints <- function(inputdata, tolerance = 0.4) {
  unique_times <- sort(unique(inputdata$wks_diet))
  groups <- list()

  # Create groups of close values
  while (length(unique_times) > 0) {
    ref <- unique_times[1]
    close_vals <- unique_times[abs(unique_times - ref) <= tolerance]
    groups[[length(groups) + 1]] <- close_vals
    unique_times <- setdiff(unique_times, close_vals)}

  # Create a lookup table for replacements
  replacements <- lapply(groups, function(g) {
    rep(mean(g), length(g))}) %>% unlist()

  value_map <- data.frame(
    original = unlist(groups),
    new = replacements)

  # Join with original data
  inputdata <- inputdata %>%
    left_join(value_map, by = c("wks_diet" = "original")) %>%
    mutate(wks_diet = ifelse(is.na(new), wks_diet, new)) %>%
    select(-new)%>%
    mutate(wks_diet = round(wks_diet,digits=1))
  return(inputdata)
}
home <- normalizePath("~") # bc Windows does not start inuser dir but in user/documents dir
parent <- dirname(home)
PATHS <- list(
  metabolomics = list(
    untargeted_xlsx = paste0(parent,"/OneDrive - Universität Salzburg/AG_Tumorimmunologie - Dokumente/Data/Freia Krause/01_Experiments/FK49_CD-HFD_13wks/FK49_untargetedLiverMetabolomics/Report_M086_untargeted_20260113.xlsx"),
    targeted_xlsx   = paste0(parent,"/OneDrive - Universität Salzburg/AG_Tumorimmunologie - Dokumente/Data/Freia Krause/01_Experiments/FK49_CD-HFD_13wks/FK49_HILIC09/Report_M086_HILIC09_20251222.xlsx"),
    rawdata         = paste0(parent,"/OneDrive - Universität Salzburg/AG_Tumorimmunologie - Dokumente/Data/Freia Krause/01_Experiments/FK49_CD-HFD_13wks/FK49_HILIC09"),
    output          = paste0(parent,"/OneDrive - Universität Salzburg/AG_Tumorimmunologie - Dokumente/Data/Freia Krause/01_Experiments/FK49_CD-HFD_13wks/FK49_Analysis/02_GeneratedData/metabolomics")

  ),
  BA = list(
    input =  paste0(parent,"/OneDrive - Universität Salzburg/AG_Tumorimmunologie - Dokumente/Data/Freia Krause/01_Experiments/FK49_CD-HFD_13wks/FK49_BA/Report_M087_BACID01_20251222_withMeta.csv"),
    output = paste0(parent,"/OneDrive - Universität Salzburg/AG_Tumorimmunologie - Dokumente/Data/Freia Krause/01_Experiments/FK49_CD-HFD_13wks/FK49_Analysis/02_GeneratedData/BA")
  ),
  microbiome = list(
    input_16S    = paste0(parent, "/OneDrive - Universität Salzburg/AG_Tumorimmunologie - Dokumente/Data/Freia Krause/01_Experiments/FK49_CD-HFD_13wks/FK49_Microbiome/16S"),
    input_meta   = paste0(parent, "/OneDrive - Universität Salzburg/AG_Tumorimmunologie - Dokumente/Data/Freia Krause/01_Experiments/FK49_CD-HFD_13wks/FK49_Microbiome/FK49_CD-HFD_13wks_Microbiome_Meta.csv"),
    output       = paste0(parent, "/OneDrive - Universität Salzburg/AG_Tumorimmunologie - Dokumente/Data/Freia Krause/01_Experiments/FK49_CD-HFD_13wks/FK49_Analysis/02_GeneratedData/FK49_Microbiome"),
    output_stats = paste0(parent, "/OneDrive - Universität Salzburg/AG_Tumorimmunologie - Dokumente/Data/Freia Krause/01_Experiments/FK49_CD-HFD_13wks/FK49_Analysis/02_GeneratedData/FK49_Microbiome/Statistics"),
    output_plots = paste0(parent, "/OneDrive - Universität Salzburg/AG_Tumorimmunologie - Dokumente/Data/Freia Krause/01_Experiments/FK49_CD-HFD_13wks/FK49_Analysis/02_GeneratedData/FK49_Microbiome/Plots")
  ),
  proteomics = list(
    input  = paste0(parent,"/OneDrive - Universität Salzburg/AG_Tumorimmunologie - Dokumente/Data/Freia Krause/01_Experiments/FK49_CD-HFD_13wks/FK49_FK46_Proteomics_Phosphoproteomics"),
    output = paste0(parent,"/OneDrive - Universität Salzburg/AG_Tumorimmunologie - Dokumente/Data/Freia Krause/01_Experiments/FK49_CD-HFD_13wks/FK49_Analysis/02_GeneratedData/FK49_Proteomics")
  ),
  TEM = list(
    input  = paste0(parent, "/OneDrive - Universität Salzburg/AG_Tumorimmunologie - Dokumente/Data/Freia Krause/01_Experiments/FK49_CD-HFD_13wks/FK49_Analysis/02_GeneratedData/FK49_TEM/QuPath"),
    output = paste0(parent, "/OneDrive - Universität Salzburg/AG_Tumorimmunologie - Dokumente/Data/Freia Krause/01_Experiments/FK49_CD-HFD_13wks/FK49_Analysis/02_GeneratedData/FK49_TEM")
  ),
  legendplex = list(
    FK49_input = paste0(parent,"/OneDrive - Universität Salzburg/AG_Tumorimmunologie - Dokumente/Data/Freia Krause/01_Experiments/FK49_CD-HFD_13wks/FK49_Legendplex/02_generated"),
    FK46_input = paste0(parent,"/OneDrive - Universität Salzburg/AG_Tumorimmunologie - Dokumente/Data/Freia Krause/01_Experiments/FK46_iALmice_high Fat diet 52 weeks 7d after injection/FK46_Legendplex/02_generated"),
    FK46_output =paste0(parent,"/OneDrive - Universität Salzburg/AG_Tumorimmunologie - Dokumente/Data/Freia Krause/01_Experiments/FK46_iALmice_high Fat diet 52 weeks 7d after injection/Analysis/02_GeneratedData/Legendplex"),
    FK49_output =paste0(parent,"/OneDrive - Universität Salzburg/AG_Tumorimmunologie - Dokumente/Data/Freia Krause/01_Experiments/FK49_CD-HFD_13wks/FK49_Analysis/02_GeneratedData/Legendplex")
  ),
  MASH = list(
    FK49_output  =paste0(parent,"/OneDrive - Universität Salzburg/AG_Tumorimmunologie - Dokumente/Data/Freia Krause/01_Experiments/FK49_CD-HFD_13wks/FK49_Analysis/02_GeneratedData/NASH_Score"),
    FK49_input = paste0(parent,"/OneDrive - Universität Salzburg/AG_Tumorimmunologie - Dokumente/Data/Freia Krause/01_Experiments/FK49_CD-HFD_13wks/FK49_Analysis/01_RawData"),
    FK46_output = paste0(parent,"/OneDrive - Universität Salzburg/AG_Tumorimmunologie - Dokumente/Data/Freia Krause/01_Experiments/FK46_iALmice_high Fat diet 52 weeks 7d after injection/Analysis/02_GeneratedData/NASH_Score"),
    FK46_input = paste0(parent,"/OneDrive - Universität Salzburg/AG_Tumorimmunologie - Dokumente/Data/Freia Krause/01_Experiments/FK46_iALmice_high Fat diet 52 weeks 7d after injection/Analysis/01_RawData")
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
   FK49_output = paste0(parent,"/OneDrive - Universität Salzburg/AG_Tumorimmunologie - Dokumente/Data/Freia Krause/01_Experiments/FK49_CD-HFD_13wks/FK49_Analysis/02_GeneratedData/Tumor_Ascites_Granuloma"),
   FK49_input = paste0(parent,"/OneDrive - Universität Salzburg/AG_Tumorimmunologie - Dokumente/Data/Freia Krause/01_Experiments/FK49_CD-HFD_13wks/FK49_Analysis/01_RawData"),
   FK46_output = paste0(parent,"/OneDrive - Universität Salzburg/AG_Tumorimmunologie - Dokumente/Data/Freia Krause/01_Experiments/FK46_iALmice_high Fat diet 52 weeks 7d after injection/Analysis/02_GeneratedData/Tumor_Ascites_Granuloma"),
   FK46_input = paste0(parent,"/OneDrive - Universität Salzburg/AG_Tumorimmunologie - Dokumente/Data/Freia Krause/01_Experiments/FK46_iALmice_high Fat diet 52 weeks 7d after injection/Analysis/01_RawData")
 ),
 Organs = list(
   FK49_output = paste0(parent,"/OneDrive - Universität Salzburg/AG_Tumorimmunologie - Dokumente/Data/Freia Krause/01_Experiments/FK49_CD-HFD_13wks/FK49_Analysis/02_GeneratedData/Weight_Organs"),
   FK49_input = paste0(parent,"/OneDrive - Universität Salzburg/AG_Tumorimmunologie - Dokumente/Data/Freia Krause/01_Experiments/FK49_CD-HFD_13wks/FK49_Analysis/01_RawData"),
   FK46_output = paste0(parent,"/OneDrive - Universität Salzburg/AG_Tumorimmunologie - Dokumente/Data/Freia Krause/01_Experiments/FK46_iALmice_high Fat diet 52 weeks 7d after injection/Analysis/02_GeneratedData/Weight_Organs"),
   FK46_input = paste0(parent,"/OneDrive - Universität Salzburg/AG_Tumorimmunologie - Dokumente/Data/Freia Krause/01_Experiments/FK46_iALmice_high Fat diet 52 weeks 7d after injection/Analysis/01_RawData")
 ),
 BH_baseline = list(
   input = paste0(parent,"/OneDrive - Universität Salzburg/AG_Tumorimmunologie - Dokumente/Data/Freia Krause/01_Experiments/FK49_CD-HFD_13wks/FK49_Analysis/01_RawData")
 )

)

# ============================================================
# process_metabolome: BA-style censored tracking via pivot
# ============================================================
# Reads a data frame with metadata columns + metabolite columns.
# Metabolite columns may contain "<LOD" strings (censored values).
# Returns a wide data frame with:
#   - metadata columns (meta_cols)
#   - per metabolite: <Metab> (numeric, LOD-imputed),
#                     <Metab>_censored (logical),
#                     <Metab>_direction (character: "<" or NA),
#                     <Metab>_raw (character: original value)
#
# LOD imputation: <LOD replaced with min_detected / 2 per metabolite.
# EtOH -> Ctrl renaming (like process_BA).

process_metabolome <- function(df, meta_cols) {
  # Create derived metadata columns and convert European decimal commas
  df <- df %>%
    mutate(
      Treatment = ifelse(Treatment == "EtOH", "Ctrl", Treatment),
      Sample = as.factor(Sample),
      T_D_S = paste0(Treatment, "_", Diet, "_", Sex),
      T_D   = paste0(Treatment, "_", Diet),
      T_S   = paste0(Treatment, "_", Sex)
    ) %>%
    mutate(across(where(is.character), ~ gsub(",", ".", .)))

  # Identify metabolite columns (everything not in meta_cols)
  metab_cols <- setdiff(colnames(df), meta_cols)

  df %>%
    # Convert all metabolite columns to character so pivot_longer
    # can combine numeric and <LOD columns into one value_raw column
    mutate(across(all_of(metab_cols), as.character)) %>%
    pivot_longer(all_of(metab_cols), names_to = "metabolite",
                 values_to = "value_raw") %>%
    mutate(
      value_raw     = as.character(value_raw),
      censored      = (value_raw == "<LOD"),
      direction     = ifelse(censored, "<", NA_character_),
      numeric_value = as.numeric(gsub("\\.(?=.*\\.)", "", value_raw, perl = TRUE))
    ) %>%
    group_by(metabolite) %>%
    mutate(
      numeric_value = if_else(
        value_raw == "<LOD",
        min(numeric_value[!censored & !is.na(numeric_value)], na.rm = TRUE) / 2,
        numeric_value
      )
    ) %>%
    ungroup() %>%
    pivot_wider(
      names_from  = metabolite,
      values_from = c(value_raw, censored, direction, numeric_value),
      names_glue  = "{metabolite}_{.value}"
    ) %>%
    mutate(
      Treatment = factor(Treatment, levels = c("Ctrl", "TAM")),
      Sex       = factor(Sex,       levels = c("female", "male")),
      Diet      = factor(Diet,      levels = c("ND", "CDHFD13"))
    ) %>%
    rename_with(~ str_replace(.x, "_numeric_value$", ""),
                ends_with("_numeric_value")) %>%
    rename_with(~ str_replace(.x, "_value_raw$", "_raw"),
                ends_with("_value_raw"))
}

# ============================================================
# subset_data: subsetting + normalisation (used by other scripts)
# ============================================================

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
    BA_secondary_con =c("GDCA","TDCA","GLCA","TLCA","TUDCA","THDCA","TMCA"),
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
PARAMETERS$metabolomics <- list(
  metabolite_abbrev = c(
    # --- targeted (65 metabolites) ---
    "alpha-Ketoglutaric acid" = "2-Oxoglutaric acid",
    "AMP" = "AMP",
    "ATP" = "ATP",
    "CDP" = "CDP",
    "Citric acid" = "Citric acid",
    "CMP" = "CMP",
    "D-Sedoheptulose 7-phosphate" = "Sedoheptulose 7-phosphate",
    "Erythrose 4-phosphate" = "Erythrose 4-phosphate",
    "Fumarate" = "Fumaric acid",
    "GDP" = "GDP",
    "Gluconolactone" = "Gluconolactone",
    "Glucose" = "Glucose",
    "Glucose 1-phosphate" = "Glucose 1-phosphate",
    "Glucose 6-phosphate" = "Glucose 6-phosphate",
    "GMP" = "GMP",
    "GTP" = "GTP",
    "Guanosine" = "Guanosine",
    "Lactic acid" = "Lactic acid",
    "Malic acid" = "Malic acid",
    "Oxaloacetic acid" = "Oxaloacetic acid",
    "Ribose" = "Ribose",
    "Ribose 5-phosphate" = "Ribose 5-phosphate",
    "Ribulose 5-phosphate" = "Ribulose 5-phosphate",
    "Succinic acid" = "Succinic acid",
    "UMP" = "UMP",
    "Uridine" = "Uridine",
    "4-Hydroxyproline" = "cis-4-Hydroxyproline",
    "Aminoadipic Acid" = "Aminoadipic acid",
    "2-Aminoisobutyric Acid" = "Aminoisobutyric acid",
    "Amino-n-butyric Acid" = "Amino-n-butyric Acid",
    "Butyrylcarnitine" = "CAR 4:0",
    "Carnitine" = "Carnitine",
    "Citrulline" = "Citrulline",
    "Creatinine" = "Creatinine",
    "Ethanolamine" = "Ethanolamine",
    "Glycine" = "Glycine",
    "Isovalerylcarnitine" = "CAR 4:0;3Me",
    "Kynurenine" = "Kynurenine",
    "L-Acetylcarnitine" = "CAR 2:0",
    "L-Alanine" = "Alanine",
    "L-Arginine" = "Arginine",
    "L-Asparagine" = "Asparagine",
    "L-Aspartic Acid" = "Aspartic acid",
    "L-Cystathionine" = "Cystathionine",
    "L-Cystine" = "Cystine",
    "L-Glutamic Acid" = "Glutamic acid",
    "L-Glutamine" = "Glutamine",
    "L-Histidine" = "Histidine",
    "L-Isoleucine" = "Isoleucine",
    "L-Leucine" = "Leucine",
    "L-Lysine" = "Lysine",
    "L-Methionine" = "Methionine",
    "L-Phenylalanine" = "Phenylalanine",
    "L-Proline" = "Proline",
    "L-Sarcosine" = "Sarcosine",
    "L-Serine" = "Serine",
    "L-Threonine" = "Threonine",
    "L-Tryptophan" = "Tryptophan",
    "L-Tyrosine" = "Tyrosine",
    "L-Valine" = "Valine",
    "Methylhistidine" = "Methylhistidine",
    "Ornithine" = "Ornithine",
    "Propionylcarnitine" = "CAR 3:0",
    "Taurine" = "Taurine",
    "Urea" = "Urea",

    # --- negative (93 metabolites) ---
    "(2R,5S)-5-Amino-1,2,6,7-tetrahydroxyoctane-3,4-dione" = "(2R,5S)-5-Amino-1,2,6,7-tetrah",
    "(2S)-4,5,7-Trihydroxy-2-(3,4,5-trihydroxyphenyl)-4H-chromen-3-one" = "(2S)-4,5,7-Trihydroxy-2-(3,4,5",
    "(2S,3R,4R,6R,7S,8R,9R)-1,2,3,4,6,7,8,9-Octahydroxydecan-5-one" = "(2S,3R,4R,6R,7S,8R,9R)-1,2,3,4",
    "(3R,5R,6R)-1,2,3,5,6,7-Hexahydroxyheptan-4-one" = "(3R,5R,6R)-1,2,3,5,6,7-Hexahyd",
    "(3S,4R,5S)-3,4-Dihydroxy-5-(2-hydroxyacetyl)oxolan-2-one" = "(3S,4R,5S)-3,4-Dihydroxy-5-(2-",
    "[(2R,3R,4R,5S)-3,4-Diacetyloxy-5-(2,4-dioxopyrimidin-1-yl)oxolan-2-yl]methyl acetate" = "[(2R,3R,4R,5S)-3,4-Diacetyloxy",
    "1-(3-Chloro-4-hydroxyphenyl)-2-((1,1-dimethylethyl)amino)-1-propanone" = "1-(3-Chloro-4-hydroxyphenyl)-2",
    "1,5,6-Trihydroxyhexane-2,3-dione" = "1,5,6-Trihydroxyhexane-2,3-dio",
    "1,5-Dinitrocyclohexa-1,4-diene" = "1,5-Dinitrocyclohexa-1,4-diene",
    "1-[(2R,4R,5R)-3,4-Dihydroxy-5-(hydroxymethyl)-2-oxolanyl]-2-pyrimidinone" = "1-[(2R,4R,5R)-3,4-Dihydroxy-5-",
    "1-Hydroxypyrrole-2,5-diol" = "1-Hydroxypyrrole-2,5-diol",
    "1-Methylpyrrole-2,3,5-triol" = "1-Methylpyrrole-2,3,5-triol",
    "1-Naphthalenesulfonic acid, 3,4-dihydro-3,4-dioxo-" = "1-Naphthalenesulfonic acid, 3,",
    "2,6-Dinitrophenol" = "2,6-Dinitrophenol",
    "2'-C-Methylcytidine" = "2'-C-Methylcytidine",
    "2'-Deoxyadenosine 3'-monophosphate" = "2'-Deoxy-3'-AMP",
    "2-Phenylethanol glucuronide" = "2-Phenylethanol glucuronide",
    "3-Hydroxyisovaleric acid" = "3-Hydroxyisovaleric acid",
    "4-Hydroxybutyric acid" = "4-Hydroxybutyric acid",
    "4-Hydroxynicotinamide" = "4-Hydroxynicotinamide",
    "4-Hydroxyquinoline" = "4-Hydroxyquinoline",
    "4-Isopropoxyphenol" = "4-Isopropoxyphenol",
    "4-Mercaptobutyric acid" = "4-Mercaptobutyric acid",
    "4-Mercapto-ethyl-pyridine" = "4-Mercapto-ethyl-pyridine",
    "4-Methylbenzenecarbothioamide" = "4-Methylbenzenecarbothioamide",
    "5-Hydroxyindoleacetic acid" = "5-Hydroxyindoleacetic acid",
    "Butanoic acid, [(diethoxyphosphinyl)oxy]methyl ester" = "Butanoic acid, [(diethoxyphosp",
    "Cinnamic acid" = "trans-Cinnamic acid",
    "Citraconic acid" = "Citraconic acid",
    "D-2-Hydroxyglutaric acid" = "D-2-Hydroxyglutaric acid",
    "D-alpha-Aminobutyric acid" = "2R-Aminobutyric acid",
    "D-Arabitol" = "Arabitol",
    "D-Fructose" = "Fructose",
    "DL-Dopa" = "DOPA",
    "D-Ribulose 5-phosphate" = "Ribulose 5-phosphate",
    "Erythronic acid" = "Erythronic acid",
    "Ethyl glucuronide" = "Ethyl glucuronide",
    "Galactose 1-phosphate" = "Galactose 1-phosphate",
    "Gluconic acid" = "Gluconic acid",
    "Glutarylglycine" = "Glutarylglycine",
    "Hexanoylglycine" = "NA-Gly 6:0",
    "Hippuric acid" = "Hippuric acid",
    "Homocitrulline" = "Homocitrulline",
    "Hydroxyisocaproic acid" = "Hydroxyisocaproic acid",
    "Hydroxyphenyllactic acid" = "Hydroxyphenyllactic acid",
    "Hydroxypropionic acid" = "Hydroxypropionic acid",
    "Indolelactic acid" = "Indolelactic acid",
    "Isodesmosine" = "Isodesmosine",
    "Isovalerylglycine" = "Isovalerylglycine",
    "Ketoleucine" = "Ketoleucine",
    "L-Arabinose" = "L-Arabinose",
    "L-Homoserine" = "Homoserine",
    "L-Rhamnulose" = "L-Rhamnulose",
    "Malonic acid" = "Malonic acid",
    "Mannitol" = "Mannitol",
    "Methylglutaric acid" = "Methylglutaric acid",
    "N(1)-Ethylchlorpropamide" = "N(1)-Ethylchlorpropamide",
    "N-(2-Amino-2-oxoethyl)-2-methylprop-2-enamide" = "N-(2-Amino-2-oxoethyl)-2-methy",
    "N-(3-Pyrene)maleimide" = "N-(3-Pyrene)maleimide",
    "N(4)-Acetylsulfisoxazole" = "N(4)-Acetylsulfisoxazole",
    "N-(4-Aminobenzoyl)-L-glutamic acid" = "4-Aminobenzoylglutamic acid",
    "N(epsilon)-(carboxyethyl)lysine" = "N(6)-(2-Carboxyethyl)-lysine",
    "N-4-Carboxybenzylglucamine dithiocarbamate" = "N-4-Carboxybenzylglucamine dit",
    "N-Acetylcyanamide" = "N-Acetylcyanamide",
    "N-Acetyl-DL-phenylalanine" = "N-Acetyl-DL-phenylalanine",
    "N-Acetyl-L-aspartic acid" = "N-Acetylaspartic acid",
    "N-Allylglycine" = "N-Allylglycine",
    "N-Butyrylglycine" = "NA-Gly 4:0",
    "n-carboxymethyllysine" = "n-carboxymethyllysine",
    "n-chloroaniline" = "n-chloroaniline",
    "Neuraminic acid" = "Neuraminic acid",
    "n-formimidoyl-glutamic acid" = "Formiminoglutamic acid",
    "N-Formylglycine" = "N-Formylglycine",
    "N-Hbgd" = "N-Hbgd",
    "N-Hydroxybenzamide" = "N-Hydroxybenzamide",
    "N-hydroxy-L-tryptophan" = "N-hydroxy-L-tryptophan",
    "N-Hydroxymethyl-N-methylformamide" = "N-Hydroxymethyl-N-methylformam",
    "N-Methyl-4-nitroaniline" = "N-Methyl-4-nitroaniline",
    "N-Methyl-DL-aspartic acid" = "N-Methyl-DL-aspartic acid",
    "N-Methylsuccinimide" = "N-Methylsuccinimide",
    "N-NITROSO-N-METHYLURETHANE" = "N-NITROSO-N-METHYLURETHANE",
    "N-Oxalylglycine" = "N-Oxalylglycine",
    "N-Oxalyl-L-alanine" = "N-Oxalyl-L-alanine",
    "N-Propyl-N-nitrosourea" = "N-Propyl-N-nitrosourea",
    "ortho-Hydroxyphenylacetic acid" = "2-Hydroxyphenylacetic acid",
    "Oxypurinol" = "Oxypurinol",
    "Phenylacetylglycine" = "Phenylacetylglycine",
    "Phenyllactic acid" = "Phenyllactic acid",
    "Propionylglycine" = "N-Propionylglycine",
    "Pseudouridine" = "Pseudouridine",
    "Suberic acid" = "Suberic acid",
    "Tyramine glucuronide" = "Coenzyme Q10",
    "Valeric acid" = "Valeric acid",

    # --- positive (484 metabolites) ---
    "(2S)-5-[Carbamimidoyl(methyl)amino]-2-(methylamino)pentanoic acid" = "(2S)-5-[Carbamimidoyl(methyl)a",
    "(2R,3R,4S,5R)-2-Amino-3,4,5,6-tetrahydroxyhexanal" = "(2R,3R,4S,5R)-2-Amino-3,4,5,6-",
    "Methyl L-argininate" = "Methyl L-argininate",
    "Ethyl lysine" = "Ethyl lysine",
    "isoleucine glutamate" = "isoleucine glutamate",
    "N-(3-Aminopropyl)-N-methylcarbamic acid tert-butyl ester" = "N-(3-Aminopropyl)-N-methylcarb",
    "Benzoyl glucuronide (Benzoic acid)" = "Benzoyl glucuronide (Benzoic a",
    "epsilon-(gamma-Glutamyl)lysine" = "Epsilon-(gamma-Glutamyl)-lysine",
    "2-Hydrazinopyridine" = "2-Hydrazinopyridine",
    "D-Histidine" = "D-Histidine",
    "Menadione bisulfite" = "Menadione bisulfite",
    "1-Methylpiperazine" = "1-Methylpiperazine",
    "Meldonium" = "Meldonium",
    "1-Piperideine" = "1-Piperideine",
    "(1r,2s)-2-Aminocyclopentanecarboxylic acid" = "(1r,2s)-2-Aminocyclopentanecar",
    "Triethyl orthoformate" = "Triethyl orthoformate",
    "Visnagin" = "Visnagin",
    "piperidine-1-carboxamide" = "piperidine-1-carboxamide",
    "Bis((hexahydro-4-methyl-1H-1,4-diazepin-1-yl)thiocarbonyl)disulfide" = "Bis((hexahydro-4-methyl-1H-1,4",
    "Z-PP-CHO" = "Z-PP-CHO",
    "N(6)-Methyllysine" = "N(6)-Methyllysine",
    "6-[Ethyl-(3-isobutoxy-4-isopropylphenyl)amino]nicotinic acid" = "6-[Ethyl-(3-isobutoxy-4-isopro",
    "6-{3-[(Pyridin-2-yl)disulfanyl]propanamido}hexanoic acid" = "6-{3-[(Pyridin-2-yl)disulfanyl",
    "Methoxsalen" = "Methoxsalen",
    "4-(Dipropylsulfamoyl)-2-nitrobenzoic acid" = "4-(Dipropylsulfamoyl)-2-nitrob",
    "D-Ornithine" = "D-Ornithine",
    "N-Mononitrosopiperazine" = "N-Mononitrosopiperazine",
    "Trimethylolpropane" = "Trimethylolpropane",
    "Pyrrole, 2-pentyl" = "Pyrrole, 2-pentyl",
    "10-Acetyl-3,7-dihydroxyphenoxazine" = "10-Acetyl-3,7-dihydroxyphenoxa",
    "1-[(5-Amino-5-carboxypentyl)amino]-1-deoxyfructose" = "Fructosyllysine",
    "Pentisomide" = "Pentisomide",
    "THURFYL NICOTINATE" = "THURFYL NICOTINATE",
    "Dibutylamine" = "N-Butylbutan-1-amine",
    "Lumichrome" = "Lumichrome",
    "Recainam" = "Recainam",
    "3h-Adrenaline" = "3h-Adrenaline",
    "Spisulosine" = "Spisulosine",
    "N-Lactoylphenylalanine" = "N-Lactoyl phenylalanine",
    "15d PGD2" = "15-Deoxy-delta-12,14-PGD2",
    "N-Deschlorobenzoyl indomethacin" = "N-Deschlorobenzoyl indomethaci",
    "2-Methylbutyroylcarnitine" = "CAR 4:0;2Me",
    "Nicotinic acid" = "Nicotinic acid",
    "N-Nitramido-N-phenylnitramide" = "N-Nitramido-N-phenylnitramide",
    "Amidephrine" = "Amidephrine",
    "Cytosine-5-carboxylic acid" = "Cytosine-5-carboxylic acid",
    "Picolinamide" = "Picolinamide",
    "Pantothenic acid" = "Pantothenic acid",
    "Oleamide" = "Oleamide",
    "Spiroxamine" = "Spiroxamine",
    "L-Threo-Sphingosine C-18" = "L-threo-Sphingosine",
    "Fluoroacetamide" = "Fluoroacetamide",
    "(2R,5S)-2-(6-Aminopurin-9-yl)-5-(methylsulfanylmethyl)oxolane-3,4-diol" = "(2R,5S)-2-(6-Aminopurin-9-yl)-",
    "5-Methyl-1H-pyrazole-3-carboxylic acid" = "5-Methyl-1H-pyrazole-3-carboxy",
    "2-Nitrophenol" = "2-Nitrophenol",
    "N1-Methyl-2-pyridone-5-carboxamide" = "N1-Methyl-2-pyridone-5-carboxamide",
    "2',3'-Didehydro-2',3'-dideoxyguanosine" = "2',3'-Didehydro-2',3'-dideoxyg",
    "Swertiamarin" = "Swertiamarin",
    "6-Hydroxydopamine" = "6-Hydroxydopamine",
    "2-(1-Phenylpropan-2-yl)-1,3,2-dioxazetidine" = "2-(1-Phenylpropan-2-yl)-1,3,2-",
    "Procainamide" = "Procainamide",
    "[(3-Methylbutyl)amino]acetic acid" = "[(3-Methylbutyl)amino]acetic a",
    "1-Octanethiol" = "1-Octanethiol",
    "2-Hexenoylcarnitine" = "CAR 6:1",
    "N-(4-Aminobenzoyl)-L-glutamic acid" = "4-Aminobenzoylglutamic acid",
    "Hexaminolevulinate" = "Hexaminolevulinate",
    "(2S)-6-Amino-2-(hexanoylamino)hexanoic Acid" = "(2S)-6-Amino-2-(hexanoylamino)",
    "Phenylalanylphenylalanine" = "Phe-Phe",
    "Butanoic acid, 3-((2-((3R)-2-oxo-3-(2-(4-piperidinyl)ethyl)-1-piperidinyl)acetyl)amino)-, (3R)-" = "Butanoic acid, 3-((2-((3R)-2-o",
    "Flunixin" = "Flunixin",
    "Amobarbital" = "Amobarbital",
    "Indicine-N-oxide" = "Indicine-N-oxide",
    "[1S-[1alpha,2alpha(Z),3alpha,4alpha]]-7-[3-[[2-[(Phenylamino)carbonyl]hydrazino]methyl]-7-oxabicyclo[2.2.1]hept-2-yl]-5-heptenoicacid" = "[1S-[1alpha,2alpha(Z),3alpha,4",
    "9H-Purine-9-ol" = "9H-Purine-9-ol",
    "N-Acetylserotonin" = "N-Acetylserotonin",
    "Carbaryl" = "Carbaryl",
    "3,5-Dihydroxyphenylglycine" = "3,5-Dihydroxyphenylglycine",
    "3-Nitrophenylhydrazine" = "3-Nitrophenylhydrazine",
    "N4-Acetylcytidine" = "N4-Acetylcytidine",
    "Pimelylcarnitine" = "CAR DC7:0",
    "Allopurinol riboside" = "Allopurinol riboside",
    "Procarbazine" = "Procarbazine",
    "11-Mercaptoundecanoic acid" = "11-Mercaptoundecanoic acid",
    "Carbidopa" = "Carbidopa",
    "N,N-Diethyl-1-methylsulfinylformamide" = "N,N-Diethyl-1-methylsulfinylfo",
    "Creatinine bicarbonate" = "Creatinine bicarbonate",
    "Riboflavin" = "Riboflavin",
    "N2,N2-Dimethylguanosine" = "N2,N2-Dimethylguanosine",
    "7-(2-Hydroxyethyl)guanine" = "7-(2-Hydroxyethyl)guanine",
    "(3-Aminophenyl) nitrite" = "(3-Aminophenyl) nitrite",
    "3-Methylglutarylcarnitine" = "CAR DC5:0;3Me",
    "2-(3,4-Dihydro-2,2-dimethyl-6-nitro-2H-1,4-benzoxazin-4-yl)pyridine N-oxide" = "2-(3,4-Dihydro-2,2-dimethyl-6-",
    "Etifoxine" = "Etifoxine",
    "3-Methylthymidine" = "3-Methylthymidine",
    "Lersivirine" = "Lersivirine",
    "Tribufos" = "Tribufos",
    "Ecgonine methyl ester" = "Ecgonine methyl ester",
    "Formycin" = "Formycin",
    "5-(Cycloocten-1-yl)-5-ethyl-barbituric acid" = "5-(Cycloocten-1-yl)-5-ethyl-ba",
    "1-[(2R,3R,4S,5R)-3,4-Dihydroxy-5-(hydroxymethyl)oxolan-2-yl]oxypyrimidin-2-one" = "1-[(2R,3R,4S,5R)-3,4-Dihydroxy",
    "Diisopropanolamine" = "1-(2-hydroxypropylamino)propan-2-ol",
    "Adenine" = "Adenine",
    "PIPERIDIN-2-IMINE" = "PIPERIDIN-2-IMINE",
    "triacsin" = "triacsin",
    "N-Lactoylleucine" = "N-Lactoyl leucine",
    "2-[(2R)-2-Aminopropyl]-5-hydroxybenzoic acid" = "2-[(2R)-2-Aminopropyl]-5-hydro",
    "4-Methoxycinnamic acid" = "4-Methoxycinnamic acid",
    "Zolmitriptan" = "Zolmitriptan",
    "Val-Met" = "Val-Met",
    "Acefylline" = "Acefylline",
    "Metioprim" = "Metioprim",
    "Xanthosine" = "Xanthosine",
    "Xanthine" = "Xanthine",
    "2-Hydroxyisovalerylcarnitine" = "2-Hydroxyisovalerylcarnitine",
    "N-Ethyl-N-(2-hydroxy-3-sulfopropyl)-3-toluidine" = "N-Ethyl-N-(2-hydroxy-3-sulfopr",
    "Butylate" = "Butylate",
    "Nelarabine" = "Nelarabine",
    "Lysophosphatidylcholine" = "Lysophosphatidylcholine",
    "Nicotinamide ribotide" = "Nicotinamide ribotide",
    "Cycloate" = "Cycloate",
    "Tert-butyl N-[2-(prop-2-enamido)ethyl]carbamate" = "Tert-butyl N-[2-(prop-2-enamid",
    "Idrapril" = "Idrapril",
    "4-(1,3-Benzothiazol-2-yl)-2-methylaniline" = "4-(1,3-Benzothiazol-2-yl)-2-me",
    "1-Dodecanethiol" = "1-Dodecanethiol",
    "Thiochrome" = "Thiochrome",
    "(2-Phenylacetyl) (2R)-2,5-diaminopentanoate" = "(2-Phenylacetyl) (2R)-2,5-diam",
    "N,N,N-Trimethyl-L-alanyl-L-proline betaine" = "N,N,N-Trimethyl-alanylproline betaine",
    "2-Methylthiazolidine-4-carboxylic acid" = "2-Methylthiazolidine-4-carboxylic acid",
    "3-Hydroxy-11-norcytisine" = "3-Hydroxy-11-norcytisine",
    "3-Hydroxybutyrylcarnitine" = "CAR 4:0;3OH",
    "XANTHOPTERIN" = "Xanthopterin",
    "Carbimazole" = "Carbimazole",
    "Alaptide" = "Alaptide",
    "2,2'-Anhydrocytidine" = "Ancitabine",
    "2-Amino-3-methyl-1-butanol" = "Valinol",
    "Biopterin" = "Biopterin",
    "Hexobarbital" = "Hexobarbital",
    "N6-Succinyl Adenosine" = "Succinyladenosine",
    "2-Methyl-2-phenylsuccinimide" = "2-Methyl-2-phenylsuccinimide",
    "Benzaldehyde" = "Benzaldehyde",
    "Isoindoline" = "Isoindoline",
    "Ethyl 3-aminobenzoate" = "Tricaine",
    "Phenylethylmalonamide" = "Phenylethylmalonamide",
    "Bunitrolol" = "Bunitrolol",
    "L-Norleucine" = "Norleucine",
    "Phenethylamine glucuronide" = "Phenethylamine glucuronide",
    "Suberylglycine" = "Suberylglycine",
    "3-(4-Hydroxy-3-methoxyphenyl)-2-methylpropionic acid" = "3-(4-Hydroxy-3-methoxyphenyl)-",
    "Phenidone" = "Phenidone",
    "Azepine" = "Azepine",
    "Allobarbital" = "Allobarbital",
    "Noroxyhydrastinine" = "Noroxyhydrastinine",
    "Cyclopentylamine" = "Cyclopentylamine",
    "4-Hydroxy-benzamidine" = "4-Hydroxy-benzamidine",
    "2-Guanidinobutanoic acid" = "2-Guanidinobutanoic acid",
    "2-Aminobicyclo[3.1.0]hexane-2,6-dicarboxylic acid" = "2-Aminobicyclo[3.1.0]hexane-2,",
    "3-Pyridylacetic acid" = "3-Pyridylacetic acid",
    "pyramid" = "pyramid",
    "Uric acid" = "Uric acid",
    "Allantoin" = "Allantoin",
    "Mepirizole" = "Mepirizole",
    "thiazoline" = "thiazoline",
    "Thiazolidine-4-carboxylic acid" = "Thioproline",
    "N-(2-Aminoethyl)maleimide" = "N-(2-Aminoethyl)maleimide",
    "beta-Hydroxy-gamma-trimethylaminobutyric acid" = "beta-Hydroxy-gamma-trimethylam",
    "2-(5-Fluoropentyl)-2-methylmalonic acid" = "2-(5-Fluoropentyl)-2-methylmal",
    "4-Hydroxyantipyrine" = "4-Hydroxyantipyrine",
    "4,4'-Bipyridine" = "4,4'-Bipyridine",
    "Indoleacrylic acid" = "Indoleacrylic acid",
    "1,5-Naphthalenediamine" = "1,5-Naphthalenediamine",
    "1-Hydroxyisoquinoline" = "1(2H)-Isoquinolinone",
    "21,22,23,24-Tetrahydroporphyrin" = "21,22,23,24-Tetrahydroporphyri",
    "Glycine, L-g-glutamyl-L-cysteinyl-, 3-methyl ester" = "Glycine, L-g-glutamyl-L-cystei",
    "Pent-2-enoic acid" = "Pent-2-enoic acid",
    "D-Erythro-7,8-dihydrobiopterin" = "D-Erythro-7,8-dihydrobiopterin",
    "1,1'-(1,8-Dioxo-1,8-octanediyl)bis-2,5-pyrrolidinedione" = "1,1'-(1,8-Dioxo-1,8-octanediyl",
    "2-Butylhydroquinone" = "2-Butylhydroquinone",
    "tetrahydrothiazine" = "tetrahydrothiazine",
    "L-Methionine" = "Methionine",
    "9-Deazaguanine" = "9-Deazaguanine",
    "1-Pentanesulfonic acid" = "1-Pentanesulfonic acid",
    "1-Hydroxybenzotriazole" = "1-Hydroxybenzotriazole",
    "Hydrazinonicotinamide" = "Hydrazinonicotinamide",
    "Tert-Butyl carbamate" = "Tert-Butyl carbamate",
    "Boc-L-cysteine" = "Boc-L-cysteine",
    "1-Deoxy-1-morpholino-D-fructose" = "1-Deoxy-1-morpholinofructose",
    "Aminopropanol" = "Aminopropanol",
    "4-Acetamidobutanoic acid" = "4-Acetamidobutanoic acid",
    "1-Ethylaziridine" = "1-Ethylaziridine",
    "Cytarabine" = "Cytarabine",
    "Bisaramil" = "Bisaramil",
    "bicyclol" = "bicyclol",
    "Succinaldehyde" = "Succinaldehyde",
    "2-Aminoisobutyric acid" = "Aminoisobutyric acid",
    "2-Pyrrolidinone" = "2-Pyrrolidinone",
    "1-Methyladenosine" = "1-Methyladenosine",
    "(1R)-1-Amino-2-sulfanylethanesulfinic acid" = "(1R)-1-Amino-2-sulfanylethanes",
    "N-(2-Hydroxyethyl)acrylamide" = "N-(2-Hydroxyethyl)acrylamide",
    "Butyronitrile" = "Butyronitrile",
    "1-Methylcytosine" = "1-Methylcytosine",
    "2-(2,4-Dimethylphenyl)indan-1,3-dione" = "2-(2,4-Dimethylphenyl)indan-1,",
    "N-(2,3-Dihydroxypropyl)valine" = "N-(2,3-Dihydroxypropyl)valine",
    "Aldicarb sulfoxide" = "Aldicarb sulfoxide",
    "beta-naphthoflavone" = "beta-naphthoflavone",
    "Indirubin-3'-monoxime" = "Indirubin-3'-monoxime",
    "N(6)-(1-Carboxyethyl)-L-lysine" = "N(6)-(1-Carboxyethyl)-lysine",
    "Formyllysine" = "Formyllysine",
    "4-Hydroxybenzaldehyde" = "4-Hydroxybenzaldehyde",
    "Isocoumarin" = "Isocoumarin",
    "4-Hydroxycinnamic acid" = "cis-p-Coumaric acid",
    "4-Aminoacetophenone" = "4-Aminoacetophenone",
    "2-Amino-2-methyl-1,3-propanediol" = "2-Amino-2-methyl-1,3-propanedi",
    "2-Oxazolidinone, 3-[[(5-nitro-2-furanyl)methylene]amino]-" = "2-Oxazolidinone, 3-[[(5-nitro-",
    "Pentanamide" = "Pentanamide",
    "Amino (2S)-2,6-diaminohexanoate" = "Amino (2S)-2,6-diaminohexanoat",
    "p-Hydroxyfelbamate" = "p-Hydroxyfelbamate",
    "Ribothymidine" = "Ribothymidine",
    "gamma-Glutamylacetamide" = "gamma-Glutamylacetamide",
    "3-Hydroxyhalazepam" = "3-Hydroxyhalazepam",
    "Temazepam" = "Temazepam",
    "Perflexane" = "Perflexane",
    "{[5-(5-Nitro-2-furyl)-1,3,4-oxadiazol-2-YL]thio}acetic acid" = "{[5-(5-Nitro-2-furyl)-1,3,4-ox",
    "1,1,1,3,3-Pentafluoropropane" = "1,1,1,3,3-Pentafluoropropane",
    "2-Butanethiol" = "2-Butanethiol",
    "4-Aminophenyl phosphate" = "4-Aminophenyl phosphate",
    "Taurocyamine" = "Taurocyamine",
    "2-Hydroxy-3-(1H-imidazol-5-yl)propanoic acid" = "2-Hydroxy-3-(1H-imidazol-5-yl)",
    "(2R,4R)-4-Aminopyrrolidine-2,4-dicarboxylic acid" = "(2R,4R)-4-Aminopyrrolidine-2,4",
    "S-Nitrosomercaptoethanol" = "S-Nitrosomercaptoethanol",
    "4-Chloroaniline" = "4-Chloroaniline",
    "2-Amino-4-chloropyridine" = "2-Amino-4-chloropyridine",
    "Naphth(1,2-b)oxirene" = "Naphth(1,2-b)oxirene",
    "Bisphenol" = "Bisphenol A",
    "Thiacloprid" = "Thiacloprid",
    "Dioxybenzone" = "Dioxybenzone",
    "Serine methyl ester" = "Serine methyl ester",
    "1-Aminocyclopropanol" = "1-Aminocyclopropanol",
    "2,3-Butanedione monoxime" = "2,3-Butanedione monoxime",
    "Morph" = "Morph",
    "Ergothioneine" = "Ergothioneine",
    "3-Aminobenzanthrone" = "3-Aminobenzanthrone",
    "Pirprofen" = "Pirprofen",
    "6-Nitrochrysene" = "6-Nitrochrysene",
    "Guanidoacetic acid" = "Guanidoacetic acid",
    "N-Methylglucosamine" = "N-Methylglucosamine",
    "Pegorgotein" = "Pegorgotein",
    "Thioacetamide-S-oxide" = "Thioacetamide-S-oxide",
    "Hypotaurine" = "Hypotaurine",
    "2-Amino-1,3-propanediol" = "2-Amino-1,3-propanediol",
    "Ala-Ala-Ala" = "Ala-Ala-Ala",
    "Ethylbenzene" = "Ethylbenzene",
    "Isochroman" = "Isochroman",
    "Cumene hydroperoxide" = "Cumene hydroperoxide",
    "Citraconic acid" = "Citraconic acid",
    "Alanine acetate" = "Alanine acetate",
    "Oxazine" = "Oxazine",
    "o-Acetylhydroxylamine" = "o-Acetylhydroxylamine",
    "2-Oxoarginine" = "2-Oxoarginine",
    "10,13-Dimethyl-7-sulfanylspiro[2,6,7,8,9,11,12,14,15,16-decahydro-1H-cyclopenta[a]phenanthrene-17,5'-oxolane]-2',3-dione" = "10,13-Dimethyl-7-sulfanylspiro",
    "TRIPHENYLMETHANE" = "TRIPHENYLMETHANE",
    "Firocoxib" = "Firocoxib",
    "Guanidinosuccinic acid" = "Guanidinosuccinic acid",
    "4-Methylumbelliferylguanidinobenzoate" = "4-Methylumbelliferylguanidinob",
    "N'-[(4-Oxo-4H-chromen-3-yl)methylene]nicotinohydrazide" = "N'-[(4-Oxo-4H-chromen-3-yl)met",
    "N-(2-Hydroxyethyl)ethylenediaminetriacetic acid" = "N-(2-Hydroxyethyl)ethylenediam",
    "Glycylsarcosine" = "Glycylsarcosine",
    "Aldehydo-ascarylose" = "Aldehydo-ascarylose",
    "Glutathione" = "Glutathione",
    "Kresoxim-Methyl" = "Kresoxim-Methyl",
    "2-Methyl-1,1,4-trioxo-N-(2-pyridinyl)-3H-1$l^{6},2-benzothiazine-3-carboxamide" = "2-Methyl-1,1,4-trioxo-N-(2-pyr",
    "Fleroxacin N-oxide" = "Fleroxacin N-oxide",
    "S-Adenosylhomocysteine" = "S-Adenosylhomocysteine",
    "5C-aglycone" = "5C-Aglycone",
    "Methyl(acetoxymethyl)nitrosamine" = "Methyl(acetoxymethyl)nitrosami",
    "methyl 4-mercaptobutyrimidate" = "methyl 4-mercaptobutyrimidate",
    "Tryptophan 2-C-mannoside" = "alpha-C-Mannosyltryptophan",
    "N(epsilon)-(Carboxymethyl)hydroxylysine" = "N(epsilon)-(Carboxymethyl)hydr",
    "N-Methyl-L-histidine" = "N-Methylhistidine",
    "Mitomycin" = "Mitomycin",
    "Bis-triazole" = "Bis-triazole",
    "2-Hydroxy-4-methylthiobutyric acid" = "2-Hydroxy-4-(methylthio)butanoic acid",
    "N-Hydroxysuccinimide" = "N-Hydroxysuccinimide",
    "D-Aspartic acid" = "D-Aspartic acid",
    "trifluoroacetyl-l-lysyl-l-alaninanilide" = "trifluoroacetyl-l-lysyl-l-alan",
    "TETRAHYDROURIDINE" = "TETRAHYDROURIDINE",
    "Homocitrulline" = "Homocitrulline",
    "rac S 33138" = "rac S 33138",
    "alpha-AMINO-3-HYDROXY-5-METHYL-4-ISOXAZOLEPROPIONIC ACID" = "Ampa",
    "(+/-)-(E)-Methyl-2-[(E)-hydroxyimino]-5-nitro-6-methoxy-3-hexeneamide" = "(+/-)-(E)-Methyl-2-[(E)-hydrox",
    "4-[(2R,5R)-5-(6-Aminopurin-9-yl)-3,4-dihydroxyoxolan-2-yl]butanamide" = "4-[(2R,5R)-5-(6-Aminopurin-9-y",
    "Carbamoyl (2R)-2,5-diaminopentanoate" = "Carbamoyl (2R)-2,5-diaminopent",
    "Oxiracetam" = "Oxiracetam",
    "(2R,4R,5S,6S)-2,4-Dihydroxy-5-[(2-hydroxyacetyl)amino]-6-[(1S,2S)-1,2,3-trihydroxypropyl]oxane-2-carboxylic acid" = "(2R,4R,5S,6S)-2,4-Dihydroxy-5-",
    "(2R,3S,4R,5R,8R,9S,10R,11R)-5-Amino-1,2,3,4,8,9,10,11-octahydroxydodecane-6,7-dione" = "(2R,3S,4R,5R,8R,9S,10R,11R)-5-",
    "Rubiadin" = "Rubiadin",
    "Cystine" = "Cystine",
    "Zipeprol" = "Zipeprol",
    "Leteprinim" = "Leteprinim",
    "2,6-Diamino-5-hydroxyhexanoic acid" = "5-Hydroxylysine",
    "Neramexane" = "Neramexane",
    "Varenicline" = "Varenicline",
    "Bethanidine" = "Bethanidine",
    "Hippuryl-L-arginine" = "Hippuryl-L-arginine",
    "2-Octenoylcarnitine" = "CAR 8:1(2E)",
    "Forchlorfenuron" = "Forchlorfenuron",
    "(2R,3R)-2-Aminooctadecane-1,3-diol" = "(2R,3R)-2-Aminooctadecane-1,3-",
    "Aminosalicylic Acid" = "Aminosalicylic acid",
    "5-Isoxazoleacetic acid, 4,5-dihydro-3-(4-hydroxyphenyl)-, methyl ester" = "5-Isoxazoleacetic acid, 4,5-di",
    "Uracil" = "Uracil",
    "Isopilocarpinic acid" = "Isopilocarpinic acid",
    "Glutarylglycine" = "Glutarylglycine",
    "3-Hydroxyhexanoylcarnitine" = "CAR 6:0;3OH",
    "8-Aminoguanosine" = "8-Aminoguanosine",
    "Celgosivir" = "Celgosivir",
    "delta-Guanidinovaleric acid" = "DGVA",
    "9,9-Dimethyl-1-(sulfinylamino)decane" = "9,9-Dimethyl-1-(sulfinylamino)",
    "4-Pyridylethylmercaptan" = "4-Pyridylethylmercaptan",
    "2-Cyanohept-2-enoic acid" = "2-Cyanohept-2-enoic acid",
    "8-Hydroxy-3-methyl-3,4-dihydrotetraphene-1,7,12(2H)-trione" = "8-Hydroxy-3-methyl-3,4-dihydro",
    "Nitracrine" = "Nitracrine",
    "Amoxapine" = "Amoxapine",
    "10H-Pyrido(3,2-b)(1,4)benzothiazine, 10-(2-piperidinoethyl)-" = "10H-Pyrido(3,2-b)(1,4)benzothi",
    "7-[(2-Benzyl-3-sulfanylpropanoyl)amino]heptanoic acid" = "7-[(2-Benzyl-3-sulfanylpropano",
    "(2R,3S)-Epoxiconazole" = "(2R,3S)-Epoxiconazole",
    "Glycyl-prolyl-glutamic acid" = "Gly-Pro-Glu",
    "1-Hydroxy-2,2,5,5-tetramethylpyrrolidine-3-carboxylic acid" = "1-Hydroxy-2,2,5,5-tetramethylp",
    "(1r,3r)-1-Aminocyclopentane-1,3-dicarboxylic acid" = "(1r,3r)-1-Aminocyclopentane-1,",
    "N2-Methylguanine" = "N2-Methylguanine",
    "RESAZURIN" = "Resazurin",
    "alpha-Hydroxymetoprolol" = "alpha-Hydroxymetoprolol",
    "(1R,7S,13S,15S)-2,15-Dihydroxy-7-methyl-6-oxabicyclo[11.3.0]hexadeca-3,11-dien-5-one" = "(1R,7S,13S,15S)-2,15-Dihydroxy",
    "(1R,4R,5S,6R)-4-Amino-2-oxabicyclo[3.1.0]hexane-4,6-dicarboxylic acid" = "(1R,4R,5S,6R)-4-Amino-2-oxabic",
    "Methyl gallate" = "Gallic acid methyl ester",
    "Azetirelin" = "Azetirelin",
    "Etimizol" = "Etimizol",
    "3[N-Morpholino]propane sulfonic acid" = "3-(N-morpholino)propanesulfonic acid",
    "(3S,8As)-7-hydroxy-3-(hydroxymethyl)-2,3,6,7,8,8a-hexahydropyrrolo[1,2-a]pyrazine-1,4-dione" = "(3S,8As)-7-hydroxy-3-(hydroxym",
    "N-[(3s)-2-Oxotetrahydrofuran-3-Yl]butanamide" = "N-Butyryl-L-homoserine lactone",
    "ISOVALERYLUREA" = "ISOVALERYLUREA",
    "TRIAZIQUONE" = "TRIAZIQUONE",
    "1-(4-Methyl-1H-pyrazol-5-yl)piperidine" = "1-(4-Methyl-1H-pyrazol-5-yl)pi",
    "1-Naphthyl isocyanate" = "1-Naphthyl isocyanate",
    "1-Nitrosopiperazine" = "1-Nitrosopiperazine",
    "2-[(3-Chlorophenyl)diazenyl]acetonitrile" = "2-[(3-Chlorophenyl)diazenyl]ac",
    "Divinyl sulfone" = "Divinyl sulfone",
    "N-(2-Amino-2-oxoethyl)-2-methylprop-2-enamide" = "N-(2-Amino-2-oxoethyl)-2-methy",
    "VASICINONE" = "Vasicinone",
    "Diethyldithiophosphate" = "Diethyldithiophosphate",
    "Beta-Guanidinopropionic acid" = "3-Guanidinopropanoic acid",
    "5-Iodo-2-pyrimidinone-2'-deoxyribose" = "5-Iodo-2-pyrimidinone-2'-deoxy",
    "Cochinchinenin" = "Cochinchinenin",
    "Diazomethyl ketone" = "Diazomethyl ketone",
    "S-Methyl N,N-Diethylthiocarbamate" = "S-Methyl N,N-Diethylthiocarbam",
    "Clobenpropit" = "Clobenpropit",
    "N(4)-Acetylsulfisoxazole" = "N(4)-Acetylsulfisoxazole",
    "2,4-Bis(2-methoxyethoxy)-1,3,5-triazine" = "2,4-Bis(2-methoxyethoxy)-1,3,5",
    "2-Iminiopropanoate" = "2-Iminiopropanoate",
    "[(2R,3S,4R,5R)-3,4,5,6-Tetrahydroxy-1-oxohexan-2-yl] (2S)-2-aminopropanoate" = "[(2R,3S,4R,5R)-3,4,5,6-Tetrahy",
    "1,2,2-Trimethylpropyl dimethylphosphinate" = "1,2,2-Trimethylpropyl dimethyl",
    "Secnidazole" = "Secnidazole",
    "2-Methylthioadenosine" = "2-Methylthioadenosine",
    "2-Furancarboximidamide" = "2-Furancarboximidamide",
    "Diuron" = "Diuron",
    "4-N-[1-(Dimethylamino)ethyl]benzene-1,4-diamine" = "4-N-[1-(Dimethylamino)ethyl]be",
    "Cystemustine" = "Cystemustine",
    "(2E)-3-[3-(sulfooxy)phenyl]prop-2-enoic acid" = "3-Hydroxycinnamate sulfate",
    "7-Hydroxyamoxapine" = "7-Hydroxyamoxapine",
    "Armillarisin" = "Armillarisin",
    "1H-Imidazole-1-ethanol, alpha-(1-aziridinylmethyl)-2-nitro-" = "1H-Imidazole-1-ethanol, alpha-",
    "N8-Acetylspermidine" = "N8-Acetylspermidine",
    "2-Methoxy-1,3-benzenediol" = "2-Methoxy-1,3-benzenediol",
    "3-Butyl-1-hydroxy-1-(4-methylphenyl)sulfonylurea" = "3-Butyl-1-hydroxy-1-(4-methylp",
    "Malonylglycine" = "Malonylglycine",
    "Milacemide" = "Milacemide",
    "Pregabalin" = "Pregabalin",
    "4,6-Bis(methylthio)hexanoic acid" = "4,6-Bis(methylthio)hexanoic ac",
    "2-(2-Aminopropanylamino)-4-carbamoylbutyric acid" = "2-(2-Aminopropanylamino)-4-car",
    "Cytosine" = "Cytosine",
    "1-Pyrrolidineethanol" = "1-Pyrrolidineethanol",
    "o-Phenetidine" = "o-Phenetidine",
    "8-Hydroxy-7-methylguanine" = "8-Hydroxy-7-methylguanine",
    "Androsta-4,16-dien-3-one" = "Androsta-4,16-dien-3-one",
    "3-Hydroxynorvaline" = "3-Hydroxynorvaline",
    "Dihydroxyethyldithiocarbamate" = "Dihydroxyethyldithiocarbamate",
    "Phenazine-2,3-diamine" = "Phenazine-2,3-diamine",
    "4-(1,3-Dihydroxy-4,4,5,5-tetramethyl-imidazolidin-2-yl)benzoic acid" = "4-(1,3-Dihydroxy-4,4,5,5-tetra",
    "Adipoylglycine" = "Adipoylglycine",
    "Kelatorphan" = "Kelatorphan",
    "2,4(1H,3H)-Pyrimidinedione, 5-fluoro-1-(tetrahydro-2-furanyl)-, (R)-" = "2,4(1H,3H)-Pyrimidinedione, 5-",
    "2-Quinolinylmethanol" = "2-Quinolinylmethanol",
    "Piroximone" = "Piroximone",
    "Aspirin eugenol ester" = "Aspirin eugenol ester",
    "Thienamycin" = "Thienamycin",
    "Diminazene" = "Diminazene",
    "(2R,3R,4R,5R)-2-Amino-4,5,6-trihydroxy-3-[(2R)-1-oxopropan-2-yl]oxyhexanal" = "(2R,3R,4R,5R)-2-Amino-4,5,6-tr",
    "Cytosine deoxyribonucleoside" = "Cytosine deoxyribonucleoside",
    "Dihydrothymine" = "Dihydrothymine",
    "2-n-Propylthiazolidine-4-carboxylic acid" = "2-n-Propylthiazolidine-4-carbo",
    "Glycylproline" = "Gly-Pro",
    "Kinetin riboside" = "Kinetin riboside",
    "1-beta-D-Arabinofuranosyl-5-fluoro-(1H,3H)-pyrimidine-2,4-dione" = "1-beta-D-Arabinofuranosyl-5-fl",
    "2,6-Diaminopurine 2',3'-dideoxyriboside" = "2,6-Diaminopurine 2',3'-dideox",
    "Tryptophol" = "Tryptophanol",
    "Diphacinone" = "Diphacinone",
    "(Z)-N-Feruloyl-5-hydroxyanthranilic acid" = "Avenanthramide 2f",
    "Tetradecylamine" = "Tetradecylamine",
    "9-Octadecenamide, N-(2-hydroxyethyl)-, (9Z)-" = "Oleoyl-EA",
    "(2S,3R)-2-Amino-4-octadecene-3-ol" = "1-Deoxysphingosine",
    "Y-27632 Dihydrochloride" = "Y-27632 Dihydrochloride",
    "S-Farnesyl cysteine" = "Farnesylcysteine",
    "O-Phosphoethanolamine" = "Phosphoethanolamine",
    "N-OMEGA-PROPYL-L-ARGININE" = "N-OMEGA-PROPYL-L-ARGININE",
    "(2R)-3-Sulfanylpropane-1,2-diol" = "(2R)-3-Sulfanylpropane-1,2-dio",
    "Aminocarb" = "Aminocarb",
    "Dibenzo(b,f)thiepin-3-methanol, 5,5-dioxide" = "Dibenzo(b,f)thiepin-3-methanol",
    "Benzoxazole" = "Benzoxazole",
    "METAZACHLOR" = "Metazachlor",
    "1-Nonanethiol" = "1-Nonanethiol",
    "2,2,6,6-Tetramethylpiperidin-1-ol" = "2,2,6,6-Tetramethylpiperidin-1",
    "2,3-Dihydro-1H-pyrrole-2-carboxylic acid" = "2,3-Dihydro-1H-pyrrole-2-carbo",
    "Betonicine" = "Betonicine",
    "5-Hydroxy-L-tryptophan" = "5-Hydroxy-tryptophan",
    "2-(4-(Methylamino)phenyl)benzo[d]thiazol-6-ol" = "2-(4-(Methylamino)phenyl)benzo",
    "N-(9H-Fluoren-9-ylacetyl)-L-phenylalanine" = "N-(9H-Fluoren-9-ylacetyl)-L-ph",
    "N-Methacryloyl-L-glutamic acid" = "N-Methacryloyl-L-glutamic acid",
    "Homocysteine" = "Homocysteine",
    "Dihydrotanshinone" = "Dihydrotanshinone",
    "(2-Chlorophenyl)diphenylmethane" = "(2-Chlorophenyl)diphenylmethane",
    "N6-Carboxymethyllysine" = "N(6)-Carboxymethyl-lysine",
    "Caffeine, 8-((3-methoxypropyl)amino)-" = "Caffeine, 8-((3-methoxypropyl)",
    "DL-Cysteine" = "Cysteine",
    "Metronidazole phosphate" = "Metronidazole phosphate",
    "(2S)-6-Oxa-1-azabicyclo[3.1.0]hexane-2-carboxylic acid" = "(2S)-6-Oxa-1-azabicyclo[3.1.0]",
    "Isoniazid" = "Isoniazid",
    "tetrahydrothiophene" = "tetrahydrothiophene",
    "Chlorambucil-tertiary butyl ester" = "Chlorambucil-tertiary butyl es",
    "n-octyl-beta-D-thioglucopyranoside" = "n-octyl-beta-D-thioglucopyrano",
    "Morantel" = "Morantel",
    "Mecamylamine" = "Mecamylamine",
    "Aceprometazine" = "Aceprometazine",
    "Linoleoyl ethanolamide" = "Linoleoyl-EA",
    "Tiglylcarnitine" = "CAR 4:1;2Me",
    "Methotrimeprazine" = "Methotrimeprazine",
    "Dapivirine" = "Dapivirine",
    "Ketamine" = "Ketamine",
    "Tyr-Gly" = "Tyr-Gly",
    "N-Desmethyllevomepromazine" = "N-Desmethyllevomepromazine",
    "N-Acetyl-2-hydroxybenzamide" = "N-Acetyl-2-hydroxybenzamide",
    "3-Chlorobiphenyl" = "3-Chlorobiphenyl",
    "Norketamine" = "Norketamine",
    "5-(Hexahydro-2-oxo-1H-thieno[3,4-D]imidazol-6-YL)pentanal" = "5-(Hexahydro-2-oxo-1H-thieno[3",
    "4-Nitrophenyl butyrate" = "4-Nitrophenyl butyrate",
    "Dihydrokainic acid" = "Dihydrokainic acid",
    "Tert-Butyl (2-(2-hydroxyethoxy)ethyl)carbamate" = "Tert-Butyl (2-(2-hydroxyethoxy",
    "1,5-Benzothiazepin-4(5H)-one, 2,3-dihydro-3-hydroxy-2-(4-methoxyphenyl)-5-(2-(methylamino)ethyl)-, (2S,3S)-" = "1,5-Benzothiazepin-4(5H)-one, ",
    "Telbivudine" = "Telbivudine",
    "2-Propylthiazolo[4,5-c]quinolin-4-amine" = "2-Propylthiazolo[4,5-c]quinoli",
    "3-(Phosphonomethyl)-5-phenyl-D-phenylalanine" = "3-(Phosphonomethyl)-5-phenyl-D",
    "Ethosuximide" = "Ethosuximide",
    "4,7-Dioxoheptanoic acid" = "4,7-Dioxoheptanoic acid",
    "Clofibrate" = "Clofibrate",
    "6-Hydroxymethylpterin" = "2-Amino-6-Hydroxymethyl-3H-Pteridin-4-one",
    "(2Z)-2-Morpholin-4-yliminoacetonitrile" = "(2Z)-2-Morpholin-4-yliminoacet",
    "1,4-Pentadien-3-one, 1,5-diphenyl-" = "1,4-Pentadien-3-one, 1,5-diphe",
    "N,N-Dimethylacrylamide" = "N,N-Dimethylacrylamide",
    "1-Methyl-4-piperidyl acetate" = "1-Methyl-4-piperidyl acetate",
    "Trazodone" = "Trazodone",
    "Proline betaine" = "Proline betaine",
    "(2R,3S,5R)-5-(4,5-Dihydroimidazo[2,1-f]purin-3-yl)-2-(hydroxymethyl)oxolan-3-ol" = "(2R,3S,5R)-5-(4,5-Dihydroimida",
    "Sulforaphane-cysteine-glycine" = "Sulforaphane-cysteine-glycine",
    "3-Ethylphenylsulfate" = "3-Ethylphenyl sulfate",
    "Progabide" = "Progabide",
    "Beclobrate" = "Beclobrate",
    "Dimethenamid" = "Dimethenamid",
    "Pirotiodecane" = "Pirotiodecane",
    "Diaminopimelic acid" = "Diaminopimelic acid",
    "Terbufos" = "Terbufos",
    "2-(N-Morpholino)-ethanesulfonic acid" = "2-(N-Morpholino)ethanesulfonic acid",
    "6-Diazo-5-oxo-L-norleucine" = "6-Diazo-5-oxo-L-norleucine",
    "Parbendazole" = "Parbendazole",
    "O-Succinyl-L-homoserine" = "o-Succinylhomoserine",
    "1-Heptanesulfonic acid" = "1-Heptanesulfonic acid",
    "Allyl-isopropyl-acetylharnstoff" = "Allyl-isopropyl-acetylharnstof",
    "1,3,7-Trimethyl-8-(3-chlorostyryl)xanthine" = "1,3,7-Trimethyl-8-(3-chlorosty"
  )
)

# ============================================================
# PARAMETERS$microbiome — FK49 16S rRNA microbiome analysis
# ============================================================
PARAMETERS$microbiome <- list(
  tax_levels           = c("Phylum", "Class", "Order", "Family", "Genus"),
  alpha_measures       = c("Observed", "Chao1", "Shannon", "Simpson", "InvSimpson"),
  prevalence_threshold = 0.05,      # keep taxa present in >= 5% of samples
  min_reads            = 2000,      # QC filter: remove samples below this
  timepoints           = c("F1", "F2", "F3", "F4"),   # F5 excluded (sequencing failure)
  da_methods           = c("ancombc2", "aldex2"),
  fdr_method           = "BH",
  fdr_threshold        = 0.05,
  feces_order          = c("F1", "F2", "F3", "F4"),
  feces_labels         = c(F1 = "-1 wks (ND)", F2 = "0 wks (ND)",
                           F3 = "3 wks (CDHFD)", F4 = "7 wks (CDHFD)"),
  diet_short_colors    = c(ND = "darkorange3", CDHFD = "darkviolet"),
  target_genera        = c("Lactobacillus")  # hypothesis-driven targeted analysis
)


# # Note Reactom enriched pathways hierachy -----
# Metabolism
# │
# ├── Biological oxidations 
# |   |
# |   |── Aflatoxin activation and detoxification
# │   ├── Phase I - Functionalization of compounds
# │   │   └── Cytochrome P450 - arranged by substrate type
# │   │       ├── Xenobiotics
# │   │       │  
# │   │       └── [weitere nicht-signifikante Reactome-Zweige]
# │   │
# │   └── Phase II - Conjugation of compounds
# │       └── Glutathione conjugation
# │
# ├── Metabolism of lipids
# │   ├── Metabolism of steroids
# │   │   └── Bile acid and bile salt metabolism
# │   │       └── Recycling of bile acids and salts
# │   │
# │   ├── Fatty acid metabolism
# │   │
# │   └── Biosynthesis of specialized proresolving mediators (SPMs)
# │       └── Biosynthesis of DHA-derived SPMs
# │
# ├── Drug ADME
# │   ├── Aspirin ADME
# │   ├── Prednisone ADME
# │   └── Paracetamol ADME
# │
# ├── Platelet activation, signaling and aggregation
# │   ├── Platelet degranulation
# │   │   └── Response to elevated platelet cytosolic Ca2+
#   │   │
# │   ├── GRB2:SOS provides linkage to MAPK signaling for Integrins
# │   │   └── [auch unter Integrin signaling / Platelet Aggregation]
# │   │
# │   └── p130Cas linkage to MAPK signaling for integrins
# │       └── [auch unter Integrin signaling / Platelet Aggregation]
# │
# ├── Coagulation pathway
# │   └── Fibrin formation
# │
# ├── Regulation of TLR by endogenous ligand
# │   ├── Innate Immune System
# │   └── Toll-like Receptor Cascades
# │
# ├── Post-translational protein modification
# │   └── Post-translational protein phosphorylation
# │
# └── Regulation of Insulin-like Growth Factor (IGF) transport
# and uptake by Insulin-like Growth Factor Binding Proteins (IGFBPs)
# 
PARAMETERS$Proteom <- list(
  Pathway_parents =c(
    "Biological oxidations",
    "Metabolism of lipids",
    "Drug ADME",
    "Platelet activation, signaling and aggregation",
    "Coagulation pathway",
    "Regulation of TLR by endogenous ligand",
    "Post-translational protein modification",
    "Regulation of Insulin-like Growth Factor (IGF) transport
      and uptake by Insulin-like Growth Factor Binding Proteins (IGFBPs)"
  ),
  Children_of_parents = c(
    # Biological oxidations
    "Phase I - Functionalization of compounds",
    "Cytochrome P450 - arranged by substrate type",
    "Xenobiotics",
    "Aflatoxin activation and detoxification",
    "Phase II - Conjugation of compounds",
    "Glutathione conjugation",
  
    # Metabolism of lipids
    "Metabolism of steroids",
    "Bile acid and bile salt metabolism",
    "Recycling of bile acids and salts",
    "Fatty acid metabolism",
    "Biosynthesis of specialized proresolving mediators (SPMs)",
    "Biosynthesis of DHA-derived SPMs",
  
    # Drug ADME
    "Aspirin ADME",
    "Prednisone ADME",
    "Paracetamol ADME",
    
    # Platelet activation, signaling and aggregation
    "Platelet degranulation",
    "Response to elevated platelet cytosolic Ca2+",
    "GRB2:SOS provides linkage to MAPK signaling for Integrins",
    "p130Cas linkage to MAPK signaling for integrins",
  
    # Coagulation pathway
    "Fibrin formation",
  
    # Regulation of TLR by endogenous ligand
    # keine signifikanten Kinder in deiner Liste
  
    # Post-translational protein modification
    "Post-translational protein phosphorylation"
  ),
Pathway_lowest = c(
    # Biological oxidations
    "Aflatoxin activation and detoxification",
    #"Cytochrome P450 - arranged by substrate type",
    "Glutathione conjugation",
  
    # Metabolism of lipids
    "Recycling of bile acids and salts",
    "Fatty acid metabolism",
    "Biosynthesis of DHA-derived SPMs",
  
    # Drug ADME
    "Aspirin ADME",
    "Prednisone ADME",
    "Paracetamol ADME",
  
    # Platelet activation, signaling and aggregation
    "Response to elevated platelet cytosolic Ca2+",
    "GRB2:SOS provides linkage to MAPK signaling for Integrins",
    "p130Cas linkage to MAPK signaling for integrins",
  
    # Coagulation pathway
    "Fibrin formation",
  
    #Regulation of TLR by endogenous ligand
      "Toll-like Receptor Cascades",
      "Innate Immune System",
  
    # Post-translational protein modification
    "Post-translational protein phosphorylation",
    # Regulation of Insulin-like Growth Factor (IGF) transport and uptake by Insulin-like Growth Factor Binding Proteins (IGFBPs)
    "Regulation of Insulin-like Growth Factor (IGF) transport
      and uptake by Insulin-like Growth Factor Binding Proteins (IGFBPs)"
  )
)