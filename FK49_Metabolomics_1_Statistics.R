gc()
rm(list = ls())
home <- normalizePath("~") # bc Windows does not start inuser dir but in user/documents dir
parent <- dirname(home)
setwd(paste0(parent,"/OneDrive - Universität Salzburg/AG_Tumorimmunologie - Dokumente/Data/Freia Krause/01_Experiments/FK46_FK49_RProject"))

library(dplyr)
library(tidyr)
library(ggplot2)
library('corrr')
library(ggcorrplot)
library("FactoMineR")
library(dplyr)
library(factoextra)
library(ggrepel)
library(emmeans)
source("FK49_Definitions.R")
home <- normalizePath("~") # bc Windows does not start inuser dir but in user/documents dir
parent <- dirname(home) # goves users dir works for me since i aways have the data from work in the onedrive at the same location and only local home dir changes
ExpId = "FK49"

process_metabolome <- function(df, meta_cols) {
  df %>%
    dplyr::mutate(
      Sample = as.factor(Sample),
      T_D_S = paste0(Treatment, "_", Diet, "_", Sex),
      T_D   = paste0(Treatment, "_", Diet),
      T_S   = paste0(Treatment, "_", Sex)
    ) %>%
    dplyr::mutate(
      across(where(is.character), ~na_if(., "<LOD")),
      across(!all_of(meta_cols), as.numeric),
      across(where(is.numeric), ~ifelse(is.na(.), min(., na.rm = TRUE)/2, .))
    ) %>%
    dplyr::mutate(
      Treatment = factor(Treatment, levels = c("EtOH","TAM")),
      Sex       = factor(Sex, levels = c("female","male")),
      Diet      = factor(Diet, levels = c("ND","CDHFD13"))
    )
}

if (ExpId=="FK49") {
  meta_cols <- c("Animal","Sample","Sex","Treatment","Diet","ExpID","T_D_S","T_D" ,"T_S")
  metabolome_positive <- read.csv(paste0(parent,"/OneDrive - Universität Salzburg/AG_Tumorimmunologie - Dokumente/Data/Freia Krause/01_Experiments/FK49_CD-HFD_13wks/FK49_untargetedLiverMetabolomics/Report_M086_untargeted_20260113_positive.csv"), 
                                  sep=";", stringsAsFactors=FALSE, check.names=FALSE) %>%
    dplyr::mutate(Sample = SampleNo) %>%
    dplyr::select(-SampleCode, -SampleNo) %>%
    process_metabolome(meta_cols)
  saveRDS(metabolome_positive, file= paste0(parent,"/OneDrive - Universität Salzburg/AG_Tumorimmunologie - Dokumente/Data/Freia Krause/01_Experiments/FK49_CD-HFD_13wks/FK49_Analysis/01_RawData/FK49_metabolome_positive.rds"))
 
   metabolome_negative <- read.csv(paste0(parent,"/OneDrive - Universität Salzburg/AG_Tumorimmunologie - Dokumente/Data/Freia Krause/01_Experiments/FK49_CD-HFD_13wks/FK49_untargetedLiverMetabolomics/Report_M086_untargeted_20260113_negative.csv"),
                                    sep=";",stringsAsFactors = FALSE,check.names=FALSE) %>%  
                          dplyr::mutate(Sample = SampleNo) %>%
                          dplyr::select(-SampleCode, -SampleNo) %>%
                          process_metabolome(meta_cols)
  saveRDS(metabolome_negative, file= paste0(parent,"/OneDrive - Universität Salzburg/AG_Tumorimmunologie - Dokumente/Data/Freia Krause/01_Experiments/FK49_CD-HFD_13wks/FK49_Analysis/01_RawData/FK49_metabolome_negative.rds"))
  
  
  metabolome_targeted <- read.csv(paste0(parent,"/OneDrive - Universität Salzburg/AG_Tumorimmunologie - Dokumente/Data/Freia Krause/01_Experiments/FK49_CD-HFD_13wks/FK49_HILIC09/Report_M086_HILIC09_20251222.csv"),
                                  sep=";" , check.names=FALSE) %>%
                        dplyr::select(-SampleID) %>%
                        dplyr::rename(Sample = SampleNo) %>%
                        process_metabolome(meta_cols)
  saveRDS(metabolome_targeted, file= paste0(parent,"/OneDrive - Universität Salzburg/AG_Tumorimmunologie - Dokumente/Data/Freia Krause/01_Experiments/FK49_CD-HFD_13wks/FK49_Analysis/01_RawData/FK49_metabolome_targeted.rds"))
  
  targeted_pwd <- paste0(parent, "/OneDrive - Universität Salzburg/AG_Tumorimmunologie - Dokumente/Data/Freia Krause/01_Experiments/FK49_CD-HFD_13wks/FK49_Analysis/02_GeneratedData/targetedLivMetabolome")
  untargeted_pwd <- paste0(parent, "/OneDrive - Universität Salzburg/AG_Tumorimmunologie - Dokumente/Data/Freia Krause/01_Experiments/FK49_CD-HFD_13wks/FK49_Analysis/02_GeneratedData/untargetedLivMetabolome")
  
  # Define the output folder names
  output_folders <- c("CDHFD_a", "CDHFD_m", "CDHFD_f", "ND_f","NDvsCDHFD")
  
  
  # Create the output folders under both targeted and untargeted paths
  create_output_folders(targeted_pwd, output_folders)
  create_output_folders(untargeted_pwd, output_folders)

}  else if (ExpId == "FK46"){
  print("You dont have data for this Experiment")
} else{
  print("You dont have data for this Experiment")
  }


no_filter         <- list(Sex_filter  = c("male", "female"),  Diet_filter = c("CDHFD13","ND"),      ExpID_filter = c("FK49","BH"))
CDHFD_a_filters   <- list(Sex_filter  = c("male", "female"),  Diet_filter = c("CDHFD13"),      ExpID_filter = c("FK49"))
CDHFD_f_filters   <- list(Sex_filter  = c("female"),          Diet_filter = c("CDHFD13"),      ExpID_filter = c("FK49"))
CDHFD_m_filters   <- list(Sex_filter  = c("male"),            Diet_filter = c("CDHFD13"),      ExpID_filter = c("FK49"))
ND_f_filters      <- list(Sex_filter  = c("female"),          Diet_filter = c("ND"),      ExpID_filter = c("BH")) # ND experiment only has females
NDvsCHFD_filter   <- list(Sex_filter  = c("female"),          Diet_filter = c("ND","CDHFD13"),      ExpID_filter = c("BH","FK49")) # 

## Subsets Dataset in subset and performs data preprocessing -----
### normalization, log transform, scaling
subset_data <- function(data,
                        Sex_filter  = c("male", "female"), 
                        Diet_filter = c("CDHFD13", "ND"), 
                        ExpID_filter = c("FK49", "BH"),
                        method = NULL
                        ) {
  
  
  
  # Subset
  dataset <- data %>%
        dplyr::filter( Sex  %in% Sex_filter,
        Diet %in% Diet_filter,
        ExpID %in% ExpID_filter    )
  
  # Metadata columns
  meta_cols <- c( "Sample", "Animal", "Sex", "Treatment", "Diet", "ExpID", "T_D_S", "T_D", "T_S")
  
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
 # I save raw data normalized data log data and scaled data
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
sub_data <-do.call(subset_data,c(list(data = metabolome_negative, method = "untargeted"), CDHFD_a_filters))
ggsave(plot= sub_data$Preprocessing, width= 12, height= 6, dpi=300, filename="BG_prepocessing_CDHFD_A_neg.png", path =paste0(untargeted_pwd,"/CDHFD_a") )

sub_data <-do.call(subset_data,c(list(data = metabolome_positive, method = "untargeted"), CDHFD_a_filters))
ggsave(plot= sub_data$Preprocessing, width= 12, height= 6, dpi=300, filename="BG_prepocessing_CDHFD_A_pos.png", path =paste0(untargeted_pwd,"/CDHFD_a"))

sub_data <-do.call(subset_data,c(list(data = metabolome_positive, method = "untargeted"), ND_f_filters))
ggsave(plot= sub_data$Preprocessing, width= 12, height= 6, dpi=300, filename="BG_prepocessing_ND_f_pos.png", path =paste0(untargeted_pwd,"/ND_f"))

sub_data <-do.call(subset_data, c(list(data = metabolome_negative, method = "untargeted"), ND_f_filters))
ggsave(plot= sub_data$Preprocessing, width= 12, height= 6, dpi=300, filename="BG_prepocessing_ND_f_neg.png", path =paste0(untargeted_pwd,"/ND_f"))


sub_data <-do.call(subset_data,c(list(data = metabolome_negative, method = "untargeted"), CDHFD_f_filters))
ggsave(plot= sub_data$Preprocessing, width= 12, height= 6, dpi=300, filename="BG_prepocessing_CDHFD_f_neg.png", path =paste0(untargeted_pwd,"/CDHFD_f"))
sub_data <-do.call(subset_data,c(list(data = metabolome_positive, method = "untargeted"), CDHFD_f_filters))
ggsave(plot= sub_data$Preprocessing, width= 12, height= 6, dpi=300, filename="BG_prepocessing_CDHFD_f_pos.png", path =paste0(untargeted_pwd,"/CDHFD_f"))

sub_data <-do.call(subset_data,c(list(data = metabolome_negative, method = "untargeted"), CDHFD_m_filters))
ggsave(plot= sub_data$Preprocessing, width= 12, height= 6, dpi=300, filename="BG_prepocessing_CDHFD_m_neg.png", path =paste0(untargeted_pwd,"/CDHFD_m"))
sub_data <-do.call(subset_data,c(list(data = metabolome_positive, method = "untargeted"), CDHFD_m_filters))
ggsave(plot= sub_data$Preprocessing, width= 12, height= 6, dpi=300, filename="BG_prepocessing_CDHFD_m_pos.png", path =paste0(untargeted_pwd,"/CDHFD_m"))

sub_data <-do.call(subset_data,c(list(data = metabolome_targeted, method = "targeted"), ND_f_filters))
ggsave(plot= sub_data$Preprocessing, width= 9, height= 6, dpi=300, filename="BG_prepocessing_ND_f.png", path =paste0(targeted_pwd,"/ND_f"))
sub_data <-do.call(subset_data,c(list(data = metabolome_targeted, method = "targeted"), CDHFD_a_filters))
ggsave(plot= sub_data$Preprocessing, width= 9, height= 6, dpi=300, filename="BG_prepocessing_CDHFD_A.png", path =paste0(targeted_pwd,"/CDHFD_a"))
sub_data <-do.call(subset_data,c(list(data = metabolome_targeted, method = "targeted"), CDHFD_m_filters))
ggsave(plot= sub_data$Preprocessing, width= 9, height= 6, dpi=300, filename="BG_prepocessing_CDHFD_m.png", path =paste0(targeted_pwd,"/CDHFD_m"))
sub_data <-do.call(subset_data,c(list(data = metabolome_targeted, method = "targeted"), CDHFD_f_filters))
ggsave(plot= sub_data$Preprocessing, width= 9, height= 6, dpi=300, filename="BG_prepocessing_CDHFD_f.png", path =paste0(targeted_pwd,"/CDHFD_f"))

plot_PCA<- function(data,
                    plot_variables= "Diet",
                    Sex_filter,
                    Diet_filter,
                    ExpID_filter, 
                    NAME= "",
                    method = NULL,
                    FOLDER= NULL){
  
  df<-subset_data(data,
                  Sex_filter =Sex_filter,
                  Diet_filter = Diet_filter,
                  ExpID_filter = ExpID_filter,
                  method = method)
  

  if(method == "targeted"){
    PATH = targeted_pwd
  }else if (method == "untargeted"){
    PATH = untargeted_pwd
  }else (message(" Select a method 'targeted' or 'untargeted' so that we can correctly save the data."))
  
  
  pca <- prcomp(df$scaled_log_values, center = FALSE, scale. = FALSE)
  
  # Summary of PCA
  pca_summary <- summary(pca)
  pc1_variance <- pca_summary$importance[2, 1]  
  pc2_variance <- pca_summary$importance[2, 2]
  
  loadings <- pca$rotation
  
  pc1_loadings <- loadings[, 1]  
  pc2_loadings <- loadings[, 2]  
  
  
  top_pc1_contributors <- sort(abs(pc1_loadings), decreasing = TRUE)
  top_pc1_variables <- names(top_pc1_contributors)
    top_pc2_contributors <- sort(abs(pc2_loadings), decreasing = TRUE)
  top_pc2_variables <- names(top_pc2_contributors)
  
  # Prepare the data for writing to CSV
  top_contributors <- data.frame(
    PC1_Variables = top_pc1_variables[1:10],
    PC2_Variables = top_pc2_variables[1:10]
  )
  
  # Write the data to a CSV file
  write.csv(top_contributors, file = paste0(PATH,"/",FOLDER,"/PC1_PC2_Top10_",NAME,".csv"), row.names = FALSE)
  cat("Top contributors saved to CSV:", paste0(PATH,"/",FOLDER,"/PC1_PC2_Top10_",NAME,".csv"), "\n")
  rm(pc1_loadings,pc2_loadings,top_pc1_contributors,top_pc1_variables,top_pc2_contributors,top_pc2_variables,top_contributors)
  gc()
  # 
  # # Visualizations
  # eigen<-fviz_eig(pca, addlabels = TRUE)
  # print(eigen)
  # 
  # variable<-fviz_pca_var(pca, col.var = "black")
  # print(variable)
  # 
  #  variable_color<-fviz_pca_var(pca, col.var = "cos2",
  #                               alpha.var="contrib",
  #                               select.var = list(cos2 = 0.8),
  #                               #gradient.cols = c("black", "orange", "green"),
  #                               repel = TRUE)
  #  print(variable_color)
  # 
  #  cos<-fviz_cos2(pca, choice = "var", axes = 1:2)
  #  print(cos)
  
  for (v in plot_variables){
   
     if (v == "Sex") {
      color_map <- Sex_colors
    } else if (v == "Treatment") {
      color_map <- Treatment_colors
    } else if (v == "T_D_S") {
      color_map <- T_D_S_colors
    } else if (v == "T_S") {
      color_map <- T_S_colors
    } else if (v == "Diet") {
      color_map <- Diet_colors
    } else {
      stop("Unknown plot variable. Please provide a valid variable name.")
    }
    
    plottiplot<-fviz_pca_ind(pca, 
                             habillage = df$metadata[[v]],
                             addEllipses = TRUE,
                             ellipse.alpha=0,
                             geom="point",
                             pointsize = 4.2, 
                             invisible="quali")+
      labs(title = "PCA", 
           x = paste("PC1 (", round(pc1_variance * 100, 2), "%)", sep = ""),
           y = paste("PC2 (", round(pc2_variance * 100, 2), "%)", sep = "")) +
      scale_color_manual(values = color_map) +
      theme_classic()
      
    print(plottiplot)
    ggsave(plot = plottiplot, width= 12, height= 12, dpi=300, filename=paste0("PCA_", NAME,"_",v, ".png"),path = paste0(PATH,"/",FOLDER ) )
    
  }
    gc()
}  


do.call(plot_PCA,c(list( data = metabolome_positive,    plot_variables = c("Sex","Treatment","T_S"),FOLDER = "CDHFD_a", NAME= "CDHFD_a_pos", method = "untargeted" ), CDHFD_a_filters))
do.call(plot_PCA,c(list( data = metabolome_negative,    plot_variables = c("Sex","Treatment","T_S"),FOLDER = "CDHFD_a", NAME= "CDHFD_a_neg" , method = "untargeted" ), CDHFD_a_filters))
do.call(plot_PCA,c(list( data = metabolome_positive,    plot_variables = c("Treatment"),FOLDER = "ND_f", NAME= "ND_f_pos" , method = "untargeted"),  ND_f_filters))
do.call(plot_PCA,c(list( data = metabolome_negative,    plot_variables = c("Treatment"),FOLDER = "ND_f", NAME= "ND_f_neg" , method = "untargeted"),  ND_f_filters))
do.call(plot_PCA,c(list( data = metabolome_positive,    plot_variables = c("Treatment"),FOLDER = "CDHFD_m", NAME= "CDHFD_m_pos", method = "untargeted" ),  CDHFD_m_filters))
do.call(plot_PCA,c(list( data = metabolome_negative,    plot_variables = c("Treatment"),FOLDER = "CDHFD_m", NAME= "CDHFD_m_neg" , method = "untargeted" ), CDHFD_m_filters))
do.call(plot_PCA,c(list( data = metabolome_positive,    plot_variables = c("Treatment"),FOLDER = "CDHFD_f", NAME= "CDHFD_f_pos", method = "untargeted" ),  CDHFD_f_filters))
do.call(plot_PCA,c(list( data = metabolome_negative,    plot_variables = c("Treatment"),FOLDER = "CDHFD_f", NAME= "CDHFD_f_neg", method = "untargeted" ),  CDHFD_f_filters))

do.call(plot_PCA,c(list( data = metabolome_negative,    plot_variables = c("Diet"),FOLDER = "NDvsCDHFD", NAME= "NDvsCDHFD_neg", method = "untargeted" ),  NDvsCHFD_filter))
do.call(plot_PCA,c(list( data = metabolome_positive,    plot_variables = c("Diet"),FOLDER = "NDvsCDHFD", NAME= "NDvsCDHFD_pos", method = "untargeted" ),  NDvsCHFD_filter))


do.call(plot_PCA,c(list( data = metabolome_targeted,    plot_variables = c("Treatment"), FOLDER = "CDHFD_a", NAME= "CDHFD_a", method = "targeted" ),  CDHFD_a_filters))
do.call(plot_PCA,c(list( data = metabolome_targeted,    plot_variables = c("Treatment"), FOLDER = "CDHFD_m", NAME= "CDHFD_m", method = "targeted" ),  CDHFD_f_filters))
do.call(plot_PCA,c(list( data = metabolome_targeted,    plot_variables = c("Treatment"), FOLDER = "CDHFD_f",NAME= "CDHFD_f", method = "targeted" ),  CDHFD_m_filters))
do.call(plot_PCA,c(list( data = metabolome_targeted,    plot_variables = c("Treatment"), FOLDER = "ND_f", NAME= "ND_f", method = "targeted" ),  ND_f_filters))

do.call(plot_PCA,c(list( data = metabolome_targeted,    plot_variables = c("Treatment"), FOLDER = "NDvsCDHFD", NAME= "NDvsCDHFD", method = "targeted" ),  NDvsCHFD_filter))



perform_ttest <- function(data,
                          comparison = "Treatment",
                          Sex_filter,
                          Diet_filter,
                          ExpID_filter,
                          top_n = 8,
                          data_type = c("log", "raw", "norm", "scaled"),
                          use_lm = FALSE,
                          XLIM=c(-9,9),
                          YLIM=c(0,7.5),
                          NAME= "",
                          method = NULL,
                          FOLDER = NULL
                          )
                          {
  
  
  data_type <- match.arg(data_type)
  as.character(data_type)
  cat("Data type received:", data_type, "\n")
  
  if(method == "targeted"){
    PATH = targeted_pwd
  }else if (method == "untargeted"){
    PATH = untargeted_pwd
  }else (message(" Select a method 'targeted' or 'untargeted' so that we can correctly save the data."))
  
  # Subset the data
  df <- subset_data(data, 
                    Sex_filter = Sex_filter, 
                    Diet_filter = Diet_filter, 
                    ExpID_filter = ExpID_filter,
                    method = method)
  
  metab_names <- colnames(df$log_values)
 
 if (data_type == "raw") {
   selected_values <- df$raw_values
  } else if (data_type == "log") {
    selected_values <- df$log_values
  } else if (data_type == "norm") {
    selected_values <- df$norm_values
  } else if (data_type == "scaled") {
    selected_values <- df$scaled_log_values
  } else {
    stop("Invalid data type")  # Error if data_type is not valid
  }
  
  for (c in comparison){
    
    if (c == "Sex") {
      color_map <- Sex_colors
    } else if (c == "Treatment") {
      color_map <- Treatment_colors
    } else if (c == "T_D_S") {
      color_map <- T_D_S_colors
    } else if (c == "T_S") {
      color_map <- T_S_colors
    } else if (c == "Diet") {
      color_map <- Diet_colors
    } else {
      stop("Unknown plot variable. Please provide a valid variable name.")
    }
  }
  
  
  if(use_lm){
    lm_metab <- function(metab_name){
      values <- selected_values[[metab_name]]
      df_metab <- data.frame(
        value = values,
        Treatment = as.factor(df$metadata[[comparison]]),
        Sex = as.factor(df$metadata$Sex) )
      
      fit <- lm(value ~ Treatment * Sex, data = df_metab)
      coefs <- summary(fit)$coefficients
      
      # log2FC estimates
      mean_TAM <- mean(df_metab$value[df_metab$Treatment=="TAM"])
      mean_EtOH <- mean(df_metab$value[df_metab$Treatment=="EtOH"])
      log2FC_overall <- mean_TAM - mean_EtOH
      log2FC_female <- mean(df_metab$value[df_metab$Treatment=="TAM" & df_metab$Sex=="female"]) -
        mean(df_metab$value[df_metab$Treatment=="EtOH" & df_metab$Sex=="female"])
      log2FC_male <- mean(df_metab$value[df_metab$Treatment=="TAM" & df_metab$Sex=="male"]) -
        mean(df_metab$value[df_metab$Treatment=="EtOH" & df_metab$Sex=="male"])
      
      res <- data.frame(
        Metabolite = metab_name,
        Intercept  = coefs["(Intercept)", "Estimate"],
        Treatment  = coefs["TreatmentTAM", "Estimate"],
        Sex = coefs["Sexmale", "Estimate"],
        Treatment_Sex = coefs["TreatmentTAM:Sexmale", "Estimate"],
        p_Intercept = coefs["(Intercept)", "Pr(>|t|)"],
        p_Treatment = coefs["TreatmentTAM", "Pr(>|t|)"],
        p_Sex = coefs["Sexmale", "Pr(>|t|)"],
        p_Treatment_Sex = coefs["TreatmentTAM:Sexmale", "Pr(>|t|)"],
        log2FC = log2FC_overall,
        log2FC_female = log2FC_female,
        log2FC_male = log2FC_male,
        mean_TAM = mean_TAM,
        mean_EtOH = mean_EtOH,
        sd_TAM = sd(df_metab$value[df_metab$Treatment == "TAM"]),
        sd_EtOH = sd(df_metab$value[df_metab$Treatment == "EtOH"]),
        stringsAsFactors = FALSE,
        Comparison = paste0("Down in TAM Up in TAM")
      )
      
      return(list(results = res, df_metab = df_metab, fit = fit))
    }
    
    # Run LM for all metabolites
    lm_list <- lapply(metab_names, lm_metab)
    results <- do.call(rbind, lapply(lm_list, `[[`, "results"))
    
    # Adjust p-values
    results$adj.p.value <- p.adjust(results$p_Treatment, method="fdr")
    results$negLog10FDR <- -log10(results$adj.p.value)
    results$significant <- results$adj.p.value < 0.05 & abs(results$log2FC) > 0.5  
    results$trend <- results$adj.p.value < 0.1 & abs(results$log2FC) > 0.3
    results$direction <-case_when(results$adj.p.value < 0.05 &results$log2FC > 0.5  ~ "UP",
                                  results$adj.p.value < 0.05 &results$log2FC < -0.5 ~ "DOWN",
                                  TRUE ~"NS")
    interesting_results <- results %>%  arrange(adj.p.value, desc(abs(log2FC)))
    significant_results <-results %>% filter(significant == TRUE)
    top_labels <- interesting_results %>% head(top_n)
    top_metabolites <- head(interesting_results$Metabolite, top_n)
    all_results_saving<-results%>% 
      arrange(log2FC) %>%
      dplyr::select(
        Metabolite, 
        p_Treatment, 
        adj.p.value, 
        log2FC, 
        log2FC_female, 
        log2FC_male, 
        mean_TAM, 
        mean_EtOH, 
        sd_TAM, 
        sd_EtOH
      )
    write.csv2(all_results_saving, file = paste0(PATH,"/",FOLDER,"/",NAME,"_LM_results.csv"), row.names = FALSE)
    
    # After computing the results and determining the trend metabolites
    trend_results <- results %>% 
      filter(trend == TRUE) %>%
      arrange(log2FC) %>%
      dplyr::select(
        Metabolite, 
        p_Treatment, 
        adj.p.value, 
        log2FC, 
        log2FC_female, 
        log2FC_male, 
        mean_TAM, 
        mean_EtOH, 
        sd_TAM, 
        sd_EtOH
      )
    
    # Save the trend metabolites data to a CSV
    write.csv2(trend_results, file = paste0(PATH,"/",FOLDER,"/",NAME,"_trend_metabolites_LM.csv"), row.names = FALSE)
    
    cat("Saved trend metabolites to ", paste0(PATH,"/",FOLDER,"/",NAME,"_trend_metabolites_LM.csv"), "\n")
    
    sig_results <- results %>% 
      filter(significant == TRUE) %>%
      arrange(log2FC) %>%
      dplyr::select(
        Metabolite, 
        p_Treatment, 
        adj.p.value, 
        log2FC, 
        log2FC_female, 
        log2FC_male, 
        mean_TAM, 
        mean_EtOH, 
        sd_TAM, 
        sd_EtOH
      )
    
    # Save the sig metabolites data to a CSV
    write.csv2(sig_results, file = paste0(PATH,"/",FOLDER,"/",NAME,"_sig_metabolites_LM.csv"), row.names = FALSE)
    
    cat("Saved sig metabolites to ", paste0(PATH,"/",FOLDER,"/",NAME,"_sig_metabolites_LM.csv"), "\n")
    
    # ==== LM Visualizations for first top metabolite ====
    # Coefficients bar plot
    first_fit <- lm_list[[which(metab_names == top_metabolites[1])]]$fit
    coefs_df <- as.data.frame(summary(first_fit)$coefficients)
    coefs_df$Term <- rownames(coefs_df)
    
    p_coefs <- ggplot(coefs_df[2:4,], aes(x=Term, y=Estimate)) +
      geom_bar(stat="identity", fill="steelblue") +
      geom_errorbar(aes(ymin=Estimate-`Std. Error`, ymax=Estimate+`Std. Error`)) +
      labs(y="Effect (log2 FC)", title=paste("LM Coefficients -", top_metabolites[1])) +
      theme_minimal()
    
    print(p_coefs)
    
    # Predicted means (emmeans)
    emm <- emmeans(first_fit, ~ Treatment * Sex)
    df_emm <- as.data.frame(emm)
    
    p_emm <- ggplot(df_emm, aes(x=Treatment, y=emmean, fill=Sex)) +
      geom_bar(stat="identity", position=position_dodge()) +
      geom_errorbar(aes(ymin=lower.CL, ymax=upper.CL), position=position_dodge(0.9), width=0.2) +
      labs(y="Predicted log2 intensity", title=paste("Predicted Means -", top_metabolites[1])) +
      theme_minimal()
    
    print(p_emm)
    
  } else {
    # T-test version
    ttest_metab <- function(metab_name){
      values <- selected_values[[metab_name]]
      df_metab <- data.frame(
        value = values,
        group = as.factor(df$metadata[[comparison]])
      )
      
      t_res <- t.test(value ~ group, data = df_metab)
      log2FC <- mean(df_metab$value[df_metab$group == levels(df_metab$group)[2]]) -
        mean(df_metab$value[df_metab$group == levels(df_metab$group)[1]])
     
       
      data.frame(Metabolite = metab_name, 
                 p.value = t_res$p.value, 
                 log2FC = log2FC,
                 mean_1 = mean(df_metab$value[df_metab$group == levels(df_metab$group)[1]]),
                 mean_2 =  mean(df_metab$value[df_metab$group == levels(df_metab$group)[2]]),
                 sd_1   =sd(df_metab$value[df_metab$group == levels(df_metab$group)[1]]),
                 sd_2  = sd(df_metab$value[df_metab$group == levels(df_metab$group)[2]]),
                 Comparison = paste0("Down in ", levels(df_metab$group)[2], " | Up in ", levels(df_metab$group)[2])
      )
    }
    
    results <- lapply(metab_names, ttest_metab) %>% do.call(rbind, .)
    results$adj.p.value <- p.adjust(results$p.value, method="fdr")
    results$negLog10FDR <- -log10(results$adj.p.value)
    results$significant <- results$adj.p.value < 0.05 & abs(results$log2FC) > 0.5
    results$trend <- results$adj.p.value < 0.1 & abs(results$log2FC) > 0.3
    results$Comparison <-results$Comparison
    results$direction <-case_when(results$adj.p.value < 0.05 &results$log2FC > 0.5  ~ "UP",
                                  results$adj.p.value < 0.05 &results$log2FC < -0.5 ~ "DOWN",
                                  TRUE ~"NS")
    interesting_results <- results %>%  arrange(adj.p.value, desc(abs(log2FC)))
    significant_results <-results %>% filter(significant == TRUE)
    top_labels <- interesting_results %>% head(top_n)
    top_metabolites <- head(interesting_results$Metabolite, top_n)
    
    all_results_saving<-results%>% 
      arrange(log2FC) %>%
      dplyr::select(
        Metabolite, 
        p.value, 
        adj.p.value, 
        log2FC, 
        mean_1, 
        mean_2, 
        sd_1, 
        sd_2
      )
    write.csv2(all_results_saving, file = paste0(PATH,"/",FOLDER,"/",NAME,"_test_results.csv"), row.names = FALSE)
    
    # After computing the results and determining the trend metabolites
    sig_results <- results %>% 
      filter(significant == TRUE) %>%
      arrange(log2FC) %>%
      dplyr::select(
        Metabolite, 
        p.value, 
        adj.p.value, 
        log2FC, 
        mean_1, 
        mean_2, 
        sd_1, 
        sd_2
      )
    
    # Save the trend metabolites data to a CSV
    write.csv2(sig_results, file = paste0(PATH,"/",FOLDER,"/",NAME, "_sig_metabolites_ttest.csv"), row.names = FALSE)
        cat("Saved sig metabolites to ", paste0(PATH,"/",FOLDER,"/",NAME, "_sig_metabolites_ttest.csv"), "\n")
        
        trend_results <- results %>% 
          filter(trend == TRUE) %>%
          arrange(log2FC) %>%
          dplyr::select(
            Metabolite, 
            p.value, 
            adj.p.value, 
            log2FC, 
            mean_1, 
            mean_2, 
            sd_1, 
            sd_2
          )
        # Save the trend metabolites data to a CSV
        write.csv2(trend_results, file = paste0(PATH,"/",FOLDER,"/",NAME, "_trend_metabolites_ttest.csv"), row.names = FALSE)
        
        cat("Saved trend metabolites to ", paste0(PATH,"/",FOLDER,"/",NAME, "_trend_metabolites.csv"), "\n")
    }
  
  # ==== Violin / Boxplots for Top Metabolites ====
  
  if(method == "untargeted"){
  df_values <- df$norm_values ## for violin visualization, I want to use the normalized values (relative abundance of metabolite per sample)
  y_label_method<-paste0("relative abundance [%]")
  } else if(method == "targeted"){
    df_values <-  df$raw_values
    y_label_method<-paste0("Abundance  [pM/mg]")
  }else(message("I need method 'untargeted' or 'targeted' to decide which data to use vor violin vizualization."))
  
  df_plot <- cbind(df$metadata, df_values[, top_metabolites, drop=FALSE])
  df_long <- df_plot %>%
    pivot_longer(cols = all_of(top_metabolites), names_to = "Metabolite", values_to = "value")
  
  df_long <- df_long %>%
    left_join(results[, c("Metabolite", "adj.p.value")], by="Metabolite")
  
  # Dynamically compute y-axis limits based on the range of values for each metabolite
  y_limits <- df_long %>%
    group_by(Metabolite) %>%
    summarise(y_min = min(value, na.rm = TRUE), y_max = max(value, na.rm = TRUE))
  
  # Create the violin plot
  p_violin <- ggplot(df_long, aes(x=.data[[comparison]], y=value, fill=.data[[comparison]])) +
    geom_violin(trim=FALSE, alpha=0.5) +
    geom_boxplot(width=0.1, outlier.shape=NA) +
    geom_jitter(width=0.1, size=1.5, alpha=0.7) +
    scale_fill_manual(values=color_map) +
    facet_wrap(~Metabolite, scales="free_y", ncol = 5, labeller = label_wrap_gen(width = 15)) +  # Wrap facet labels
    labs(y=y_label_method,
         x=comparison, title=paste0("Top ", top_n, " Metabolites")) +
    theme_classic() +
    theme(strip.text = element_text(size = 10, hjust = 0.5))
  
  # Add dynamic y-limits and p-value labels
  p_violin <- p_violin +
    geom_blank(data = y_limits, aes(y = y_min, yend = y_max, group = Metabolite), inherit.aes = FALSE) +  # Set the y-axis range dynamically
    geom_text(data = y_limits, 
              aes(x = 1.5, y = y_max + 0.05 * (y_max - y_min),  # Slightly above y_max for the label
                  label = paste0("adj.p = ", signif(results$adj.p.value[match(Metabolite, results$Metabolite)], 3))),
              inherit.aes = FALSE, size = 3, vjust = -0.5)  # Adjust vertical position as necessary
  
  # Print and save the plot
  print(p_violin)
  ggsave(plot= p_violin, filename= paste0(NAME,"_violins.png"), width= 18, height = 9, dpi = 300, path =paste0(PATH,"/",FOLDER,"/" ))
  # Volcano plot
  p_volcano <- ggplot(results,aes(x = log2FC, y = negLog10FDR)) +
    geom_point(aes(fill = direction), alpha = 0.5, size = 3,stroke = 0.5,
                 shape=21,color= "black")+
  scale_fill_manual(values = c("DOWN"="blue","NS"="grey60", "UP"="firebrick")) +
    geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed", color = "grey80") +
    geom_hline(yintercept = 0, linetype = "longdash", color = "grey50") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey80") +
    labs(
      title = paste(results$Comparison[1]),
      x = expression(paste("FC [", log[2], "]")),       
      y = expression(paste("-log"[10], "(adj.p.value)")))+
    theme_classic()+  
    theme(panel.grid= element_line(color ="grey90", linewidth = 0.1))+
    geom_text_repel(data = results %>% filter(significant == TRUE), 
                    aes(label = Metabolite),
                    size = 3.5,
                    max.overlaps = 25) +
    xlim(XLIM)+
    ylim(YLIM)
    
  
  print(p_volcano)
  ggsave(plot= p_volcano, filename= paste0(NAME,"_volcano.png"), width= 6, height = 9, dpi = 300, path =paste0(PATH,"/",FOLDER,"/" ))
  
  ### Do HeatMap
  
  interesting_metabolites <- results %>%
    filter(trend == TRUE) %>%
    pull(Metabolite )
  
  if(length(interesting_metabolites)>2){
    
   annotation_colors <- list(
    Sex = Sex_colors,
    Treatment = Treatment_colors,
    Diet = Diet_colors )
    
  heatmap_data <- t(selected_values[, interesting_metabolites, drop = FALSE])
  
  colnames(heatmap_data) <- df$metadata$Animal  # Assuming SampleID is in metadata
    ann <- data.frame(Treatment = df$metadata$Treatment,Sex=df$metadata$Sex, Diet= df$metadata$Diet)
  rownames(ann) <- df$metadata$Animal
  heatmap_height <-nrow(heatmap_data)/8+3
  
  p_heat<- pheatmap::pheatmap(heatmap_data,
           scale = "row",
           cluster_rows = TRUE, 
           cluster_cols = TRUE,
           annotation_col  = ann,
           annotation_colors = annotation_colors)
  
  print(p_heat)
  ggsave(plot= p_heat, filename= paste0(NAME,"_heatmap_clutered.png"),limitsize = FALSE, 
         width= 20, height = heatmap_height, dpi = 300,bg = "white", path =paste0(PATH,"/",FOLDER,"/" ))
  p_heat<- pheatmap::pheatmap(heatmap_data,
                              scale = "row",
                              cluster_rows = FALSE, 
                              cluster_cols = FALSE,
                              annotation_col  = ann,
                              annotation_colors = annotation_colors)
  
  print(p_heat)
  ggsave(plot= p_heat, filename= paste0(NAME,"_heatmap.png"),limitsize = FALSE, 
         width= 20, height = heatmap_height, dpi = 300,bg = "white", path =paste0(PATH,"/",FOLDER,"/" ))
  } else {
    print("Not enough significant changes for heatmap.")
}
  
  
        
  return(results)
}

# Helper function to run t-test and save result
run_and_save <- function(data, NAME, comparison, top_n = 5, data_type = "log",
                         use_lm = FALSE, XLIM = c(-5,5), YLIM = c(0,5), 
                         method = "untargeted", filters, FOLDER = NULL) {
  
  # Run the perform_ttest function
  result <- do.call(perform_ttest, c(list(
    data = data,
    NAME = NAME,
    comparison = comparison,
    top_n = top_n,
    data_type = data_type,
    use_lm = use_lm,
    XLIM = XLIM,
    YLIM = YLIM,
    method = method,
    FOLDER=FOLDER
  ), filters))
  
  # Choose folder based on method
  out_folder <- ifelse(method == "targeted", targeted_pwd, untargeted_pwd)
  out_folder =paste0(out_folder,"/",FOLDER,"/")
  # Make folder if it doesn't exist
  if (!dir.exists(out_folder)) dir.create(out_folder, recursive = TRUE)
  # Save as RDS
  saveRDS(result, file = file.path(out_folder, paste0(NAME, ".rds")))
  return(result)
}

run_and_save(data = metabolome_positive, FOLDER = "NDvsCDHFD",  NAME ="NDvsCDHFD_p", comparison = "Diet",      top_n = 5, data_type = "log", use_lm = FALSE, XLIM= c(-25,25), YLIM= c(0,35),  method = "untargeted", NDvsCHFD_filter)
run_and_save(data = metabolome_negative, FOLDER = "NDvsCDHFD",  NAME ="NDvsCDHFD_n", comparison = "Diet",      top_n = 5, data_type = "log", use_lm = FALSE, XLIM= c(-10,10), YLIM= c(0,5),   method = "untargeted", NDvsCHFD_filter)
run_and_save(data = metabolome_positive, FOLDER = "CDHFD_a",    NAME ="CDHFD_a_p",  comparison = "Treatment", top_n = 5, data_type = "log", use_lm = TRUE,  XLIM= c(-5,5),   YLIM= c(0,6.5), method = "untargeted", CDHFD_a_filters)
run_and_save(data = metabolome_negative, FOLDER = "CDHFD_a",    NAME ="CDHFD_a_n",  comparison = "Treatment", top_n = 5, data_type = "log", use_lm = TRUE,  XLIM= c(-5,5),   YLIM= c(0,6.5), method = "untargeted", CDHFD_a_filters)
run_and_save(data = metabolome_positive, FOLDER = "ND_f",       NAME ="ND_f_p",     comparison = "Treatment", top_n = 5, data_type = "log", use_lm = FALSE, XLIM= c(-11,11), YLIM= c(0,6.5), method = "untargeted", ND_f_filters)
run_and_save(data = metabolome_negative, FOLDER = "ND_f",       NAME ="ND_f_n",     comparison = "Treatment", top_n = 5, data_type = "log", use_lm = FALSE, XLIM= c(-5,5),   YLIM= c(0,6.5), method = "untargeted", ND_f_filters)
run_and_save(data = metabolome_positive, FOLDER = "CDHFD_m",    NAME ="CDHFD_m_p",  comparison = "Treatment", top_n = 5, data_type = "log", use_lm = FALSE, XLIM= c(-5,5),   YLIM= c(0,6.5), method = "untargeted", CDHFD_m_filters)
run_and_save(data = metabolome_negative, FOLDER = "CDHFD_m",    NAME ="CDHFD_m_n",  comparison = "Treatment", top_n = 5, data_type = "log", use_lm = FALSE, XLIM= c(-5,5),   YLIM= c(0,6.5), method = "untargeted", CDHFD_m_filters)
run_and_save(data = metabolome_positive, FOLDER = "CDHFD_f",    NAME ="CDHFD_f_p",  comparison = "Treatment", top_n = 5, data_type = "log", use_lm = FALSE, XLIM= c(-5,5),   YLIM= c(0,6.5), method = "untargeted", CDHFD_f_filters)
run_and_save(data = metabolome_negative ,FOLDER = "CDHFD_f",    NAME ="CDHFD_f_n",  comparison = "Treatment", top_n = 5, data_type = "log", use_lm = FALSE, XLIM= c(-5,5),   YLIM= c(0,6.5), method = "untargeted", CDHFD_f_filters)
### resulting trend and sig metabolites were saved into csv files. tried matching names in metaboanalyst but did not work completely, after various attempt of aoutmated matching, started to google metabolite names. Most of them I found in HMDB and filled in information (HMBD;KEGG,PubChem etc). Some were not found there thatn I used wahtever database i found it in.
run_and_save(data = metabolome_targeted, FOLDER = "CDHFD_a",  NAME ="CDHFD_a",  comparison = "Treatment",  top_n = 5, data_type = "log", use_lm = TRUE,  XLIM= c(-5,5),   YLIM= c(0,6.5), method = "targeted", CDHFD_a_filters)
run_and_save(data = metabolome_targeted, FOLDER = "CDHFD_m",  NAME ="CDHFD_m",  comparison = "Treatment",  top_n = 5, data_type = "log", use_lm = FALSE, XLIM= c(-5,5),   YLIM= c(0,6.5), method = "targeted", CDHFD_m_filters)
run_and_save(data = metabolome_targeted, FOLDER = "CDHFD_f",  NAME ="CDHFD_f",  comparison = "Treatment",  top_n = 5, data_type = "log", use_lm = FALSE, XLIM= c(-5,5),   YLIM= c(0,6.5), method = "targeted", CDHFD_f_filters)
run_and_save(data = metabolome_targeted, FOLDER = "ND_f",     NAME ="ND_f",     comparison = "Treatment",  top_n = 5, data_type = "log", use_lm = FALSE, XLIM= c(-5,5),   YLIM= c(0,6.5), method = "targeted", ND_f_filters)
run_and_save(data = metabolome_targeted, FOLDER = "NDvsCDHFD",NAME ="NDvsCDHFD",comparison = "Diet",       top_n = 5, data_type = "log", use_lm = FALSE, XLIM= c(-15,15), YLIM= c(0,12),  method = "targeted", NDvsCHFD_filter)

rm(list = ls())
gc()
