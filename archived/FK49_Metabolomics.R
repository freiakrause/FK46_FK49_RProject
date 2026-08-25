#on new machine ceck where home directory is
#normalizePath("~")
# [1] "C:\\Users\\b1084855\\Documents"#
# there create plain text file named .Renviron
# Cononte of file:
# RENV_PATHS_LIBRARY=C:/renv/library
# RENV_PATHS_CACHE=C:/renv/cache
# TMPDIR=C:/rtemp
# TEMP=C:/rtemp
# TMP=C:/rtemp
# R will now read this files upont startup. In R once crete the following dirs:
# code in R
# dir.create("C:/renv/library", recursive = TRUE, showWarnings = FALSE)
# dir.create("C:/renv/cache", recursive = TRUE, showWarnings = FALSE)
# dir.create("C:/rtemp", showWarnings = FALSE)
# then restore renv: Now you should have the librarys that are saves as info on ondreive installed on yosr local comoputer
# renv::restore()
# check if paths are correctly read from R:
# Sys.getenv("RENV_PATHS_LIBRARY")
# Sys.getenv("TMPDIR")
# Now renv should saf info on librarys in project ononedrive but librarys should be locally installed. so there should not be issues with lirbary installation due to onedribe
# From https://pmc.ncbi.nlm.nih.gov/articles/PMC9032224/
# 1.Compound Detecion MS or NMR
# 2.Data pre-processing 
# 3.Data processing data normalization and compund indentification 
# 4.Statistical ANalyis  ( I guwess i do this and MetaboNET did the prevoius steps)
# 5.Function Analysis Enrichment, Pathway
# 6.Omic Dasta integration Transpritomis, protemoics microbiome
rm(list = ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library('corrr')
library(ggcorrplot)
library("FactoMineR")
library(dplyr)
library(factoextra)
library(MetaboAnalystR)
home <- normalizePath("~") # bc Windows does not start inuser dir but in user/documents dir
parent <- dirname(home) # goves users dir works for me since i aways have the data from work in the onedrive at the same location and only local home dir changes
ExpId = "FK49"
meta_cols <- c("Animal","Sample","Sex","Treatment","Diet","ExpID","T_D_S","T_D" ,"T_S")# cols containing metadata

if (ExpId=="FK49") {
  metabolome <- read.csv(paste0(parent,"/OneDrive - Universität Salzburg/Freigegebene Dokumente - AG_Tumorimmunologie/Data/Freia Krause/01_Experiments/FK49_CD-HFD_13wks/FK49_HILIC09/Report_M086_HILIC09_20251222.csv"),sep=";")
  setwd("C:/Users/b1084855/OneDrive - Universität Salzburg/Freigegebene Dokumente - AG_Tumorimmunologie/Data/Freia Krause/01_Experiments/FK49_CD-HFD_13wks/FK49_Analysis")
  load(paste0("01_RawData/",ExpId,"_Data_prepared.Rda"))
  }  else if (ExpId == "FK46"){
  print("You dont have data for this Experiment")
} else{
  print("You dont have data for this Experiment")}

metabolome<- metabolome %>%select(-Sample.ID)%>%       #is only the id the facility used for samples
  mutate(T_D_S= paste0(Treatment,"_",Diet,"_",Sex))%>% # i want these combined variable to play around and find out
  mutate(T_D= paste0(Treatment,"_",Diet))%>%
  mutate(T_S= paste0(Treatment,"_",Sex))%>%
  rename(Sample = Sample.No.)%>% 
  mutate(Sample = as.factor(Sample))%>%
  mutate(across(where(is.character),~na_if(., "<LOD")))%>%# data contains "<LOD" this is anoying and intrudces characters 
  mutate(across(!all_of(meta_cols), as.numeric))%>% #is supposed to be change to minimium / 2 in literature need to check if this is good practice
  mutate(across(where(is.numeric), ~ifelse(is.na(.), min(., na.rm = TRUE)/2, .))) # LOD is now na and is changed to min/2


## Define some usefull functions -----
subset_data <- function(data,
                        Sex_filter  = c("male", "female"),
                        Diet_filter = c("CDHFD13", "ND"),
                        ExpID_filter = c("FK49", "BH")
                        ) {
  
  # Subset
  dataset <- data %>%
    filter(
      Sex  %in% Sex_filter,
      Diet %in% Diet_filter,
      ExpID %in% ExpID_filter    )
  
  # Metadata columns
  meta_cols <- c(
    "Sample", "Animal", "Sex", "Treatment",
    "Diet", "ExpID", "T_D_S", "T_D", "T_S")
  
  # Extract metabolite matrix
  numerical_data <- dataset %>%
    select(-all_of(meta_cols)) %>%
    select(where(is.numeric))
  
  # Remove zero-variance features
  nzv <- apply(numerical_data, 2, var, na.rm = TRUE) > 0
  numerical_data <- numerical_data[, nzv, drop = FALSE]
  
  # Optional scaling
  eps <- min(numerical_data[numerical_data > 0], na.rm = TRUE) / 2
  #eps hier noch addieren innerhalb des log()? Macht das einen unterschied?
  data_log <- log2(numerical_data+eps)
  data_scaled <- as.data.frame(scale(data_log))

  
  return(list(
    raw_values  = numerical_data,
    log_values = data_log,
    scaled_log_values = data_scaled,
    metadata   = dataset[, meta_cols, drop = FALSE]
  ))
}
  
plot_PCA<- function(data,
                    plot_variables= "Diet",
                    Sex_filter,
                    Diet_filter,
                    ExpID_filter){
  
  df<-subset_data(data,
              Sex_filter =Sex_filter,
              Diet_filter = Diet_filter,
              ExpID_filter = ExpID_filter)

  

  pca <- prcomp(df$scaled_log_values, center = FALSE, scale. = FALSE)
  
  # Summary of PCA
  print(summary(pca))
  
  # Visualizations
  eigen<-fviz_eig(pca, addlabels = TRUE)
  print(eigen)
  
  #variable<-fviz_pca_var(pca, col.var = "black")
  #print(variable)
  
  variable_color<-fviz_pca_var(pca, col.var = "cos2",
                               alpha.var="contrib",
                               select.var = list(cos2 = 0.96),
               #gradient.cols = c("black", "orange", "green"),
               repel = TRUE)
  print(variable_color)
  
  cos<-fviz_cos2(pca, choice = "var", axes = 1:2)
  print(cos)
  
  for (v in plot_variables){
    
    plottiplot<-fviz_pca_ind(pca, habillage = df$metadata[[v]],addEllipses = TRUE,
                             #ellipse.level=0.95,
                             geom="point",
                             pointsize = 4)+
                labs(title ="PCA", x = "PC1", y = "PC2")+
                scale_color_brewer(palette="Dark2") +
                theme_minimal()
    print(plottiplot)
    
  }
}  
plot_PCA(data = metabolome, 
         Diet_filter = c("CDHFD13", "ND"),
        Sex_filter  = c("male", "female"),
        ExpID_filter = c("FK49", "BH"),
        plot_variables = c("Sex","Treatment","Diet", "T_D_S","T_D","T_S", "ExpID"))

plot_PCA(data = metabolome, Diet_filter = c("CDHFD13"),     Sex_filter= c("male","female"), ExpID_filter = c("FK49"), plot_variables = c("Sex","Treatment", "T_D_S","T_D","T_S"))
plot_PCA(data = metabolome,  Diet_filter = c("ND","CDHFD13"),Sex_filter  = c("female"),   ExpID_filter = c("FK49","BH"),plot_variables = c( "Diet"), )

plot_PCA(data = metabolome, Sex_filter= c("male"),Diet_filter = c("CDHFD13"),plot_variables = c("Treatment", "T_D_S"))
plot_PCA(data = metabolome, Sex_filter= c("female"),Diet_filter = c("CDHFD13"),plot_variables = c("Treatment", "T_D_S"))
  
perform_ttest <- function(data,
                          comparison = "Treatment",
                          Sex_filter,
                          Diet_filter,
                          ExpID_filter,
                          top_n = 8,
                          data_type = c("log", "raw"),
                          use_lm = FALSE) {
  
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(emmeans)
  
  data_type <- match.arg(data_type)
  
  # Subset the data
  df <- subset_data(data, Sex_filter = Sex_filter, Diet_filter = Diet_filter, ExpID_filter = ExpID_filter)
  metab_names <- colnames(df$log_values)
  
  if(use_lm){
    lm_metab <- function(metab_name){
      values <- if(data_type == "log") df$log_values[[metab_name]] else df$raw_values[[metab_name]]
      df_metab <- data.frame(
        value = values,
        Treatment = as.factor(df$metadata[[comparison]]),
        Sex = as.factor(df$metadata$Sex)
      )
      
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
        Intercept = coefs["(Intercept)", "Estimate"],
        Treatment = coefs["TreatmentTAM", "Estimate"],
        Sex = coefs["Sexmale", "Estimate"],
        Treatment_Sex = coefs["TreatmentTAM:Sexmale", "Estimate"],
        p_Intercept = coefs["(Intercept)", "Pr(>|t|)"],
        p_Treatment = coefs["TreatmentTAM", "Pr(>|t|)"],
        p_Sex = coefs["Sexmale", "Pr(>|t|)"],
        p_Treatment_Sex = coefs["TreatmentTAM:Sexmale", "Pr(>|t|)"],
        log2FC_overall = log2FC_overall,
        log2FC_female = log2FC_female,
        log2FC_male = log2FC_male,
        stringsAsFactors = FALSE
      )
      
      return(list(results = res, df_metab = df_metab, fit = fit))
    }
    
    # Run LM for all metabolites
    lm_list <- lapply(metab_names, lm_metab)
    results <- do.call(rbind, lapply(lm_list, `[[`, "results"))
    
    # Adjust p-values
    results$adj.p_Treatment <- p.adjust(results$p_Treatment, method="fdr")
    
    # Select top metabolites
    top_metabolites <- head(results$Metabolite[order(results$adj.p_Treatment)], top_n)
    
    # ==== LM Visualizations for first top metabolite ====
    # Coefficients bar plot
    first_fit <- lm_list[[which(metab_names == top_metabolites[1])]]$fit
    coefs_df <- as.data.frame(summary(first_fit)$coefficients)
    coefs_df$Term <- rownames(coefs_df)
    
    p_coefs <- ggplot(coefs_df[2:4,], aes(x=Term, y=Estimate)) +
      geom_bar(stat="identity", fill="steelblue") +
      geom_errorbar(aes(ymin=Estimate-`Std. Error`, ymax=Estimate+`Std. Error`), width=0.2) +
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
      values <- if(data_type == "log") df$log_values[[metab_name]] else df$raw_values[[metab_name]]
      df_metab <- data.frame(
        value = values,
        group = as.factor(df$metadata[[comparison]])
      )
      
      t_res <- t.test(value ~ group, data = df_metab)
      log2FC <- mean(df_metab$value[df_metab$group == levels(df_metab$group)[2]]) -
        mean(df_metab$value[df_metab$group == levels(df_metab$group)[1]])
      
      data.frame(Metabolite = metab_name, p.value = t_res$p.value, log2FC = log2FC)
    }
    
    results <- lapply(metab_names, ttest_metab) %>% do.call(rbind, .)
    results$adj.p.value <- p.adjust(results$p.value, method="fdr")
    
    top_metabolites <- head(results$Metabolite[order(results$adj.p.value)], top_n)
  }
  
  # ==== Violin / Boxplots for Top Metabolites ====
  df_values <- if(data_type=="log") df$log_values else df$raw_values
  df_plot <- cbind(df$metadata, df_values[, top_metabolites, drop=FALSE])
  df_long <- df_plot %>%
    pivot_longer(cols = all_of(top_metabolites), names_to = "Metabolite", values_to = "value")
  
  if(use_lm){
    df_long <- df_long %>%
      left_join(results[, c("Metabolite","adj.p_Treatment")], by="Metabolite")
    p_violin <- ggplot(df_long, aes(x=.data[[comparison]], y=value, fill=.data[[comparison]])) +
      geom_violin(trim=FALSE, alpha=0.5) +
      geom_boxplot(width=0.1, outlier.shape=NA) +
      geom_jitter(width=0.1, size=1.5, alpha=0.7) +
      facet_wrap(~Metabolite, scales="free_y") +
      geom_text(aes(x=1.5, y=max(value, na.rm=TRUE),
                    label=paste0("adj.p = ", signif(adj.p_Treatment,3))),
                inherit.aes=FALSE, size=3, vjust=-0.5) +
      labs(y=ifelse(data_type=="log", expression("log"[2]*" intensity"), "Concentration"),
           x=comparison,
           title=paste0("Top ", top_n, " Metabolites (LM)")) +
      theme_minimal() +
      scale_fill_brewer(palette="Set2")
    
  } else {
    df_long <- df_long %>%
      left_join(results[, c("Metabolite","adj.p.value")], by="Metabolite")
    p_violin <- ggplot(df_long, aes(x=.data[[comparison]], y=value, fill=.data[[comparison]])) +
      geom_violin(trim=FALSE, alpha=0.5) +
      geom_boxplot(width=0.1, outlier.shape=NA) +
      geom_jitter(width=0.1, size=1.5, alpha=0.7) +
      facet_wrap(~Metabolite, scales="free_y") +
      geom_text(aes(x=1.5, y=max(value, na.rm=TRUE),
                    label=paste0("adj.p = ", signif(adj.p.value,3))),
                inherit.aes=FALSE, size=3, vjust=-0.5) +
      labs(y=ifelse(data_type=="log", expression("log"[2]*" intensity"), "Concentration"),
           x=comparison,
           title=paste0("Top ", top_n, " Metabolites (T-Test)")) +
      theme_minimal() +
      scale_fill_brewer(palette="Set2")
  }
  
  print(p_violin)
  
  return(results)
}


CDHFD_all <- perform_ttest(  data = metabolome,  comparison = "Treatment", 
                             Sex_filter  = c("female","male"),  Diet_filter = c("CDHFD13"),  
                             ExpID_filter = c("FK49"),  top_n = 5,  data_type = "log",
                             use_lm = TRUE)
ND_females <- perform_ttest( data = metabolome,comparison = "Treatment", 
                             Sex_filter  = c("female"),
                             Diet_filter = c("ND"),ExpID_filter = c("BH"), 
                             top_n = 5,  data_type = "log",use_lm = FALSE)
CDHFD_males <- perform_ttest(  data = metabolome,  comparison = "Treatment",  
                               Sex_filter  = c("male"),  Diet_filter = c("CDHFD13"),  
                               ExpID_filter = c("FK49"),  top_n = 5,  data_type = "log",use_lm = FALSE)
CDHFD_females <- perform_ttest(  data = metabolome,  comparison = "Treatment",  
                                 Sex_filter  = c("female"),  Diet_filter = c("CDHFD13"),  
                                 ExpID_filter = c("FK49"),  top_n = 5,  data_type = "log",use_lm = FALSE)




plot_violins <- function(data,
                         metabolite = "Glucose.6.phosphate",
                         comparison = "Treatment",
                         fill = "Treatment",
                         Sex_filter,
                         Diet_filter,
                         ExpID_filter,
                         data_type = c("log", "raw")) {
  
  data_type <- match.arg(data_type)  # wählt zwischen "log" und "raw"
  
  # Subset data
  df <- subset_data(data,
                    Sex_filter  = Sex_filter,
                    Diet_filter = Diet_filter,
                    ExpID_filter = ExpID_filter)
  
  # Wähle die Datenquelle
  values <- if (data_type == "log") {
    df$log_values[[metabolite]]
  } else {
    df$raw_values[[metabolite]]
  }
  
  # Erstelle Plot-Datenframe
  df_plot <- data.frame(
    value = values,
    group = df$metadata[[comparison]],
    fill_col = df$metadata[[fill]]
  )
  
  # Plot
  ggplot(df_plot, aes(x = group, y = value, fill = fill_col)) +
    geom_violin(trim = FALSE, alpha = 0.5) +
    geom_boxplot(width = 0.1, outlier.shape = NA) +
    geom_jitter(width = 0.1, size = 1.5, alpha = 0.7) +
    labs(title = paste("Metabolite:", metabolite, "(", data_type, "values )"),
         x = comparison,
         y = ifelse(data_type == "log",
                    expression("log"[2]*" intensity"),
                    "pM/mg")) +
    theme_minimal() +
    scale_fill_brewer(palette = "Set2")
}
plot_violins(
  data = metabolome,
  metabolite = "Guanosine",
  comparison = "Diet",
  fill = "Diet",
  Sex_filter  = c("female"),
  Diet_filter = c("ND","CDHFD13"),
  ExpID_filter = c("FK49","BH"),
  data_type = "log"
)
