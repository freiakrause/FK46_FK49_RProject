gc()
rm(list = ls())
home <- normalizePath("~") # bc Windows does not start inuser dir but in user/documents dir
parent <- dirname(home)
setwd(paste0(parent,"/OneDrive - Universität Salzburg/AG_Tumorimmunologie - Freia Krause/01_Experiments/FK46_FK49_RProject"))

library(dplyr)
library(tidyr)
library(ggplot2)
library('corrr')
library(ggcorrplot)
library("FactoMineR")
library(factoextra)
library(ggrepel)
library(emmeans)
source("FK49_Definitions.R")
home <- normalizePath("~") # bc Windows does not start inuser dir but in user/documents dir
parent <- dirname(home) # gives users dir works for me since i aways have the data from work in the onedrive at the same location and only local home dir changes
ExpId = "FK49"

#BA contains all animal but not 199 since it was massive outlier. technical problem, quality problem or Cholemia
#BA filtered does not contain 199, 186 and 187 since 186 and 187 were also a little bnit ouliters but not sure if this shows Cholemia in them. So I  want to analys both datasets and disucc with PI/Peers which to include.
BA         <-readRDS(file= paste0(parent,"/OneDrive - Universität Salzburg/AG_Tumorimmunologie - Freia Krause/01_Experiments/FK49_CD-HFD_13wks/FK49_Analysis/01_RawData/FK49_BA_preprocessed.rds"))
BA_filtered<-readRDS(file= paste0(parent,"/OneDrive - Universität Salzburg/AG_Tumorimmunologie - Freia Krause/01_Experiments/FK49_CD-HFD_13wks/FK49_Analysis/01_RawData/FK49_BAfiltered_preprocessed.rds"))

#Define functions important for Stats and Viz ---

plot_PCA<- function(data,
                    plot_variables= "Diet",
                    Sex_filter= NULL,
                    Diet_filter= NULL,
                    ExpID_filter= NULL,
                    filters=list(),
                    NAME= "",
                    method = NULL,
                    FOLDER= NULL){
  
  df<-subset_data(data,
                  Sex_filter =Sex_filter,
                  Diet_filter = Diet_filter,
                  ExpID_filter = ExpID_filter,
                  filters=filters,
                  method = method)
  
  
  
  PATH = BApwd
  
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
  write.csv(top_contributors, file = paste0(PATH,"/",FOLDER,"/Stats/PC1_PC2_Top10_",NAME,".csv"), row.names = FALSE)
  cat("Top contributors saved to CSV:", paste0(PATH,"/",FOLDER,"/Stats/PC1_PC2_Top10_",NAME,".csv"), "\n")
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
    } else if (v == "Timepoint") {
      color_map <- c("red","blue")
    } else if (v == "Time_Treat") {
      color_map <- c("red","blue","purple","royalblue")
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
    ggsave(plot = plottiplot, width= 12, height= 12, dpi=300, 
           filename=paste0("PCA_", NAME,"_",v, ".png"),
           path = paste0(PATH,"/",FOLDER ) )
    
  }
  gc()
}  
perform_ttest <- function(data,
                          comparison = "Treatment",
                          Sex_filter= NULL,
                          Diet_filter= NULL,
                          ExpID_filter= NULL,
                          filters=list(),
                          top_n = 8,
                          data_type = c("log", "raw", "norm", "scaled"),
                          use_lm = FALSE,
                          XLIM=c(-9,9),
                          YLIM=c(0,7.5),
                          NAME= "",
                          method = NULL,
                          FOLDER = NULL
){
  # Load libraries
  library(dplyr)
  library(ggplot2)
  library(ggrepel)
  library(tidyr)
  library(pheatmap)
  library(lme4)
  library(lmerTest)
  library(emmeans)
  
  data_type <- match.arg(data_type)
  cat("Data type received:", data_type, "\n")
  
  
  PATH = BApwd
  
  
  # Subset data
  df<-subset_data(data,
                  Sex_filter =Sex_filter,
                  Diet_filter = Diet_filter,
                  ExpID_filter = ExpID_filter,
                  filters=filters,
                  method = method)
  
  metab_names <- colnames(df$log_values)
  
  selected_values <- switch(data_type,
                            raw = df$raw_values,
                            log = df$log_values,
                            norm = df$norm_values,
                            scaled = df$scaled_log_values)
  
  # Determine which comparison to plot (first if multiple)
  plot_comp <- comparison[1]
  
  # Color mapping for plotting
  color_map <- switch(plot_comp,
                      Sex = Sex_colors,
                      Treatment = Treatment_colors,
                      T_D_S = T_D_S_colors,
                      T_S = T_S_colors,
                      Diet = Diet_colors,
                      Timepoint = c("red","blue"),
                      Time_Treat = c("red","blue","purple","royalblue"),
                      stop("Unknown comparison for color mapping."))
  
  # ---- Linear model branch ----
  if(use_lm){
    
    lm_metab <- function(metab_name){
      
      values <- selected_values[[metab_name]]
      
      df_metab <- data.frame(
        value = values,
        Treatment = factor(df$metadata$Treatment,levels=c("TAM","Ctrl")),
        Sex = as.factor(df$metadata$Sex),
        Timepoint = as.factor(df$metadata$Timepoint),
        Animal = as.factor(df$metadata$Animal),
        Diet = as.factor(df$metadata$Diet)
      )
      
      # keep only animals with repeated measures
      df_metab <- df_metab %>%
        dplyr::group_by(Animal) %>%
        dplyr::filter(dplyr::n_distinct(Timepoint) > 1) %>%
        dplyr::ungroup()
        
      print(df_metab)
      fit <- lmer(value ~ Treatment * Timepoint * Sex + (1 | Animal),
                  data = df_metab)
      
      anova_res <- as.data.frame(anova(fit))
      
      # =========================================================
      # CONTRAST 1: Treatment effect at each Timepoint (use T2)
      # =========================================================
      emm_trt <- emmeans(fit, ~ Treatment | Timepoint)
      contr_treat <- as.data.frame(contrast(emm_trt, "pairwise"))

      # =========================================================
      # CONTRAST 2: Trajectory effect (difference-in-differences)
      # =========================================================
      emm_traj <- emmeans(fit, ~ Timepoint * Treatment)
      contr_traj <- as.data.frame(
        contrast(emm_traj, interaction = "pairwise")
      )
      
      # =========================================================
      # CONTRAST 3: Sex modifies treatment (at T2)
      # =========================================================
      emm_sex <- emmeans(fit, ~ Treatment | Timepoint * Sex)
      contr_sex <- as.data.frame(contrast(emm_sex, "pairwise"))
      
      # keep only T2 for sex effect
      contr_sex <- subset(contr_sex, Timepoint == levels(df_metab$Timepoint)[2])
      
      # =========================================================
      # EXTRACT KEY VALUES SAFELY
      # =========================================================
      tp2 <- levels(df_metab$Timepoint)[2]
      tp1 <- levels(df_metab$Timepoint)[1]
      
      get_val <- function(df, cond){
        out <- df %>% dplyr::filter(!!rlang::parse_expr(cond))
        if(nrow(out) == 0) return(c(NA, NA))
        c(out$estimate[1], out$p.value[1])
      }
      
      # Treatment at T2
      c_treat_tp2 <- get_val(contr_treat, paste0("Timepoint=='", tp2, "'"))
      c_treat_tp1 <- get_val(contr_treat, paste0("Timepoint=='", tp1, "'"))
      write.csv2(c_treat_tp2, file = paste0(PATH,"/",FOLDER,"/Stats/",NAME,"_LM_resultsc.csv"), row.names = FALSE)
      
      # Trajectory effects
      # (emmeans interaction output already structured)
      traj_est <- contr_traj$estimate[1]
      traj_p   <- contr_traj$p.value[1]
      
      # Sex effect at T2 (collapsed across sexes within treatment)
      c_sex_tp2 <- get_val(contr_sex, paste0("Timepoint=='", tp2, "'"))
      
      # =========================================================
      # OUTPUT TABLE
      # =========================================================
      data.frame(
        Metabolite = metab_name,
        # ANOVA (global effects)
        p_Treatment = anova_res["Treatment", "Pr(>F)"],
        p_Sex = anova_res["Sex", "Pr(>F)"],
        p_Timepoint = anova_res["Timepoint", "Pr(>F)"],
        p_Treat_Time = anova_res["Treatment:Timepoint", "Pr(>F)"],
        p_3way = anova_res["Treatment:Timepoint:Sex", "Pr(>F)"],
        
        # CONTRAST 1: endpoint effect
        logFC_Treat_T2 = c_treat_tp2[1],
        p_Treat_T2 = c_treat_tp2[2],
        logFC_Treat_T1 = c_treat_tp1[1],
        p_Treat_T1 = c_treat_tp1[2],
        
        # CONTRAST 2: trajectory (MAIN RESULT)
        logFC_Trajectory = traj_est,
        p_Trajectory = traj_p,
        
        # CONTRAST 3: sex modulation at T2
        logFC_Sex_T2 = c_sex_tp2[1],
        p_Sex_T2 = c_sex_tp2[2]
        
        
      )
    }
    
    results <- do.call(rbind, lapply(metab_names, lm_metab))
    
    results <- results %>%
      mutate(
        ANOVA_p_treat = p_Treatment,
        adj_p_Treat_T2 = p.adjust(p_Treat_T2, method = "fdr"),
        adj_p_Treat_T1 = p.adjust(p_Treat_T1, method = "fdr"),
        adj_p_Trajectory = p.adjust(p_Trajectory, method = "fdr"),
        adj_p_Sex_T2 = p.adjust(p_Sex_T2, method = "fdr")
      ) %>%
      mutate(
        p.value = p_Treat_T2,
        log2FC = logFC_Treat_T2,
        adj.p.value = adj_p_Treat_T2,
        negLog10FDR = -log10(adj.p.value),
        significant = adj.p.value < 0.05 & abs(log2FC) > 0.5,
        trend = adj.p.value < 0.1 & abs(log2FC) > 0.3,
        direction=case_when(adj_p_Treat_T2 < 0.05 & logFC_Treat_T2 > 0 ~ "UP",
                            adj_p_Treat_T2 < 0.05 & logFC_Treat_T2 < 0 ~ "DOWN",
                            TRUE             ~"NS")
        
        )
    
    all_results_saving<-results%>% 
      arrange(log2FC) %>%
      dplyr::select(
        Metabolite, 
        ANOVA_p_treat,
        adj_p_Treat_T2,
        adj_p_Treat_T1,
        logFC_Treat_T1,
        adj_p_Trajectory,
        adj_p_Sex_T2,
        p.value, 
        adj.p.value, 
        log2FC, 
        negLog10FDR,
        direction,
        significant, 
        trend      )
    write.csv2(all_results_saving, file = paste0(PATH,"/",FOLDER,"/Stats/",NAME,"_LM_results.csv"), row.names = FALSE)
    
    # After computing the results and determining the trend metabolites
    sig_results <- results %>% 
      filter(significant == TRUE) %>%
      arrange(log2FC) %>%
      dplyr::select(
        Metabolite, 
        p.value, 
        adj.p.value, 
       log2FC,
       negLog10FDR,
       direction,
        significant, 
        trend
      )
    
    # Save the trend metabolites data to a CSV
    write.csv2(sig_results, file = paste0(PATH,"/",FOLDER,"/Stats/",NAME, "_sig_BA_LM.csv"), row.names = FALSE)
    cat("Saved sig metabolites to ", paste0(PATH,"/",FOLDER,"/Stats/",NAME, "_sig_BA_LM.csv"), "\n")
    
    trend_results <- results %>% 
      filter(trend == TRUE) %>%
      arrange(log2FC) %>%
      dplyr::select(
        Metabolite, 
        p.value, 
        adj.p.value, 
        log2FC,
        negLog10FDR,
        direction,
        significant, 
        trend
      )
    # Save the trend metabolites data to a CSV
    write.csv2(trend_results, file = paste0(PATH,"/",FOLDER,"/Stats/",NAME, "_trend_BA_LM.csv"), row.names = FALSE)
    
    cat("Saved trend metabolites to ", paste0(PATH,"/",FOLDER,"/Stats/",NAME, "_trend_BA_LM.csv"), "\n")
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
    write.csv2(all_results_saving, file = paste0(PATH,"/",FOLDER,"/Stats/",NAME,"_test_results.csv"), row.names = FALSE)
    
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
    write.csv2(sig_results, file = paste0(PATH,"/",FOLDER,"/Stats/",NAME, "_sig_metabolites_ttest.csv"), row.names = FALSE)
    cat("Saved sig metabolites to ", paste0(PATH,"/",FOLDER,"/Stats/",NAME, "_sig_metabolites_ttest.csv"), "\n")
    
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
    write.csv2(trend_results, file = paste0(PATH,"/",FOLDER,"/Stats/",NAME, "_trend_metabolites_ttest.csv"), row.names = FALSE)
    
    cat("Saved trend metabolites to ", paste0(PATH,"/",FOLDER,"/Stats/",NAME, "_trend_metabolites_ttest.csv"), "\n")
  }
  # ---- Dynamic plotting ----
  df_values <- if(method == "untargeted") df$norm_values else df$raw_values
  y_label_method <- if(method == "untargeted") "Relative abundance [%]" else "Concentration [umol/L]"
  
  top_metabolites <- head(results %>% arrange(adj.p.value) %>% pull(Metabolite), top_n)
  df_plot <- cbind(df$metadata, df_values[, top_metabolites, drop=FALSE])
  print(colnames(results))
  print(colnames(df_plot))

  df_long <- df_plot %>%
    pivot_longer(cols = all_of(top_metabolites), names_to = "Metabolite", values_to = "value") %>%
    left_join(results[, c("Metabolite", "adj.p.value","adj_p_Trajectory","adj_p_Treat_T1","adj_p_Sex_T2","direction")], by="Metabolite")
  
  df_long$Metabolite<-factor(df_long$Metabolite, levels = BA_sort)
  saveRDS(df_long,file =paste0(PATH,"/",FOLDER,"/",NAME,".rds" ))
}
 

# Helper function to run t-test and save result
run_and_save <- function(data, NAME, comparison, top_n = 5, data_type = "log",
                         use_lm = FALSE, XLIM = c(-5,5), YLIM = c(0,5), 
                         method = "untargeted",  filters=list(), FOLDER = NULL) {
  
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
    filters=filters,
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
## Subsets Dataset in subset and performs data preprocessing -----
### normalization, log transform, scaling
BA_a_filter =list(Sex= c("male", "female"),Diet= c("ND","CDHFD13"),ExpID= c("FK49"))
BA_f_filter   = list(Sex= c("female"),Diet = c("ND","CDHFD13"), ExpID = c("FK49"))
BA_m_filter   = list(Sex= c("male"),  Diet = c("ND","CDHFD13"), ExpID= c("FK49"))
BA_tp1_filter = list(Sex= c("male","female"),  Diet = c("ND"), ExpID = c("FK49"))
BA_tp2_filter = list(Sex= c("male","female"),  Diet = c("CDHFD13"), ExpID = c("FK49"))

BA_f_filter2   = list(Sex= c("female"),Diet = c("CDHFD13"), ExpID = c("FK49"))
BA_m_filter2   = list(Sex= c("male"),  Diet = c("CDHFD13"), ExpID = c("FK49"))

sub_data <-do.call(subset_data,c(list(data = BA, method = "targeted"), BA_a_filter))
sub_data<-subset_data(data=BA, method= "targeted", filters=list(Sex=c("female")))
ggsave(plot= sub_data$Preprocessing, width= 12, height= 6, dpi=300,
       filename="BG_prepocessing_BA_a.png", path =paste0(BApwd,"/BA_a") )


# GDCA und GCDCA haben gesamt >75% LOD, dabei 100%LOD in T2, >50% T1 für Stats din, für PCA raus
BA_PCA<-BA_filtered%>%select(-GDCA,-GCDCA)
plot_PCA(data = BA_PCA, plot_variables = c("Sex","Treatment","Timepoint"),
         FOLDER = "BA_a", NAME= "TP-1_11_BA_a_filtered", method="targeted" )

plot_PCA(data = BA_PCA,  plot_variables = c("Sex", "Treatment"),  filters = list(Timepoint = "-1"),
         NAME = "TP-1_filtered",  FOLDER = "BA_a",  method= "targeted")

plot_PCA( data = BA_PCA,plot_variables = c("Sex", "Treatment"),
          filters = list(Timepoint = "11"), NAME = "TP11_filtered",FOLDER = "BA_a",method= "targeted")

plot_PCA( data = BA_PCA,plot_variables = c("Treatment"),  filters = list(Timepoint = "11",Sex= "male"),
          NAME = "TP11_m_filtered",  FOLDER = "BA_a",  method= "targeted")

plot_PCA(  data = BA_PCA,  plot_variables = c("Treatment"),  filters = list(Timepoint = "11",Sex= "female"),
           NAME = "TP11_f_filtered",  FOLDER = "BA_a",  method= "targeted")


#29.04.26 
# Right now I deceided the best analsisy wold be: 
#complete data set (only 199 excluded bc xx fold more BA than all other
#unsure if technical erroer or cholemia
#linear model: lmer(value ~ Treatment * Timepoint * Sex + (1 | Animal)
#Does TAM change BA at TP2?
#emm_trt <- emmeans(fit, ~ Treatment | Timepoint) use TP2
#Does TAM change trajectory from non treated to CDHFD differntly than Ctrl?
#emm_traj <- emmeans(fit, ~ Timepoint * Treatment)
#Does Sex modifie treatment at t2
#emm_sex <- emmeans(fit, ~ Treatment | Timepoint * Sex)
#Heatmap prints all included animals at all timepoints
#Violin/Bar plots prints included animals at tp2 (bc thats the number 1)
#Trajector is not change for any BA, need to think about how to print
# Sex modifies few BA, need to think about how to print
# maybe need to change function (make it smaller to only account for these specific questions.)
run_and_save(data = BA_filtered, FOLDER = "BA_a",  NAME ="BA_a_filtered", comparison = c("Treatment"),   
             top_n = 21, data_type = "log", use_lm = TRUE, XLIM= c(-10,10), YLIM= c(0,5), 
             method = "targeted", BA_a_filter)

run_and_save(data = BA_filtered, FOLDER = "BA_f",  NAME ="BA_f_filtered", comparison = c("Treatment"),   
             top_n = 21, data_type = "log", use_lm = FALSE, XLIM= c(-25,25), YLIM= c(0,35), 
             method = "targeted", BA_f_filter2)

run_and_save(data = BA_filtered, FOLDER = "BA_m",  NAME ="BA_m_filtered", comparison = c("Treatment"),   
             top_n = 21, data_type = "log", use_lm = FALSE, XLIM= c(-25,25), YLIM= c(0,35), 
             method = "targeted", BA_m_filter2)

