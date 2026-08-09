rm(list=ls())
gc()
library(lme4)
library(lmerTest)
library(emmeans)
library(dplyr)
library(ggplot2)
source("FK49_Definitions.R")

# Read Dataframe to analyse ------
Mito_level <- readRDS(file.path(PATHS$TEM$output,"CleanData/Mito_Dataframe.rds"))%>%
  filter(Cell_Type != "nonHepatocyte")
Cell_level <- readRDS(file.path(PATHS$TEM$output,"CleanData/Cell_Dataframe.rds"))%>%
  filter(Cell_Type != "nonHepatocyte") 




# QC ob meine Analyse sich im laufe der zeit verändert hat
## for values on Mito level
Mito_values <- Mito_level %>%
  select(-Sex, -Treatment, -Animal, -Image_ID, -Mito_ID,-Cell_ID, -Pedigree_ID,
         -Shape,-Perimeter, -Analysis_ID,-Cell_Type) %>%
  colnames()
Mito_level <- Mito_level %>%mutate(Analysis_ID = as.numeric(as.character(Analysis_ID)))

for(value in Mito_values){
  model <- lm(reformulate("Analysis_ID", response = value), data = Mito_level)
  slope <- coef(model)["Analysis_ID"]
  r2    <- summary(model)$r.squared
  drift <- slope * (max(Mito_level$Analysis_ID, na.rm = TRUE) -
                    min(Mito_level$Analysis_ID, na.rm = TRUE))
  
  plot <- ggplot(Mito_level,aes(x = Analysis_ID, y = .data[[value]])) +
    geom_point(size = 0.5, alpha = 0.1) +
    geom_smooth(method = "lm", se = TRUE) +
    labs(title = value,
         subtitle = paste0(
        "Slope = ", signif(slope, 3),
        " | Drift = ", signif(drift, 3),
        " | R² = ", round(r2, 3)))
  ggsave(plot=plot, filename = paste0("2_QC_AnalysisID_",value,".png"),
  path= file.path(PATHS$TEM$output,"/Background/"),
  dpi= 300, width=7, height=4 )
  print(plot)
}


Cell_values <- Cell_level %>%
  select(-Sex, -Treatment, -Animal, -Image_ID, -Cell_ID, -Cell_Type,-Analysis_ID,-Cell_Area,-n_Mito,-Freq_Donut,-Freq_LostCristae,-Total_Mito_Area,
         -starts_with("SD_") ,-starts_with("Mean")) %>% #nimmt die werte raus, die auf mito ebene schongetestet wurden
  colnames()

Cell_level <- Cell_level %>%mutate(Analysis_ID = as.numeric(as.character(Analysis_ID)))
# For values on cell level
for(value in Cell_values){
  model <- lm(reformulate("Analysis_ID", response = value), data = Cell_level)
  slope <- coef(model)["Analysis_ID"]
  r2    <- summary(model)$r.squared
  drift <- slope * (max(Cell_level$Analysis_ID, na.rm = TRUE) -
                    min(Cell_level$Analysis_ID, na.rm = TRUE))
  
  plot <- ggplot(Cell_level,aes(x = Analysis_ID, y = .data[[value]])) +
    geom_point(size = 0.5, alpha = 0.1) +
    geom_smooth(method = "lm", se = TRUE) +
    labs(title = value,
         subtitle = paste0(
           "Slope = ", signif(slope, 3),
           " | Drift = ", signif(drift, 3),
           " | R² = ", round(r2, 3)))
  ggsave(plot=plot, filename = paste0("2_QC_AnalysisID_",value,".png"),
         path= file.path(PATHS$TEM$output,"/Background/"),
         dpi= 300, width=7, height=4 )
  print(plot)
}


# I am not exactly sure what to make of these slopes. 
# It looks not like massive changes. Si i will ignore analysis id. Still, some paramters drift. 
# But i dont feel like i affected them thorugh my leraning

# Statistics on Parameters directly derived from Mitos -----

 ## lme instead of lm because it takes into account that mitos from same cells are more similar 
 ## thatn mitos from not same cell

Mito_models <- list()
Mito_results<-data.frame()

for (value in Mito_values) {
  formula <- as.formula(paste0(value," ~Treatment * Sex + (1|Animal) + (1|Animal:Image_ID) + (1|Animal:Image_ID:Cell_ID)"))
  model <- lmer(formula, data = Mito_level)
  Mito_models[[value]] <- model
  # Treatment effect separated by Sex
  emm_TS <- emmeans(model, pairwise ~ Treatment | Sex,adjust = "none",pbkrtest.limit = 10000)
  contrast_table_TS <- as.data.frame(emm_TS$contrasts)
  # Treatment effect 
  emm_T <- emmeans(model, pairwise ~ Treatment,adjust = "none",pbkrtest.limit = 10000) #pbkrtest da durch die vielen mitos system automatisch eine freiheitsgrad berechnung runter macht da sons viel rechenleistung. ich will aber die mitos drin haben
  contrast_table_T <- as.data.frame(emm_T$contrasts)
  
  # Interaction
  interaction <- summary(model)$coefficients["TreatmentTAM:Sexmale",]
  
  Mito_results <- rbind(
    Mito_results,
    data.frame(
      Variable = value,
      Treatment_Effect =  contrast_table_T$estimate[contrast_table_T$contrast=="EtOH - TAM"],
      Treatment_p =  contrast_table_T$p.value[contrast_table_T$contrast=="EtOH - TAM"],
      Female_Treatment_Effect =  contrast_table_TS$estimate[contrast_table_TS$Sex=="female"],
      Female_p = contrast_table_TS$p.value[contrast_table_TS$Sex=="female"],
      Male_Treatment_Effect =  contrast_table_TS$estimate[contrast_table_TS$Sex=="male"],
      Male_p =contrast_table_TS$p.value[contrast_table_TS$Sex=="male"],
      Sex_Treatment_Interaction =interaction["Estimate"],
      Sex_Treatment_Interaction_p =  interaction["Pr(>|t|)"],
      row.names = NULL
    )
  )
  
}

write.csv2(Mito_results,
           file=file.path(PATHS$TEM$output,"Statistics/01_Mito_level_Paramters.csv"),
           row.names = FALSE)

# Statistics on Paramters derived from Cell Lvel -----
## lme instead of lm because it takes into account that mitos from same cells are more similar 
## thatn mitos from not same cell


Cell_models <- list()
Cell_results<-data.frame()

for (value in Cell_values) {
  formula <- as.formula(paste0(value," ~Treatment * Sex + (1|Animal) + (1|Animal:Image_ID)"))
  model <- lmer(formula, data = Cell_level)
  Cell_models[[value]] <- model
  
  # Treatment effect separated by Sex
  emm_TS <- emmeans(model, pairwise ~ Treatment | Sex,adjust = "none",pbkrtest.limit = 10000)
  contrast_table_TS <- as.data.frame(emm_TS$contrasts)
  # Treatment effect 
  emm_T <- emmeans(model, pairwise ~ Treatment,adjust = "none",pbkrtest.limit = 10000) #pbkrtest da durch die vielen mitos system automatisch eine freiheitsgrad berechnung runter macht da sons viel rechenleistung. ich will aber die mitos drin haben
  contrast_table_T <- as.data.frame(emm_T$contrasts)
  
  # Interaction
  interaction <- summary(model)$coefficients["TreatmentTAM:Sexmale",]
  
  Cell_results <- rbind(
    Cell_results,
    data.frame(
      Variable = value,
      Treatment_Effect =  contrast_table_T$estimate[contrast_table_T$contrast=="EtOH - TAM"],
      Treatment_p =  contrast_table_T$p.value[contrast_table_T$contrast=="EtOH - TAM"],
      Female_Treatment_Effect =  contrast_table_TS$estimate[contrast_table_TS$Sex=="female"],
      Female_p = contrast_table_TS$p.value[contrast_table_TS$Sex=="female"],
      Male_Treatment_Effect =  contrast_table_TS$estimate[contrast_table_TS$Sex=="male"],
      Male_p =contrast_table_TS$p.value[contrast_table_TS$Sex=="male"],
      Sex_Treatment_Interaction =interaction["Estimate"],
      Sex_Treatment_Interaction_p =  interaction["Pr(>|t|)"],
      row.names = NULL
    )
  )
  
}
write.csv2(Cell_results, file=file.path(PATHS$TEM$output,"Statistics/01_Cells_level_Paramters.csv"),row.names = FALSE)
#09.08.26 Not sure if I should adjust within Mito_level and Cell Level or adjust in combined level.
#No I adjust on combined level. In the end everything is not significant no mater how I adjust or dont adjust.
# so i guess doesnt matter. would still be interesting to know.

Cellresults<-read.csv2( file=file.path(PATHS$TEM$output,"Statistics/01_Cells_level_Paramters.csv"))
Mitoresults <-read.csv2( file=file.path(PATHS$TEM$output,"Statistics/01_Mito_level_Paramters.csv"))
Total_Stats<-rbind(Cellresults,Mitoresults)%>% 
  mutate(
  adj_Treatment_p = p.adjust(Treatment_p, method = "fdr"), #für plot ctrl vs TAm ohne geschlechter
  adj_Female_p = p.adjust(Female_p, method = "fdr"), #für rein female plot
  adj_Male_p = p.adjust(Male_p, method = "fdr"),   #für rein male plot
  adj_Sex_Treatment_Interaction_p = p.adjust(Sex_Treatment_Interaction_p, method = "fdr")) # ist sex treatment interaction signifikan ist reiner ctrl tam plot nicht die end aussage, dann kommen noch geschlechter. ist nicht sig: Es wurde kein Sex treatment effekt festgestellt
write.csv2(Total_Stats, file=file.path(PATHS$TEM$output,"Statistics/01_Total_Stats.csv"),row.names = FALSE)



# ============================================================
# Model diagnostics / QC
# ============================================================

get_model_QC <- function(models, level_name) {
  
  results <- lapply(names(models), function(variable) {
    
    model <- models[[variable]]
    
    # -----------------------------
    # Basic model information
    # -----------------------------
    
    n_obs <- nobs(model)
    
    # -----------------------------
    # Singularity
    # -----------------------------
    
    singular <- isSingular(model, tol = 1e-4)
    
    # -----------------------------
    # Convergence
    # -----------------------------
    
    conv_messages <- model@optinfo$conv$lme4$messages
    
    convergence_ok <- is.null(conv_messages)
    
    convergence_message <- if (convergence_ok) {
      NA_character_
    } else {
      paste(conv_messages, collapse = " | ")
    }
    
    # -----------------------------
    # Random effects
    # -----------------------------
    
    vc <- as.data.frame(VarCorr(model))
    
    # Random-effect variances
    animal_var <- vc$vcov[
      vc$grp == "Animal"
    ]
    
    image_var <- vc$vcov[
      vc$grp == "Animal:Image_ID"
    ]
    
    cell_var <- vc$vcov[
      vc$grp == "Animal:Image_ID:Cell_ID"
    ]
    
    residual_var <- vc$vcov[
      vc$grp == "Residual"
    ]
    
    # Number of levels for each random effect
    re_levels <- lapply(ranef(model), nrow)
    
    n_animal <- if ("Animal" %in% names(re_levels)) {
      re_levels$Animal
    } else NA
    
    n_image <- if ("Animal:Image_ID" %in% names(re_levels)) {
      re_levels$`Animal:Image_ID`
    } else NA
    
    n_cell <- if ("Animal:Image_ID:Cell_ID" %in% names(re_levels)) {
      re_levels$`Animal:Image_ID:Cell_ID`
    } else NA
    
    # -----------------------------
    # ICC
    # -----------------------------
    
    total_var <- sum(vc$vcov)
    ICC_animal <- ifelse(length(animal_var) == 1,animal_var / total_var,  NA)
    ICC_image <- ifelse(length(image_var) == 1,image_var / total_var, NA)
    ICC_cell <- ifelse(length(cell_var) == 1,  cell_var / total_var,NA)
    
    # -----------------------------
    # Residual diagnostics
    # -----------------------------
    
    residuals_model <- residuals(model)
    
    residual_mean <- mean(residuals_model, na.rm = TRUE)
    residual_sd <- sd(residuals_model, na.rm = TRUE)
    
    residual_min <- min(residuals_model, na.rm = TRUE)
    residual_max <- max(residuals_model, na.rm = TRUE)
    
    # Standardized residuals
    standardized_residuals <- scale(residuals_model)[,1]
    
    n_extreme_residuals <- sum(abs(standardized_residuals) > 3, na.rm = TRUE)
    
    # -----------------------------
    # Fixed effects
    # -----------------------------
    
    coef_table <- summary(model)$coefficients
    
    treatment_p <- if ("TreatmentTAM" %in% rownames(coef_table)) {
      coef_table["TreatmentTAM", "Pr(>|t|)"]
    } else NA
    
    sex_p <- if ("Sexmale" %in% rownames(coef_table)) {
      coef_table["Sexmale", "Pr(>|t|)"]
    } else NA
    
    interaction_p <- if ("TreatmentTAM:Sexmale" %in% rownames(coef_table)) {
      coef_table["TreatmentTAM:Sexmale", "Pr(>|t|)"]
    } else NA
    
    # -----------------------------
    # Return one row
    # -----------------------------
    
    data.frame(
      Level = level_name,
      Variable = variable,
      
      N = n_obs,
      
      Singular = singular,
      Convergence_OK = convergence_ok,
      Convergence_Message = convergence_message,
      
      N_Animals = n_animal,
      N_Images = n_image,
      N_Cells = n_cell,
      Var_Animal = ifelse(length(animal_var) == 1,animal_var, NA),
      Var_Image = ifelse(length(image_var) == 1, image_var, NA),
      Var_Cell = ifelse(length(cell_var) == 1,  cell_var, NA),
      Var_Residual = ifelse(length(residual_var) == 1,residual_var, NA),
      
      ICC_Animal = ICC_animal,
      ICC_Image = ICC_image,
      ICC_Cell = ICC_cell,
      
      Residual_Mean = residual_mean,
      Residual_SD = residual_sd,
      Residual_Min = residual_min,
      Residual_Max = residual_max,
      N_Extreme_Residuals = n_extreme_residuals,
      
      AIC = AIC(model),
      BIC = BIC(model),
      
      Treatment_p = treatment_p,
      Sex_p = sex_p,
      Treatment_Sex_Interaction_p = interaction_p,
      
      stringsAsFactors = FALSE
    )
  })
  
  bind_rows(results)
}


# ============================================================
# Run QC
# ============================================================

Mito_model_QC <- get_model_QC(Mito_models,level_name = "Mito")

Cell_model_QC <- get_model_QC(Cell_models,level_name = "Cell")


# Combine
Model_QC <- bind_rows(Mito_model_QC,Cell_model_QC)




write.csv2( Model_QC, file = file.path(PATHS$TEM$output,"Statistics/01_Model_QC.csv"),row.names = FALSE)
Model_QC
Model_QC %>%
  select(
    Level,
    Variable,
    N,
    Singular,
    Convergence_OK,
    N_Animals,
    N_Images,
    N_Cells,
    Var_Animal,
    Var_Image,
    Var_Cell,
    ICC_Animal,
    ICC_Image,
    ICC_Cell,
    N_Extreme_Residuals
  )
Model_QC %>%filter( Singular == TRUE |  Convergence_OK == FALSE )

# Model diagnostics for all Mito-level models ----------------------------

for (value in names(Mito_models)) {
  
  model <- Mito_models[[value]]
  
  png(
    filename = file.path(
      PATHS$TEM$output,
      "Background",
      paste0("QC_Model_", value, ".png")
    ),
    width = 1800,
    height = 1800,
    res = 300
  )
  
  par(mfrow = c(2, 2))
  
  # 1. Residuals vs fitted
  plot(
    model,
    main = paste(value, "- Residuals vs Fitted")
  )
  
  # 2. Q-Q plot
  qqnorm(
    resid(model),
    main = paste(value, "- Q-Q plot")
  )
  qqline(resid(model))
  
  # 3. Histogram of residuals
  hist(  resid(model),
    main = paste(value, "- Residuals"),
    xlab = "Residuals")
  
  # 4. Residuals vs observation order
  plot(  resid(model),
    main = paste(value, "- Residuals vs Observation"),
    xlab = "Observation",
    ylab = "Residual")
  abline(h = 0, lty = 2)
  
  dev.off()}

for (value in names(Cell_models)) {
  
  model <- Cell_models[[value]]
  
  png(filename = file.path(PATHS$TEM$output,"Background",paste0("QC_Model_", value, ".png") ),
    width = 1800, height = 1800,  res = 300 )
  
  par(mfrow = c(2, 2))
  
  # 1. Residuals vs fitted
  plot(model, main = paste(value, "- Residuals vs Fitted"))
  
  # 2. Q-Q plot
  qqnorm(
    resid(model),
    main = paste(value, "- Q-Q plot") )
  qqline(resid(model))
  
  # 3. Histogram of residuals
  hist(  resid(model), main = paste(value, "- Residuals"), xlab = "Residuals")
  
  # 4. Residuals vs observation order
  plot(  resid(model),   main = paste(value, "- Residuals vs Observation"), xlab = "Observation", ylab = "Residual")
  abline(h = 0, lty = 2)
  dev.off()
  }
