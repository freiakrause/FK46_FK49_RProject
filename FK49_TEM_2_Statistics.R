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
  select(-Sex, -Treatment, -Animal, -Image_ID, -Mito_ID,-Cell_ID, -Pedigree_ID,-Shape, -Analysis_ID,-Cell_Type) %>%
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
  select(-Sex, -Treatment, -Animal, -Image_ID, -Cell_ID, -Cell_Type,-Analysis_ID,-Cell_Area,
         -starts_with("SD_") ,-starts_with("Mean")) %>% #nimmt die werte raus, die auf mito ebene schongetestet wurden
  colnames()
Cell_values
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

summary(Mito_level)
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

Cellresults<-read.csv2( file=file.path(PATHS$TEM$output,"Statistics/01_Cells_level_Paramters.csv"))
Mitoresults <-read.csv2( file=file.path(PATHS$TEM$output,"Statistics/01_Mito_level_Paramters.csv"))
Total_Stats<-rbind(Cellresults,Mitoresults)%>% 
  mutate(
  adj_Treatment_p = p.adjust(Treatment_p, method = "fdr"), #für plot ctrl vs TAm ohne geschlechter
  adj_Female_p = p.adjust(Female_p, method = "fdr"), #für rein female plot
  adj_Male_p = p.adjust(Male_p, method = "fdr"),   #für rein male plot
  adj_Sex_Treatment_Interaction_p = p.adjust(Sex_Treatment_Interaction_p, method = "fdr")) # ist sex treatment interaction signifikan ist reiner ctrl tam plot nicht die end aussage, dann kommen noch geschlechter. ist nicht sig: Es wurde kein Sex treatment effekt festgestellt
write.csv2(Total_Stats, file=file.path(PATHS$TEM$output,"Statistics/01_Total_Stats.csv"),row.names = FALSE)
