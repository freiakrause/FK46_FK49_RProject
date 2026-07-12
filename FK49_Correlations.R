
# libraries --------------------------------------------------------------------
library(tidyverse)
library(pheatmap) 
library(svDialogs)
# library(rstudioapi)
library(ggfortify)
library(factoextra)
library(broom)
library(purrr)
library('corrr')
library(ggcorrplot)
library("FactoMineR")
library(survminer)
library(patchwork)
library(superb)
library(ggbreak)
library(tibble)
library(waffle)
library(rstatix)
library(lmerTest)
library(emmeans)
library(grid)
library(ggnewscale)
library(NADA2)
library(effsize)

ExpId="FK49"

if (ExpId=="FK49") {
  setwd("C:/Users/b1084855/OneDrive - Universität Salzburg/AG_Tumorimmunologie - Dokumente/Data/Freia Krause/01_Experiments/FK49_CD-HFD_13wks/FK49_Analysis/")
  
}  else if (ExpId == "FK46"){
  setwd("C:/Users/b1084855/OneDrive - Universität Salzburg/AG_Tumorimmunologie - Dokumente/Data/Freia Krause/01_Experiments/FK46_iALmice_high Fat diet 52 weeks 7d after injection/Analysis")
} else{
  print("Let me set a folder Path and define ExpId and Folderpath")}
Exigo_cols= c("ALB","TP","GLOB","A.G","TB","GGT","AST","ALT","ALP","AMY","Crea","UA","BUN","GLU","TC","TG")
cytokines <-c("IL23","IL1a" ,"IFNy" ,"TNFa" , "MCP1" , "IL12p70" ,"IL1ß", "IL10" ,"IL6" ,"IL27" ,"IL17A" , "IFNß" ,"GMCSF") 
cytokine_pattern <- paste(cytokines, collapse = "|")
to_correlate<-c(cytokines,"SAA_last",
                Exigo_cols,
                "NASH_I","NASH_S","NASH_B","NASH_SAF",
                "Liver","Weight","Spleen")
load(paste0("01_RawData/",ExpId,"_Data_prepared.Rda"))

SAA_last <- data %>%
  filter(!is.na(SAA)) %>%
  arrange(Animal, DOW) %>%
  group_by(Animal) %>%
  slice_tail(n = 1) %>%
  ungroup() %>%
  select(Animal, SAA_last = SAA)
endpoint <- data %>%
  filter(!is.na(Liver)) %>%
  left_join(SAA_last, by = "Animal")

endpoint <- endpoint %>%
  mutate(across(c(all_of(to_correlate)), ~ as.numeric(.)))
endpoint <- endpoint %>%
  mutate(
    Sex_num = case_when(
      Sex %in% c("male", "M") ~ 1,
      Sex %in% c("female", "F") ~ 0,
      TRUE ~ NA_real_
    ),
    Treatment_num = case_when(
      Treatment == "TAM" ~ 1,
      Treatment == "ctrl" ~ 0,
      TRUE ~ NA_real_
    )
  )


d_sex_M<-endpoint%>%filter(Sex=="male")
d_sex_F<-endpoint%>%filter(Sex=="female")
d_Treat_T<-endpoint%>%filter(Treatment=="TAM")
d_Treat_C<-endpoint%>%filter(Treatment=="ctrl")
cor_all <- endpoint %>%
  select(all_of(to_correlate),"SAA_last","Sex_num", "Treatment_num") %>%
  cor(use = "pairwise.complete.obs", method = "spearman")
cor_sex_M <- d_sex_M %>%
  select(all_of(to_correlate), "SAA_last","Treatment_num") %>%
  cor(use = "pairwise.complete.obs", method = "spearman")
cor_sex_F <- d_sex_F %>%
  select(all_of(to_correlate),"SAA_last", "Treatment_num") %>%
  cor(use = "pairwise.complete.obs", method = "spearman")
cor_all_plot   <- cor_all
cor_sex_M_plot <- cor_sex_M
cor_sex_F_plot <- cor_sex_F

cor_all_plot[is.na(cor_all_plot)]     <- 0
cor_sex_M_plot[is.na(cor_sex_M_plot)] <- 0
cor_sex_F_plot[is.na(cor_sex_F_plot)] <- 0
pheatmap(
  cor_sex_M_plot,
  color=colorRampPalette(c("blue", "grey90", "grey95", "red"))(100),
  breaks = seq(-1, 1, length.out = 101),  # 100 intervals  clustering_distance_rows = "correlation",
  clustering_distance_cols = "correlation",
  clustering_distance_rows = "correlation",
  main = "Correlation in Males",
  display_numbers = TRUE,
  number_format = "%.1f"
)
pheatmap(
  cor_sex_F_plot,
  color=colorRampPalette(c("blue", "grey90", "grey95", "red"))(100),
  breaks = seq(-1, 1, length.out = 101),  # 100 intervals  clustering_distance_rows = "correlation",
  clustering_distance_cols = "correlation",
  clustering_distance_rows = "correlation",
  main = "Correlation in Females",
  display_numbers = TRUE,
  number_format = "%.1f"
)
cor_all_heat<-pheatmap(
  cor_all_plot,
  color=colorRampPalette(c("blue", "grey90", "grey95", "red"))(100),
  breaks = seq(-1, 1, length.out = 101),  # 100 intervals  clustering_distance_rows = "correlation",
  clustering_distance_cols = "correlation",
  clustering_distance_rows = "correlation",
  main = "Correlation",
  display_numbers = FALSE,
  number_format = "%.1f"
)
ggsave(cor_all_heat, file= "Correlationall.png",width= 8, height= 8,
       path="C:/Users/b1084855/OneDrive - Universität Salzburg/AG_Tumorimmunologie - Dokumente/Data/Freia Krause/01_Experiments/FK49_CD-HFD_13wks/FK49_Analysis/02_GeneratedData/FK49_Correlation")

relevant_cor <- cor_all %>%
  as.data.frame() %>%
  mutate(Var1 = rownames(cor_all)) %>%
  tidyr::pivot_longer(-Var1, names_to = "Var2", values_to = "r") %>%
  filter(abs(r) >= 0.7)



# Convert selected variables to numeric
data_cor <- endpoint %>%
  mutate(across(all_of(to_correlate), ~ as.numeric(.))) %>%
  mutate(
    Sex_num = ifelse(Sex == "male", 1, ifelse(Sex == "female", 0, NA)),
    Treatment_num = ifelse(Treatment == "TAM", 1, ifelse(Treatment == "ctrl", 0, NA))
  )

to_correlate_extended <- c(
  setdiff(to_correlate, "SAA_last"),
  "Sex_num",
  "Treatment_num"
)
var_pairs <- expand.grid(var1 = to_correlate_extended,
                         var2 = to_correlate_extended,
                         stringsAsFactors = FALSE) %>%
  filter(var1 != var2)  # remove self-correlations

cor_results <- var_pairs %>%
  mutate(
    test = map2(var1, var2, ~ cor.test(data_cor[[.x]], data_cor[[.y]], method = "spearman")),
    rho = map_dbl(test, ~ .x$estimate),      # Spearman rho
    p_value = map_dbl(test, ~ .x$p.value)    # p-value
  ) %>%
  select(var1, var2, rho, p_value)%>%
  mutate(
    p_adj = p.adjust(p_value, method = "fdr")  # adjust p-values for multiple testing
  )

