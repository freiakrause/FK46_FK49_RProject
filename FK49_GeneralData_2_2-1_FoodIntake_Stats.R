rm(list=ls())
gc()
library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)
library(lmerTest)
library(lme4)      
library(emmeans)
library(rstatix)
library(patchwork)
library(readr)     
source("FK49_Definitions.R")
# Read Raw Inputdata after general Data manipulation ------------------------------------------------------
# ### 10.08.26 hier morgen weiter machen
# #statisti is noch wild. ich möchte nur ein linear model mit allen nötigen effetcs drin haben.
#  dann im interaction sehen, was important ist und von da entscheidne, ob einsummary plot oder verschiedenen plots übre zeit.
#  Food_consumed ~ days_diet * Treatment * Sex + Batch + (1 + days_diet| Cage)
#  #oder wks diet?


output_pwd<-file.path(PATHS$FoodIntake$output)
load(file.path(PATHS$FoodIntake$input,"FK49_Data_prepared.Rda"))
exigo <- c("ALB", "TP", "GLOB","A.G", "TB", "GGT", "AST", "ALT", "ALP", "AMY","Crea","UA","BUN","GLU","TC","TG")

#Plot Food and Water ------------------
#Check weather Batch1 and Batch 2 behave the same in food and water consumption. In Batch 1 animals were single housed, in Batch 2 group housed 1-5 animals/group. Consumption is cagewise
#before i can combine batches for analysis I nedd to test weather the change in methodology changed the data.
d_for_food<-data%>%select(-all_of(exigo),-matches("Exigo"),-matches("direction"),-matches("value"),-matches("censored"),-matches("NASH"),-"DFactor")%>%
  filter(!Animal %in% c("EC1", "EC2","EC3"))%>% #,#"161",#"166","160"))%>% #excluded 161,166,160 bc they alwyas destroyed food and left it in theri cages - affecting food measurment
  filter(!(Animal == "164"  & days_diet == 24))%>% #excluded this measurment bc from d4 to d5 weight incresead and this is not possible-> measurment error
  filter(!(Animal == "164"  & days_diet == 25))%>%##excluded this measurment bc from d4 to d5 weight incresead and this is not possible-> measurment error; weight shift from d3 to d4 seems to be an outlier, maybe on d4 was wrong measurement
  #filter(!(Cage=="17"))%>%
  #filter(!(Cage=="15"))%>%
  filter(!(Block=="0"))%>% # batch 2 measurements in first weeks were not complete
  filter(!(Animal=="203" & DOW=="2025-07-30"))%>%
  filter(!(is.na(Food_consumed)))%>%
  filter((Block%in% c("1","4","8","12")))

### I summarize measurements as cagewise since I/animal caretakes only measrued cageweise. So we can not give values vor individual mice(in bacth 2)
d_food_cagewise <- d_for_food %>%
  group_by(Cage,Sex,days_diet, wks_diet,Treatment, BATCH, Block) %>%
  summarise(Food_consumed = unique(Food_consumed),  # all animals have same value
            n_animals = first(n_animals),           # optional
            Water_consumed = unique(Water_consumed),
            rel.weight=mean(rel.weight),
            .groups = "drop",)  # all animals have same value)
#rm(d_for_food)
### Removing outliers due which potentially occured due to spilling water and Food Krümelmonster mice -------
outliers_food <- d_food_cagewise %>%group_by(BATCH, Treatment) %>%identify_outliers(Food_consumed) %>%
  ungroup() %>%
  filter(is.outlier == TRUE) %>%
  select(Cage, days_diet, BATCH, Treatment, Food_consumed) %>%
  mutate(Food_consumed_outlier = TRUE)

# Identify extreme outliers for Water ---
outliers_water <- d_food_cagewise %>%
  group_by(BATCH, Treatment) %>%
  identify_outliers(Water_consumed) %>%
  ungroup() %>%
  filter(is.outlier == TRUE) %>%
  select(Cage, days_diet, BATCH, Treatment, Water_consumed) %>%
  mutate(Water_consumed_outlier = TRUE)

# Replace only those extreme outlier values with NA ---
d_food_cagewise <- d_food_cagewise %>%
  left_join(outliers_food, by = c("Cage", "days_diet", "BATCH", "Treatment", "Food_consumed")) %>%
  left_join(outliers_water, by = c("Cage", "days_diet", "BATCH", "Treatment", "Water_consumed")) %>%
  mutate( Food_consumed = if_else(!is.na(Food_consumed_outlier), NA_real_, Food_consumed),
          Water_consumed = if_else(!is.na(Water_consumed_outlier), NA_real_, Water_consumed)) %>%
  select(-Food_consumed_outlier, -Water_consumed_outlier)

d_food_analysis <- d_food_cagewise %>%
  filter(!is.na(Food_consumed))
saveRDS(d_food_analysis,file = file.path(PATHS$FoodIntake$input,"FoodIntake_Cagewise.rds"))

d_water_analysis <- d_food_cagewise %>%
  filter(!is.na(Water_consumed))
saveRDS(d_water_analysis,file = file.path(PATHS$FoodIntake$input,"WaterIntake_Cagewise.rds"))

### Check if unbalanced or missing data ----
summary(d_food_analysis)
any(is.na(d_food_analysis$Food_consumed)) #FALSE
table(d_food_analysis$BATCH, d_food_analysis$Treatment) 
summary(d_water_analysis)
any(is.na(d_water_analysis$Water_consumed))#FALSE
table(d_water_analysis$BATCH, d_water_analysis$Treatment) 
#ctrl1 80 TAM1 78 Ctrl2 16, TAM2 47 might be unbalance in batch 2

# Food Statistics  -----
# ### Fit linear mied effects model -----
#model_F <- lmer(Food_consumed ~ Treatment * Sex *days_diet + BATCH +(1 + days_diet| Cage), data = d_food_analysis)
#ranef(model_F)
## Cages bascially dont have random individual slope and random individual baseline
formula_F_interaction <- "Food_consumed ~ Treatment * Sex * days_diet + BATCH + (1 | Cage)"
model_F_Interaction <- lmer(formula_F_interaction, data = d_food_analysis)
ranef(model_F_Interaction)
VarCorr(model_F_Interaction)
summary(model_F_Interaction)
anova_table_F_interaction<-anova(model_F_Interaction, type = 3)
formula_F_additive    <- "Food_consumed ~ Treatment + Sex + days_diet + BATCH + (1 | Cage)"
model_F_additive <- lmer(formula_F_additive, data = d_food_analysis)
ranef(model_F_additive)
VarCorr(model_F_additive)
summary(model_F_additive)
anova_table_F_additive<-anova(model_F_additive, type = 3)
anova_F_interaction <- as.data.frame(anova_table_F_interaction) %>%  rownames_to_column("Effect")
anova_F_additive <- as.data.frame(anova_table_F_additive) %>%  rownames_to_column("Effect")
make_interpretation <- function(effect, p) {
  case_when(
    effect == "Treatment" & p > 0.05 ~   "No evidence for a Treatment effect",
    effect == "Treatment" & p <= 0.05 ~   "Evidence for a Treatment effect",
    effect == "Sex" & p <= 0.05 ~"Evidence for a Sex effect",
    effect == "Sex" & p > 0.05 ~"No Evidence for a Sex effect",
    effect == "days_diet" & p <= 0.05 ~  "Food consumption changed significantly over diet duration",
    effect == "days_diet" & p >0.05 ~  "No Evidence for a Time effect",
    effect == "BATCH" & p  <= 0.05 ~  "Significant difference between experimental batches",
    effect == "BATCH" & p > 0.05 ~  "No significant difference between experimental batches",
    effect == "Treatment:Sex" & p > 0.05 ~  "No evidence for a Treatment * Sex interaction",
    effect == "Treatment:Sex" & p <= 0.05 ~  "Evidence for a Treatment * Sex interaction",
    effect == "Treatment:days_diet" & p > 0.05 ~  "No evidence for a Treatment * time interaction",
    effect == "Treatment:days_diet" & p <= 0.05 ~  "Evidence for a Treatment * time interaction",
    effect == "Sex:days_diet" & p > 0.05 ~ "No evidence for a Sex * time interaction",
    effect == "Sex:days_diet" & p <= 0.05 ~ "Evidence for a Sex * time interaction",
    effect == "Treatment:Sex:days_diet" & p > 0.05 ~   "No evidence for a Treatment * Sex * time interaction",
    effect == "Treatment:Sex:days_diet" & p <= 0.05 ~   "Evidence for a Treatment * Sex * time interaction",
    TRUE ~ ""
  )
}

results_F_interaction <- anova_F_interaction %>%
  mutate(Model = "Interaction",  Formula = formula_F_interaction,  Interpretation = make_interpretation(Effect, `Pr(>F)`) )
results_F_additive <- anova_F_additive %>%
  mutate(Model = "Additive",  Formula = formula_F_additive,  Interpretation = make_interpretation(Effect, `Pr(>F)`) )

food_model_results <- bind_rows(results_F_interaction,results_F_additive)
write_csv2(food_model_results, file.path(output_pwd,"Statistics/Food_consumption_LMM_results.csv"))
# 
# 
# anova_label_F <- paste0("ANOVA over linear mixed-effects model\n",
#                         "Treatment: F = ", round(food_model_results$F[1], 2), ", p = ", signif(food_model_results$`Pr(>F)`[1], 3), "\n",
#                         "Time: F = ", round(food_model_results$F[2], 2), ", p < ", format.pval(food_model_results$`Pr(>F)`[2], digits = 1), "\n",
#                         "Interaction: F = ", round(food_model_results$F[3], 2), ", p = ", signif(food_model_results$`Pr(>F)`[3], 3))

block_days <- d_food_analysis %>%
  group_by(Block) %>% summarise(days_diet = median(days_diet, na.rm = TRUE),.groups = "drop"  )

# Calculate the estimated marginal means for both Treatment and wks_diet
emm_F_additive <- emmeans( model_F_additive,~ Treatment | days_diet,at = list(days_diet = block_days$days_diet))
emm_df_F_additive <- as.data.frame(emm_F_additive) %>%  left_join(block_days, by = "days_diet")
write_csv2(emm_df_F_additive,file.path(output_pwd, "Statistics/Food_consumption_LMM_EMMs_by_Block.csv"))
# Perform pairwise contrasts between Treatment levels (ctrl vs TAM) at all levels of wks_diet
pwc_F <- contrast(emm_F_additive,method = "pairwise",adjust = "bonferroni")
pwc_df_F <- as.data.frame(pwc_F)
posthoc_label_F<- "Post Hoc: Pairwise with Bonferroni correction"

pwc_df_F <- as.data.frame(pwc_F)%>% left_join(  block_days,  by = "days_diet")
pwc_df_F <- pwc_df_F %>%
  mutate(rounded_p_value = ifelse( is.na(p.value), "NA", format.pval(p.value, digits = 3, eps = 0.001)),
         significance = case_when(
           is.na(p.value) ~ "NA",
           p.value <= 0.001 ~ "***",
           p.value <= 0.01 ~ "**",
           p.value <= 0.05 ~ "*",
           TRUE ~ "NS" ))
write_csv2( pwc_df_F, file.path(output_pwd,"Statistics/Food_consumption_LMM_posthoc_by_Block.csv"))


# Water Statisticsc -----
# ## Fit linear Mixed Effects model
formula_W_interaction <- "Water_consumed ~ Treatment * Sex * days_diet + BATCH + (1 | Cage)"
model_W_Interaction <- lmer(formula_W_interaction, data = d_water_analysis)
ranef(model_W_Interaction)
VarCorr(model_W_Interaction)
summary(model_W_Interaction)
anova_table_W_interaction<-anova(model_W_Interaction, type = 3)
# in interaction model noo  interaction effect of Treatment Sex effect, no evidence of sex diet effect, no evidence of t s time effect
#don t use result from additive model because addive assumes treatment is effect doe snot differ ver tim. but interacton model shows taht
# formula_W_additive    <- "Water_consumed ~ Treatment + Sex + days_diet + BATCH + (1 | Cage)"
# model_W_additive <- lmer(formula_W_additive, data = d_water_analysis)
# ranef(model_W_additive)
# VarCorr(model_W_additive)
# summary(model_W_additive)
# anova_table_W_additive<-anova(model_W_additive, type = 3)

formula_W_simple    <- "Water_consumed ~ Treatment * days_diet+ Sex + BATCH + (1 | Cage)"
model_W_simple <- lmer(formula_W_simple, data = d_water_analysis)

ranef(model_W_simple)
VarCorr(model_W_simple)
summary(model_W_simple)
anova_table_W_simple<-anova(model_W_simple, type = 3)

# Model formulas
anova_W_interaction <- as.data.frame(anova_table_W_interaction) %>%  rownames_to_column("Effect")
# anova_W_additive <- as.data.frame(anova_table_W_additive) %>%  rownames_to_column("Effect")
anova_W_simple <- as.data.frame(anova_table_W_simple) %>%  rownames_to_column("Effect")
# result from simple Water consuption differed between treamten and magnitude of treatment effect changed over course of diet.

results_W_interaction <- anova_W_interaction %>%
  mutate(Model = "Interaction",  Formula = formula_W_interaction,  Interpretation = make_interpretation(Effect, `Pr(>F)`) )
# results_W_additive <- anova_W_additive %>%
  # mutate(Model = "Additive",  Formula = formula_W_additive,  Interpretation = make_interpretation(Effect, `Pr(>F)`) )
results_W_simple <- anova_W_simple %>%
  mutate(Model = "Simple",  Formula = formula_W_simple,  Interpretation = make_interpretation(Effect, `Pr(>F)`) )
# interaction model showd 
water_model_results <- bind_rows(results_W_interaction,results_W_simple)
write_csv2(water_model_results, file.path(output_pwd,"Statistics/Water_consumption_LMM_results.csv"))
# 
# 
# anova_label_F <- paste0("ANOVA over linear mixed-effects model\n",
#                         "Treatment: F = ", round(food_model_results$F[1], 2), ", p = ", signif(food_model_results$`Pr(>F)`[1], 3), "\n",
#                         "Time: F = ", round(food_model_results$F[2], 2), ", p < ", format.pval(food_model_results$`Pr(>F)`[2], digits = 1), "\n",
#                         "Interaction: F = ", round(food_model_results$F[3], 2), ", p = ", signif(food_model_results$`Pr(>F)`[3], 3))

block_days <- d_water_analysis %>%
  group_by(Block) %>% summarise(days_diet = median(days_diet, na.rm = TRUE),.groups = "drop"  )
# Calculate the estimated marginal means for both Treatment and wks_diet
emm_W_simple <- emmeans( model_W_simple,~ Treatment | days_diet,at = list(days_diet = block_days$days_diet))
emm_df_W_simple <- as.data.frame(emm_W_simple) %>%  left_join(block_days, by = "days_diet")
write_csv2(emm_df_W_simple,file.path(output_pwd, "Statistics/Water_consumption_LMM_EMMs_by_Block.csv"))
# Perform pairwise contrasts between Treatment levels (ctrl vs TAM) at all levels of wks_diet
pwc_W <- contrast(emm_W_simple,method = "pairwise",adjust = "bonferroni")
pwc_df_W <- as.data.frame(pwc_W)
posthoc_label_W<- "Post Hoc: Pairwise with Bonferroni correction"

pwc_df_W <- as.data.frame(pwc_W)%>% left_join(  block_days,  by = "days_diet")
pwc_df_W <- pwc_df_W %>%
  mutate(rounded_p_value = ifelse( is.na(p.value), "NA", format.pval(p.value, digits = 3, eps = 0.001)),
         significance = case_when(
           is.na(p.value) ~ "NA",
           p.value <= 0.001 ~ "***",
           p.value <= 0.01 ~ "**",
           p.value <= 0.05 ~ "*",
           TRUE ~ "NS" ))
write_csv2( pwc_df_W, file.path(output_pwd,"Statistics/Water_consumption_LMM_posthoc_by_Block.csv"))

## Statistic Analysis is now reasonable. Now i need to adjust plottings and incude stats in plots

rm(list=ls())
gc