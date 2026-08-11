rm(list=ls())
gc()
library(knitr)
library(car)
library(dplyr)
library(broom)  # for tidy() function
library(tidyr)
library(dplyr)
library(survival)
library(lubridate)
library(tidyverse)
library(ggsurvfit)
library(survminer)
library(patchwork)
library(superb)
library(ggbreak)
library(tibble)
#library(waffle)
library(rstatix)
library(lmerTest)
library(emmeans)
library(grid)
library(ggnewscale)
library(NADA2)
library(effsize)
source("FK49_Definitions.R")
# Read Raw Inputdata after general Data manipulation ------------------------------------------------------
output_pwd <-file.path(PATHS$general_data$FK49_output)
load(file.path(output_pwd,"01_RawData/FK49_Data_prepared.Rda"))
exigo <- c("ALB", "TP", "GLOB","A.G", "TB", "GGT", "AST", "ALT", "ALP", "AMY","Crea","UA","BUN","GLU","TC","TG")
output_pwd <-file.path(PATHS$general_data$FK49_output)
#

  
  # Dotplot Organ Weights --------------------------------------------------------------
## Function for Organ Weights --------------------------------------------------------------
do_organ_weight <- function(inputdata, value, batch = "ALL", sex = "both", y_title, path_images,colors = c("black", "darkred")) {
  ### Data manipulation -----------------------------------------------------------------
  d <- inputdata %>%
    #select(Animal, Sex, Treatment, Weight, Liver, Fat, Spleen, Ascites.no.yes, Tumor.no.yes, wks_diet, BATCH) %>%
    filter(complete.cases(.data[[value]])) %>%
    mutate(Tumor.no.yes = as.factor(Tumor.no.yes),
      Ascites.no.yes = as.factor(Ascites.no.yes),
      event_status = case_when(
        Tumor.no.yes == "0" & Ascites.no.yes == "0" ~ "normal",
        Tumor.no.yes == "1" & Ascites.no.yes == "0" ~ "tumor",
        Tumor.no.yes == "0" & Ascites.no.yes == "1" ~ "ascites",
        Tumor.no.yes == "1" & Ascites.no.yes == "1" ~ "both",
        TRUE ~ "unknown")) %>%
    mutate(event_status = factor(event_status, levels = c("normal", "ascites", "tumor", "both", "unknown")))
  
  # Optional: Filter by batch and sex
  if (batch != "ALL") d <- d %>% filter(BATCH == batch)
  if (sex != "both") d <- d %>% filter(Sex == sex)
  
  # Convert selected value to numeric if not already
  d[[value]] <- as.numeric(d[[value]])
  
  ### Summary stats ----------------------------------------------------------------------
  stats <- d %>%
    group_by(Treatment) %>%
    summarise(
      Mean = mean(.data[[value]], na.rm = TRUE),
      SD = sd(.data[[value]], na.rm = TRUE),
      .groups = "drop")
  
  ### Normality test ---------------------------------------------------------------------
  shapiro_test_ctrl <- shapiro.test(d[[value]][d$Treatment == "ctrl"])
  shapiro_test_tam <- shapiro.test(d[[value]][d$Treatment == "TAM"])
  print(shapiro_test_tam)
  print(shapiro_test_ctrl)
  ### T-test ------------------------------------------------------------------------------
  t_test_result <- t.test(as.formula(paste(value, "~ Treatment")), data = d)
  p_value <- t_test_result$p.value
  p_value_label <- paste("p =", format(p_value, digits = 3))
  
  ### Effect size -------------------------------------------------------------------------
  library(effsize)
  cohen_d_result <- cohen.d(as.formula(paste(value, "~ Treatment")), data = d)
  print(cohen_d_result)
  
  # Jitter to handle ties
  d[[value]] <- jitter(d[[value]], amount = 0.001)
  d$wks_diet <- jitter(d$wks_diet, amount = 0.001)
  
 
  ### Plot 1 ------------------------------------------------------------------------------
  p1 <- ggplot(stats, aes(x = Treatment, y = Mean, fill = Treatment)) +
    geom_bar(stat = "identity", color = "black", alpha = 0.5, width = 0.75, position = "dodge") +
    geom_point(data = d, fill = "lightgrey", color = "black",
               aes(y = .data[[value]], shape = event_status),
               position =  position_jitterdodge(0.1, dodge.width = 0.75), alpha = 0.8, size = 5.3, stroke = 1.8)+
    geom_errorbar(aes(ymin = Mean - SD, ymax = Mean + SD),  position = position_dodge(width = 0.75), width = 0.2) +
    scale_fill_manual(values = colors, labels = c("Ctrl", "TAM")) +
    scale_shape_manual(values = c(21, 22, 24, 25, 26)) +
    scale_y_continuous(name = y_title) +
    labs(x = "Treatment") +
    annotate("text", x = 1.5, y = max(d[[value]], na.rm = TRUE) * 1.1, label = p_value_label, size = 6, fontface = "italic") +
    theme_minimal() +
    theme(
      legend.position = "bottom",
      axis.line = element_line(color = "black", linewidth = 0.5),
      axis.ticks = element_line(color = "black", linewidth = 0.5),
      axis.title = element_text(size = 20, face = "bold"),
      axis.title.x = element_blank(),
      axis.text = element_text(size = 19, face = "bold"),
      plot.title = element_blank(),
      legend.title = element_text(size = 10),
      legend.text = element_text(size = 10),
      panel.grid = element_blank()) +
  guides(
    shape = guide_legend(title = "Status", order = 2 ,nrow = 2, byrow = TRUE),
    fill = guide_legend(title = "Treatment", order = 3,override.aes = list(shape = 21), nrow =2, byrow = TRUE),
    color = "none")
  
  
  ### Save Plots -------------------------------------------------------------------------
  value_clean <- gsub("[^[:alnum:]_]", "_", value)
  filename1 <- paste0("FK49_", value_clean, "_Treatment_Batch", batch, "_", sex, ".png")
  #filename3 <- paste0("FK49_", value_clean, "_Treatment_Batch", batch, "_", sex, ".svg")
 
  ggsave(filename = filename1, plot = p1, path = path_images, width = 4, height = 11, dpi = 300)
  #ggsave(filename = filename3, plot = p1, path = file.path(path_images, "background"), width = 4.5, height = 9)
 }

## Call the function with desired arguments --------------------------------------------------------------
path_for_saving_images<-"02_GeneratedData/Weight_Organs"
do_organ_weight(data, value = "Liver_rel", batch = "ALL", sex = "both", y_title = "Liver/BW [%]",path_for_saving_images)
do_organ_weight(data, value = "Liver_rel", batch = "ALL", sex = "male", y_title = "Liver/BW [%]",path_for_saving_images)
do_organ_weight(data, value = "Liver_rel", batch = "ALL", sex = "female", y_title = "Liver/BW [%]",path_for_saving_images)

do_organ_weight(data, value = "Spleen_rel",batch = "ALL", sex = "both",y_title= "Spleen/BW [%]",path_for_saving_images)
do_organ_weight(data, value = "Spleen_rel",batch = "ALL", sex = "male",y_title= "Spleen/BW [%]",path_for_saving_images)
do_organ_weight(data, value = "Spleen_rel",batch = "ALL", sex = "female",y_title= "Spleen/BW [%]",path_for_saving_images)

do_organ_weight(data, value = "Spleen",batch = "ALL", sex = "both",y_title= "Spleen [mg]",path_for_saving_images)
do_organ_weight(data, value = "Spleen",batch = "ALL", sex = "male",y_title= "Spleen [mg]",path_for_saving_images)
do_organ_weight(data, value = "Spleen",batch = "ALL", sex = "female",y_title= "Spleen [mg]",path_for_saving_images)

do_organ_weight(data, value = "Liver",batch = "ALL", sex = "both",y_title= "Liver [g]",path_for_saving_images)
do_organ_weight(data, value = "Liver",batch = "ALL", sex = "male",y_title= "Liver [g]",path_for_saving_images)
do_organ_weight(data, value = "Liver",batch = "ALL", sex = "female",y_title= "Liver [g]",path_for_saving_images)

do_organ_weight(data, value = "Fat",batch = "ALL", sex = "male",y_title= "Fat [g]",path_for_saving_images)
do_organ_weight(data, value = "Fat",batch = "ALL", sex = "female",y_title= "Fat [g]",path_for_saving_images)
do_organ_weight(data, value = "Fat",batch = "ALL", sex = "both",y_title= "Fat [g]",path_for_saving_images)

do_organ_weight(data, value = "Fat_rel",batch = "ALL", sex = "both",y_title= "Fat/BW [%]",path_for_saving_images)
do_organ_weight(data, value = "Fat_rel",batch = "ALL", sex = "male",y_title= "Fat/BW [%]",path_for_saving_images)
do_organ_weight(data, value = "Fat_rel",batch = "ALL", sex = "female",y_title= "Fat/BW [%]",path_for_saving_images)

dd<-data%>%group_by(Animal,Treatment,BATCH,Tumor.no.yes,Ascites.no.yes, Sex,EP_weight,rel_EP_weight)%>%summarize(wks_diet=max(wks_diet))

do_organ_weight(dd, value = "EP_weight", batch = "ALL", sex = "both", y_title = "Weight at Endpoint [g]",path_for_saving_images)
do_organ_weight(dd, value = "EP_weight", batch = "ALL", sex = "female", y_title = "Weight at Endpoint [g]",path_for_saving_images)
do_organ_weight(dd, value = "EP_weight", batch = "ALL", sex = "male", y_title = "Weight at Endpoint [g]",path_for_saving_images)

do_organ_weight(dd, value = "rel_EP_weight", batch = "ALL", sex = "both", y_title = "rel.Weight at Endpoint [%]",path_for_saving_images)
do_organ_weight(dd, value = "rel_EP_weight", batch = "ALL", sex = "female", y_title = "rel.Weight at Endpoint [%]",path_for_saving_images)
do_organ_weight(dd, value = "rel_EP_weight", batch = "ALL", sex = "male", y_title = "rel.Weight at Endpoint [%]",path_for_saving_images)



# Plot Serum SAA ----
## Summarize for Plotting -----
stats_SAA <- data %>%
  group_by(Treatment,wks_diet) %>%
  summarise(mean_SAA = mean(SAA, na.rm = TRUE),
            sd_SAA   = sd(SAA, na.rm = TRUE),
            n=sum(!is.na(SAA)),
            .groups = "drop")%>%filter(!is.na(mean_SAA))
## Statistical Analysis -----
#Fit the linear mixed-effects model
model_SAA <- lmer(SAA ~ Treatment * as.factor(wks_diet) + (1 | Animal), data = data) #wks diet needed to be facotr so that it gives analysis for all timpioint and not mean
anova_table_SAA<-anova(model_SAA, type = 3)
anova_label_SAA <- paste0("ANOVA over linear mixed-effects model\n",
                        "Treatment: F = ", round(anova_table_SAA$F[1], 2), ", p = ", signif(anova_table_SAA$`Pr(>F)`[1], 3), "\n",
                        "Time: F = ", round(anova_table_SAA$F[2], 2), ", p < ", format.pval(anova_table_SAA$`Pr(>F)`[2], digits = 1), "\n",
                        "Interaction: F = ", round(anova_table_SAA$F[3], 2), ", p = ", signif(anova_table_SAA$`Pr(>F)`[3], 3))

emm_SAA <- emmeans(model_SAA, ~ Treatment | wks_diet) # Calculate the estimated marginal means for both Treatment and wks_diet
pwc_SAA <- contrast(emm_SAA, method = "pairwise", adjust = "bonferroni") # Perform pairwise contrasts between Treatment levels (ctrl vs TAM) at all levels of wks_diet
posthoc_label_SAA<- "Post Hoc: Pairwise with Bonferroni correction"

pwc_df_SAA <- as.data.frame(pwc_SAA)
pwc_df_rounded_SAA <- pwc_df_SAA %>%
  mutate(rounded_p_value = ifelse(is.na(p.value), "NA", round(p.value, 3))) %>%
  mutate(significance = case_when(
    is.na(p.value) ~ "NA",                     # For NA p-values
    p.value < 0.001 ~ "***",                   # p < 0.001 is highly significant
    p.value >= 0.001 & p.value < 0.01 ~ "**",  # 0.001 ≤ p < 0.01 is significant
    p.value >= 0.01 & p.value < 0.05 ~ "*",    # 0.01 ≤ p < 0.05 is moderately significant
    p.value >= 0.05 ~ "NS",                    # p ≥ 0.05 is not significant
    TRUE ~ "NA"   ))  %>%                      # Default case
  select(wks_diet, rounded_p_value, significance)

## Plot Serum SAA overtime -----
plot <- ggplot(data = stats_SAA,aes(x = wks_diet, y = mean_SAA, color = Treatment, fill = Treatment)) +
  geom_ribbon(data = stats_SAA, aes(x = wks_diet, ymin = mean_SAA-sd_SAA,  ymax = mean_SAA + sd_SAA,  fill = Treatment,  group = Treatment), alpha = 0.1, linetype = 0)+
  geom_line(data = stats_SAA,   aes(x = wks_diet, y = mean_SAA, color = Treatment, group = Treatment), linewidth = 1) +
  geom_point(data = stats_SAA,  aes(x = wks_diet, y = mean_SAA, color = Treatment),  size = 3, stroke = 1.1) +
  geom_text(aes(y = c(mean_SAA[4]-20,mean_SAA[5:6]+20 ,mean_SAA[4]-20,mean_SAA[5:6]+20), label = n), position = position_dodge(width = 0.4), hjust = 0.5, size = 3, show.legend = FALSE)+
  scale_color_manual(values = c("grey30", "palevioletred4")) +
  scale_fill_manual (values = c("grey30", "palevioletred4")) +
  scale_y_continuous(limits = c(-15,550), breaks = seq(0, 550, by = 100)) +
  scale_x_continuous(name = "Time on CD-HFD [wks]", limits = c(-1.2,12), breaks =  seq(-1, 12, by = 1)) +
  ylab("Serum SAA [pg/mL]") +
  theme_bw() +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.ticks.length = unit(4, "pt"),
        legend.position = "right")+
  ggtitle("Serum SAA levels") +
  guides(x = guide_axis(cap = "upper", minor.ticks = FALSE), y = guide_axis(cap = "upper")) +
  theme_bw() +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.ticks.length = unit(4, "pt"))+
  annotate("text", x = pwc_df_rounded_SAA$wks_diet,
           y = c(stats_SAA$mean_SAA[4]-35,stats_SAA$mean_SAA[5:6]+30),   label = pwc_df_rounded_SAA$significance,
           size = 2.5, color = "black",    fontface = "italic") +
  annotate("text", x=0.3, y = 200,   label = anova_label_SAA,     size = 2,   hjust = 0,    color = "black",   fontface = "italic") +
  annotate("text", x=0.3, y = 150, label = paste0(posthoc_label_SAA,"\n*** = p 0.001"),   size = 2,   hjust = 0,    color = "black",   fontface = "italic")

# Saving Plot --- 
ggsave(filename = "FK49_SAA_ELISA.png", plot = plot,  path = "02_GeneratedData/", width = 9, height = 6,dpi = 300)

rm(emm_SAA,model_SAA,pwc_df_SAA,pwc_df_rounded_SAA,pwc_SAA,anova_table_SAA,anova_label_SAA,posthoc_label_SAA,stats_SAA)

