rm(list=ls())
gc()
#### hier weiter machen, um SAA Statsiticen als csv zu speichern
library(tidyverse)
library(lmerTest)
library(emmeans)
source("FK49_Definitions.R")
# Read Raw Inputdata after general Data manipulation ------------------------------------------------------
output_pwd <-file.path(PATHS$Organs$FK49_output)
load(file.path(PATHS$Organs$FK49_input,"FK49_Data_prepared.Rda"))
#plot Serum SAA ----
## Summarize for Plotting -----
stats_SAA <- data %>% 
  group_by(Treatment,wks_diet) %>%
  summarise(mean_SAA = mean(SAA, na.rm = TRUE),
            sd_SAA   = sd(SAA, na.rm = TRUE),
            n=sum(!is.na(SAA)),
            .groups = "drop")%>%filter(!is.na(mean_SAA))%>%
  mutate(Method= "Summary")
## Statistical Analysis -----
#Fit the linear mixed-effects model
data<-data%>%select(Animal,Treatment,Sex,SAA,wks_diet)%>%mutate(wks_diet_factor=as.factor(wks_diet)) #wks diet needed to be facotr so that it gives analysis for all timpioint and not mean
#model_SAA <- lmer(SAA ~ Treatment * wks_diet_factor + (1 +wks_diet_factor |Animal), data = data) # rando-effects parameters and residual variance are probably unidentifiable
model_SAA <- lmer(SAA ~ Treatment * wks_diet_factor + (1 |Animal), data = data) #wks diet needed to be facotr so that it gives analysis for all timpioint and not mean

anova_table_SAA <- anova(model_SAA, type = 3)%>%mutate(Method= "ANOVA")%>%rownames_to_column("Effects")

anova_label_SAA <- paste0("ANOVA over linear mixed-effects model\n",
                        "Treatment: F = ", round(anova_table_SAA$F[1], 2), ", p = ", signif(anova_table_SAA$`Pr(>F)`[1], 3), "\n",
                        "Time: F = ", round(anova_table_SAA$F[2], 2), ", p < ", format.pval(anova_table_SAA$`Pr(>F)`[2], digits = 1), "\n",
                        "Interaction: F = ", round(anova_table_SAA$F[3], 2), ", p = ", signif(anova_table_SAA$`Pr(>F)`[3], 3))

emm_SAA <- emmeans(model_SAA, ~ Treatment | wks_diet_factor) # Calculate the estimated marginal means for both Treatment and wks_diet
pwc_SAA <- contrast(emm_SAA, method = "pairwise", adjust = "bonferroni") # Perform pairwise contrasts between Treatment levels (ctrl vs TAM) at all levels of wks_diet
posthoc_label_SAA<- "Post Hoc: Pairwise with Bonferroni correction"

pwc_df_SAA <- as.data.frame(pwc_SAA)%>%mutate(Method= "PostHoc")
pwc_df_rounded_SAA <- pwc_df_SAA %>%
  mutate(rounded_p_value = ifelse(is.na(p.value), "NA", round(p.value, 3))) %>%
  mutate(significance = case_when(
    is.na(p.value) ~ "NA",                     # For NA p-values
    p.value < 0.001 ~ "***",                   # p < 0.001 is highly significant
    p.value >= 0.001 & p.value < 0.01 ~ "**",  # 0.001 ≤ p < 0.01 is significant
    p.value >= 0.01 & p.value < 0.05 ~ "*",    # 0.01 ≤ p < 0.05 is moderately significant
    p.value >= 0.05 ~ "NS",                    # p ≥ 0.05 is not significant
    TRUE ~ "NA"   ))  %>%                      # Default case
  select(wks_diet_factor, rounded_p_value, significance)%>%
  mutate(wks_diet = as.numeric(as.character(wks_diet_factor)))

## Plot Serum SAA overtime -----

plot <- ggplot(data = stats_SAA,aes(x = wks_diet, y = mean_SAA, color = Treatment, fill = Treatment)) +
  geom_ribbon(aes(x = wks_diet, ymin = mean_SAA-sd_SAA,  ymax = mean_SAA + sd_SAA,  group = Treatment), alpha = 0.1, linetype = 0)+
  geom_line(aes(x = wks_diet, y = mean_SAA, group = Treatment), linewidth = 1) +
  geom_point(aes(x = wks_diet, y = mean_SAA, ),  size = 3, stroke = 0.2) +
  geom_text(aes(y = c(mean_SAA[4]-20,mean_SAA[5:6]+20 ,mean_SAA[4]-20,mean_SAA[5:6]+20), label = n), position = position_dodge(width = 0.4), hjust = 0.5, size = 3, show.legend = FALSE)+
  scale_color_manual(values = Treatment_colors[c("Ctrl","TAM")]) +
  scale_fill_manual (values = Treatment_colors[c("Ctrl","TAM")]) +
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
  annotate("text", x=0.3, y = 150, label = paste0(posthoc_label_SAA,"\n*** < p 0.001"),   size = 2,   hjust = 0,    color = "black",   fontface = "italic")

# Saving Plot --- 
ggsave(filename = "FK49_SAA_ELISA.png", plot = plot,  path = output_pwd, width = 9, height = 6,dpi = 300)
str(stats_SAA)
str(data)

stats_SAA<-stats_SAA%>%mutate(wks_diet=as.factor(wks_diet))
ALL_Statistics<-bind_rows(stats_SAA,anova_table_SAA,pwc_df_SAA)%>%select(Method,colnames(stats_SAA),colnames(anova_table_SAA),colnames(pwc_df_SAA))
write.csv2(ALL_Statistics,file=file.path(output_pwd,"Statistics/FK49_SAA_Stats.csv"))


rm(list=ls())
gc()

