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
# ### 10.08.26 hier morgen weiter machen
# #statisti is noch wild. ich möchte nur ein linear model mit allen nötigen effetcs drin haben.
#  dann im interaction sehen, was important ist und von da entscheidne, ob einsummary plot oder verschiedenen plots übre zeit.
#  Food_consumed ~ days_diet * Treatment * Sex + Batch + (1 + days_diet| Cage)
#  #oder wks diet?


output_pwd <-file.path(PATHS$general_data$FK49_output)
load(file.path(output_pwd,"01_RawData/FK49_Data_prepared.Rda"))
exigo <- c("ALB", "TP", "GLOB","A.G", "TB", "GGT", "AST", "ALT", "ALP", "AMY","Crea","UA","BUN","GLU","TC","TG")
output_pwd <-file.path(PATHS$general_data$FK49_output)
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
rm(d_for_food)
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


## Can I combine BAtches 1 and 2 ------
###  -------
path_for_output=file.path(output_pwd,"02_GeneratedData/FoodIntake/Background/Test_Batches")
### Visual normality check -------
ShapiroF<-shapiro.test(d_food_cagewise$Food_consumed)#normalverteilt
ShapiroW<-shapiro.test(d_food_cagewise$Water_consumed)#normalverteilt
p1<-ggplot(d_food_cagewise, aes(sample = Food_consumed)) +   stat_qq() + stat_qq_line() + facet_wrap(~BATCH)
p2<-ggplot(d_food_cagewise, aes(sample = Water_consumed)) +   stat_qq() + stat_qq_line() + facet_wrap(~BATCH)
ggsave(p1,file=file.path(path_for_output,"QQ_Food.png"),dpi= 300, width=5, height=5)
ggsave(p2,file=file.path(path_for_output,"QQ_Water.png"),dpi= 300, width=5, height=5)
rm(p1,p2)
### Check if unbalanced or missing data ----
summary(d_food_cagewise)
any(is.na(d_food_cagewise$Food_consumed)) #FALSE
any(is.na(d_food_cagewise$Water_consumed))#FALSE
table(d_food_cagewise$BATCH, d_food_cagewise$Treatment) 
#ctrl1 80 TAM1 78 Ctrl2 16, TAM2 47 might be unbalance in batch 2

### Test for variance -------
LeveneF<-leveneTest(Food_consumed ~ BATCH, data = d_food_cagewise) #Variance is the same between batches
LeveneW<-leveneTest(Water_consumed ~ BATCH, data = d_food_cagewise)#Variance is the same between batches, but smaller p


###Is food consumed different across batches (independent of sex, treatment, time)  -------
TF_B <- t.test( Food_consumed~BATCH,var.equal= TRUE,data=d_food_cagewise,alternative = "two.sided") #on this level batch has no effect
WF_B <- wilcox.test( Food_consumed~BATCH,var.equal= TRUE,data=d_food_cagewise,exact = FALSE, correct = FALSE, conf.int = FALSE) #on this level batch has no effect
TW_B <- t.test( Water_consumed~BATCH,var.equal= TRUE,data=d_food_cagewise,alternative = "two.sided") #on this level batch might have effect

###Is food consumed different across sexes (independent of batch, treatment, time)  -------
TF_S <- t.test( Food_consumed~Sex,var.equal= TRUE,data=d_food_cagewise,alternative = "less")#I assume males should eat more than females,on this level sex has no effect
WF_S <- wilcox.test( Food_consumed~Sex,var.equal= TRUE,data=d_food_cagewise,exact = FALSE, correct = FALSE, conf.int = FALSE) #on this level batch has no effect
TW_S <- t.test( Water_consumed~Sex,var.equal= TRUE,data=d_food_cagewise,alternative = "less") #I assume males shoul eat more than females,on this level sex has no effect
pF<-ggplot(d_food_cagewise, aes(x = Treatment, y = Food_consumed, fill = BATCH)) +  geom_boxplot(position = position_dodge(0.8)) +theme_bw()
pW<-ggplot(d_food_cagewise, aes(x = Treatment, y = Water_consumed, fill = BATCH)) +  geom_boxplot(position = position_dodge(0.8)) +  theme_bw()
ggsave(filename = "FK49_Food_Batch_Effect.png", plot = pF, path = path_for_output, width = 4.5, height = 9, dpi = 300)
ggsave(filename = "FK49_Water_Batch_Effect.png", plot = pW, path = path_for_output, width = 4.5, height = 9, dpi = 300)


# Convert htest objects to tidy data frames
TF_B <- broom::tidy(TF_B) %>% mutate(table = "t.test for Batch")
WF_B <- broom::tidy(WF_B) %>% mutate(table = "Wilcoxon for Batch")
WF_S <- broom::tidy(WF_S) %>% mutate(table = "Wilcoxon for Sex")
TF_S <- broom::tidy(TF_S) %>% mutate(table = "t.test for Sex")
TW_B <- broom::tidy(TW_B) %>% mutate(table = "t.test for Batch")
TW_S <- broom::tidy(TW_S) %>% mutate(table = "t.test for Sex")

# Convert Levene and outlier results to data frames if needed
LeveneF <- as.data.frame(LeveneF) %>% mutate(table = "Levene")
LeveneW <- as.data.frame(LeveneW) %>% mutate(table = "Levene")
outliersF <- as.data.frame(outliers_food) %>% mutate(table = "Outliers")
outliersW <- as.data.frame(outliers_water) %>% mutate(table = "Outliers")
ShapiroF<-broom::tidy(ShapiroF) %>% mutate(table = "Shapiro for Normality")
ShapiroW<-broom::tidy(ShapiroW) %>% mutate(table = "Shapiro for Normality")

# Combine results
StatsOutputF <- bind_rows(outliersF, LeveneF, TF_B, WF_B,WF_S, TF_S) %>% relocate(table)
StatsOutputW <- bind_rows(outliersW, LeveneW, TW_B, TW_S) %>% relocate(table)

# Save to CSV
write.csv2( StatsOutputF, file = paste0(path_for_output,"/FK49_Statistics_Food_BatchTest.csv"),
            row.names = FALSE,na = "",fileEncoding = "UTF-8")

write.csv2(StatsOutputW,file = paste0(path_for_output,"/FK49_Statistics_Water_BatchTest.csv"),
           row.names = FALSE,na = "",fileEncoding = "UTF-8")

rm(pF,pW,outliersF,outliersW,outliers_water,outliers_food,LeveneW,LeveneF, TF_B, WF_B, TF_S,TW_B, TW_S,WF_S,StatsOutputF,StatsOutputW,ShapiroF,ShapiroW)
###I decided to analyse batch 1 and 2 together
#although methodology was changed (1 animal per cage vs 1-5 animals per cage)

## Food / Water Summarized Data -------
###Test Normalverteilung statistisch und visuell -------
#I summarize repeated measures for same cages over time to not inflate my statistical power with this type of visualisation and test
d_for_boxplot <- d_food_cagewise %>%
  group_by(Cage, Treatment,Sex) %>%  # group by cage AND treatment
  summarize(mean_Food = mean(Food_consumed, na.rm = TRUE),
            mean_Water = mean(Water_consumed, na.rm = TRUE),.groups = "drop")

stats_FW <- d_for_boxplot %>%
  group_by(Treatment) %>%
  summarise(Mean_F = mean(d_for_boxplot$mean_Food, na.rm = TRUE),SD_F = sd(d_for_boxplot$mean_Food, na.rm = TRUE),
            Mean_W = mean(d_for_boxplot$mean_Water, na.rm = TRUE),  SD_W = sd(d_for_boxplot$mean_Water, na.rm = TRUE),.groups = "drop" )

### FOOD data  -------
#### Normality  -------
ShapiroF<- shapiro.test(d_for_boxplot$mean_Food) #nicht normalverteilt
hist(d_for_boxplot$mean_Food)
qqline(d_for_boxplot$mean_Food)

#### Difference  -------
#t.test( mean_Food~Treatment,var.equal= TRUE,data=d_for_boxplot,alternative = "two.sided") #
WF<-wilcox.test(mean_Food~Treatment, data = d_for_boxplot,exact = FALSE, correct = FALSE, conf.int = FALSE) #
#### Plot and Save summarized data -------

pF<- ggplot(stats_FW, aes(x = Treatment, y = Mean_F, fill = Treatment)) +
  geom_bar(stat = "identity", color = "black", alpha = 0.5, width = 0.75, position = "dodge") +
  geom_point(data = d_for_boxplot, fill = "lightgrey", color = "black",shape=21,aes(y = mean_Food),   position = position_jitterdodge(0.1, dodge.width = 0.75), alpha = 0.8, size = 5.3, stroke = 1.8) +
  geom_errorbar(aes(ymin = Mean_F - SD_F, ymax = Mean_F + SD_F),  position = position_dodge(width = 0.75), width = 0.2) +
  scale_fill_manual(values = c("black","darkred"), labels = c("Ctrl", "TAM")) +
  scale_y_continuous(name = paste0("Mean Daily Food Intake per animal [g]"),limits = c(0,20), breaks = seq(0, 20, by = 2)) +
  labs(x = "Treatment") +
  annotate("text", x = 1.5, y = max(d_for_boxplot$mean_Food, na.rm = TRUE) * 1.1, label = "ns", size = 6, fontface = "italic") +
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
  guides( shape = guide_legend(title = "Status", nrow = 1), fill = guide_legend(title = "Treatment", override.aes = list(shape = 21), nrow = 1), color = "none")
ggsave(filename = "FK49_Food Consumption_Summary_ALL.png", plot = pF, path = "02_GeneratedData/FoodIntake", width = 4, height = 11, dpi = 300)


### Water data  -------
#### Normality  -------
ShapiroW<-shapiro.test(d_for_boxplot$mean_Water)#normalverteilt
hist(d_for_boxplot$mean_Water)
qqline(d_for_boxplot$mean_Water)
#### Difference  -------  
TW<- t.test( mean_Water~Treatment,var.equal= TRUE,data=d_for_boxplot,alternative = "two.sided") # 
#### Plot and Save summarized data -------  
pW<- ggplot(stats_FW, aes(x = Treatment, y = Mean_W, fill = Treatment)) +
  geom_bar(stat = "identity", color = "black", alpha = 0.5, width = 0.75, position = "dodge") +
  geom_point(data = d_for_boxplot, fill = "lightgrey", color = "black",shape=21,aes(y = mean_Water),   position = position_jitterdodge(0.1, dodge.width = 0.75), alpha = 0.8, size = 5.3, stroke = 1.8) +
  geom_errorbar(aes(ymin = Mean_W - SD_W, ymax = Mean_W + SD_W),  position = position_dodge(width = 0.75), width = 0.2) +
  scale_fill_manual(values = c("black","darkred"), labels = c("Ctrl", "TAM")) +
  scale_y_continuous(name = paste0("Mean Daily Water Intake per animal [g]"), limits = c(0,16), breaks = seq(0, 16, by = 2)) +
  labs(x = "Treatment") +
  annotate("text", x = 1.5, y = max(d_for_boxplot$mean_Food, na.rm = TRUE) * 1.1, label = "ns", size = 6, fontface = "italic") +
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
  guides( shape = guide_legend(title = "Status", nrow = 1), fill = guide_legend(title = "Treatment", override.aes = list(shape = 21), nrow = 1),            color = "none")
ggsave(filename = "FK49_Water Consumption_Summary_ALL.png", plot = pW, path = "02_GeneratedData/FoodIntake", width = 4, height = 11, dpi = 300)

rm(d_for_boxplot,stats_FW,pF,pW) 
## Food and water consumption per block ----
d_for_Block <- d_food_cagewise %>%
  group_by(Cage, Treatment, Sex, Block, BATCH) %>%  
  summarise(mean_Food = mean(Food_consumed, na.rm = TRUE),
            mean_Water = mean(Water_consumed, na.rm = TRUE),
            n = n(),  .groups = "drop") %>%
  mutate( wks_diet = case_when(
    Block == "1"  ~ 0.3575,
    Block == "4"  ~ 3.3575,
    Block == "8"  ~ 7.3575,
    Block == "12" ~ 11.3575,
    TRUE ~ NA_real_ ) )

stats2_FW <- d_for_Block %>%
  group_by(Treatment,wks_diet) %>%
  summarise(Mean_F = mean(mean_Food, na.rm = TRUE),
            SD_F   = sd(mean_Food, na.rm = TRUE),
            Mean_W = mean(mean_Water, na.rm = TRUE),
            SD_W   = sd(mean_Water, na.rm = TRUE), n=n(),
            .groups = "drop")
### Food ----
# Fit the linear mixed-effects model 
model_F <- lmer(mean_Food ~ Treatment * Block + (1 | Cage), data = d_for_Block)
anova_table_F<-anova(model_F, type = 3)
anova_label_F <- paste0("ANOVA over linear mixed-effects model\n",
                        "Treatment: F = ", round(anova_table_F$F[1], 2), ", p = ", signif(anova_table_F$`Pr(>F)`[1], 3), "\n",
                        "Time: F = ", round(anova_table_F$F[2], 2), ", p < ", format.pval(anova_table_F$`Pr(>F)`[2], digits = 1), "\n",
                        "Interaction: F = ", round(anova_table_F$F[3], 2), ", p = ", signif(anova_table_F$`Pr(>F)`[3], 3))


# Calculate the estimated marginal means for both Treatment and wks_diet
emm_F <- emmeans(model_F, ~ Treatment | Block)

# Perform pairwise contrasts between Treatment levels (ctrl vs TAM) at all levels of wks_diet
pwc_F <- contrast(emm_F, method = "pairwise", adjust = "bonferroni")
posthoc_label_F<- "Post Hoc: Pairwise with Bonferroni correction"

pwc_df_F <- as.data.frame(pwc_F)
pwc_df_rounded_F <- pwc_df_F %>%
  mutate(rounded_p_value = ifelse(is.na(p.value), "NA", round(p.value, 3))) %>%
  mutate(significance = case_when(
    is.na(p.value) ~ "NA",                     # For NA p-values
    p.value < 0.001 ~ "***",                   # p < 0.001 is highly significant
    p.value >= 0.001 & p.value < 0.01 ~ "**",  # 0.001 ≤ p < 0.01 is significant
    p.value >= 0.01 & p.value < 0.05 ~ "*",    # 0.01 ≤ p < 0.05 is moderately significant
    p.value >= 0.05 ~ "NS",                    # p ≥ 0.05 is not significant
    TRUE ~ "NA"   ))  %>%                      # Default case
  select(Block, rounded_p_value, significance)


plot <- ggplot(data = stats2_FW, aes(x = wks_diet, y = Mean_F, color = Treatment, fill = Treatment)) +
  geom_ribbon(data = stats2_FW,  aes(x = wks_diet, ymin = Mean_F,  ymax = Mean_F + SD_F,  fill = Treatment,  group = Treatment), alpha = 0.1, linetype = 0)+
  geom_line(data = stats2_FW,   aes(x = wks_diet, y = Mean_F, color = Treatment, group = Treatment), linewidth = 1) +
  geom_point(data = stats2_FW,aes(x = wks_diet, y = Mean_F, color = Treatment),  size = 3, stroke = 1.1) +
  geom_text(aes(y = c(Mean_F[1:4]+1.8,Mean_F[1:4]+1.8) , label = n), position = position_dodge(width = 0.4), hjust = 0.5, size = 3, show.legend = FALSE)+
  scale_color_manual(values = c("grey30", "palevioletred4")) +
  scale_fill_manual (values = c("grey30", "palevioletred4")) +
  scale_y_continuous(limits = c(0,20), breaks = seq(0, 20, by = 2)) +
  scale_x_continuous(name = "Time on CD-HFD [wks]", limits = c(0,12), breaks =  seq(0, 12, by = 1)) +
  ylab("Mean Daily Food Intake per animal [g]") +
  theme_bw() +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.ticks.length = unit(4, "pt"),
        legend.position = "right")+
  ggtitle("Food Consumption Batch1+2 all sexes, per cage animal mean measurments per block") +
  guides(x = guide_axis(cap = "upper", minor.ticks = FALSE), y = guide_axis(cap = "upper")) +
  theme_bw() +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.ticks.length = unit(4, "pt"))+
  annotate("text", x = as.numeric(as.character(pwc_df_rounded_F$Block))-0.6425,
           y = stats2_FW$Mean_F[1:4]+1.3,   label = pwc_df_rounded_F$significance,
           size = 2.5, color = "black",    fontface = "italic") +
  annotate("text", x=0.3, y = 4,   label = anova_label_F,     size = 2,   hjust = 0,    color = "black",   fontface = "italic") +
  annotate("text", x=0.3, y = 2.8, label = posthoc_label_F,   size = 2,   hjust = 0,    color = "black",   fontface = "italic")

# Saving Plot ---
ggsave(filename = "FK49_Food_Consumption_Summary_Block.png", plot = plot,  path = "02_GeneratedData/FoodIntake", width = 9, height = 6,dpi = 300)

rm(emm_F,model_F,pwc_df_F,pwc_df_rounded_F,pwc_F,anova_table_F,anova_label_F,posthoc_label_F)

### Water ----
# Fit the linear mixed-effects model 
model_W <- lmer(mean_Water ~ Treatment * Block + (1 | Cage), data = d_for_Block)
anova_table_W<-anova(model_W, type = 3)
anova_label_W <- paste0("ANOVA over linear mixed-effects model\n",
                        "Treatment: F = ", round(anova_table_W$F[1], 2), ", p = ", signif(anova_table_W$`Pr(>F)`[1], 3), "\n",
                        "Time: F = ", round(anova_table_W$F[2], 2), ", p < ", format.pval(anova_table_W$`Pr(>F)`[2], digits = 1), "\n",
                        "Interaction: F = ", round(anova_table_W$F[3], 2), ", p = ", signif(anova_table_W$`Pr(>F)`[3], 3))


# Calculate the estimated marginal means for both Treatment and wks_diet
emm_W <- emmeans(model_W, ~ Treatment | Block)

# Perform pairwise contrasts between Treatment levels (ctrl vs TAM) at all levels of wks_diet
pwc_W <- contrast(emm_W, method = "pairwise", adjust = "bonferroni")
posthoc_label_W<- "Post Hoc: Pairwise with Bonferroni correction"

pwc_df_W <- as.data.frame(pwc_W)
pwc_df_rounded_W <- pwc_df_W %>%
  mutate(rounded_p_value = ifelse(is.na(p.value), "NA", round(p.value, 3))) %>%
  mutate(significance = case_when(
    is.na(p.value) ~ "NA",                     # For NA p-values
    p.value < 0.001 ~ "***",                   # p < 0.001 is highly significant
    p.value >= 0.001 & p.value < 0.01 ~ "**",  # 0.001 ≤ p < 0.01 is significant
    p.value >= 0.01 & p.value < 0.05 ~ "*",    # 0.01 ≤ p < 0.05 is moderately significant
    p.value >= 0.05 ~ "NS",                    # p ≥ 0.05 is not significant
    TRUE ~ "NA"   ))  %>%                      # Default case
  select(Block, rounded_p_value, significance)


plot <- ggplot(data = stats2_FW, aes(x = wks_diet, y = Mean_W, color = Treatment, fill = Treatment)) +
  geom_ribbon(data = stats2_FW,  aes(x = wks_diet, ymin = Mean_W,  ymax = Mean_W + SD_W,  fill = Treatment,  group = Treatment), alpha = 0.1, linetype = 0)+
  geom_line(data = stats2_FW,   aes(x = wks_diet, y = Mean_W, color = Treatment, group = Treatment), linewidth = 1) +
  geom_point(data = stats2_FW,aes(x = wks_diet, y = Mean_W, color = Treatment),  size = 3, stroke = 1.1) +
  
  geom_text(aes(y = c(Mean_W[1:4]+0.75,Mean_W[1:4]+0.75) , label = n), position = position_dodge(width = 0.4), hjust = 0.5, size = 3, show.legend = FALSE)+
  scale_color_manual(values = c("grey30", "palevioletred4")) +
  scale_fill_manual(values = c("grey30", "palevioletred4")) +
  scale_y_continuous(limits = c(0,16), breaks = seq(0,16, by= 2)) +
  scale_x_continuous(name = "Time on CD-HFD [wks]", limits = c(0,12), breaks = seq(0,12, by= 1)) +
  ylab("Mean Daily Water Intake per animal [g]") +
  theme_bw() +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.ticks.length = unit(4, "pt"),
        legend.position = "right")+
  ggtitle("Water Consumption Batch1+2 all sexes, per cage animal mean measurments per block") +
  guides(x = guide_axis(cap = "upper", minor.ticks = FALSE), y = guide_axis(cap = "upper")) +
  theme_bw() +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.ticks.length = unit(4, "pt"))+
  annotate("text", x = as.numeric(as.character(pwc_df_rounded_W$Block))-0.6425,
           y = stats2_FW$Mean_W[1:4]+0.5,   label = pwc_df_rounded_W$significance,
           size = 2.5, color = "black",    fontface = "italic") +
  annotate("text", x=0.3, y = 4,   label = anova_label_W,     size = 2,   hjust = 0,    color = "black",   fontface = "italic") +
  annotate("text", x=0.3, y = 2.8, label = posthoc_label_W,   size = 2,   hjust = 0,    color = "black",   fontface = "italic")


# ---  Saving Plot --- 
ggsave(filename = "FK49_Water_Consumption_Summary_Block.png", plot = plot,  path = "02_GeneratedData/FoodIntake", width = 9, height = 6,dpi = 300)

rm(emm_W,model_W,pwc_df_W,pwc_df_rounded_W,pwc_W,anova_table_W,anova_label_W,posthoc_label_W)
rm(d_for_Block,stats2_FW)
# Function for Food Curves ----------------------------------------------
do_food_curve <- function(inputdata, value, value_label = NULL, unit = "g",
                          batch = "ALL", sex = "both",Block, N, path_images, 
                          ymin = NULL, ymax = NULL,offset=1, savestats = "NO"){
  # --- Setup ---
  value_label_final <- if (is.null(value_label)) deparse(substitute(value)) else value_label
  file_base <- paste0("FK49_", value_label_final, "_Batch", batch, "_", sex,"_B",Block,"_n", N)
  
  # --- Filter by sex ---
  filtered <- inputdata %>%
    filter(!is.na({{value}})) %>%
    filter(case_when(
      sex == "female" ~ Sex == "female",
      sex == "male" ~ Sex == "male",
      sex == "EC1" ~ Sex == "EC1",
      sex == "EC2" ~ Sex == "EC2",
      sex == "EC" ~ Sex == "EC",
      sex == "both" ~ TRUE))
  
  # --- BATCH filtering ---
  if (batch == "ALL") {
    common_timepoints <- filtered %>%       # Find common timepoints across both batches
      filter(BATCH %in% c(1, 2)) %>%
      group_by(days_diet, BATCH) %>%
      summarise(n = n(), .groups = "drop") %>%
      group_by(days_diet) %>%
      summarise(n_batches = n_distinct(BATCH)) %>%
      filter(n_batches == 2) %>%
      pull(days_diet)
    filtered <- filtered %>%
      filter(BATCH %in% c(1, 2)) %>%
      filter(days_diet %in% common_timepoints) } 
  else {filtered <- filtered %>% filter(BATCH == batch) }
  print(filtered) 
  
  changed <- filtered %>%  filter(Block==!!Block) %>%
    group_by(Treatment, days_diet) %>%
    summarise(value_used = mean({{value}}, na.rm = TRUE),n = n(), sd = sd({{value}}, na.rm = TRUE)) %>%
    filter(n > N)
  
  ## Statistical Tests of Weight Curves --------------------------------------------------------
  value_str<-deparse(substitute(value))
  stat_tests <- filtered %>% convert_as_factor(Cage, days_diet)
  
  STATS <- tryCatch({ 
    stat_tests %>%group_by(Treatment, days_diet) %>%  get_summary_stats({{value}}, type = "mean_sd")}, error = function(e) 
    {  message("⚠️ Fehler bei get_summary_stats(): ", e$message)
      return(NULL)})
  
  outliers <- tryCatch({stat_tests %>%
      group_by(Treatment, days_diet) %>%
      identify_outliers({{value}})}, error = function(e) {message("⚠️ Fehler bei identify_outliers(): ", e$message)
        return(NULL)})
  
  shapiro_results <- tryCatch({ stat_tests %>%
      group_by(Treatment, days_diet) %>%
      summarise(
        n = sum(!is.na({{value}})),
        all_same = all({{value}} == {{value}}[1], na.rm = TRUE),
        shapiro_p = ifelse(n >= 3 && !all_same,
                           tryCatch({shapiro.test({{value}})$p.value }, error = function(e) NA_real_), NA_real_),
        result = case_when(
          is.na(shapiro_p) ~ "not tested",
          shapiro_p >= 0.05 ~ "normal",
          TRUE ~ "non-normal"),
        note = ifelse(n >= 3, "OK", "Too few samples"), .groups = "drop" )}, error = function(e) {
          message("⚠️ Fehler bei shapiro.test(): ", e$message)
          return(NULL)})
  
  #print(shapiro_results,n=40)
  filtered$days_diet <- factor(filtered$days_diet)
  
  # --- Fit the linear mixed-effects model --- 
  model <- lmer(as.formula(paste(deparse(substitute(value)), "~ Treatment * days_diet + (1 | Cage)")), data = filtered)
  model <- lmer(mean_Food ~ Treatment * wks_diet + Sex + (1 + wks_diet | Cage),data = d_food_week, REML = TRUE )
  anova_table<-anova(model, type = 3)
  
  anova_label <- paste0("ANOVA over linear mixed-effects model\n",
                        "Treatment: F = ", round(anova_table$F[1], 2), ", p = ", signif(anova_table$`Pr(>F)`[1], 3), "\n",
                        "Time: F = ", round(anova_table$F[2], 2), ", p < ", format.pval(anova_table$`Pr(>F)`[2], digits = 1), "\n",
                        "Interaction: F = ", round(anova_table$F[3], 2), ", p = ", signif(anova_table$`Pr(>F)`[3], 3))
  
  
  # --- Calculate the estimated marginal means for both Treatment and days_diet --- 
  emm <- emmeans(model, ~ Treatment | days_diet)
  
  # --- Perform pairwise contrasts between Treatment levels (ctrl vs TAM) at all levels of days_diet --- 
  pwc <- contrast(emm, method = "pairwise", adjust = "bonferroni")
  posthoc_label<- "Post Hoc: Pairwise with Bonferroni correction"
  
  pwc_df <- as.data.frame(pwc)
  pwc_df_rounded <- pwc_df %>%
    mutate(rounded_p_value = ifelse(is.na(p.value), "NA", round(p.value, 3))) %>%
    mutate(days_diet = as.character(days_diet)) %>%
    mutate(days_diet = as.numeric(days_diet)) %>%
    mutate(significance = case_when(
      is.na(p.value) ~ "NA",                     # For NA p-values
      p.value < 0.001 ~ "***",                   # p < 0.001 is highly significant
      p.value >= 0.001 & p.value < 0.01 ~ "**", # 0.001 ≤ p < 0.01 is significant
      p.value >= 0.01 & p.value < 0.05 ~ "*",    # 0.01 ≤ p < 0.05 is moderately significant
      p.value >= 0.05 ~ "NS",                    # p ≥ 0.05 is not significant
      TRUE ~ "NA"   ))  %>%                           # Default case
    select(days_diet, rounded_p_value, significance)
  
  
  # --- Get max weight per days_diet from 'changed' (across both Treatment groups) --- 
  max_values <- changed %>%
    group_by(days_diet) %>%
    summarise(max_value = max(value_used, na.rm = TRUE)) %>%
    ungroup() %>%
    mutate(y_position = max_value + max_value * 0.05)
  
  # --- Join y_position back into pwc_df_rounded --- 
  pwc_df_annotated <- pwc_df_rounded %>%
    left_join(max_values, by = "days_diet")
  changed <- changed %>%left_join(pwc_df_annotated %>% select(days_diet, y_position), by = "days_diet")
  # Variables for Plot Setup
  # Calculate default y limits only if not provided
  if (is.null(ymin) || is.null(ymax)) {
    mean_value <- mean(changed$value_used, na.rm = TRUE)
    sd_value <- sd(changed$value_used, na.rm = TRUE)
    min_value_default <- round(mean_value - 3 * sd_value)
    max_value_default <- round(mean_value + 3 * sd_value) + 4
    min_value <- ifelse(is.null(ymin), min_value_default, ymin)
    max_value <- ifelse(is.null(ymax), max_value_default, ymax)} 
  else {min_value <- ymin
  max_value <- ymax }
  
  range_value <- max_value - min_value
  step_size <- ceiling((range_value * 0.2) / 2) *2
  breaks_value <- seq(min_value, max_value, by = step_size)
  breaks_value <- round(breaks_value /2) * 2
  
  min_x <- round(min(changed$days_diet, na.rm = TRUE))
  max_x <- round(max(changed$days_diet, na.rm = TRUE))
  breaks_x <- seq(min_x, max_x, by = 1)
  
  unit_label <- unit
  
  # Plot
  plot <- ggplot(data = changed, aes(x = days_diet, y = value_used, color = Treatment, fill = Treatment)) +
    geom_ribbon(aes(y = value_used, ymin = value_used, ymax = value_used + sd), alpha = 0.1, linetype = 0) +
    geom_point(data = d, fill = "lightgrey", color = "black",
               aes(y = .data[[value]], shape = event_status),
               position =  position_jitterdodge(0.1, dodge.width = 0.75), alpha = 0.8, size = 3, stroke = 1.1)+
    geom_line(linewidth = 1) +
    geom_text(aes(y = y_position + offset, label = n), position = position_dodge(width = 0.4), hjust = 0.5, size = 3, show.legend = FALSE)+
    scale_color_manual(values = c("grey30", "palevioletred4")) +
    scale_fill_manual(values = c("grey30","palevioletred4")) +
    scale_x_continuous(name = "Time on CD-HFD [days]",
                       limits = c(min_x-0.1, max_x+0.1),
                       breaks = breaks_x,
                       minor_breaks = seq(min_x, max_x, by = 1)) +
    scale_y_continuous(name = sprintf("%s [%s]", deparse(substitute(value)), unit_label),
                       limits = c(min_value, max_value),
                       breaks = breaks_value) +
    xlab("Time on CD-HFD [d]") +
    ylab(sprintf("%s [%s]", deparse(substitute(value)), unit_label)) +
    theme_bw() +
    #ggtitle(sprintf("%s of %ss \n from batch %s (n > %d)", value_label_final, sex, batch, N)) +
    guides(x = guide_axis(cap = "upper", minor.ticks = TRUE),
           y = guide_axis(cap = "upper")) +
    theme(axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank(),
          axis.ticks.length = unit(4, "pt"))+
    annotate("text", 
             x = pwc_df_annotated$days_diet, 
             y = pwc_df_annotated$y_position+offset/2,
             label = pwc_df_annotated$significance, 
             size = 2.5, 
             color = "black", 
             fontface = "italic")
  
  # ---  Saving Plot --- 
  ggsave(filename = paste0(file_base, ".png"), plot = plot,  path = path_images, width = 3, height = 4,dpi = 300)
  #ggsave(filename = paste0(file_base, ".pdf"), plot = plot, path = path_images, width = 9, height = 6, dpi = 300, device = cairo_pdf)
  #dev.off()
  # Optionally save stats tables
  if (savestats == "YES") {
    outliers <- mutate(outliers, table = "Outliers")
    shapiro_results <-mutate(shapiro_results, table = "Shapiro")
    anova_table <- mutate(anova_table, table = "ANOVA")
    pwc_df <- mutate(pwc_df, table = "Pairwise Comparison")
    StatsOutput <- bind_rows(outliers,shapiro_results, anova_table,pwc_df)%>% relocate( table)
    write.csv2(
      StatsOutput,
      file = file.path(paste0(path_images,"/Statistics"), paste0(file_base, "_StatsOutput.csv")),
      row.names = FALSE,
      na = "",
      fileEncoding = "UTF-8")
  }
  
  # --- Return output --- 
  return(list(
    outliers = outliers,
    shapiro_results = shapiro_results,
    anova_table = anova_table,
    posthoc = pwc_df,
    model= model,
    plot = plot
  ))
  
}

##
path_for_saving_images<-"02_GeneratedData/FoodIntake"
### Plot Food Consumption Block 1-12 Both Batches
do_food_curve(d_food_cagewise,value = Food_consumed,  value_label = "Food Consumption", unit = "g",  batch = "ALL",  sex = "both", Block= 1,  N = 0,path_images = path_for_saving_images,ymin=0, ymax=20,offset=1)
do_food_curve(d_food_cagewise,value = Food_consumed,  value_label = "Food Consumption", unit = "g",  batch = "ALL",  sex = "both", Block= 4,  N = 0,path_images = path_for_saving_images,ymin=0, ymax=20,offset=1)
do_food_curve(d_food_cagewise,value = Food_consumed,  value_label = "Food Consumption", unit = "g",  batch = "ALL",  sex = "both", Block= 8,  N = 0,path_images = path_for_saving_images,ymin=0, ymax=20,offset=1)
do_food_curve(d_food_cagewise,value = Food_consumed,  value_label = "Food Consumption", unit = "g",batch = "ALL",  sex = "both", Block= 12,  N = 0, path_images = path_for_saving_images,ymin=0, ymax=20,offset=1)
### Water Consumption Block 1-12 Both Batches##
#Water Consumption was different between batches but not between groups. Checked for single batches and did not find treatment depended differences.
#since no differences in groups decided to to visualization with both batches together despite change in methodology
do_food_curve(d_food_cagewise,value = Water_consumed,  value_label = "Water Consumption", unit = "g",batch = "ALL",  sex = "both", Block=1,  N = 0,path_images = path_for_saving_images,ymin=0, ymax=16,offset=1)
do_food_curve(d_food_cagewise,value = Water_consumed,  value_label = "Water Consumption", unit = "g",batch = "ALL",  sex = "both", Block=4,  N = 0,path_images = path_for_saving_images,ymin=0, ymax=16,offset=1)
do_food_curve(d_food_cagewise,value = Water_consumed,  value_label = "Water Consumption", unit = "g",batch = "ALL",  sex = "both", Block=8,  N = 0,path_images = path_for_saving_images,ymin=0, ymax=16,offset=1)
do_food_curve(d_food_cagewise,value = Water_consumed,  value_label = "Water Consumption", unit = "g",batch = "ALL",  sex = "both", Block=12,  N = 0,path_images = path_for_saving_images,ymin=0, ymax=16,offset=1)
gc()

path_for_saving_images<-"02_GeneratedData/FoodIntake/Background"
#1,4,8,12 relevant Blocks

### Food Consumption Block 1-12 Batch1
do_food_curve(d_food_cagewise,value = Food_consumed,  value_label = "Food Consumption", unit = "g",  batch = "1",  sex = "both", Block= 1,  N = 0,path_images = path_for_saving_images,ymin=0, ymax=16,offset=1)
do_food_curve(d_food_cagewise,value = Food_consumed,  value_label = "Food Consumption", unit = "g",  batch = "1",  sex = "both", Block= 4,  N = 0,path_images = path_for_saving_images,ymin=0, ymax=34,offset=2)
do_food_curve(d_food_cagewise,value = Food_consumed,  value_label = "Food Consumption", unit = "g",  batch = "1",  sex = "both", Block= 8,  N = 0,path_images = path_for_saving_images,ymin=0, ymax=16,offset=1)
do_food_curve(d_food_cagewise,value = Food_consumed,  value_label = "Food Consumption", unit = "g",batch = "1",  sex = "both", Block= 12,  N = 0, path_images = path_for_saving_images,ymin=0, ymax=15,offset=1)

### Food Consumption Block 1-8 Batch2 0 1 4 8 12
do_food_curve(d_food_cagewise,value = Food_consumed,  value_label = "Food Consumption", unit = "g",  batch = "2",  sex = "both", Block= 1,  N = 0,path_images = path_for_saving_images,ymin=0, ymax=35,offset=2)
do_food_curve(d_food_cagewise,value = Food_consumed,  value_label = "Food Consumption", unit = "g",  batch = "2",  sex = "both", Block= 4,  N = 0,path_images = path_for_saving_images,ymin=0, ymax=35,offset=2)
do_food_curve(d_food_cagewise,value = Food_consumed,  value_label = "Food Consumption", unit = "g",  batch = "2",  sex = "both", Block= 8,  N = 0,path_images = path_for_saving_images,ymin=0, ymax=30,offset=2)
do_food_curve(d_food_cagewise,value = Food_consumed,  value_label = "Food Consumption", unit = "g",batch = "2",  sex = "both", Block= 12,  N = 0, path_images = path_for_saving_images,ymin=0, ymax=15,offset=1)

gc()
### Food Consumption Block 1-12 Both female
do_food_curve(d_food_cagewise,value = Food_consumed,  value_label = "Food Consumption", unit = "g",  batch = "ALL",  sex = "female", Block= 1,  N = 0,path_images = path_for_saving_images,ymin=0, ymax=35,offset=2)
do_food_curve(d_food_cagewise,value = Food_consumed,  value_label = "Food Consumption", unit = "g",  batch = "ALL",  sex = "female", Block= 4,  N = 0,path_images = path_for_saving_images,ymin=0, ymax=35,offset=2)
do_food_curve(d_food_cagewise,value = Food_consumed,  value_label = "Food Consumption", unit = "g",  batch = "ALL",  sex = "female", Block= 8,  N = 0,path_images = path_for_saving_images,ymin=0, ymax=30,offset=2)
do_food_curve(d_food_cagewise,value = Food_consumed,  value_label = "Food Consumption", unit = "g",batch = "ALL",  sex = "female", Block= 12,  N = 0, path_images = path_for_saving_images,ymin=0, ymax=15,offset=1)

### Food Consumption Block 1-12 Both male
do_food_curve(d_food_cagewise,value = Food_consumed,  value_label = "Food Consumption", unit = "g",  batch = "ALL",  sex = "male", Block= 1,  N = 0,path_images = path_for_saving_images,ymin=0, ymax=35,offset=2)
do_food_curve(d_food_cagewise,value = Food_consumed,  value_label = "Food Consumption", unit = "g",  batch = "ALL",  sex = "male", Block= 4,  N = 0,path_images = path_for_saving_images,ymin=0, ymax=35,offset=2)
do_food_curve(d_food_cagewise,value = Food_consumed,  value_label = "Food Consumption", unit = "g",  batch = "ALL",  sex = "male", Block= 8,  N = 0,path_images = path_for_saving_images,ymin=0, ymax=30,offset=2)
do_food_curve(d_food_cagewise,value = Food_consumed,  value_label = "Food Consumption", unit = "g",batch = "ALL",  sex = "male", Block= 12,  N = 0, path_images = path_for_saving_images,ymin=0, ymax=15,offset=1)


### Water Consumption Block 1-12 Batch1

do_food_curve(d_food_cagewise,value = Water_consumed,  value_label = "Water Consumption", unit = "g",batch = "1",  sex = "both", Block=1,  N = 0,path_images = path_for_saving_images,ymin=0, ymax=45,offset=2)
do_food_curve(d_food_cagewise,value = Water_consumed,  value_label = "Water Consumption", unit = "g",batch = "1",  sex = "both", Block=4,  N = 0,path_images = path_for_saving_images,ymin=0, ymax=45,offset=2)
do_food_curve(d_food_cagewise,value = Water_consumed,  value_label = "Water Consumption", unit = "g",batch = "1",  sex = "both", Block=8,  N = 0,path_images = path_for_saving_images,ymin=0, ymax=45,offset=2)
do_food_curve(d_food_cagewise,value = Water_consumed,  value_label = "Water Consumption", unit = "g",batch = "1",  sex = "both", Block=12,  N = 0,path_images = path_for_saving_images,ymin=0, ymax=45,offset=2)

### Water Consumption Block 1-12 Batch2
do_food_curve(d_food_cagewise,value = Water_consumed,  value_label = "Water Consumption", unit = "g",batch = "2",  sex = "both", Block=1,  N = 0,path_images = path_for_saving_images,ymin=0, ymax=45,offset=2)
do_food_curve(d_food_cagewise,value = Water_consumed,  value_label = "Water Consumption", unit = "g",batch = "2",  sex = "both", Block=4,  N = 0,path_images = path_for_saving_images,ymin=0, ymax=45,offset=2)
do_food_curve(d_food_cagewise,value = Water_consumed,  value_label = "Water Consumption", unit = "g",batch = "2",  sex = "both", Block=8,  N = 0,path_images = path_for_saving_images,ymin=0, ymax=45,offset=2)
do_food_curve(d_food_cagewise,value = Water_consumed,  value_label = "Water Consumption", unit = "g",batch = "2",  sex = "both", Block=12,  N = 0,path_images = path_for_saving_images,ymin=0, ymax=45,offset=2)
rm(d_food_cagewise,do_food_curve)

# # Function for Food Curves ----------------------------------------------
# do_EC_curve <- function(inputdata, value, value_label = NULL, unit = "g",
#                           batch = "ALL", sex = "both",block_value = NULL, N, path_images, 
#                           ymin = NULL, ymax = NULL){
#   # --- Setup ---
#   value_label_final <- if (is.null(value_label)) deparse(substitute(value)) else value_label
#   file_base <- paste0("FK46_", value_label_final, "_Batch", batch, "_", sex,"_B",block_value,"_n", N)
#   
#   # --- Filter by sex ---
#   filtered <- inputdata %>%
#     filter(!is.na({{value}})) %>%
#     filter(case_when(
#       sex == "female" ~ Sex == "female",
#       sex == "male" ~ Sex == "male",
#       sex == "EC1" ~ Sex == "EC1",
#       sex == "EC2" ~ Sex == "EC2",
#       sex == "both" ~ TRUE))
#   
#   # --- BATCH filtering ---
#   if (batch == "ALL") {
#     common_timepoints <- filtered %>%       # Find common timepoints across both batches
#       filter(BATCH %in% c(1, 2)) %>%
#       group_by(DOW_in_block, BATCH) %>%
#       summarise(n = n(), .groups = "drop") %>%
#       group_by(days_diet) %>%
#       summarise(n_batches = n_distinct(BATCH)) %>%
#       filter(n_batches == 2) %>%
#       pull(DOW_in_block)
#     
#     filtered <- filtered %>%
#       filter(BATCH %in% c(1, 2)) %>%
#       filter(DOW_in_block %in% common_timepoints)
#   } 
#   else {
#     filtered <- filtered %>% filter(BATCH == batch)
#   }
#   print(filtered) 
#   changed <- filtered %>%
#     filter(if (is.null(block_value)) Block %in% c("0", "1", "2", "3") else Block == block_value) %>%
#     group_by(Diet, DOW_in_block) %>%
#     summarise(
#       weight = mean({{value}}, na.rm = TRUE),
#       n = n(),
#       sd = sd({{value}}, na.rm = TRUE),
#       .groups = "drop"
#     ) %>%
#     filter(n > N)
#   
#   ## Statistical Tests of Weight Curves --------------------------------------------------------
#   ### Stat Tests --------------------------------------------------------
#   value_str<-deparse(substitute(value))
# 
#   
#   if (nrow(changed) == 0) {
#     warning("No data available after filtering for plotting. Check filters (e.g., Block, N).")
#     return(NULL)
#   }
#   # Variables for Plot Setup
#   # Calculate default y limits only if not provided
#   if (is.null(ymin) || is.null(ymax)) {
#     mean_value <- mean(changed$weight, na.rm = TRUE)
#     sd_value <- sd(changed$weight, na.rm = TRUE)
#     min_value_default <- round(mean_value - 3 * sd_value)
#     max_value_default <- round(mean_value + 3 * sd_value) + 4
#     min_value <- ifelse(is.null(ymin), min_value_default, ymin)
#     max_value <- ifelse(is.null(ymax), max_value_default, ymax)} 
#   else {min_value <- ymin
#     max_value <- ymax }
#   
#   range_value <- max_value - min_value
#   step_size <- ceiling((range_value * 0.2) / 2) *2
#   breaks_value <- seq(min_value, max_value, by = step_size)
#   breaks_value <- round(breaks_value /2) * 2
#   
#   min_x <- 1
#   max_x <- 5
#   breaks_x <- seq(min_x, max_x, by = 1)
#   
#   unit_label <- unit
#   
#   # Plot
#   plot <- ggplot(data = changed, aes(x = DOW_in_block, y = weight, color = Diet, fill = Diet)) +
#     geom_ribbon(aes(y = weight, ymin = weight, ymax = weight + sd), alpha = 0.2, linetype = 0) +
#     geom_point(size = 3) +
#     geom_line(linewidth = 1) +
#     geom_text(aes(label = n), hjust = 0, vjust = -1, size = 3, show.legend = FALSE) +
#     scale_color_manual(values = c("purple2", "peachpuff3")) +
#     scale_fill_manual(values = c("purple2","peachpuff3")) +
#     scale_x_continuous(name = "Measurment days",
#                        limits = c(min_x, max_x),
#                        breaks = breaks_x,
#                        minor_breaks = seq(min_x, max_x, by = 1)) +
#     scale_y_continuous(name = sprintf("%s [%s]", deparse(substitute(value)), unit_label),
#                        limits = c(min_value, max_value),
#                        breaks = breaks_value) +
#     xlab("Time on CD-HFD [d]") +
#     ylab(sprintf("%s [%s]", deparse(substitute(value)), unit_label)) +
#     theme_bw() +
#     #ggtitle(sprintf("%s of %ss \n from batch %s (n > %d)", value_label_final, sex, batch, N)) +
#     guides(x = guide_axis(cap = "upper", minor.ticks = TRUE),
#            y = guide_axis(cap = "upper")) +
#     theme(axis.line = element_line(colour = "black"),
#           panel.grid.major = element_blank(),
#           panel.grid.minor = element_blank(),
#           panel.border = element_blank(),
#           panel.background = element_blank(),
#           axis.ticks.length = unit(4, "pt"))
#     
#   print(plot)
#   
#   # ---  Saving Plot --- 
#   ggsave(filename = paste0(file_base, ".png"), plot = plot,  path = path_images, width = 3, height = 4,dpi = 300)
#   #ggsave(filename = paste0(file_base, ".pdf"), plot = plot, path = path_images, width = 9, height = 6, dpi = 300, device = cairo_pdf)
#   dev.off()
#   # Optionally save stats tables
#   # write.csv(anova_table, file = file.path(path_images, paste0(file_base, "_anova.csv")), row.names = FALSE)
#   # write.csv(pwc, file = file.path(path_images, paste0(file_base, "_pairwise.csv")), row.names = FALSE)
#   
#   # --- Return output --- 
#   return(list(
#     #plot = plot
#   ))
#   
# }
# d<-data%>%filter(Animal %in% c("EC1", "EC2"))
# do_EC_curve(d,value = Food_consumed,  value_label = "Food Loss", unit = "g",  
#               batch = "1",  sex = "both", block_value= NULL, N = 0,
#               path_images = path_for_saving_images,ymin=0, ymax=16)
# do_EC_curve(d,value = Food_consumed,  value_label = "Food Loss", unit = "g",  
#             batch = "1",  sex = "EC1", block_value= 1, N = 0,
#             path_images = path_for_saving_images,ymin=0, ymax=16)
# do_EC_curve(d,value = Food_consumed,  value_label = "Food Loss", unit = "g",  
#             batch = "1",  sex = "EC1", block_value= 0, N = 0,
#             path_images = path_for_saving_images,ymin=0, ymax=16)
# do_EC_curve(d,value = Food_consumed,  value_label = "Food Loss", unit = "g",  
#             batch = "1",  sex = "EC1", block_value= 2, N = 0,
#             path_images = path_for_saving_images,ymin=0, ymax=16)
# do_EC_curve(d,value = Water_consumed,  value_label = "Water Loss", unit = "g",  
#             batch = "1",  sex = "both", block_value=NULL, N = 0,
#             path_images = path_for_saving_images,ymin=0, ymax=45)
# do_EC_curve(d,value = Water_consumed,  value_label = "Water Loss", unit = "g",  
#             batch = "1",  sex = "both", block_value=0, N = 0,
#             path_images = path_for_saving_images,ymin=0, ymax=45)
# do_EC_curve(d,value = Water_consumed,  value_label = "Water Loss", unit = "g",  
#             batch = "1",  sex = "both", block_value=1, N = 0,
#             path_images = path_for_saving_images,ymin=0, ymax=45)
# do_EC_curve(d,value = Water_consumed,  value_label = "Water Loss", unit = "g",  
#             batch = "1",  sex = "both", block_value=2, N = 0,
#             path_images = path_for_saving_images,ymin=0, ymax=45)
# 
# 

