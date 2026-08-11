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

# Hie weiter machen und Plotting sortieren und korrekte statistic hinzufügen
d_food <-readRDS(file = file.path(PATHS$FoodIntake$input,"FoodIntake_Cagewise.rds"))
d_water<-readRDS(file = file.path(PATHS$FoodIntake$input,"WaterIntake_Cagewise.rds"))
LMM_Food <-read.csv2(file = file.path(PATHS$FoodIntake$output,"Statistics/Food_consumption_LMM_results.csv"))
LMM_Water <-read.csv2(file = file.path(PATHS$FoodIntake$output,"Statistics/Water_consumption_LMM_results.csv"))
Posthoc_Food <-read.csv2(file = file.path(PATHS$FoodIntake$output,"Statistics/Food_consumption_LMM_posthoc_by_Block.csv"))
Posthoc_Water<-read.csv2(file = file.path(PATHS$FoodIntake$output,"Statistics/Water_consumption_LMM_posthoc_by_Block.csv"))

output_pwd <-file.path(PATHS$FoodIntake$output)
# Food / Water Plotting Data -------
##Food Barplot summed over all -----
d_for_barplot <- d_food %>%
  group_by(Cage, Treatment,Sex) %>%  # group by cage AND treatment
  summarize(mean_Food = mean(Food_consumed, na.rm = TRUE),.groups = "drop")

stats_F <- d_food %>%
  group_by(Treatment) %>%
  summarise(Mean_F = mean(d_for_barplot$mean_Food, na.rm = TRUE),SD_F = sd(d_for_barplot$mean_Food, na.rm = TRUE),.groups = "drop" )
#hier später korrekte statistic/ pvalue einführen
pF<- ggplot(stats_F, aes(x = Treatment, y = Mean_F, fill = Treatment)) +
  geom_bar(stat = "identity", color = "black", alpha = 0.5, width = 0.75, position = "dodge") +
  geom_point(data = d_for_barplot, fill = "lightgrey", color = "black",shape=21,aes(y = mean_Food),   position = position_jitterdodge(0.1, dodge.width = 0.75), alpha = 0.8, size = 5.3, stroke = 1.8) +
  geom_errorbar(aes(ymin = Mean_F - SD_F, ymax = Mean_F + SD_F),  position = position_dodge(width = 0.75), width = 0.2) +
  scale_fill_manual(values = Treatment_colors[c("Ctrl","TAM")], labels = c("Ctrl", "TAM")) +
  scale_y_continuous(name = paste0("Mean Daily Food Intake per animal [g]"),limits = c(0,20), breaks = seq(0, 20, by = 2)) +
  labs(x = "Treatment") +
  annotate("text", x = 1.5, y = max(d_for_barplot$mean_Food, na.rm = TRUE) * 1.1, 
           label = "ns", size = 6, fontface = "italic") +
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
  guides( shape = guide_legend(title = "Status", nrow = 1), 
          fill = guide_legend(title = "Treatment", override.aes = list(shape = 21), nrow = 1), color = "none")
ggsave(filename = "FK49_Food Consumption_Summary_ALL.png", plot = pF, path = file.path(output_pwd,"/Plots"), width = 4, height = 11, dpi = 300)


# ##Water Barplot summed over all -----
# pW<- ggplot(stats_FW, aes(x = Treatment, y = Mean_W, fill = Treatment)) +
#   geom_bar(stat = "identity", color = "black", alpha = 0.5, width = 0.75, position = "dodge") +
#   geom_point(data = d_for_barplot, fill = "lightgrey", color = "black",shape=21,aes(y = mean_Water),   position = position_jitterdodge(0.1, dodge.width = 0.75), alpha = 0.8, size = 5.3, stroke = 1.8) +
#   geom_errorbar(aes(ymin = Mean_W - SD_W, ymax = Mean_W + SD_W),  position = position_dodge(width = 0.75), width = 0.2) +
#   scale_fill_manual(values = Treatment_colors[c("Ctrl","TAM")], labels = c("Ctrl", "TAM")) +
#   scale_y_continuous(name = paste0("Mean Daily Water Intake per animal [g]"), limits = c(0,16), breaks = seq(0, 16, by = 2)) +
#   labs(x = "Treatment") +
#   annotate("text", x = 1.5, y = max(d_for_barplot$mean_Food, na.rm = TRUE) * 1.1, label = "ns", size = 6, fontface = "italic") +
#   theme_minimal() +
#   theme(
#     legend.position = "bottom",
#     axis.line = element_line(color = "black", linewidth = 0.5),
#     axis.ticks = element_line(color = "black", linewidth = 0.5),
#     axis.title = element_text(size = 20, face = "bold"),
#     axis.title.x = element_blank(),
#     axis.text = element_text(size = 19, face = "bold"),
#     plot.title = element_blank(),
#     legend.title = element_text(size = 10),
#     legend.text = element_text(size = 10),
#     panel.grid = element_blank()) +
#   guides( shape = guide_legend(title = "Status", nrow = 1), fill = guide_legend(title = "Treatment", override.aes = list(shape = 21), nrow = 1),            color = "none")
# ggsave(filename = "FK49_Water Consumption_Summary_ALL.png", plot = pW, path =  file.path(output_pwd,"/Plots"), width = 4, height = 11, dpi = 300)

rm(d_for_barplot,stats_FW,pF,pW) 
## Plot Food and water consumption per block ----
d_for_Block <- d_food %>%
  group_by(Cage, Treatment, Sex, Block, BATCH) %>%  
  summarise(mean_Food = mean(Food_consumed, na.rm = TRUE),
            n = n(),  .groups = "drop") %>%
  mutate( wks_diet = case_when(
    Block == "1"  ~ 0.3575,
    Block == "4"  ~ 3.3575,
    Block == "8"  ~ 7.3575,
    Block == "12" ~ 11.3575,
    TRUE ~ NA_real_ ) )

stats2_F <- d_for_Block %>%
  group_by(Treatment,wks_diet) %>%
  summarise(Mean_F = mean(mean_Food, na.rm = TRUE),
            SD_F   = sd(mean_Food, na.rm = TRUE),
         n=n(),
            .groups = "drop")
### Plotting Food per Block----
plot <- ggplot(data = stats2_F, aes(x = wks_diet, y = Mean_F, color = Treatment, fill = Treatment)) +
  geom_ribbon(data = stats2_F,  aes(x = wks_diet, ymin = Mean_F,  ymax = Mean_F + SD_F,  fill = Treatment,  group = Treatment), alpha = 0.1, linetype = 0)+
  geom_line(data = stats2_F,   aes(x = wks_diet, y = Mean_F, color = Treatment, group = Treatment), linewidth = 1) +
  geom_point(data = stats2_F,aes(x = wks_diet, y = Mean_F, color = Treatment),  size = 3, stroke = 1.1) +
  geom_text(aes(y = c(Mean_F[1:4]+1.8,Mean_F[1:4]+1.8) , label = n), position = position_dodge(width = 0.4), hjust = 0.5, size = 3, show.legend = FALSE)+
  scale_color_manual(values = Treatment_colors[c("Ctrl","TAM")]) +
  scale_fill_manual (values = Treatment_colors[c("Ctrl","TAM")]) +
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
        axis.ticks.length = unit(4, "pt"))
# +
#   annotate("text", x = as.numeric(as.character(pwc_df_rounded_F$Block))-0.6425,
#            y = stats2_FW$Mean_F[1:4]+1.3,   label = pwc_df_rounded_F$significance,
#            size = 2.5, color = "black",    fontface = "italic") +
#   annotate("text", x=0.3, y = 4,   label = anova_label_F,     size = 2,   hjust = 0,    color = "black",   fontface = "italic") +
#   annotate("text", x=0.3, y = 2.8, label = posthoc_label_F,   size = 2,   hjust = 0,    color = "black",   fontface = "italic")

# Saving Plot ---
ggsave(filename = "FK49_Food_Consumption_Summary_Block.png", plot = plot,  path = file.path(output_pwd,"/Plots"), width = 9, height = 6,dpi = 300)

rm(emm_F,model_F,pwc_df_F,pwc_df_rounded_F,pwc_F,anova_table_F,anova_label_F,posthoc_label_F)

### Plotting Water per Block ----
plot <- ggplot(data = stats2_FW, aes(x = wks_diet, y = Mean_W, color = Treatment, fill = Treatment)) +
  geom_ribbon(data = stats2_FW,  aes(x = wks_diet, ymin = Mean_W,  ymax = Mean_W + SD_W,  fill = Treatment,  group = Treatment), alpha = 0.1, linetype = 0)+
  geom_line(data = stats2_FW,   aes(x = wks_diet, y = Mean_W, color = Treatment, group = Treatment), linewidth = 1) +
  geom_point(data = stats2_FW,aes(x = wks_diet, y = Mean_W, color = Treatment),  size = 3, stroke = 1.1) +
  geom_text(aes(y = c(Mean_W[1:4]+0.75,Mean_W[1:4]+0.75) , label = n), position = position_dodge(width = 0.4), hjust = 0.5, size = 3, show.legend = FALSE)+
  scale_color_manual(values = Treatment_colors[c("Ctrl","TAM")]) +
  scale_fill_manual(values = Treatment_colors[c("Ctrl","TAM")]) +
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

ggsave(filename = "FK49_Water_Consumption_Summary_Block.png", plot = plot,  path =  file.path(output_pwd,"/Plots"), width = 9, height = 6,dpi = 300)

rm(emm_W,model_W,pwc_df_W,pwc_df_rounded_W,pwc_W,anova_table_W,anova_label_W,posthoc_label_W)
rm(d_for_Block,stats2_FW)

#Function for Food Curves ----------------------------------------------
# do_food_curve <- function(inputdata, value, value_label = NULL, unit = "g",
#                           batch = "ALL", sex = "both",Block, N, path_images, 
#                           ymin = NULL, ymax = NULL,offset=1, savestats = "NO"){
#   # --- Setup ---
#   value_label_final <- if (is.null(value_label)) deparse(substitute(value)) else value_label
#   file_base <- paste0("FK49_", value_label_final, "_Batch", batch, "_", sex,"_B",Block,"_n", N)
#   
#   # --- Filter by sex ---
#   filtered <- inputdata %>%
#     filter(!is.na({{value}})) %>%
#     filter(case_when(
#       sex == "female" ~ Sex == "female",
#       sex == "male" ~ Sex == "male",
#       sex == "EC1" ~ Sex == "EC1",
#       sex == "EC2" ~ Sex == "EC2",
#       sex == "EC" ~ Sex == "EC",
#       sex == "both" ~ TRUE))
#   
#   # --- BATCH filtering ---
#   if (batch == "ALL") {
#     common_timepoints <- filtered %>%       # Find common timepoints across both batches
#       filter(BATCH %in% c(1, 2)) %>%
#       group_by(days_diet, BATCH) %>%
#       summarise(n = n(), .groups = "drop") %>%
#       group_by(days_diet) %>%
#       summarise(n_batches = n_distinct(BATCH)) %>%
#       filter(n_batches == 2) %>%
#       pull(days_diet)
#     filtered <- filtered %>%
#       filter(BATCH %in% c(1, 2)) %>%
#       filter(days_diet %in% common_timepoints) } 
#   else {filtered <- filtered %>% filter(BATCH == batch) }
#   print(filtered) 
#   
#   changed <- filtered %>%  filter(Block==!!Block) %>%
#     group_by(Treatment, days_diet) %>%
#     summarise(value_used = mean({{value}}, na.rm = TRUE),n = n(), sd = sd({{value}}, na.rm = TRUE)) %>%
#     filter(n > N)
#   
#   ## Statistical Tests of Weight Curves --------------------------------------------------------
#   value_str<-deparse(substitute(value))
#   stat_tests <- filtered %>% convert_as_factor(Cage, days_diet)
#   
#   STATS <- tryCatch({ 
#     stat_tests %>%group_by(Treatment, days_diet) %>%  get_summary_stats({{value}}, type = "mean_sd")}, error = function(e) 
#     {  message("⚠️ Fehler bei get_summary_stats(): ", e$message)
#       return(NULL)})
#   
#   outliers <- tryCatch({stat_tests %>%
#       group_by(Treatment, days_diet) %>%
#       identify_outliers({{value}})}, error = function(e) {message("⚠️ Fehler bei identify_outliers(): ", e$message)
#         return(NULL)})
#   
#   
#   filtered$days_diet <- factor(filtered$days_diet)
#   
#   # --- Fit the linear mixed-effects model --- 
#   model <- lmer(as.formula(paste(deparse(substitute(value)), "~ Treatment * days_diet + (1 | Cage)")), data = filtered)
#   model <- lmer(mean_Food ~ Treatment * wks_diet + Sex + (1 + wks_diet | Cage),data = d_food_week, REML = TRUE )
#   anova_table<-anova(model, type = 3)
#   
#   anova_label <- paste0("ANOVA over linear mixed-effects model\n",
#                         "Treatment: F = ", round(anova_table$F[1], 2), ", p = ", signif(anova_table$`Pr(>F)`[1], 3), "\n",
#                         "Time: F = ", round(anova_table$F[2], 2), ", p < ", format.pval(anova_table$`Pr(>F)`[2], digits = 1), "\n",
#                         "Interaction: F = ", round(anova_table$F[3], 2), ", p = ", signif(anova_table$`Pr(>F)`[3], 3))
#   
#   
#   # --- Calculate the estimated marginal means for both Treatment and days_diet --- 
#   emm <- emmeans(model, ~ Treatment | days_diet)
#   
#   # --- Perform pairwise contrasts between Treatment levels (ctrl vs TAM) at all levels of days_diet --- 
#   pwc <- contrast(emm, method = "pairwise", adjust = "bonferroni")
#   posthoc_label<- "Post Hoc: Pairwise with Bonferroni correction"
#   
#   pwc_df <- as.data.frame(pwc)
#   pwc_df_rounded <- pwc_df %>%
#     mutate(rounded_p_value = ifelse(is.na(p.value), "NA", round(p.value, 3))) %>%
#     mutate(days_diet = as.character(days_diet)) %>%
#     mutate(days_diet = as.numeric(days_diet)) %>%
#     mutate(significance = case_when(
#       is.na(p.value) ~ "NA",                     # For NA p-values
#       p.value < 0.001 ~ "***",                   # p < 0.001 is highly significant
#       p.value >= 0.001 & p.value < 0.01 ~ "**", # 0.001 ≤ p < 0.01 is significant
#       p.value >= 0.01 & p.value < 0.05 ~ "*",    # 0.01 ≤ p < 0.05 is moderately significant
#       p.value >= 0.05 ~ "NS",                    # p ≥ 0.05 is not significant
#       TRUE ~ "NA"   ))  %>%                           # Default case
#     select(days_diet, rounded_p_value, significance)
#   
#   
#   # --- Get max weight per days_diet from 'changed' (across both Treatment groups) --- 
#   max_values <- changed %>%
#     group_by(days_diet) %>%
#     summarise(max_value = max(value_used, na.rm = TRUE)) %>%
#     ungroup() %>%
#     mutate(y_position = max_value + max_value * 0.05)
#   
#   # --- Join y_position back into pwc_df_rounded --- 
#   pwc_df_annotated <- pwc_df_rounded %>%
#     left_join(max_values, by = "days_diet")
#   changed <- changed %>%left_join(pwc_df_annotated %>% select(days_diet, y_position), by = "days_diet")
#   # Variables for Plot Setup
#   # Calculate default y limits only if not provided
#   if (is.null(ymin) || is.null(ymax)) {
#     mean_value <- mean(changed$value_used, na.rm = TRUE)
#     sd_value <- sd(changed$value_used, na.rm = TRUE)
#     min_value_default <- round(mean_value - 3 * sd_value)
#     max_value_default <- round(mean_value + 3 * sd_value) + 4
#     min_value <- ifelse(is.null(ymin), min_value_default, ymin)
#     max_value <- ifelse(is.null(ymax), max_value_default, ymax)} 
#   else {min_value <- ymin
#   max_value <- ymax }
#   
#   range_value <- max_value - min_value
#   step_size <- ceiling((range_value * 0.2) / 2) *2
#   breaks_value <- seq(min_value, max_value, by = step_size)
#   breaks_value <- round(breaks_value /2) * 2
#   
#   min_x <- round(min(changed$days_diet, na.rm = TRUE))
#   max_x <- round(max(changed$days_diet, na.rm = TRUE))
#   breaks_x <- seq(min_x, max_x, by = 1)
#   
#   unit_label <- unit
#   
#   # Plot
#   plot <- ggplot(data = changed, aes(x = days_diet, y = value_used, color = Treatment, fill = Treatment)) +
#     geom_ribbon(aes(y = value_used, ymin = value_used, ymax = value_used + sd), alpha = 0.1, linetype = 0) +
#     geom_point(data = d, fill = "lightgrey", color = "black",
#                aes(y = .data[[value]], shape = event_status),
#                position =  position_jitterdodge(0.1, dodge.width = 0.75), alpha = 0.8, size = 3, stroke = 1.1)+
#     geom_line(linewidth = 1) +
#     geom_text(aes(y = y_position + offset, label = n), position = position_dodge(width = 0.4), hjust = 0.5, size = 3, show.legend = FALSE)+
#     scale_color_manual(values = c("grey30", "palevioletred4")) +
#     scale_fill_manual(values = c("grey30","palevioletred4")) +
#     scale_x_continuous(name = "Time on CD-HFD [days]",
#                        limits = c(min_x-0.1, max_x+0.1),
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
#           axis.ticks.length = unit(4, "pt"))+
#     annotate("text", 
#              x = pwc_df_annotated$days_diet, 
#              y = pwc_df_annotated$y_position+offset/2,
#              label = pwc_df_annotated$significance, 
#              size = 2.5, 
#              color = "black", 
#              fontface = "italic")
#   
#   # ---  Saving Plot --- 
#   ggsave(filename = paste0(file_base, ".png"), plot = plot,  path = path_images, width = 3, height = 4,dpi = 300)
#   #ggsave(filename = paste0(file_base, ".pdf"), plot = plot, path = path_images, width = 9, height = 6, dpi = 300, device = cairo_pdf)
#   #dev.off()
#   # Optionally save stats tables
#   if (savestats == "YES") {
#     outliers <- mutate(outliers, table = "Outliers")
#     
#     anova_table <- mutate(anova_table, table = "ANOVA")
#     pwc_df <- mutate(pwc_df, table = "Pairwise Comparison")
#     StatsOutput <- bind_rows(outliers, anova_table,pwc_df)%>% relocate( table)
#     write.csv2(
#       StatsOutput,
#       file = file.path(paste0(path_images,"/Statistics"), paste0(file_base, "_StatsOutput.csv")),
#       row.names = FALSE,
#       na = "",
#       fileEncoding = "UTF-8")
#   }
#   
#   # --- Return output --- 
#   return(list(
#     outliers = outliers,
#     anova_table = anova_table,
#     posthoc = pwc_df,
#     model= model,
#     plot = plot
#   ))
#   
# }

##



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

