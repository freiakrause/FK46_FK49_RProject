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
library(waffle)
library(rstatix)
library(lmerTest)
library(emmeans)
library(grid)
library(ggnewscale)
library(NADA2)
library(effsize)
# Read Raw Inputdata after general Data manipulation ------------------------------------------------------
setwd("C:/Users/b1084855/OneDrive - Universität Salzburg/Freigegebene Dokumente - AG_Tumorimmunologie/Data/Freia Krause/01_Experiments/FK46_iALmice_high Fat diet 52 weeks 7d after injection/Analysis/")
load("01_RawData/FK46_Data_prepared.Rda")
# Function for Weight Curves ----------------------------------------------
do_weight_curve <- function(inputdata, value, value_label = NULL, unit = "g",
                            batch = "ALL", sex = "both", N, path_images){
  # --- Setup ---
  value_label_final <- if (is.null(value_label)) deparse(substitute(value)) else value_label
  file_base <- paste0("FK46_", value_label_final, "_Batch", batch, "_", sex, "_n", N)
  
  # --- Filter by sex ---
  filtered <- inputdata %>%
    filter(!is.na({{value}})) %>%
    filter(case_when(
      sex == "female" ~ Sex == "female",
      sex == "male" ~ Sex == "male",
      sex == "both" ~ TRUE))
  
  # --- BATCH filtering ---
  if (batch == "ALL") {
    common_timepoints <- filtered %>%       # Find common timepoints across both batches
      filter(BATCH %in% c(1, 2)) %>%
      group_by(wks_diet, BATCH) %>%
      summarise(n = n(), .groups = "drop") %>%
      group_by(wks_diet) %>%
      summarise(n_batches = n_distinct(BATCH)) %>%
      filter(n_batches == 2) %>%
      pull(wks_diet)
    
    filtered <- filtered %>%
      filter(BATCH %in% c(1, 2)) %>%
      filter(wks_diet %in% common_timepoints)
  } 
  else {
    filtered <- filtered %>% filter(BATCH == batch)
  }
  print(filtered) 
  changed <- filtered %>%
    group_by(Treatment, wks_diet) %>%
    summarise(weight = mean({{value}}, na.rm = TRUE),n = n(), sd = sd({{value}}, na.rm = TRUE)) %>%
    filter(n > N)
  
  ## Statistical Tests of Weight Curves --------------------------------------------------------
  ### Stat Tests --------------------------------------------------------
  value_str<-deparse(substitute(value))
  
  stat_tests <- filtered %>%convert_as_factor(Animal,wks_diet)
  
  stat_tests %>%
    group_by(Treatment, wks_diet) %>%
    get_summary_stats({{value}}, type = "mean_sd")
  
  bxp <- ggboxplot( stat_tests, x = "wks_diet", y = value_str,color = "Treatment", palette = "jco")
  print(bxp)
  
  outliers<-stat_tests %>%
    group_by(Treatment, wks_diet) %>%
    identify_outliers({{value}})
  
  n<-stat_tests %>%group_by(Treatment, wks_diet) %>%summarize(n=n())%>%filter(n>3)
  
  # shapiro_results <- stat_tests %>%
  #   group_by(Treatment, wks_diet) %>%
  #   summarise(n = sum(!is.na({{value}})),
  #             all_same = all({{value}} == 100, na.rm = TRUE),      # Check if all rel.weight values are the same for this timepoint
  #             shapiro_p = ifelse(all_same, NA_real_,
  #                                ifelse(n >= 3, shapiro.test(pick({{value}})[[1]])$p.value, NA_real_)),
  #             result = ifelse(shapiro_p >= 0.05, "normal", "non-normal"),
  #             note = ifelse(n >= 3, "OK", "Too few samples"), .groups = "drop")
  # 
  #print(shapiro_results,n=40)
  filtered$wks_diet <- factor(filtered$wks_diet)
  
  # --- Fit the linear mixed-effects model --- 
  model <- lmer(as.formula(paste(deparse(substitute(value)), "~ Treatment * wks_diet + (1 | Animal)")), data = filtered)
  anova_table<-anova(model, type = 3)
  
  anova_label <- paste0("ANOVA over linear mixed-effects model\n",
                        "Treatment: F = ", round(anova_table$F[1], 2), ", p = ", signif(anova_table$`Pr(>F)`[1], 3), "\n",
                        "Time: F = ", round(anova_table$F[2], 2), ", p < ", format.pval(anova_table$`Pr(>F)`[2], digits = 1), "\n",
                        "Interaction: F = ", round(anova_table$F[3], 2), ", p = ", signif(anova_table$`Pr(>F)`[3], 3))
  
  
  # --- Calculate the estimated marginal means for both Treatment and wks_diet --- 
  emm <- emmeans(model, ~ Treatment | wks_diet)
  
  # --- Perform pairwise contrasts between Treatment levels (Ctrl vs TAM) at all levels of wks_diet --- 
  pwc <- contrast(emm, method = "pairwise", adjust = "bonferroni")
  posthoc_label<- "Post Hoc: Pairwise with Bonferroni correction"
  
  pwc_df <- as.data.frame(pwc)
  pwc_df_rounded <- pwc_df %>%
    mutate(rounded_p_value = ifelse(is.na(p.value), "NA", round(p.value, 3))) %>%
    mutate(wks_diet = as.character(wks_diet)) %>%
    mutate(wks_diet = as.numeric(wks_diet)) %>%
    mutate(significance = case_when(
      is.na(p.value) ~ "NA",                     # For NA p-values
      p.value < 0.001 ~ "***",                   # p < 0.001 is highly significant
      p.value >= 0.001 & p.value < 0.01 ~ "**", # 0.001 ≤ p < 0.01 is significant
      p.value >= 0.01 & p.value < 0.05 ~ "*",    # 0.01 ≤ p < 0.05 is moderately significant
      p.value >= 0.05 ~ "NS",                    # p ≥ 0.05 is not significant
      TRUE ~ "NA"   ))  %>%                           # Default case
    select(wks_diet, rounded_p_value, significance)
  
  # --- Get max weight per wks_diet from 'changed' (across both Treatment groups) --- 
  max_weights <- changed %>%
    group_by(wks_diet) %>%
    summarise(max_weight = max(weight, na.rm = TRUE)) %>%
    ungroup() %>%
    mutate(y_position = max_weight + max_weight * 0.1)
  
  # --- Join y_position back into pwc_df_rounded --- 
  pwc_df_annotated <- pwc_df_rounded %>%
    left_join(max_weights, by = "wks_diet")
  
  # Variables for Plot Setup
  mean_value <- mean(changed$weight, na.rm = TRUE)
  sd_value <- sd(changed$weight, na.rm = TRUE)
  min_value <- round(mean_value - 3 * sd_value)
  max_value <- round(mean_value + 3 * sd_value)
  range_value <- max_value - min_value
  step_size <- ceiling((range_value * 0.2) / 5) * 5
  breaks_value <- seq(min_value, max_value, by = step_size)
  breaks_value <- round(breaks_value / 5) * 5
  
  min_x <- round(min(changed$wks_diet, na.rm = TRUE))
  max_x <- round(max(changed$wks_diet, na.rm = TRUE) + 1)
  breaks_x <- c(-1, 0, seq(4, max_x, by = 4))
  unit_label <- unit
  
  # Plot
  plot <- ggplot(data = changed, aes(x = wks_diet, y = weight, color = Treatment, fill = Treatment)) +
    geom_ribbon(aes(y = weight, ymin = weight - sd, ymax = weight + sd), alpha = 0.1, linetype = 0) +
    geom_point(size = 3) +
    geom_line(linewidth = 1) +
    geom_text(aes(label = n), hjust = 0, vjust = -1, size = 3, show.legend = FALSE) +
    scale_color_manual(values = c("grey30", "palevioletred4")) +
    scale_fill_manual(values = c("grey30", "palevioletred4")) +
    scale_x_continuous(name = "Time on CD-HFD [wks]",
                       limits = c(min_x, max_x),
                       breaks = breaks_x,
                       minor_breaks = seq(min_x, max_x, by = 1)) +
    scale_y_continuous(name = sprintf("%s [%s]", deparse(substitute(value)), unit_label),
                       limits = c(min_value, max_value),
                       breaks = breaks_value) +
    xlab("Time on CD-HFD [wks]") +
    ylab(sprintf("%s [%s]", deparse(substitute(value)), unit_label)) +
    theme_bw() +
    ggtitle(sprintf("%s of %ss from batch %s (n > %d)", value_label_final, sex, batch, N)) +
    guides(x = guide_axis(cap = "upper", minor.ticks = TRUE),
           y = guide_axis(cap = "upper")) +
    theme(axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank(),
          axis.ticks.length = unit(4, "pt"))+
    annotate("text", 
             x = pwc_df_annotated$wks_diet, 
             y = pwc_df_annotated$y_position,
             label = pwc_df_annotated$significance, 
             size = 2.5, 
             color = "black", 
             fontface = "italic")+
    annotate("text", x = min_x + 2, y = min_value + (max_value - min_value) * 0.05,  # 5% above the bottom
             label = anova_label,
             size = 2,
             hjust = 0,
             color = "black", 
             fontface = "italic")+
    annotate("text", x = min_x + 12, y = min_value + (max_value - min_value) * 0.05,  # 5% above the bottom
             label = posthoc_label,
             size = 2,
             hjust = 0,
             color = "black", 
             fontface = "italic")
  print(plot)
  
  # ---  Saving Plot --- 
  ggsave(filename = paste0(file_base, ".png"), plot = plot,  path = path_images, width = 9, height = 6,dpi = 300)
  #ggsave(filename = paste0(file_base, ".pdf"), plot = plot, path = path_images, width = 9, height = 6, dpi = 300, device = cairo_pdf)
  
  # Optionally save stats tables
  # write.csv(anova_table, file = file.path(path_images, paste0(file_base, "_anova.csv")), row.names = FALSE)
  # write.csv(pwc, file = file.path(path_images, paste0(file_base, "_pairwise.csv")), row.names = FALSE)
  
  # --- Return output --- 
  return(list(
    outliers = outliers,
    #shapiro_results = shapiro_results,
    anova_table = anova_table,
    posthoc = pwc_df,
    model= model
    #plot = plot
  ))
  
}

#
path_for_saving_images<-"02_GeneratedData/background"
##Run Function for Weight Curves and save output --------------------------------------------------------
do_weight_curve(data, value=Weight,value_label = "Body Weight", unit = "g", batch=2, sex="male",N=0,path_for_saving_images)
do_weight_curve(data, value=Weight,value_label = "Body Weight", unit = "g", batch=1, sex="male",N=0,path_for_saving_images)
do_weight_curve(data, value=Weight,value_label = "Body Weight", unit = "g", batch=2, sex="female",N=0,path_for_saving_images)
do_weight_curve(data, value=Weight,value_label = "Body Weight", unit = "g", batch=1, sex="female",N=0,path_for_saving_images)
do_weight_curve(data, value=rel.weight, value_label = "rel. BW",unit = "perc", batch="1", sex="male",N=0,path_for_saving_images)
do_weight_curve(data, value=rel.weight, value_label = "rel. BW",unit = "perc", batch="1", sex="female",N=0,path_for_saving_images)
do_weight_curve(data, value=rel.weight, value_label = "rel. BW",unit = "perc", batch="2", sex="female",N=0,path_for_saving_images)
do_weight_curve(data, value=rel.weight, value_label = "rel. BW",unit = "perc", batch="2", sex="male",N=0,path_for_saving_images)
path_for_saving_images<-"02_GeneratedData"
do_weight_curve(data, value=Score, value_label = "Score",unit = "", batch="ALL", sex="male",N=0,path_for_saving_images)
do_weight_curve(data, value=Score, value_label = "Score",unit = "", batch="ALL", sex="female",N=0,path_for_saving_images)
do_weight_curve(data, value=Score, value_label = "Score",unit = "", batch="ALL", sex="both",N=0,path_for_saving_images)

do_weight_curve(data, value=rel.weight, value_label = "rel. BW",unit = "perc", batch="ALL", sex="male",N=0,path_for_saving_images)
do_weight_curve(data, value=rel.weight, value_label = "rel. BW",unit = "perc", batch="ALL", sex="female",N=0,path_for_saving_images)
gc()