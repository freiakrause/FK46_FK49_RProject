rm(list=ls())
gc()
library(tidyverse)
library(rstatix)
library(lmerTest)
library(emmeans)
library(ggpubr)
source("FK49_Definitions.R")
# Hier weiter machen Comments einfügen und stats rausschreiben
# Read Raw Inputdata after general Data manipulation ------------------------------------------------------
load(file.path(PATHS$general_data$FK49_output,"01_RawData/FK49_Data_prepared.Rda"))

# Function for Weight Curves ----------------------------------------------
#Function assumes if you say BATCH== "ALL that you have batch 1 and 2. 
#If this is not right, change it in function to either take all numerical batches or to the numbers you have
#
do_weight_curve <- function(inputdata, value, value_label = NULL, unit = "g",
                            batch = "ALL", sex = "both", N, path_images,savestats = "NO"){
  
  value_label_final <- if (is.null(value_label)) deparse(substitute(value)) else value_label #label for plotting and saving
  file_base <- paste0("FK49_", value_label_final, "_Batch", batch, "_", sex, "_n", N)
  
  # filter for specified sex -----
  filtered <- inputdata %>%
    select(Sex, BATCH,Treatment,DOW,wks_diet,Animal,Block,days_diet, Cage ,{{value}})%>%
    filter(!is.na({{value}})) %>%
    filter(case_when(
      sex == "female" ~ Sex == "female",
      sex == "male" ~ Sex == "male",
      sex == "both" ~ TRUE))
  
  # BATCH filtering -----
  if (batch == "ALL") {
    common_timepoints <- filtered %>%       # Find common time points across both batches
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
  # Filter out food days -----
   #they should not appear in this overall plot. 
   #Here I only want weekly measurements not the Food/Water Intake Daily weights
    filtered <- filtered %>%
    group_by(Animal, Block) %>%
    filter((Block %in% c("0") & days_diet == -7) | # In block 0 and 1 all batches were weight on the Monday,
           (Block %in% c("1") & days_diet == 0)|    #DOW_in_Block 1 so thats the day i want to represent the week
           (!Block %in% c("0", "1")  & days_diet == as.numeric(as.character(Block))*7-4) )%>%
    ungroup()
  # I summarize data from "filtered" dataset in "Mean_SD_data" to be able to plot the mean and sd later on  
  Mean_SD_data <- filtered %>%
    group_by(Treatment, wks_diet) %>%
    summarise(weight = mean({{value}}, na.rm = TRUE),n = n(), sd = sd({{value}}, na.rm = TRUE)) %>%
    filter(n > N)
  
  ## Statistical Tests of Weight Curves --------------------------------------------------------
  # ChatGPT did and helped a lot here. So i don't know everything exactly.
  
  value_string<-deparse(substitute(value))
  stat_tests <- filtered %>%convert_as_factor(Animal,wks_diet)
  stat_tests %>%group_by(Treatment, wks_diet) %>% get_summary_stats({{value}}, type = "mean_sd")
  ggboxplot(stat_tests, x = "wks_diet", y = value_string, color = "Treatment", palette = "jco")# could save this but dont do it here
  outliers<-stat_tests%>%group_by(Treatment, wks_diet) %>%identify_outliers({{value}})
  
  # linear mixed effects model for overall test of the data ----
   # does time have effect? 
   # does treatment have effect, 
   # do time and treatment have interaction effet 
   # 

  model_data <- filtered %>%
    mutate( Animal = factor(Animal),  
            Treatment = factor(Treatment),  
            wks_diet = as.numeric(as.character(wks_diet)) )  #wks_diet needs to be numeric for random slope
  #fit linear mixed-effects model 
  formula<-as.formula(paste(deparse(substitute(value)),
                            "~ Treatment * wks_diet + (1 + wks_diet | Animal)")) 
  #vfixed effects within interaction of treatment and time
  # (1+wks_diet|Animal) random effects
  #1 is random interept - each animal kann have individual basleine value at first measuremtn
  # wks diet means random slope each animal can have individual slope
  # |Animal data has repreated measures clustered within individual animals
                           
  model <- lmer(formula, REML = TRUE, data = model_data)
  anova_table<-anova(model, type = 3)
  anova_label <- paste0("ANOVA over linear mixed-effects model\n",
                        "Treatment: F = ", round(anova_table$F[1], 2), ", p = ", signif(anova_table$`Pr(>F)`[1], 3), "\n",
                        "Time: F = ", round(anova_table$F[2], 2), ", p < ", format.pval(anova_table$`Pr(>F)`[2], digits = 1), "\n",
                        "Interaction: F = ", round(anova_table$F[3], 2), ", p = ", signif(anova_table$`Pr(>F)`[3], 3))
  
  # To test at each specific timepoint and not the overall dataset -----
  
  # estimated marginal means for both Treatment and wks_diet
  time_points <- sort(unique(model_data$wks_diet))
  emm <- emmeans(model, ~ Treatment | wks_diet, at = list(wks_diet = time_points))
  
  # pairwise contrasts between Treatment levels (ctrl vs TAM) at all levels of wks_diet
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
      p.value >= 0.001 & p.value < 0.01 ~ "**",  # 0.001 ≤ p < 0.01 is significant
      p.value >= 0.01 & p.value < 0.05 ~ "*",    # 0.01 ≤ p < 0.05 is moderately significant
      p.value >= 0.05 ~ "NS",                    # p ≥ 0.05 is not significant
      TRUE ~ "NA"   ))  %>%                      # Default case
    select(wks_diet, rounded_p_value, significance)
  

  # Variables for Plot Setup
  mean_value <- mean(Mean_SD_data$weight, na.rm = TRUE)
  sd_value <-  sd(Mean_SD_data$weight, na.rm = TRUE)
  min_value <- round(mean_value - 4 * sd_value)
  max_value <- round(mean_value+5*sd_value)
  
  range_value <- max_value - min_value
  step_size <- ceiling((range_value * 0.2) / 5) * 5
  breaks_value <- seq(min_value, max_value, by = step_size)
  breaks_value <- round(breaks_value / 5) * 5
  
  min_x <- round(min(Mean_SD_data$wks_diet, na.rm = TRUE))
  max_x <- round(max(Mean_SD_data$wks_diet, na.rm = TRUE) + 1)
  breaks_x <- seq(-1, max_x, by = 1)
  
  unit_label <- unit
  
  max_weights <- Mean_SD_data %>%
    group_by(wks_diet) %>%
    summarise(max_weight = max(weight, na.rm = TRUE)) %>%
    ungroup() %>%
    mutate(y_position = max_weight + max_weight * 0.1)
  
  # Join y_position back into pwc_df_rounded
  pwc_df_annotated <- pwc_df_rounded %>%
    left_join(max_weights, by = "wks_diet")
  
  # Plot
  plot <- ggplot(data = Mean_SD_data, aes(x = wks_diet, y = weight, color = Treatment, fill = Treatment)) +
    geom_ribbon(aes(y = weight, ymin = weight - sd, ymax = weight + sd), alpha = 0.1, linetype = 0) +
    geom_point(size = 3) +
    geom_line(linewidth = 1) +
    geom_text(aes(label = n), hjust = 0, vjust = -1, size = 3, show.legend = FALSE) +
    scale_color_manual(values = c(Treatment_colors[c("Ctrl","TAM")],"black","pink")) +
    scale_fill_manual(values = c(Treatment_colors[c("Ctrl","TAM")],"black","pink")) +
    scale_x_continuous(name = "Time on CD-HFD [wks]", limits = c(min_x, max_x),
                       breaks = breaks_x, minor_breaks = seq(min_x, max_x, by = 1)) +
    scale_y_continuous(name = sprintf("%s [%s]", deparse(substitute(value)), unit_label),
                       limits = c(min_value, max_value), breaks = breaks_value) +
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
             size = 2.5, color = "black", fontface = "italic")+
    annotate("text", x = min_x + 2, y = min_value + range_value * 0.05,  # 5% above the bottom
             label = anova_label, size = 2,   hjust = 0,  color = "black", 
             fontface = "italic")+
    annotate("text", x = min_x + 6, y = min_value + range_value * 0.05,  # 5% above the bottom
             label = posthoc_label,size = 2, hjust = 0,  color = "black",      fontface = "italic")
  
  
  # ---  Saving Plot --- 
  ggsave(filename = paste0(file_base, ".png"), plot = plot,  path = path_images, width = 9, height = 6,dpi = 300)
  #ggsave(filename = paste0(file_base, ".pdf"), plot = plot, path = path_images, width = 9, height = 6, dpi = 300, device = cairo_pdf)
  #dev.off()
  
  # Optionally save stats tables
  if (savestats == "YES") {
    outliers <- outliers %>%  mutate(wks_diet = factor(wks_diet), table = "Outliers")
    anova_table <- anova(model, type = 3) %>% as.data.frame() %>% tibble::rownames_to_column(var = "Term")    
    pwc_df <- pwc_df %>%mutate(wks_diet = factor(wks_diet),table = "Pairwise Comparison")
    StatsOutput <- bind_rows(outliers,anova_table,pwc_df)%>% relocate( table)
    write.csv2( StatsOutput,
                file = file.path(paste0(path_images,"/Statistics"), paste0(file_base, "_StatsOutput.csv")),
                row.names = FALSE,  na = "",  fileEncoding = "UTF-8"  )
  }
  
  #--- Return output ---
  return(list(
    outliers = outliers,
    anova_table = anova_table,
    posthoc = pwc_df,
    model= model,
    plot = plot
  ))
  
}


# Run Function for Weight Curves and save output --------------------------------------------------------
path_for_saving_images<-file.path(PATHS$general_data$FK49_output,"02_GeneratedData/Weight_Organs")

## Relative Weight as summary from Batch 1 and 2 together --------------------------------------------------------
do_weight_curve(data, value=rel.weight, value_label = "rel. BW",
                unit = "perc", batch="ALL", sex="male", N=0,path_for_saving_images,savestats="YES")
do_weight_curve(data, value=rel.weight, value_label = "rel. BW",unit = "perc", batch="ALL", sex="female",N=0,path_for_saving_images,savestats="YES")
do_weight_curve(data, value=Weight,value_label = "Body Weight", unit = "g", batch="ALL", sex="female",N=0,path_for_saving_images,savestats="YES")
do_weight_curve(data, value=Weight,value_label = "Body Weight", unit = "g", batch="ALL", sex="male",N=0,path_for_saving_images,savestats="YES")

gc()

path_for_saving_images<-file.path(PATHS$general_data$FK49_output,"02_GeneratedData/Weight_Organs/background")
## Absolute Body Weight Single and combined Batches  --------------------------------------------------------
do_weight_curve(data, value=Weight,value_label = "Body Weight", unit = "g", batch=2, sex="male",N=0,path_for_saving_images)
do_weight_curve(data, value=Weight,value_label = "Body Weight", unit = "g", batch=1, sex="male",N=0,path_for_saving_images)
do_weight_curve(data, value=Weight,value_label = "Body Weight", unit = "g", batch=2, sex="female",N=0,path_for_saving_images)
do_weight_curve(data, value=Weight,value_label = "Body Weight", unit = "g", batch=1, sex="female",N=0,path_for_saving_images)


## Relative Body Weight Single Batches  --------------------------------------------------------
do_weight_curve(data, value=rel.weight, value_label = "rel. BW",unit = "perc", batch="1", sex="male",N=0,path_for_saving_images)
do_weight_curve(data, value=rel.weight, value_label = "rel. BW",unit = "perc", batch="1", sex="female",N=0,path_for_saving_images)
do_weight_curve(data, value=rel.weight, value_label = "rel. BW",unit = "perc", batch="2", sex="female",N=0,path_for_saving_images)
do_weight_curve(data, value=rel.weight, value_label = "rel. BW",unit = "perc", batch="2", sex="male",N=0,path_for_saving_images)
gc()
