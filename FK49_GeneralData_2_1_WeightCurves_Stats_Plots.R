rm(list=ls())
gc()
library(tidyverse)
library(openxlsx)
library(rstatix)
library(lmerTest)
library(emmeans)
library(ggpubr)
source("FK49_Definitions.R")

ExpId = "FK49" # Set to "FK49" or "FK46"

# Read Raw Inputdata after general Data manipulation ------------------------------------------------------
load(file.path(PATHS$general_data[[paste0(ExpId,"_output")]],"01_RawData",paste0(ExpId,"_Data_prepared.Rda")))

# Function for Weight Curves ----------------------------------------------
#Function assumes if you say BATCH== "ALL that you have batch 1 and 2. 
#If this is not right, change it in function to either take all numerical batches or to the numbers you have

do_weight_curve <- function(inputdata, value, value_label = NULL, unit = "g",batch = "ALL", sex = "both", N, path_images,savestats = "NO"){
  
  value_label_final <- if (is.null(value_label)) deparse(substitute(value)) else value_label #label for plotting and saving
  file_base <- paste0(ExpId, "_", value_label_final, "_Batch", batch, "_", sex, "_n", N)
  
  # filter for specified sex -----
  filtered <- inputdata %>%
    dplyr::select(any_of(c("Sex", "BATCH", "Treatment", "DOW", "wks_diet", "Animal", "Block", "days_diet", "Cage")), {{value}})%>%
    dplyr::filter(!is.na({{value}})) %>%
    dplyr::filter(case_when(
      sex == "female" ~ Sex == "female",
      sex == "male" ~ Sex == "male",
      sex == "both" ~ TRUE))
  
  # BATCH filtering -----
  if (batch == "ALL") {
    common_timepoints <- filtered %>%       # Find common time points across both batches
      dplyr::filter(BATCH %in% c(1, 2)) %>%
      group_by(wks_diet, BATCH) %>%
      summarise(n = n(), .groups = "drop") %>%
      group_by(wks_diet) %>%
      summarise(n_batches = n_distinct(BATCH)) %>%
     # dplyr::filter(n_batches ==2) %>%
      pull(wks_diet)
    
    filtered <- filtered %>%
      dplyr::filter(BATCH %in% c(1, 2)) %>%
      dplyr::filter(wks_diet %in% common_timepoints)
  } 
  else {
    filtered <- filtered %>%  dplyr::filter(BATCH == batch)
  }
  # Filter out food days (FK49 only — FK46 has no food/water daily weighing) -----
   #they should not appear in this overall plot. 
   #Here I only want weekly measurements not the Food/Water Intake Daily weights
  if ("Block" %in% names(filtered)) {
    filtered <- filtered %>%
    group_by(Animal, Block) %>%
      dplyr::filter((Block %in% c("0") & days_diet == -7) | # In block 0 and 1 all batches were weight on the Monday,
           (Block %in% c("1") & days_diet == 0)|    #DOW_in_Block 1 so thats the day i want to represent the week
           (!Block %in% c("0", "1")  & days_diet == as.numeric(as.character(Block))*7-4) )%>%
    ungroup()
  }
  # I summarize data from "filtered" dataset in "Mean_SD_data" to be able to plot the mean and sd later on  
  Mean_SD_data <- filtered %>%
    group_by(Treatment, wks_diet) %>%
    summarise(weight = mean({{value}}, na.rm = TRUE),n = n(), sd = sd({{value}}, na.rm = TRUE)) %>%
    dplyr::filter(n > N)
  
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
    dplyr::mutate( Animal = factor(Animal),  
            Treatment = factor(Treatment),  
            wks_diet = as.numeric(as.character(wks_diet)),   #wks_diet needs to be numeric for random slope
            Sex = factor(Sex))
  #fit linear mixed-effects model 
  if( sex == "both"){ 
  formula<-as.formula(paste(deparse(substitute(value)),
  "~ Treatment * wks_diet*Sex + (1 + wks_diet | Animal)")) }else {
  
  formula<-as.formula(paste(deparse(substitute(value)),
                            "~ Treatment * wks_diet + (1 + wks_diet | Animal)")) }
  
  #vfixed effects within interaction of treatment and time
  # (1+wks_diet|Animal) random effects
  #1 is random interept - each animal kann have individual basleine value at first measuremtn
  # wks diet means random slope each animal can have individual slope
  # |Animal data has repreated measures clustered within individual animals
                           
  model <- lmer(formula, REML = TRUE, data = model_data)
  anova_table<-anova(model, type = 3)
  get_anova_result <- function(term) {
    if (term %in% rownames(anova_table)) {
      paste0( term,   ": F = ", round(anova_table[term, "F value"], 2),   ", p = ", format.pval(anova_table[term, "Pr(>F)"], digits = 3)  )
    } else {NULL}
  }
  
  if (sex == "both") {
    
    anova_label <- paste(
      "ANOVA over linear mixed-effects model",
      get_anova_result("Treatment"),
      get_anova_result("wks_diet"),
      get_anova_result("Sex"),
      get_anova_result("Treatment:wks_diet"),
      get_anova_result("Treatment:Sex"),
      get_anova_result("wks_diet:Sex"),
      get_anova_result("Treatment:wks_diet:Sex"),
      sep = "\n"
    )
    
  } else {
    
    anova_label <- paste(
      "ANOVA over linear mixed-effects model",
      get_anova_result("Treatment"),
      get_anova_result("wks_diet"),
      get_anova_result("Treatment:wks_diet"),
      sep = "\n"
    )
  }
  # To test at each specific timepoint and not the overall dataset -----
  
  # estimated marginal means for both Treatment and wks_diet
  time_points <- sort(unique(model_data$wks_diet))
  emm <- emmeans(model, ~ Treatment | wks_diet, at = list(wks_diet = time_points))
  
  # pairwise contrasts between Treatment levels (ctrl vs TAM) at all levels of wks_diet
  pwc <- contrast(emm, method = "pairwise", adjust = "bonferroni")
  posthoc_label<- "Post Hoc: Pairwise with Bonferroni correction"
  
  pwc_df <- as.data.frame(pwc)
  pwc_df_rounded <- pwc_df %>%
    dplyr::mutate(rounded_p_value = ifelse(is.na(p.value), "NA", round(p.value, 3))) %>%
    dplyr::mutate(wks_diet = as.character(wks_diet)) %>%
    dplyr::mutate(wks_diet = as.numeric(wks_diet)) %>%
    dplyr::mutate(significance = case_when(
      is.na(p.value) ~ "NA",                     # For NA p-values
      p.value < 0.001 ~ "***",                   # p < 0.001 is highly significant
      p.value >= 0.001 & p.value < 0.01 ~ "**",  # 0.001 ≤ p < 0.01 is significant
      p.value >= 0.01 & p.value < 0.05 ~ "*",    # 0.01 ≤ p < 0.05 is moderately significant
      p.value >= 0.05 ~ "NS",                    # p ≥ 0.05 is not significant
      TRUE ~ "NA"   ))  %>%                      # Default case
    dplyr::select(wks_diet, rounded_p_value, significance)
  

  # Variables for Plot Setup
  mean_value <- mean(Mean_SD_data$weight, na.rm = TRUE)
  sd_value <-  sd(Mean_SD_data$weight, na.rm = TRUE)
  min_value <- round(mean_value - 3 * sd_value)
  max_value <- round(mean_value+4*sd_value)
  
  range_value <- max_value - min_value
  step_size <- ceiling((range_value * 0.2) / 5) * 5
  breaks_value <- seq(min_value, max_value, by = step_size)
  breaks_value <- round(breaks_value / 5) * 5
  min_x <- round(min(Mean_SD_data$wks_diet, na.rm = TRUE))
  max_x <- round(max(Mean_SD_data$wks_diet, na.rm = TRUE) )
  x_break_step <- if (max_x > 20) 4 else 1
  # Große Ticks
  breaks_x <- seq(0, max_x, by = x_break_step)
  # 41 als Tick behalten, aber nicht labeln
  breaks_x <- sort(unique(c(breaks_x, max_x)))
  # Labels nur für die regulären 4er-Schritte
  labels_x <- ifelse(breaks_x %% x_break_step == 0, breaks_x, "")
  x_break_step <- if (max_x > 20) 4 else 1
  # breaks_x <- seq(0, max_x, by = x_break_step)
  # breaks_x <- sort(unique(c(breaks_x, max_x)))
  unit_label <- unit
  
  max_weights <- Mean_SD_data %>%
    group_by(wks_diet) %>%
    summarise(max_weight = max(weight, na.rm = TRUE)) %>%
    ungroup() %>%
    dplyr::mutate(y_position = max_weight + max_weight * 0.15)
  
  n_annotations <- Mean_SD_data %>%
    group_by(wks_diet) %>%
    mutate( mean_diff = abs(weight - mean(weight)),  close = diff(range(weight)) < 0.07 * range_value) %>%
    ungroup() %>%
    group_by(wks_diet) %>%
    mutate( y_position = ifelse(close,max(weight + ifelse(is.na(sd), 0.02 * range_value, 0.5 * sd),  na.rm = TRUE) + 0.006 * range_value,
        weight + ifelse(is.na(sd), 0.02 * range_value, 0.5 * sd)),
        x_position = ifelse(  close & Treatment == "Ctrl",wks_diet - 0.29, ifelse( close & Treatment == "TAM", wks_diet + 0.29, wks_diet))) %>%
    ungroup()    
  # Join y_position back into pwc_df_rounded
  pwc_df_annotated <- pwc_df_rounded %>%
    left_join(max_weights, by = "wks_diet")
  
  # Plot
  plot <- ggplot(data = Mean_SD_data, aes(x = wks_diet, y = weight, color = Treatment, fill = Treatment)) +
    geom_ribbon(aes(y = weight, ymin = weight - sd, ymax = weight + sd), alpha = 0.1, linetype = 0) +
    geom_point(size = 3) +
    geom_line(linewidth = 1) +
    #geom_text(aes(label = n), hjust = 0, vjust = -1, size = 3, show.legend = FALSE) +
    scale_color_manual(values = c(Treatment_colors[c("Ctrl","TAM")],"black","pink")) +
    scale_fill_manual(values = c(Treatment_colors[c("Ctrl","TAM")],"black","pink")) +
    scale_x_continuous(
      name = "Time on CD-HFD [wks]",
      limits = c(min_x-1, max_x+1),
      breaks = breaks_x,
      labels = labels_x,
      minor_breaks = seq(min_x, max_x+1, by = 1),
      expand = expansion(mult = c(0.01, 0.01))
    )+
    scale_y_continuous(name = sprintf("%s [%s]", substitute(value_label_final), unit_label),
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
          axis.ticks.length = unit(4, "pt"),
          axis.title.x = element_text(size = 12, face = "bold", colour = "black"),
            axis.title.y = element_text(size = 12, face = "bold", colour = "black"),
            axis.text.x  = element_text(size = 10, face = "plain", colour = "black"),
            axis.text.y  = element_text(size = 10, face = "plain", colour = "black"),
            plot.title   = element_text(size = 12, face = "bold", colour = "black") ,
          legend.position = "top")+
    annotate("text",
             x = pwc_df_annotated$wks_diet,
             y = pwc_df_annotated$y_position,
             label = pwc_df_annotated$significance,
             size = 2.5, color = "black", fontface = "italic")+
    geom_text(data = n_annotations,   aes(x = x_position, y = y_position,   label = n  , color = Treatment),
      hjust = 0.5,  vjust = 0,  size = 2.5,  show.legend = FALSE )#+
    # annotate("text", x = min_x + 2, y = min_value + range_value * 0.05,  # 5% above the bottom
    #          label = anova_label, size = 2,   hjust = 0,  color = "black", 
    #          fontface = "italic")+
    # annotate("text", x = min_x + 10, y = min_value + range_value * 0.05,  # 5% above the bottom
    #          label = posthoc_label,size = 2, hjust = 0,  color = "black",      fontface = "italic")
    # 
  
  # ---  Saving Plot --- 
  ggsave(filename = paste0(file_base, ".png"), plot = plot,  path = path_images, width = 8, height = 5,dpi = 300)
  #ggsave(filename = paste0(file_base, ".pdf"), plot = plot, path = path_images, width = 9, height = 6, dpi = 300, device = cairo_pdf)
  #dev.off()
  
  # Optionally save stats tables
  if (savestats == "YES") {
    outliers <- outliers %>%  dplyr::mutate(wks_diet = factor(wks_diet), table = "Outliers")
    anova_table <- anova(model, type = 3) %>% as.data.frame() %>% tibble::rownames_to_column(var = "Term")    
    pwc_df <- pwc_df %>%dplyr::mutate(wks_diet = factor(wks_diet),table = "Pairwise Comparison")
    StatsOutput <- bind_rows(outliers,anova_table,pwc_df)%>% relocate( table)
    write.csv2( StatsOutput,
                file = file.path(paste0(path_images), paste0(file_base, "_StatsOutput.csv")),
                row.names = FALSE,  na = "",  fileEncoding = "UTF-8"  )
    openxlsx::write.xlsx( StatsOutput,  file = file.path(path_images, paste0(file_base, "_StatsOutput.xlsx")),  overwrite = TRUE)
    pdf(file.path(path_images, paste0(file_base, "_StatsOutput.pdf")),  width = 12,   height = 8)
    gridExtra::grid.table(StatsOutput)
    dev.off()
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
path_for_saving_images<-file.path(PATHS$general_data[[paste0(ExpId,"_output")]],"02_GeneratedData/Weight_Organs")

## Relative Weight as summary from Batch 1 and 2 together --------------------------------------------------------
do_weight_curve(data, value=rel.weight, value_label = "rel.BW", unit = "perc", batch="ALL", sex="male", N=0,path_for_saving_images,savestats="YES")
do_weight_curve(data, value=rel.weight, value_label = "rel.BW",unit = "perc", batch="ALL", sex="female",N=0,path_for_saving_images,savestats="YES")
do_weight_curve(data, value=Weight,value_label = "Body Weight", unit = "g", batch="ALL", sex="female",N=0,path_for_saving_images,savestats="YES")
do_weight_curve(data, value=Weight,value_label = "Body Weight", unit = "g", batch="ALL", sex="male",N=0,path_for_saving_images,savestats="YES")
do_weight_curve(data, value=rel.weight, value_label = "rel.BW",unit = "perc", batch="ALL", sex="both",N=0,path_for_saving_images,savestats="YES")

gc()

path_for_saving_images<-file.path(PATHS$general_data[[paste0(ExpId,"_output")]],"02_GeneratedData/Weight_Organs/background")
## Absolute Body Weight Single and combined Batches  --------------------------------------------------------
do_weight_curve(data, value=Weight,value_label = "Body Weight", unit = "g", batch=2, sex="male",N=0,path_for_saving_images)
do_weight_curve(data, value=Weight,value_label = "Body Weight", unit = "g", batch=1, sex="male",N=0,path_for_saving_images)
do_weight_curve(data, value=Weight,value_label = "Body Weight", unit = "g", batch=2, sex="female",N=0,path_for_saving_images)
do_weight_curve(data, value=Weight,value_label = "Body Weight", unit = "g", batch=1, sex="female",N=0,path_for_saving_images)


## Relative Body Weight Single Batches  --------------------------------------------------------
do_weight_curve(data, value=rel.weight, value_label = "rel. body weight",unit = "perc", batch="1", sex="male",N=0,path_for_saving_images)
do_weight_curve(data, value=rel.weight, value_label = "rel. body weight",unit = "perc", batch="1", sex="female",N=0,path_for_saving_images)
do_weight_curve(data, value=rel.weight, value_label = "rel. body weight",unit = "perc", batch="2", sex="female",N=0,path_for_saving_images)
do_weight_curve(data, value=rel.weight, value_label = "rel. body weight",unit = "perc", batch="2", sex="male",N=0,path_for_saving_images)
gc()

