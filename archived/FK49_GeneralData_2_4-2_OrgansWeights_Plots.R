rm(list=ls())
gc()
library(tidyverse)
source("FK49_Definitions.R")


output_pwd <-file.path(PATHS$Organs$output)
load(file.path(PATHS$Organs$input,"FK49_Data_prepared.Rda"))
data<- data%>%
 select(Animal, Sex, Treatment, Weight, Liver, Fat, Spleen, Liver_rel, Fat_rel, Spleen_rel,
        Ascites.no.yes, Tumor.no.yes, wks_diet, BATCH)%>%
 filter(!Animal=="EC1",!Animal=="EC2",!Animal=="EC2") 
stats<-read.csv2(file = file.path(output_pwd,"Statistics/FK49_Organ_Weight_Statistics.csv"))


# Dotplot Organ Weights --------------------------------------------------------------

do_organ_weight <- function(inputdata, statistics = NULL,value, 
                            y_title, path_images,colors = Treatment_colors[c("Ctrl","TAM")]) {
  ### Data manipulation -----------------------------------------------------------------
  d <- inputdata %>%
    filter(complete.cases(.data[[value]]))
  
  # Convert selected value to numeric if not already
  d[[value]] <- as.numeric(d[[value]])
  
  ### Summary stats ----------------------------------------------------------------------
  plotting_stats <- d %>%group_by(Treatment) %>%summarise(Mean = mean(.data[[value]], na.rm = TRUE),SD = sd(.data[[value]], na.rm = TRUE),.groups = "drop")
  
  p_value_label <- statistics %>% filter(Variable == value) %>%  pull(p_Treatment_adj)%>%as.numeric()
  p_interaction <- statistics %>% filter(Variable == value) %>%  pull(p_Interaction_adj)%>%as.numeric()
  
  p_value_label <-ifelse(p_value_label<0.05,paste0("p = ",round(p_value_label,3)),"ns") # if not sig just print ns
  print(p_value_label)
  print(p_interaction)
  
  ### Plot 1 -----
  p1 <- ggplot(plotting_stats, aes(x = Treatment, y = Mean, fill = Treatment)) +
    geom_bar(stat = "identity", color = "black", alpha = 0.5, width = 0.75, position = "dodge") +
    geom_point(data = d,  aes(y = .data[[value]], shape = Sex), fill = "lightgrey", color = "black",    
               position =  position_jitterdodge(0.1, dodge.width = 0.75), alpha = 0.8, size = 5.3, stroke = 1.8)+
    geom_errorbar(aes(ymin = Mean - SD, ymax = Mean + SD),  position = position_dodge(width = 0.75), width = 0.2) +
    scale_fill_manual(values = colors, labels = c("Ctrl", "TAM")) +
    scale_shape_manual(values = Sex_shape) +
    scale_y_continuous(name = y_title) +
    labs(x = "Treatment") +
    annotate("text", x = 1.5, y = max(d[[value]], na.rm = TRUE) * 1.1, label = p_value_label, size = 6, fontface = "italic") +
    theme_minimal() +
    theme(legend.position = "bottom",
      axis.line = element_line(color = "black", linewidth = 0.5),
      axis.ticks = element_line(color = "black", linewidth = 0.5),
      axis.title = element_text(size = 20, face = "bold"),
      axis.title.x = element_blank(),
      axis.text = element_text(size = 19, face = "bold"),
      plot.title = element_blank(),
      legend.title = element_text(size = 10),
      legend.text = element_text(size = 10),
      panel.grid = element_blank()) +
    guides( shape = guide_legend(title = "Sex", order = 2 ,nrow = 2, byrow = TRUE),
      fill = guide_legend(title = "Treatment", order = 3,override.aes = list(shape = 21), nrow =2, byrow = TRUE),
      color = "none")
  
  ### Save Plots
  value_clean <- gsub("[^[:alnum:]_]", "_", value)
  filename1 <- paste0("FK49_", value_clean,".png")
  ggsave(filename = filename1, plot = p1, path = path_images, width = 4, height = 11, dpi = 300)
  
# If there is Treatment Sex Interaction, then print plots for sexes also individually
# #since i dont have that outcome i did not bother to include the automatic ptinting of the correct p value per sex
  if (!is.na(p_interaction) &p_interaction < 0.05) {
     plotting_stats <- d %>%group_by(Treatment,Sex) %>%
      summarise(Mean = mean(.data[[value]], na.rm = TRUE),
                SD = sd(.data[[value]], na.rm = TRUE),.groups = "drop")
     #p_value_label <- statistics %>% filter(Variable == value) %>%  pull(p_Treatment_adj)%>%as.numeric(round(3))
     p1 <- ggplot(plotting_stats, aes(x = Treatment, y = Mean, fill = Treatment)) +
      facet_wrap(~Sex)+
      geom_bar(stat = "identity", color = "black", alpha = 0.5, width = 0.75, position = "dodge") +
      geom_point(data = d, fill = "lightgrey", color = "black",
                 aes(y = .data[[value]], shape =Sex),
                 position =  position_jitterdodge(0.1, dodge.width = 0.75), alpha = 0.8, size = 5.3, stroke = 1.8)+
      geom_errorbar(aes(ymin = Mean - SD, ymax = Mean + SD),  position = position_dodge(width = 0.75), width = 0.2) +
      scale_fill_manual(values = colors, labels = c("Ctrl", "TAM")) +
      scale_shape_manual(values = Sex_shape) +
      scale_y_continuous(name = y_title) +
      labs(x = "Treatment") +
     # annotate("text", x = 1.5, y = max(d[[value]], na.rm = TRUE) * 1.1, label = p_value_label, size = 6, fontface = "italic") +
      theme_minimal() +
      theme(legend.position = "bottom",
            axis.line = element_line(color = "black", linewidth = 0.5),
            axis.ticks = element_line(color = "black", linewidth = 0.5),
            axis.title = element_text(size = 20, face = "bold"),
            axis.title.x = element_blank(),
            axis.text = element_text(size = 19, face = "bold"),
            plot.title = element_blank(),
            legend.title = element_text(size = 10),
            legend.text = element_text(size = 10),
            panel.grid = element_blank()) +
      guides( shape = guide_legend(title = "Sex", order = 2 ,nrow = 2, byrow = TRUE),
              fill = guide_legend(title = "Treatment", order = 3,override.aes = list(shape = 21), nrow =2, byrow = TRUE),
              color = "none")
    filename2 <- paste0("FK49_", value_clean, "_Treatment_by_sex.png")
    ggsave(filename = filename2, plot = p1, path = path_images, width = 4, height = 11, dpi = 300)
  }else{}

}

## Call the function with desired arguments --------------------------------------------------------------
do_organ_weight(data, value = "Liver_rel", statistics=stats,  y_title = "Liver/BW [%]",output_pwd)
do_organ_weight(data, value = "Spleen_rel",statistics=stats,y_title= "Spleen/BW [%]",output_pwd)
do_organ_weight(data, value = "Spleen",statistics=stats,y_title= "Spleen [mg]",output_pwd)
do_organ_weight(data, value = "Liver",statistics=stats,y_title= "Liver [g]",output_pwd)
do_organ_weight(data, value = "Fat",statistics=stats,y_title= "Fat [g]",output_pwd)
do_organ_weight(data, value = "Fat_rel",statistics=stats,y_title= "Fat/BW [%]",output_pwd)


