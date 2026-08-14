rm(list=ls())
gc()
library(tidyr)
library(dplyr)
library(pheatmap)
library(tidyverse)
library(tibble)
library(rstatix)
library(lmerTest)
library(emmeans)
library(grid)
library(ggnewscale)
library(NADA2)
library(effsize)
source("FK49_Definitions.R")
# Read Raw Inputdata and general Data manipulation ------------------------------------------------------
ExpID="FK49"
rm=(list=ls())
gc()
library(tidyverse)
library(NADA2)
library(emmeans)
source("FK49_Definitions.R")
ExpID= "FK49"   # Decide if you want to load data from FK46 or FK49

if(ExpID == "FK49"){
  load(file = file.path(PATHS$exigo$FK49_input,  "FK49_Exigo_prepared.Rda"))
  output_pwd = file.path(PATHS$exigo$FK49_output)
  param_list = PARAMETERS$EXIGO$FK49_Exigo_Comprehensive_Panel
  }else if(ExpID == "FK46") {
  load(file = file.path(PATHS$exigo$FK46_input,  "FK46_Exigo_prepared.Rda"))
    output_pwd = file.path(PATHS$exigo$FK46_output)
    param_list = PARAMETERS$EXIGO$FK46_Exigo_Liver_Panel
  }else{print("Give me an exisiting Experiment ID to load the correct data from the correct path.")
  }
# Dotplot Exigo Panel --------------------------------------------------------------


do_Exigo <- function(inputdata, value, batch = "2", sex = "both",
                     y_title, path_images,
                     normal_range = NULL, lowlimit = NULL, hilimit = NULL,
                     p_value_override = NULL) {
  
  library(dplyr)
  library(ggplot2)
  library(NADA2)
  
  d <- inputdata %>%filter(complete.cases(.data[[value]]))
  
  # if you want to have a visual descriptive reference batch (bc i dont have a real one) 
  # write inthe function call reference_batch = "BH15, other wise r´write nothing or NULL
  censored_col  <- paste0(value, "_censored")
  direction_col <- paste0(value, "_direction")
  
  d <- d %>%
    mutate(
      value_numeric = as.numeric(.data[[value]]),
      censored = as.character(.data[[censored_col]]),
      direction = .data[[direction_col]],
      cens_logical = censored == "TRUE",
      censor_status_combined = case_when(
        censored == "TRUE" & direction == "<" ~ "Below LOD",
        censored == "TRUE" & direction == ">" ~ "Above ULOQ",
        TRUE ~ "Detected"
      )
    )
  

  
  #Generate limits for plotting from parameter values
  y_vals <- d$value_numeric
  y_max  <- max(y_vals, na.rm = TRUE)
  y_min  <- min(y_vals, na.rm = TRUE)
  y_pos <- y_max + 0.15 * (y_max - y_min)
  y_pos <- if (is.finite(y_max) & y_max > y_min) {y_max + 0.15 * (y_max - y_min)}
           else {y_max * 1.05}
  
  p1 <- ggplot(d, aes(x = Treatment, y = .data[[value]])) +
    stat_summary(fun = mean, geom = "bar", aes(fill = Treatment,color = Treatment),alpha = 0.5, width = 0.75) +
    scale_color_manual( values = c(Treatment_colors[c("Ctrl","TAM")]),
                        labels = c("Ctrl" = "Ctrl","TAM" = "TAM" ))+
    guides(color = guide_legend(order = 1, nrow = 2, byrow = TRUE)) +
    ggnewscale::new_scale_color() +
    stat_summary(fun.data = mean_sdl, fun.args = list(mult = 1), geom = "errorbar", width = 0.2, color = "black") +
    scale_fill_manual(name = "Treatment",
                      values = c(Treatment_colors[c("Ctrl","TAM")]),
                      labels = c("Ctrl" = "Ctrl","TAM" = "TAM" )) +
    guides(fill = guide_legend(order = 1, nrow = 2, byrow = TRUE)) +
    ggnewscale::new_scale_fill() +
    geom_point(aes(shape = Sex,  color = censor_status_combined), fill = "lightgrey",alpha=0.8,
               position = position_jitter(width = 0.15), size = 5.3, stroke = 1.8) +
    scale_shape_manual(name = "Sex", values = Sex_shape) +
    scale_color_manual(name = "Censoring", values = c("Below LOD" = "red",
                                                      "Above ULOQ" = "darkblue",
                                                      "Detected" = "black")) +
    scale_y_continuous(name = y_title, expand = expansion(mult = c(0.05, 0.15)))+
    # scale_x_discrete(labels = c("ctrl_analysis"  = "Ctrl",
    #   "TAM_analysis"   = "TAM",
    #   "ctrl_reference" = "Ctrl \n(ND)",
    #   "TAM_reference"  = "TAM \n(ND)"   )) +
    theme_minimal() +
    theme(legend.position = "bottom",
          legend.box = "vertical",
          legend.box.just = "left",
          legend.title = element_text(size = 11, face = "bold"),
          legend.text = element_text(size = 9),
          axis.line = element_line(color = "black", linewidth = 0.5),
          axis.ticks = element_line(color = "black", linewidth = 0.5),
          axis.title = element_text(size = 20, face = "bold"),
          axis.title.x = element_blank(),
          axis.text = element_text(size = 19, face = "bold"),
          panel.grid = element_blank() ) +
    guides(shape = guide_legend(title = "Sex", order = 2, nrow = 1, byrow = TRUE),
           color = guide_legend(title = "Censoring", order = 4, nrow = 1, byrow = TRUE),
            fill =guide_legend(title = "Group", order = 3, nrow = 1, byrow = TRUE, ))

  # SAFE reference lines
  if (!is.null(normal_range) && all(is.finite(normal_range)))
    p1 <- p1 + geom_hline(yintercept = normal_range, linetype = "dotted")
  
 # if (!is.null(lowlimit) && length(lowlimit) > 0 && is.finite(lowlimit[1]))
 #   p1 <- p1 + geom_hline(yintercept = lowlimit[1], linetype = "dashed", color = "red")
  
 # if (!is.null(hilimit) && length(hilimit) > 0 && is.finite(hilimit[1]))
 #   p1 <- p1 + geom_hline(yintercept = hilimit[1], linetype = "dashed", color = "blue")
  
  
  return(list(
    plot_raw = p1,
    y_pos = y_pos
  ))
}
stats <- read.csv2(file.path(output_pwd, "FK49_Exigo_Statistics.csv"))
plots <- lapply(param_list, function(p) {
  # corresponding statistical result
  stat_row <- stats %>%
    filter(parameter == p$value)
  # create plot
  res <- do_Exigo(
    inputdata = d1,
    value = p$value,
    batch = "ALL",
    sex = "both",
    y_title = p$y_title,
    path_images = output_pwd,
    normal_range = p$normal_range %||% NULL,
    lowlimit = p$lowlimit %||% NULL,
    hilimit = p$hilimit %||% NULL
  )
  
  # add adjusted p-value
  p_final <- res$plot_raw +
    annotate(
      "text",
      x = 1.5,
      y = res$y_pos,
      label = paste0(
        "adj p = ",
        format.pval(stat_row$p_adj, digits = 3)
      ),
      size = 5.5,
      fontface = "italic"
    )
  
  # filename
  fname_val <- gsub("[^[:alnum:]_]", "_", p$value)
  filename <- paste0(ExpID, "_", fname_val, "_Treatment.png")
  
  # save final plot
  ggsave(
    filename = filename,
    plot = p_final,
    path = output_pwd,
    width = 4,
    height = 11,
    dpi = 300
  )
  
  # return plot
  p_final
})
stats <- read.csv2(file.path(output_pwd, "BH15_Exigo_Statistics.csv"))

plots <- lapply(param_list, function(p) {
  # corresponding statistical result
  stat_row <- stats %>%
    filter(parameter == p$value)
  ExpID="BH15"
  # create plot
  res <- do_Exigo(
    inputdata = baseline_data,
    value = p$value,
    batch = "ALL",
    sex = "both",
    y_title = p$y_title,
    path_images = output_pwd,
    normal_range = p$normal_range %||% NULL,
    lowlimit = p$lowlimit %||% NULL,
    hilimit = p$hilimit %||% NULL
  )
  
  # add adjusted p-value
  p_final <- res$plot_raw +
    annotate(
      "text",
      x = 1.5,
      y = res$y_pos,
      label = paste0(
        "adj p = ",
        format.pval(stat_row$p_adj, digits = 3)
      ),
      size = 5.5,
      fontface = "italic"
    )
  
  # filename
  fname_val <- gsub("[^[:alnum:]_]", "_", p$value)
  filename <- paste0(ExpID, "_", fname_val, "_Treatment.png")
  
  # save final plot
  ggsave(
    filename = filename,
    plot = p_final,
    path = output_pwd,
    width = 4,
    height = 11,
    dpi = 300
  )
  
  # return plot
  p_final
})



# 
# 
# 
# gc()
# d_mat <- exp_and_BaseLine %>%select(where(is.numeric)) %>%
#   as.matrix()
# rownames(d_mat)<-exp_and_BaseLine$Animal
# d_mat<-t(d_mat)
# ann <- data.frame(Treatment = exp_and_BaseLine$Treatment,Sex=exp_and_BaseLine$Sex,Batch=exp_and_BaseLine$BATCH)
# rownames(ann) <- exp_and_BaseLine$Animal
# 
# #rotates transposes the matrix
# d_scaled <- t(scale(t(d_mat)))
# #row_annot <- data.frame(adj.p = round(padjusted[c("NASH_I","NASH_B","NASH_S",  "NASH_SAF")],5))
# #rownames(row_annot) <- c("NASH_I", "NASH_B","NASH_S", "NASH_SAF")
# cor_spearman <- cor(d_scaled, method = "spearman",use = "pairwise.complete.obs")
# p_cor<-pheatmap(cor_spearman,
#                 cluster_rows = T,
#                 cluster_cols =T,
#                 fontsize        = 10,
#                 fontsize_main   = 7,
#                 fontsize_row    = 9,
#                 fontsize_col    = 9,
#                 fontsize_number = 8,
#                 annotation_col  = ann,
#                 #annotation_row  = row_annot,
#                 color = colorRampPalette(c("white", "orange", "red"))(50),
#                 # annotation_colors = list(Treatment = c("Ctrl" = "#4D4D4DBF", "TAM" = "#8B0000BF")),
#                 main = paste0("Correlation HeatMap BaseLine vs ",ExpID),
#                 #labels_col = "",
#                 #labels_row = c("Inflammation","Ballooning","Steatosis",  "Total Score"),
#                 border_color = "black",
#                 cellwidth = 10,
#                 cellheight = 10,
#                 angle_col = 90,
#                 gaps_row = 3,
#                 display_numbers = FALSE,
#                 number_format = "%.3f"#,
#                 #legend_breaks = c(0,2,4),
#                 #legend_labels = c("0","2 ", "4")
# )
# p_cor
# ggsave(filename = "FK49_BH15_Exigo_Correlation.png", plot = p_cor, path = file.path(output_pwd,"background/"), 
#        width = 10, height = 10, dpi = 300,bg="white")
# 
# d_s<-pheatmap(d_scaled,
#               cluster_rows = T,
#               cluster_cols = T,
#               fontsize        = 10,
#               fontsize_main   = 7,
#               fontsize_row    = 9,
#               fontsize_col    = 9,
#               fontsize_number = 8,
#               annotation_col  = ann,
#               #annotation_row  = row_annot,
#               color = colorRampPalette(c("white", "orange", "red"))(50),
#               #annotation_colors = list(Treatment = c("Ctrl" = "#4D4D4DBF", "TAM" = "#8B0000BF")),
#               main = paste0("Scaled HeatMap BaseLine vs ",ExpID),
#               #labels_col = "",
#               #labels_row = c("Inflammation","Ballooning","Steatosis",  "Total Score"),
#               border_color = "black",
#               cellwidth = 10,
#               cellheight = 10,
#               angle_col = 90,
#               gaps_row = 3,
#               display_numbers = FALSE,
#               number_format = "%.3f",
#               legend_breaks = c(0,2,4),
#               legend_labels = c("0","2 ", "4"))
# ggsave(filename = paste0(ExpID,"_BH15_Exigo_scaled.png"), plot = d_s, path = paste0(output_pwd,"/background/"), 
#        width = 10, height = 10, dpi = 300,bg="white")

# if(ExpID=="FK46"){
#   d1 <- d1 %>% select(Treatment,  c("ALB","TP","GLOB","A.G","TB","GGT","AST","ALT","ALP","TC"))
#   baseline_data  <- baseline_data %>% select(Treatment,  c("ALB","TP","GLOB","A.G","TB","GGT","AST","ALT","ALP","TC"))
#   d  <- d %>% select(BATCH,  c("ALB","TP","GLOB","A.G","TB","GGT","AST","ALT","ALP","TC"))
#   params <- names(d %>% select( c("ALB","TP","GLOB","A.G","TB","GGT","AST","ALT","ALP","TC")))
#   
#   } else if(ExpID=="FK49"){
#     d1<-d1%>%select(Animal, Sex,Treatment,BATCH,-TV,Exigo_cols,Tumor.no.yes,Ascites.no.yes)
#     baseline_data<-baseline_data%>%select(Treatment,Exigo_cols)
#     d<-d%>%select(Animal, Sex,Treatment,BATCH,-TV,Exigo_cols,Tumor.no.yes,Ascites.no.yes)
#     params <- names(d %>% select( Exigo_cols))
#   } else
#     { print("Do you need to exclude some columns for perfect merging and data presentation?")}
# 
# effect_mean_d1 <- d1  %>%
#   group_by(Treatment) %>%  
#   summarise(across(everything(), mean, na.rm=TRUE)) %>%
#   pivot_longer(-Treatment) %>%
#   pivot_wider(names_from = Treatment, values_from = value) %>%
#   mutate(effect = TAM - Ctrl)
# 
# effect_fc_d1 <- effect_mean_d1 %>%
#   mutate(log2FC = log2(TAM / Ctrl))
# 
# effect_mean_baseline_data <- baseline_data %>%
#   group_by(Treatment) %>%  summarise(across(everything(), mean, na.rm=TRUE)) %>%
#   pivot_longer(-Treatment) %>%
#   pivot_wider(names_from = Treatment, values_from = value) %>%
#   mutate(effect = TAM - Ctrl)
# 
# effect_fc_baseline_data <- effect_mean_baseline_data %>%
#   mutate(log2FC = log2(TAM / Ctrl))
# 
# effect_mean_d <- d %>%
#   group_by(BATCH) %>%  
#   summarise(across(everything(), mean, na.rm=TRUE)) %>%
#   pivot_longer(-BATCH) %>%
#   pivot_wider(names_from = BATCH, values_from = value) %>%
#   mutate(effect = 2 - BH15)
# 
# effect_fc_d <- effect_mean_d %>%
#   mutate(log2FC = log2(2 / BH15))
# 
# results <- data.frame(
#   parameter = params,
#   log2FC_FK = effect_fc_d1$log2FC,
#   log2FC_BH = effect_fc_baseline_data$log2FC,
#   log2FC_BATCH = effect_fc_d$log2FC)
# 
# 
# mat <- cbind(results$log2FC_FK, results$log2FC_BH,results$log2FC_BATCH)
# 
# rownames(mat) <- results$parameter
# colnames(mat) <- c("log2FC_FK","log2FC_BH","log2FC_BATCH")
# breaks_lower  <- seq(-9, -0.3, length.out = 50)
# breaks_middle <- seq(-0.2, 0.2, length.out = 20)  
# breaks_upper  <- seq(0.3, 1.5, length.out = 30)  
# 
# my_breaks <- c(breaks_lower, breaks_middle[-1], breaks_upper[-1])
# colors_lower  <- colorRampPalette(c("navy","skyblue","#E6F0FA"))(length(breaks_lower)-1)
# colors_middle <- rep("white", length(breaks_middle)-1)
# colors_upper  <- colorRampPalette(c("#FFEDE6","#FFD1B0","#FFA07A"))(length(breaks_upper)-1)
# my_colors <- c(colors_lower, colors_middle, colors_upper)
# h<-pheatmap(
#   mat,
#   cluster_cols = FALSE,
#   cluster_rows = TRUE,
#   angle_col = 0,
#   color = my_colors,
#   breaks = my_breaks,
#   gaps_col = 2,
#   border_color = "black",
#   labels_col = c("CD-HFD","Baseline","CD-HFD vs Base"),
#   main = "TAM vs EtOH (log2 Fold Change)",
#   legend_breaks = c(-9, -6, -3, -1.5, -0.5, 0, 0.5, 1.0, 1.5),
#   legend_labels = c("-9","-6","-3","-1.5","-0.5","0","0.5","1.0","1.5")
# )
# ggsave(filename = paste0(ExpID,"_BH15_Exigo_FC_HEatMap.png"), plot = h, path = "02_GeneratedData/Exigo/background/", 
#        width = 5, height = 10, dpi = 300)
# 
# 
