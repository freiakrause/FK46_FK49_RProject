library(knitr)
library(tidyr)
library(dplyr)
library(pheatmap)
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
# Read Raw Inputdata and general Data manipulation ------------------------------------------------------
ExpId="FK49"

if (ExpId=="FK49") {
  setwd("C:/Users/b1084855/OneDrive - Universität Salzburg/Freigegebene Dokumente - AG_Tumorimmunologie/Data/Freia Krause/01_Experiments/FK49_CD-HFD_13wks/FK49_Analysis")
  Exigo_cols= c("ALB","TP","GLOB","A.G","TB","GGT","AST","ALT","ALP","AMY","Crea","UA","BUN","GLU","TC","TG")
  }  else if (ExpId == "FK46"){
  setwd("C:/Users/b1084855/OneDrive - Universität Salzburg/Freigegebene Dokumente - AG_Tumorimmunologie/Data/Freia Krause/01_Experiments/FK46_iALmice_high Fat diet 52 weeks 7d after injection/Analysis")
    Exigo_cols=c("ALB","TP","GLOB","A.G","TB","GGT","AST","ALT","ALP","TBA","TC")
    } else{
  print("Let me set a folder Path and define ExpId and Folderpath")}
load(paste0("01_RawData/",ExpId,"_Data_prepared.Rda"))

d1<-data%>%#%>%select(Animal, Sex,Treatment,BATCH,-TV,Exigo_cols,Tumor.no.yes,Ascites.no.yes) %>%
  filter(!is.na(ALB))%>%mutate(Animal=as.character(Animal)) #%>%
  # arrange(Treatment)

#rm(data)
load("C:/Users/b1084855/OneDrive - Universität Salzburg/Freigegebene Dokumente - AG_Tumorimmunologie/Data/Freia Krause/01_Experiments/FK49_CD-HFD_13wks/FK49_Analysis/01_RawData/BH15_Data_prepared.Rda")
baseline_data<-data%>%select(Animal, Sex,Treatment,BATCH,-TV,colnames(data)[131:146]) %>%
  filter(!is.na(ALB))%>%mutate(Animal=as.character(Animal)) #%>%
  # arrange(Treatment)
d1<-d1%>%mutate(BATCH=as.character(BATCH))
exp_and_BaseLine<-bind_rows(d1,baseline_data)
#   arrange(Treatment)
if(ExpId=="FK46"){exp_and_BaseLine<-exp_and_BaseLine%>%select(-c(TBA,AMY,Crea,UA,BUN,GLU,TG))
  } else if(ExpId=="FK49"){exp_and_BaseLine=exp_and_BaseLine
    } else
    { print("Do you need to exclude some columns for perfect merging and data presentation?")}

gc()
d_mat <- exp_and_BaseLine %>%select(where(is.numeric)) %>%
  as.matrix()
rownames(d_mat)<-exp_and_BaseLine$Animal
d_mat<-t(d_mat)
ann <- data.frame(Treatment = exp_and_BaseLine$Treatment,Sex=exp_and_BaseLine$Sex,Batch=exp_and_BaseLine$BATCH)
rownames(ann) <- exp_and_BaseLine$Animal

 #rotates transposes the matrix
d_scaled <- t(scale(t(d_mat)))
#row_annot <- data.frame(adj.p = round(padjusted[c("NASH_I","NASH_B","NASH_S",  "NASH_SAF")],5))
#rownames(row_annot) <- c("NASH_I", "NASH_B","NASH_S", "NASH_SAF")
cor_spearman <- cor(d_scaled, method = "spearman",use = "pairwise.complete.obs")
p_cor<-pheatmap(cor_spearman,
         cluster_rows = T,
         cluster_cols =T,
         fontsize        = 10,
         fontsize_main   = 7,
         fontsize_row    = 9,
         fontsize_col    = 9,
         fontsize_number = 8,
         annotation_col  = ann,
         #annotation_row  = row_annot,
         color = colorRampPalette(c("white", "orange", "red"))(50),
         annotation_colors = list(Treatment = c("ctrl" = "#4D4D4DBF", "TAM" = "#8B0000BF")),
         main = paste0("Correlation HeatMap BaseLine vs ",ExpId),
         #labels_col = "",
         #labels_row = c("Inflammation","Ballooning","Steatosis",  "Total Score"),
         border_color = "black",
         cellwidth = 10,
         cellheight = 10,
         angle_col = 90,
         gaps_row = 3,
         display_numbers = FALSE,
         number_format = "%.3f",
         #legend_breaks = c(0,2,4),
         #legend_labels = c("0","2 ", "4")
         )
ggsave(filename = paste0(ExpId,"_BH15_Exigo_Correlation.png"), plot = p_cor, path = "02_GeneratedData/Exigo/background/", 
       width = 5, height = 10, dpi = 300)

d_s<-pheatmap(d_scaled,
              cluster_rows = F,
              cluster_cols = F,
              fontsize        = 10,
              fontsize_main   = 7,
              fontsize_row    = 9,
              fontsize_col    = 9,
              fontsize_number = 8,
              annotation_col  = ann,
              #annotation_row  = row_annot,
              color = colorRampPalette(c("white", "orange", "red"))(50),
              annotation_colors = list(Treatment = c("ctrl" = "#4D4D4DBF", "TAM" = "#8B0000BF")),
              main = paste0("Scaled HeatMap BaseLine vs ",ExpId),
              #labels_col = "",
              #labels_row = c("Inflammation","Ballooning","Steatosis",  "Total Score"),
              border_color = "black",
              cellwidth = 10,
              cellheight = 10,
              angle_col = 90,
              gaps_row = 3,
              display_numbers = FALSE,
              number_format = "%.3f",
              legend_breaks = c(0,2,4),
              legend_labels = c("0","2 ", "4"))
ggsave(filename = paste0(ExpId,"_BH15_Exigo_scaled.png"), plot = d_s, path = "02_GeneratedData/Exigo/background/", 
       width = 5, height = 10, dpi = 300)

# if(ExpId=="FK46"){
#   d1 <- d1 %>% select(Treatment,  c("ALB","TP","GLOB","A.G","TB","GGT","AST","ALT","ALP","TC"))
#   baseline_data  <- baseline_data %>% select(Treatment,  c("ALB","TP","GLOB","A.G","TB","GGT","AST","ALT","ALP","TC"))
#   d  <- d %>% select(BATCH,  c("ALB","TP","GLOB","A.G","TB","GGT","AST","ALT","ALP","TC"))
#   params <- names(d %>% select( c("ALB","TP","GLOB","A.G","TB","GGT","AST","ALT","ALP","TC")))
#   
#   } else if(ExpId=="FK49"){
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
#   mutate(effect = TAM - ctrl)
# 
# effect_fc_d1 <- effect_mean_d1 %>%
#   mutate(log2FC = log2(TAM / ctrl))
# 
# effect_mean_baseline_data <- baseline_data %>%
#   group_by(Treatment) %>%  summarise(across(everything(), mean, na.rm=TRUE)) %>%
#   pivot_longer(-Treatment) %>%
#   pivot_wider(names_from = Treatment, values_from = value) %>%
#   mutate(effect = TAM - ctrl)
# 
# effect_fc_baseline_data <- effect_mean_baseline_data %>%
#   mutate(log2FC = log2(TAM / ctrl))
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
# ggsave(filename = paste0(ExpId,"_BH15_Exigo_FC_HEatMap.png"), plot = h, path = "02_GeneratedData/Exigo/background/", 
#        width = 5, height = 10, dpi = 300)
# 
# 

Exigo_Comprehensive_Panel <- list(
  list(value="ALB",  y_title="Alb [g/L]",   normal_range=c(20,48), lowlimit=2),
  list(value="TP",   y_title="TP [g/L]",    normal_range=c(36,66)),
  list(value="GLOB", y_title="GLOB [g/L]"),
  list(value="A.G",  y_title="A/G"),
  list(value="TB",   y_title="TB [µmol/L]", normal_range=c(1,15), lowlimit=0.1),
  list(value="GGT",  y_title="GGT [U/L]",   lowlimit=2),
  list(value="AST",  y_title="AST [U/L]",   normal_range=c(59,247), hilimit=650, lowlimit=5),
  list(value="ALT",  y_title="ALT [U/L]",   normal_range=c(28,132)),
  list(value="ALP",  y_title="ALP [U/L]",   normal_range=c(62,209), lowlimit=5),
  list(value="AMY",  y_title="AMY [U/L]",   normal_range=c(1691,3615)),
  list(value="Crea", y_title="Crea [?]",    normal_range=c(12,71)),
  list(value="UA",   y_title="UA [µmol/L]", normal_range=c(101,321), lowlimit=10),
  list(value="BUN",  y_title="BUN [mmol/L]",normal_range=c(4,11.8)),
  list(value="GLU",  y_title="GLU [mmol/L]",normal_range=c(5,10.67)),
  list(value="TC",   y_title="TC [mmol/L]", normal_range=c(0.93,4.04)),
  list(value="TG",   y_title="TG [mmol/L]", normal_range=c(0.62,1.63)))

Exigo_Liver_Panel <- list(
  list(value="ALB",  y_title="Alb [g/L]",     normal_range=c(20,48), lowlimit=2),
  list(value="TP",   y_title="TP [g/L]",      normal_range=c(36,66)),
  list(value="GLOB", y_title="GLOB [g/L]"),
  list(value="A.G",  y_title="A/G"),
  list(value="TB",   y_title="TB [µmol/L]",   normal_range=c(0,15),  lowlimit=0.1),
  list(value="GGT",  y_title="GGT [U/L]",     lowlimit=2),
  list(value="AST",  y_title="AST [U/L]",     normal_range=c(59,247), hilimit=650, lowlimit=5),
  list(value="ALT",  y_title="ALT [U/L]",     normal_range=c(28,132)),
  list(value="ALP",  y_title="ALP [U/L]",     normal_range=c(62,209), lowlimit=5),
  list(value="TBA",  y_title="TBA [µmol/L]",  lowlimit=1),
  list(value="TC",   y_title="TC [mmol/L]",   normal_range=c(0.93,4.04))
)


if (ExpId=="FK49") {
  param_list=Exigo_Comprehensive_Panel
}  else if (ExpId == "FK46"){
  param_list=Exigo_Liver_Panel} else{
  print("I don't know which Exigo you used.")}
# Dotplot Exigo Panel --------------------------------------------------------------


do_Exigo <- function(inputdata, value, batch = "2", sex = "both",
                     y_title, path_images,
                     normal_range = NULL, lowlimit = NULL, hilimit = NULL,
                     p_value_override = NULL,
                     reference_batch = NULL) {
  
  library(dplyr)
  library(ggplot2)
  library(NADA2)
  
  # -------------------------
  # Base filtering
  # -------------------------
  
  d <- inputdata %>%
    filter(complete.cases(.data[[value]]))
  
  if (!identical(batch, "ALL")) {
    d <- d %>% filter(BATCH %in% batch)
  }
  if (sex != "both") {
    d <- d %>% filter(Sex == sex)
  }
  
  # -------------------------
  # Status variables
  # -------------------------
  
  d <- d %>%
    mutate(
      Tumor.no.yes   = as.factor(Tumor.no.yes),
      Ascites.no.yes = as.factor(Ascites.no.yes),
      event_status = case_when(
        Tumor.no.yes == "0" & Ascites.no.yes == "0" ~ "normal",
        Tumor.no.yes == "1" & Ascites.no.yes == "0" ~ "tumor",
        Tumor.no.yes == "0" & Ascites.no.yes == "1" ~ "ascites",
        Tumor.no.yes == "1" & Ascites.no.yes == "1" ~ "both",
        TRUE ~ "unknown"
      ),
      event_status = factor(
        event_status,
        levels = c("normal","ascites","tumor","both","unknown")
      )
    )
  
  # -------------------------
  # Batch roles
  # -------------------------
  
  d <- d %>%
    mutate(
      BatchRole = case_when(
        !is.null(reference_batch) & BATCH %in% reference_batch ~ "reference",
        TRUE ~ "analysis"
      ),
      TreatmentBatch = interaction(Treatment, BatchRole, sep = "_")
    )
  
  # -------------------------
  # Censoring variables
  # -------------------------
  
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
  
  d_stat <- d %>% filter(BatchRole == "analysis")
  d_plot <- d
  
  # -------------------------
  # Helper: check cen2means safety
  # -------------------------
  
  check_data_sufficient <- function(data) {
    unc <- data %>% filter(!cens_logical)
    if (nrow(unc) < 4) return(FALSE)
    unc %>%
      group_by(Treatment) %>%
      summarise(n_dist = n_distinct(value_numeric), .groups = "drop") %>%
      pull(n_dist) %>%
      all(. >= 2)
  }
  
  # -------------------------
  # Statistical test
  # -------------------------
  
  p_value <- NA
  
  if (nrow(d_stat) >= 2) {
    
    has_censored <- any(d_stat$cens_logical, na.rm = TRUE)
    
    if (has_censored) {
      
      if (!all(d_stat$cens_logical, na.rm = TRUE) &&
          check_data_sufficient(d_stat)) {
        
        cen_result <- suppressWarnings(
          with(d_stat,
               cen2means(value_numeric, cens_logical, group = Treatment))
        )
        p_value <- cen_result$pval
        
      } else {
        p_value <- NA
      }
      
    } else {
      
      x_ctrl <- d_stat$value_numeric[d_stat$Treatment == "ctrl"]
      x_tam  <- d_stat$value_numeric[d_stat$Treatment == "TAM"]
      
      s1 <- if (length(x_ctrl) > 3) shapiro.test(x_ctrl)$p.value else 0
      s2 <- if (length(x_tam)  > 3) shapiro.test(x_tam)$p.value  else 0
      
      test_result <- if (s1 < 0.05 | s2 < 0.05) {
        suppressWarnings(
          wilcox.test(value_numeric ~ Treatment,
                      data = d_stat, exact = FALSE)
        )
      } else {
        t.test(value_numeric ~ Treatment, data = d_stat)
      }
      
      p_value <- test_result$p.value
    }
  }
  
  if (!is.null(p_value_override)) {
    p_value <- p_value_override
  }
  
  # -------------------------
  # Plot
  # -------------------------
  
  y_vals <- d_stat$value_numeric
  y_max  <- max(y_vals, na.rm = TRUE)
  y_min  <- min(y_vals, na.rm = TRUE)
    y_pos <- y_max + 0.15 * (y_max - y_min)
  y_pos <- if (is.finite(y_max) && y_max > y_min) {
    y_max + 0.15 * (y_max - y_min)
  } else {
    y_max * 1.05
  }
  
  p1 <- ggplot(d_plot,   aes(x = TreatmentBatch, y = .data[[value]])) +
    stat_summary(fun = mean, geom = "bar", aes(fill = TreatmentBatch,color = TreatmentBatch),alpha = 0.5, width = 0.75) +
    scale_color_manual( values = c( "ctrl_analysis"  = "black",   "TAM_analysis"   = "black",    "ctrl_reference" = "grey", "TAM_reference"  = "grey"),
                        labels = c(   "ctrl_analysis"  = "Control",  "TAM_analysis"   = "TAM","ctrl_reference" = "Control (ND)",  "TAM_reference"  = "TAM (ND)" ))+
    guides(color = guide_legend(order = 1, nrow = 2, byrow = TRUE)) +
    ggnewscale::new_scale_color() +
    stat_summary(fun.data = mean_sdl, fun.args = list(mult = 1), geom = "errorbar", width = 0.2, color = "black") +
    scale_fill_manual(name = "Treatment",
      values = c(  "ctrl_analysis"  = "#4D4D4DBF", "TAM_analysis"   = "#8B0000BF",  "ctrl_reference" = "grey90",  "TAM_reference"  = "#ffcccc"      ),
      labels = c(   "ctrl_analysis"  = "Control",  "TAM_analysis"   = "TAM","ctrl_reference" = "Control (ND)",  "TAM_reference"  = "TAM (ND)" )) +
    guides(fill = guide_legend(order = 1, nrow = 2, byrow = TRUE)) +
    ggnewscale::new_scale_fill() +
    geom_point(aes(shape = event_status,  color = censor_status_combined), fill = "lightgrey",alpha=0.8,
               position = position_jitter(width = 0.15), size = 5.3, stroke = 1.8) +
    scale_shape_manual(name = "Status", values = c(21, 22, 24, 25, 26)) +
    scale_color_manual(name = "Censoring", values = c("Below LOD" = "red",
                                                      "Above ULOQ" = "darkblue",
                                                      "Detected" = "black")) +
    scale_y_continuous(name = y_title, expand = expansion(mult = c(0.05, 0.15)))+
    scale_x_discrete(labels = c("ctrl_analysis"  = "Ctrl",
      "TAM_analysis"   = "TAM",
      "ctrl_reference" = "Ctrl \n(ND)",
      "TAM_reference"  = "TAM \n(ND)"   )) +
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
    guides(shape = guide_legend(title = "Status", order = 2, nrow = 1, byrow = TRUE),
           color = guide_legend(title = "Censoring", order = 4, nrow = 1, byrow = TRUE),
            fill =guide_legend(title = "Group", order = 3, nrow = 1, byrow = TRUE, ))

  # SAFE reference lines
  if (!is.null(normal_range) && all(is.finite(normal_range)))
    p1 <- p1 + geom_hline(yintercept = normal_range, linetype = "dotted")
  
 # if (!is.null(lowlimit) && length(lowlimit) > 0 && is.finite(lowlimit[1]))
 #   p1 <- p1 + geom_hline(yintercept = lowlimit[1], linetype = "dashed", color = "red")
  
 # if (!is.null(hilimit) && length(hilimit) > 0 && is.finite(hilimit[1]))
 #   p1 <- p1 + geom_hline(yintercept = hilimit[1], linetype = "dashed", color = "blue")
  
  fname_val <- gsub("[^[:alnum:]_]", "_", value)
  filename  <- paste0(ExpId, "_", fname_val, "_Treatment.png")
  if (is.null(reference_batch)) {
    ggsave(filename, plot = p1, path = path_images, width = 4, height = 11, dpi = 300)
  } else {
    p1 <-p1+
      geom_vline(xintercept = 2.5,linetype = "dashed",colour = "lightgrey")
    ggsave(filename, plot = p1, path = path_images, width = 6, height = 11, dpi = 300)
  }
  return(list(
    p_value_raw = p_value,
    plot_raw = p1
  ))
}
run_all_LP_onecall <- function(inputdata, param_list,
                               batch = "2",
                               sex = "both",
                               path_images,
                               p_adjust_method = "fdr",
                               reference_batch = NULL) {
  
  library(dplyr)
  
  n <- length(param_list)
  
  raw_pvals <- numeric(n)
  raw_plots <- vector("list", n)
  
  message("STEP 1: running Exigo analysis")
  
  for (i in seq_along(param_list)) {
    
    p <- param_list[[i]]
    
    if (!(p$value %in% names(inputdata))) {
      raw_pvals[i]  <- NA
      raw_plots[[i]] <- NULL
      next
    }
    
    res <- do_Exigo(
      inputdata       = inputdata,
      value           = p$value,
      batch           = batch,
      sex             = sex,
      y_title         = p$y_title,
      path_images     = path_images,
      normal_range    = p$normal_range %||% NULL,
      lowlimit        = p$lowlimit %||% NULL,
      hilimit         = p$hilimit %||% NULL,
      reference_batch = reference_batch
    )
    
    raw_pvals[i]  <- res$p_value_raw
    raw_plots[[i]] <- res$plot_raw
  }
  
  message("STEP 2: adjusting p-values")
  
  adj_pvals <- p.adjust(raw_pvals, method = p_adjust_method)
  
  message("STEP 3: annotating plots")
  
  for (i in seq_along(param_list)) {
    
    p <- param_list[[i]]
    
    y_vals <- inputdata[[p$value]]
    y_max  <- max(y_vals, na.rm = TRUE)
    y_min  <- min(y_vals, na.rm = TRUE)
    
    y_pos <- if (is.finite(y_max) && y_max > y_min) {
      y_max + 0.15 * (y_max - y_min)
    } else {
      y_max * 1.05
    }
    
    raw_plots[[i]] <- raw_plots[[i]] +
      annotate("text",
               x = 1.5,
               y = y_pos,
               label = paste0("adj p = ",
                              format.pval(adj_pvals[i], digits = 3)),
               size = 5.5,
               fontface = "italic")
    
    fname_val <- gsub("[^[:alnum:]_]", "_", p$value)
    filename  <- paste0(ExpId, "_", fname_val, "_Treatment.png")
    
    if (is.null(reference_batch)) {
      ggsave(filename, plot = raw_plots[[i]], path = path_images, width = 4, height = 11, dpi = 300)
    } else {
      ggsave(filename, plot = raw_plots[[i]], path = path_images, width = 6, height = 11, dpi = 300)
    }
  }
  
  StatsOutput <- tibble(
    parameter = sapply(param_list, `[[`, "value"),
    p_raw = raw_pvals,
    p_adj = adj_pvals
  )
  
  write.csv2(
    StatsOutput,
    file = file.path(path_images, "Exigo_StatsOutput.csv"),
    row.names = FALSE
  )
  
  return(StatsOutput)
}
run_all_LP_onecall(
  inputdata = exp_and_BaseLine,
  param_list = param_list,
  batch = c("1", "2", "BH15"),
  sex = "both",
  path_images = "02_GeneratedData/Exigo/FK49_and_BH15",
  p_adjust_method = "fdr",
  reference_batch = "BH15"
)
run_all_LP_onecall(
  inputdata = d1,
  param_list = param_list,
  batch ="ALL",
  sex = "both",
  path_images = "02_GeneratedData/Exigo",
  p_adjust_method = "fdr",
  reference_batch = NULL
)