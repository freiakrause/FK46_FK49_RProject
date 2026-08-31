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
             TRUE ~ "unknown")
    ) %>%
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
      .groups = "drop"
    )
  
  ### Normality test ---------------------------------------------------------------------
  shapiro_test_ctrl <- shapiro.test(d[[value]][d$Treatment == "ctrl"])
  shapiro_test_tam <- shapiro.test(d[[value]][d$Treatment == "TAM"])
  
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
  
  ### Correlation -------------------------------------------------------------------------
  cor_all <- cor.test(d[[value]], d$wks_diet, method = "spearman", na.action = na.omit)
  print(paste0("Cor test all: ", cor_all))
  
  cor_normal <- cor.test(d[[value]], d$wks_diet, method = "spearman", subset = event_status == "normal", na.action = na.omit)
  print(paste0("Cor test normal: ", cor_normal))
  
  cor_tumor <- cor.test(d[[value]], d$wks_diet, method = "spearman", subset = event_status == "tumor", na.action = na.omit)
  print(paste0("Cor test tumor: ", cor_tumor))
  
  cor_ascites <- cor.test(d[[value]], d$wks_diet, method = "spearman", subset = event_status == "ascites", na.action = na.omit)
  print(paste0("Cor test ascites: ", cor_ascites))
  
  ### Plot 1 ------------------------------------------------------------------------------
  p1 <- ggplot(stats, aes(x = Treatment, y = Mean, fill = Treatment)) +
    geom_bar(stat = "identity", color = "black", alpha = 0.5, width = 0.75, position = "dodge") +
    geom_point(data = d, fill = "lightgrey", color = "black",
               aes(y = .data[[value]], shape = event_status),
               position =  position_jitterdodge(0.1, dodge.width = 0.75), alpha = 0.8, size = 5.3, stroke = 1.8)+
    geom_errorbar(aes(ymin = Mean - SD, ymax = Mean + SD),
                  position = position_dodge(width = 0.75), width = 0.2) +
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
      panel.grid = element_blank()
    ) +
    guides(
      shape = guide_legend(title = "Status", order = 2 ,nrow = 3, byrow = TRUE),
      fill = guide_legend(title = "Treatment", order = 3,override.aes = list(shape = 21), nrow =2, byrow = TRUE),
      color = "none"
    )
  print(p1)
  
  ### Plot 2 ------------------------------------------------------------------------------
  p2 <- ggplot(d, aes(x = wks_diet, y = .data[[value]])) +
    geom_point(aes(shape = event_status, fill = Treatment), position = position_dodge(width = 0.15), size = 3.5, stroke = 1, color = "black") +
    scale_fill_manual(values = colors) +
    scale_shape_manual(values = c(21, 22, 24, 25, 26)) +
    scale_y_continuous(name = y_title) +
    scale_x_continuous(name = "Time on CD-HFD [wks]", limits = c(0, 41), breaks = seq(0, 41, 4), minor_breaks = seq(0, 40, by = 1)) +
    theme_minimal() +
    theme(
      legend.position = "bottom",
      axis.line = element_line(color = "black", linewidth = 0.5),
      axis.ticks = element_line(color = "black", linewidth = 0.5),
      axis.title = element_text(size = 12, face = "bold"),
      axis.text = element_text(size = 10, face = "bold"),
      plot.title = element_blank(),
      legend.title = element_text(size = 10),
      legend.text = element_text(size = 8),
      panel.grid = element_blank()
    ) +
    guides(
      shape = guide_legend(title = "Status",order = 3, nrow = 1, byrow = TRUE),
      fill = guide_legend(title = "Treatment",order = 2, nrow = 1, byrow = TRUE, override.aes = list(shape = 21))
    )
  print(p2)
 
  ### Save Plots -------------------------------------------------------------------------
  value_clean <- gsub("[^[:alnum:]_]", "_", value)
  filename1 <- paste0("FK46_", value_clean, "_Treat", batch, "_", sex, ".png")
  filename2 <- paste0("FK46_", value_clean, "_Time_", batch, "_", sex, ".png")
  filename3 <- paste0("FK46_", value_clean, "_Treat_", batch, "_", sex, ".svg")
  filename4 <- paste0("FK46_", value_clean, "_Time_", batch, "_", sex, ".svg")
  
  ggsave(filename = filename1, plot = p1, path = path_images, width = 4, height = 11, dpi = 300)
  ggsave(filename = filename2, plot = p2, path = file.path(path_images, "background"), width = 9, height = 6, dpi = 300)
  #ggsave(filename = filename3, plot = p1, path = file.path(path_images, "background"), width = 4.5, height = 9)
  #ggsave(filename = filename4, plot = p2, path = file.path(path_images, "background"), width = 9, height = 6)
}


## Call the function with desired arguments --------------------------------------------------------------
path_for_saving_images<-"02_GeneratedData/Organ"
do_organ_weight(data, value = "Liver_rel", batch = "ALL", sex = "both", y_title = "Liver/BW [%]",path_for_saving_images)
do_organ_weight(data, value = "Liver_rel", batch = "ALL", sex = "female", y_title = "Liver/BW [%]",path_for_saving_images)
do_organ_weight(data, value = "Liver_rel", batch = "ALL", sex = "male", y_title = "Liver/BW [%]",path_for_saving_images)

do_organ_weight(data, value = "Spleen_rel",batch = "ALL", sex = "both",y_title= "Spleen/BW [%]",path_for_saving_images)
do_organ_weight(data, value = "Spleen_rel",batch = "ALL", sex = "female",y_title= "Spleen/BW [%]",path_for_saving_images)
do_organ_weight(data, value = "Spleen_rel",batch = "ALL", sex = "male",y_title= "Spleen/BW [%]",path_for_saving_images)

do_organ_weight(data, value = "Spleen",batch = "ALL", sex = "both",y_title= "Spleen [mg]",path_for_saving_images)
do_organ_weight(data, value = "Spleen",batch = "ALL", sex = "female",y_title= "Spleen [mg]",path_for_saving_images)
do_organ_weight(data, value = "Spleen",batch = "ALL", sex = "male",y_title= "Spleen [mg]",path_for_saving_images)

do_organ_weight(data, value = "Liver",batch = "ALL", sex = "both",y_title= "Liver [g]",path_for_saving_images)
do_organ_weight(data, value = "Liver",batch = "ALL", sex = "female",y_title= "Liver [g]",path_for_saving_images)
do_organ_weight(data, value = "Liver",batch = "ALL", sex = "male",y_title= "Liver [g]",path_for_saving_images)

dd<-data%>%group_by(Animal,Treatment,BATCH,Tumor.no.yes,Ascites.no.yes,Granuloma, Sex,EP_weight,rel_EP_weight)%>%summarize(wks_diet=max(wks_diet))

do_organ_weight(dd, value = "EP_weight", batch = "ALL", sex = "both", y_title = "Weight at Endpoint [g]",path_for_saving_images)
do_organ_weight(dd, value = "EP_weight", batch = "ALL", sex = "female", y_title = "Weight at Endpoint [g]",path_for_saving_images)
do_organ_weight(dd, value = "EP_weight", batch = "ALL", sex = "male", y_title = "Weight at Endpoint [g]",path_for_saving_images)

do_organ_weight(dd, value = "rel_EP_weight", batch = "ALL", sex = "both", y_title = "rel.Weight at Endpoint [%]",path_for_saving_images)
do_organ_weight(dd, value = "rel_EP_weight", batch = "ALL", sex = "female", y_title = "rel.Weight at Endpoint [%]",path_for_saving_images)
do_organ_weight(dd, value = "rel_EP_weight", batch = "ALL", sex = "male", y_title = "rel.Weight at Endpoint [%]",path_for_saving_images)

# Dotplot Exigo Liver Panel --------------------------------------------------------------
do_LP <- function(inputdata, value, batch = "ALL", sex = "both", y_title, path_images,
                  normal_range = NULL, lowlimit = NULL, hilimit = NULL) {
  # --- Data manipulation ---
  d <- inputdata %>%
    filter(complete.cases(Liver)) %>%
    mutate(across(c(Tumor.no.yes, Ascites.no.yes), as.factor),
           event_status = case_when(
             Tumor.no.yes == "0" & Ascites.no.yes == "0" ~ "normal",
             Tumor.no.yes == "1" & Ascites.no.yes == "0" ~ "tumor",
             Tumor.no.yes == "0" & Ascites.no.yes == "1" ~ "ascites",
             Tumor.no.yes == "1" & Ascites.no.yes == "1" ~ "both",
             TRUE ~ "unknown"),
           event_status = factor(event_status, levels = c("normal", "ascites", "tumor", "both", "unknown")))
  
  if (batch != "ALL") d <- d %>% filter(BATCH == batch)
  if (sex != "both") d <- d %>% filter(Sex == sex)
  
  censored_col <- paste0(value, "_censored")
  direction_col <- paste0(value, "_direction")
  
  d <- d %>%
    mutate(censored = .data[[censored_col]],
           direction = .data[[direction_col]],
           Sex = as.factor(trimws(Sex)),
           Treatment = as.factor(trimws(Treatment)),
           censored = as.character(censored),
           censor_status_combined = case_when(
             censored == "TRUE" & direction == "<" ~ "Below LOD",
             censored == "TRUE" & direction == ">" ~ "Above ULOQ",
             TRUE ~ "Detected"))
  
  # --- Prepare data for censored analysis ---
  cens_data <- d %>%
    mutate(value_numeric = as.numeric(.data[[value]]),
           cens_logical = censored == "TRUE",
           direction_factor = case_when(
             direction == "<" ~ "left",
             direction == ">" ~ "right",
             TRUE ~ "none"))
  
  # --- Tests ---
  d[[value]] <- jitter(as.numeric(d[[value]]), amount = 0.001)
  d$wks_diet <- jitter(d$wks_diet, amount = 0.001)
  
  # Shapiro-Wilk tests
  print(shapiro.test(d[[value]][d$Treatment == "ctrl"]))
  print(shapiro.test(d[[value]][d$Treatment == "TAM"]))
  
  # Helper function to check sufficient uncensored distinct values per group
  check_data_sufficient <- function(data, value_col, cens_col, group_col) {
    data %>%filter(!.data[[cens_col]]) %>%  # keep only uncensored
        group_by(.data[[group_col]]) %>%
        summarise(n_distinct_values = n_distinct(.data[[value_col]]), n_non_missing = sum(!is.na(.data[[value_col]]))) %>%
        pull(n_distinct_values) -> distinct_counts
    
    all(distinct_counts >= 2)}
  
  # Censored group comparison using cen2means, with safety check
  has_censored <- any(cens_data$cens_logical, na.rm = TRUE)
  
  if (has_censored) {
    message("Using censored-data method (cen2means)")
    if (check_data_sufficient(cens_data, "value_numeric", "cens_logical", "Treatment")) {
      #Test with censored data for 2 groups
      cen_result <- with(cens_data, cen2means(value_numeric, cens_logical, group = Treatment))
      p_value <- cen_result$pval
      message("generated P value")
      print(cen_result$pval)
      cor_p<- NULL
      correlation<-NULL
      cor_labels<- NULL
      cor_title<-NULL}
    else {message("No sufficient uncensored distinct values for censored-data test — skipping test")
      p_value <- NA
      cor_p<- NULL
      correlation<-NULL
      cor_labels<- NULL
      cor_title<-NULL}} 
  else {
    message("No censored values detected — using t.test")
    #Test for uncensored data for 2 groups
    t_test_result <- t.test(value_numeric ~ Treatment, data = cens_data)
    p_value <- t_test_result$p.value
    # --- Correlation per sex ---
    cor_results <- d %>%group_by(Sex) %>%
      summarise(cor = cor.test(.data[[value]], wks_diet, method = "spearman")$estimate,
                p_val = cor.test(.data[[value]], wks_diet, method = "spearman")$p.value,.groups = "drop")
    print(cor_results)
    
    #  Create correlation label for plot title ---
    cor_labels <- paste0("Sex: ", cor_results$Sex, " — Spearman r = ", round(cor_results$cor, 2), ", p = ", signif(cor_results$p_val, 3))
    cor_title <- paste(cor_labels, collapse = "\n")}
  
  
  p_value_label <- ifelse(is.na(p_value),"No test only censored data",paste("p =", format.pval(p_value, digits = 3)))

  # Effect size
  cohen_d_result <- cohen.d(d[[value]] ~ d$Treatment)
 
  
  
  y_max <- max(d[[value]], na.rm = TRUE)
  y_min <- min(d[[value]], na.rm = TRUE)
  y_pos <- y_max + 0.15 * (y_max - y_min)
  # Plot 1: Treatment comparison
  p1 <- ggplot(d, aes(x = Treatment, y = .data[[value]])) +
    stat_summary(fun = mean, geom = "bar", aes(fill = Treatment),alpha = 0.5, width = 0.6, color = "black") +
    stat_summary(fun.data = mean_sdl, fun.args = list(mult = 1), geom = "errorbar", width = 0.2, color = "black") +
    scale_fill_manual(name = "Treatment", values = c("black", "darkred")) +
    guides(fill = guide_legend(order = 1, nrow = 1, byrow = TRUE)) +
    ggnewscale::new_scale_fill() +
    geom_point(aes(shape = event_status, fill = Sex, color = censor_status_combined),
               position = position_jitter(width = 0.15), size = 5.3, stroke = 1.8) +
    scale_fill_manual(name = "Sex", values = c("male" = "lightblue", "female" = "pink")) +
    scale_shape_manual(name = "Status", values = c(21, 22, 24, 25, 26)) +
    scale_color_manual(name = "Censoring", values = c("Below LOD" = "darkred",
                                                      "Above ULOQ" = "darkblue",
                                                      "Detected" = "black")) +
    scale_y_continuous(name = y_title, expand = expansion(mult = c(0.05, 0.15)))+
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
        fill = guide_legend(title = "Sex", order = 3, nrow = 1, byrow = TRUE,
                          override.aes = list(shape = 21, color = NA, size = 5, fill = c("pink","lightblue"))))+
    annotate("text", x = 1.5, y = y_pos,label = p_value_label, size = 6, color = "black", fontface = "italic")
  
  if (!is.null(normal_range)) {p1 <- p1 + geom_hline(yintercept = normal_range, linetype = "dashed",color = "darkgrey", linewidth = 0.7)}
  if (!is.null(lowlimit)) { p1 <- p1 + geom_hline(yintercept = lowlimit[1], linetype = "dashed",color = "darkred", linewidth = 0.7)}
  if (!is.null(hilimit)) {p1 <- p1 + geom_hline(yintercept = hilimit[1], linetype = "dashed", color = "darkblue", linewidth = 0.7)}
  
  # Plot 2: Time correlation
  p2 <- ggplot(d, aes(x = wks_diet, y = .data[[value]])) +
    geom_smooth(method = "lm", se = TRUE, aes(color = Sex,fill= Sex), linewidth = 1,alpha=0.3) +
    scale_color_manual(name = "Sex2", values = c("male" = "lightblue", "female" = "pink"))+
    guides(color = guide_legend(order = 1, nrow = 1, byrow = TRUE)) +
    ggnewscale::new_scale_color() +
    geom_point(aes(shape = event_status, fill = Sex, color = censor_status_combined),
               position = position_dodge(width = 0.15), size = 3.5, stroke = 1) +
    scale_shape_manual(name = "Status", values = c(21, 22, 24, 25, 26)) +
    scale_color_manual(name = "Censoring", values = c("Below LOD" = "darkred","Above ULOQ" = "darkblue","Detected" = "black")) +
    scale_fill_manual(name = "Sex", values = c("male" = "lightblue", "female" = "pink")) +
    scale_y_continuous(name = y_title) +
    scale_x_continuous(name = "Time on CD-HFD [wks]",limits = c(0, 41), breaks = seq(0, 41, 4), minor_breaks = seq(0, 40, 1)) +
    labs(title = cor_title, x = "Weeks on Diet", y = value)+
    theme_minimal() +
    theme(
      legend.position = "bottom",
      axis.line = element_line(color = "black", linewidth = 0.5),
      axis.ticks = element_line(color = "black", linewidth = 0.5),
      axis.title = element_text(size = 12, face = "bold"),
      axis.text = element_text(size = 10, face = "bold"),
      legend.title = element_text(size = 11, face = "bold"),
      legend.text = element_text(size = 8),
      panel.grid = element_blank())+
    guides(
    shape = guide_legend(title = "Status", order = 2, nrow = 1, byrow = TRUE),
    color = guide_legend(title = "Censoring", order = 4, nrow = 1, byrow = TRUE),
    fill = guide_legend(title = "Sex", order = 3, nrow = 1, byrow = TRUE,
                        override.aes = list(shape = 21, color = NA, size = 5, fill = c("pink","lightblue"))))
  
  if (!is.null(normal_range)) { p2 <- p2 + geom_hline(yintercept = normal_range, linetype = "dashed",color = "darkgrey", linewidth = 0.7)}
  if (!is.null(lowlimit)) {p2 <- p2 + geom_hline(yintercept = lowlimit[1], linetype = "dashed",color = "darkred", linewidth = 0.7)}
  if (!is.null(hilimit)) {p2 <- p2 + geom_hline(yintercept = hilimit[1], linetype = "dashed", color = "darkblue", linewidth = 0.7)}
  
  # --- Save plots ---
  fname_val <- gsub("[^[:alnum:]_]", "_", value)
  filename1 <- paste0("FK46_", fname_val, "_Treatment_Batch", batch, "_", sex, ".png")
  filename2 <- paste0("/Background/FK46_", fname_val, "_Time_Batch", batch, "_", sex, ".png")
  # filename3 <- sub(".png$", ".svg", filename1)
  # filename4 <- sub(".png$", ".svg", filename2)
  
  ggsave(filename = filename1, plot = p1, path = path_images, width = 4, height = 11, dpi = 300)
  ggsave(filename = filename2, plot = p2, path = path_images, width = 9, height = 6, dpi = 300)
  # ggsave(filename = filename3, plot = p1, path = path_images, width = 4, height = 11)
  # ggsave(filename = filename4, plot = p2, path = path_images, width = 9, height = 6)
}


path_for_saving_images<-"02_GeneratedData"

do_LP(data, value = "ALB", batch = "ALL", sex = "both", y_title = "Alb [g/L]",path_for_saving_images,normal_range=c(20, 48),lowlimit=2)
do_LP(data, value = "TP", batch = "ALL", sex = "both", y_title = "TP [g/L]",path_for_saving_images,normal_range = c(36,66))
do_LP(data, value = "GLOB", batch = "ALL", sex = "both", y_title = "GLOB [g/L]",path_for_saving_images,normal_range = NULL)
do_LP(data, value = "A.G", batch = "ALL", sex = "both", y_title = "A/G",path_for_saving_images,normal_range = NULL)
do_LP(data, value = "TB", batch = "ALL", sex = "both", y_title = "TB [µmol/L]",path_for_saving_images,normal_range = c(0,15), lowlimit=0.1)
do_LP(data, value = "GGT", batch = "ALL", sex = "both", y_title = "GGT [U/L]",path_for_saving_images,normal_range = NULL,lowlimit= 2)
do_LP(data, value = "AST", batch = "ALL", sex = "both", y_title = "AST [U/L]",path_for_saving_images,normal_range = c(59,247),hilimit = 650,lowlimit=5)
do_LP(data, value = "ALT", batch = "ALL", sex = "both", y_title = "ALT [U/L]",path_for_saving_images,normal_range = c(28,132))
do_LP(data, value = "ALP", batch = "ALL", sex = "both", y_title = "ALP [U/L]",path_for_saving_images,normal_range = c(62,209),lowlimit = 5)
do_LP(data, value = "TBA", batch = "ALL", sex = "both", y_title = "TBA [µmol/L]",path_for_saving_images,normal_range = NULL,lowlimit=1)
do_LP(data, value = "TC", batch = "ALL", sex = "both", y_title = "TC [mmol/L]",path_for_saving_images,normal_range = c(0.93,4.04))


