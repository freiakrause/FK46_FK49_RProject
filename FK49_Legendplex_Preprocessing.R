#From Birgit
# protocol ---------------------------------------------------------------------

#1) After measuring the experiment, log into the Biorender Legendplex software using my Uni log in data (birgit.halwachs@plus.ac.at)
#2) Create a new experiment and enter the legendplex kit of interest and select "Save"
#3) Chose the default settings and upload your 16 standard .fcs files (labelled C0_1-C7_2)
#4) Select FL10-A for PE and FL3-A for the APC reporter channel - select "Upload"
#5) Select your sample .fcs files and "Set Replicates" to 2 with the pattern AABBCC - "Upload"
#7) The software will automatially gate your files - while this in general works very well, there can happen mistakes in individual files
#8) "Review" Gating - check all peaks are selected and optionally already remove samples with too little events - this is your first round of QC
#9) Check your results - for each analyte check the standard curve 
#10) Select "Report" to download your results - select everything but "PDF Report" in order to receive a .xlsx file - save in 02_generated of your experiment folder
#11) Open the .xlsx file locally and go to the "predicted_conc" tab
#12) Select File>Save Copy> and select CSV-UTF 8
#13) Provide your experiment metadata either in 00_meta or in a location of your choice
#14) You are ready to go :)

# libraries --------------------------------------------------------------------
library(tidyverse)
library(pheatmap) 
 library(svDialogs)
# library(rstudioapi)
library(ggfortify)
 library(factoextra)
 library(broom)
 library(purrr)

# define functions and annotation colors ---------------------------------------
annotation_colors <- list( #adjust this to your individual metadata
  Induction = c(
    "TAM" = "darkred",
    "EtOH"  = "#666666"
  ),
  Injection = c(
    "MC38" = "#669933",
    "PBS" = "#CCC",
    "na" = "darkgrey"
  ),
  Sex = c(
    "male"   = "#3366FF",
    "female" = "#FF99FF"
  ),
  batch_lysis = c(
    "1" = "lightblue",
    "2" = "darkblue"
  )
)



# load data -------------------------------------------------------------
#set working directory. If there is a new (unknown to me) experiment, path needs to be included here
ExpId = "FK49"
if (ExpId=="FK49") {
  setwd("C:/Users/b1084855/OneDrive - Universität Salzburg/Freigegebene Dokumente - AG_Tumorimmunologie/Data/Freia Krause/01_Experiments/FK49_CD-HFD_13wks/FK49_Legendplex")
  data_raw <-read.csv("02_generated/BH21.3 FK49 Legendplex_report_2025-11-19.csv", sep = ";") 
  }  else if (ExpId == "FK46"){
  setwd("C:/Users/b1084855/OneDrive - Universität Salzburg/Freigegebene Dokumente - AG_Tumorimmunologie/Data/Freia Krause/01_Experiments/FK46_iALmice_high Fat diet 52 weeks 7d after injection/FK46_Legendplex")
  data_raw <-read.csv("02_generated/BH21.5_FK46 Legendplex_report_2025-11-19.csv", sep = ";") 
    
    } else{
  print("Let me set a folder Path and define ExpId and Folderpath")}

cytokines <-c("IL23","IL1a" , "IFNy", "TNFa", "MCP1","IL12p70", "IL1ß","IL10", "IL6", "IL27", "IL17A","IFNß","GMCSF")

# data wrangling (clean up the legendplex output) ------------------------------
data_raw$well <- sub("_.*", "", data_raw$well) 
data_raw <- data_raw %>% rename(IDs = well)

data_analytes <- data_raw %>% 
  filter(sample_type != "Standard") %>%       # remove the standard values
  mutate(  Animal = gsub("#", "", IDs),
           sample_label = paste(Animal, sample, sep = "_"),
           replicate = factor(replicate)) %>%
  select(-experiment, -IDs, -sample) %>%
  rename_with(~ gsub("[A-Z][0-9]+$", "", gsub("\\.", "", .x))) %>%
  rename_with(~ gsub("α", "a", .x)) %>%
  rename_with(~ gsub("β", "ß", .x)) %>%
  rename_with(~ gsub("γ", "y", .x)) %>%
  pivot_longer(all_of(cytokines),  names_to = "parameter", values_to = "value") %>% # since some values in the assay were belo detection limit iwant to have censoring informtaion
  mutate( censored = str_detect(value, "[<>]"),                                     
          direction = case_when(
            str_detect(value, "^>") ~ ">",
            str_detect(value, "^<") ~ "<",
            TRUE ~ NA_character_),
          numeric_value = as.numeric(str_remove_all(value, "[^0-9\\.]")) )%>%      # numeric_value_diluted = numeric_value * DFactor) 
  pivot_wider(names_from  = parameter,
              values_from = c(value, censored, direction, numeric_value), #, numeric_value_diluted # cytokine_value is value or value with direction e.g <30.8
              names_glue  = "{parameter}_{.value}") %>%
      rename_with(~ str_replace(., "_numeric_value$", ""))                      #the colmn with only the cytokine name is the numeric value (wo censoring info)

 
data_long <- data_analytes %>%
  pivot_longer(
    cols = all_of(cytokines),    # your vector of cytokine column names
    names_to = "Analyte",
    values_to = "Value"
  )
# Perform a QC to check whether the technical replicates are similar -----------

# check visually if you can detect any outliers
ggplot(data_long, aes(x = replicate, y = sample_label, fill = log10(as.numeric(Value)))) +
  geom_tile(color = "white") +
  scale_fill_gradient(
    low = "white",
    high = "steelblue",
    na.value = "red",      # NAs show up clearly
    name = "log10(Value)"
  ) +
  facet_grid(~ Analyte, scales = "free") +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank()
  ) +
  labs(
    title = "Absolute Values for Selected Samples (geom_tile)",
    x = "Replicate",
    y = "Sample"
  )


# This function calculates a value called the CV (Coefficient of Variantion), which essentially just calculates
# SD/Mean --> the larger the value, the worse the quality of the replicates 
qc_replicates <- function(df, value_col = "Value", group_cols = c( "Animal", "sample_label", "Analyte"),
                          cv_threshold = 0.35, #can be adjusted but even for bad assays should be <0.5
                          range_threshold = 50,
                          min_reps = 1) {
  
  df %>%
    group_by(across(all_of(group_cols))) %>%
    summarise(
      n_reps = sum(!is.na(.data[[value_col]])),
      values = list(.data[[value_col]]),
      mean_value = mean(.data[[value_col]], na.rm = TRUE),
      sd_value = sd(.data[[value_col]], na.rm = TRUE),
      cv = ifelse(mean_value == 0 | is.na(mean_value),
                  NA_real_,
                  sd_value / mean_value),
      range = if (all(is.na(.data[[value_col]]))) {
        NA_real_
      } else {
        max(.data[[value_col]], na.rm = TRUE) -
          min(.data[[value_col]], na.rm = TRUE)
      },
      replicate_values = paste(round(.data[[value_col]], 3), collapse = ", "),
      .groups = "drop"
    ) %>%
    mutate(
      all_na = n_reps == 0,
      fail_cv = !is.na(cv_threshold) & !is.na(cv) & cv > cv_threshold,
      fail_range = !is.na(range_threshold) & !is.na(range) & range > range_threshold,
      fail_min_reps = n_reps < min_reps,
      flag_qc_fail = ifelse(all_na, FALSE, fail_cv | fail_range | fail_min_reps)
    )
}

qc_summary <- qc_replicates(data_long)

ggplot(qc_summary, aes(x = cv)) + # plot CV distribution per analyte
  geom_histogram(bins = 40) +
  facet_wrap(~ Analyte, scales = "free_y") +
  theme_minimal()

qc_failures <- qc_summary %>% filter(flag_qc_fail == T)

failed_raw <- data_long %>%
  semi_join(qc_failures, by = c( "Animal", "sample_label", "Analyte"))

qc_failures <- qc_failures %>%
  mutate(
    fail_reason = case_when(
      all_na ~ "All NA",
      fail_min_reps ~ "Too few reps",
      fail_cv ~ "High CV",
      fail_range ~ "High range",
      TRUE ~ "Unknown"
    )
  )

# Plot QC failures to check whether there are simple samples which are problematic and which analytes are probably hard to interpret
ggplot(qc_failures, 
       aes(x = Analyte,
           y = paste(Animal, sample_label, sep = "_"),
           fill = fail_reason)) +
  geom_tile(color = "white") +
  scale_fill_manual(
    values = c(
      "High CV" = "#d73027",
      "High range" = "#fc8d59",
      "Too few reps" = "#fee090",
      "All NA" = "#91bfdb"
    )
  ) +
  theme_bw(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  labs(
    title = "QC Failures Heatmap",
    x = "Analyte",
    y = "Sample",
    fill = "Failure Reason"
  )

# Optionally remove obviously erroneous replicates -----------------------------
# Based on these results I decided to kick out #152_2 (BH21.4 specific)
data_analytes <- data_analytes %>% filter(!sample_label =="202_Sample28",!Animal == "146")
# data analysis ----------------------------------------------------------------
data_sum <- data_analytes %>%
  group_by(Animal) %>%
  summarize(
    across(all_of(cytokines), ~ mean(.x, na.rm = TRUE), .names = "{.col}"),
    across(ends_with("_censored"), ~ any(.x, na.rm = TRUE), .names = "{.col}"),
    #across(ends_with("_direction"), ~ if (length(unique(.x)) == 1) unique(.x) else "xxxx", .names = "{.col}"),
    across(ends_with("_direction"), ~ paste(unique(.x),collapse = ","), .names = "{.col}"),
    
    .groups = "drop"
              )
write.csv(data_sum,file= paste0(ExpId,"_Legenplex_clean.csv"))
# OPTIONAL data subsetting -----------------------------------------------------
#meta_liver <- meta %>% filter(Tissue == "liver")
#ids_liver <- meta_liver %>% pull("IDs")
#data_liver <- data_sum %>% filter(IDs %in% ids_liver)

# plot correlation heatmaps ----------------------------------------------------
#this function will take your data to calculate a pearson correlation between all samples
make_cor_matrix <- function(data, method = "spearman") {
  d<-data%>%select(Animal,all_of(cytokines))
  d<-t(d)
  data_matrix <- as.matrix(d[, -1])  # rownames are preserved
  colnames(d) <- as.character(d[1, ])
  data_matrix <- as.matrix(d[-1,])  # colnames are preserved
  
  storage.mode(data_matrix) <- "numeric"
  
  cor_matrix <- cor(data_matrix,
                    use = "pairwise.complete.obs",
                    method = method)
  
  return(list(cor_matrix = cor_matrix, data_matrix = data_matrix, data_cor = d))
} 
#this function takes the results of make_cor_matrix() res$cor_matrix to plot a correlation heatmap
plot_cor_heatmap <- function(cor_matrix, meta, annotation_colors = NULL) {
  
  # ---- Prepare annotation ----
  annotation <- meta %>%
    dplyr::select(Induction, Sex, IDs) %>%  #batch_lysis Injection,
    mutate(across(everything(), as.factor)) %>%
    as.data.frame()
  
  rownames(annotation) <- colnames(cor_matrix)
  annotation$IDs <- NULL
  
  # ---- Plot heatmap ----
  heatmap_plot <- pheatmap::pheatmap(cluster_cols = T,cluster_rows = T,
                                     cor_matrix,
                                     annotation_col = annotation,
                                     annotation_colors = annotation_colors
  )
  
  return(heatmap_plot)
} 
cor_data <- make_cor_matrix(data = data_sum, method = "spearman")
(p_corHeatmap_ALL <- plot_cor_heatmap(cor_matrix = cor_data$cor_matrix, 
                                      meta = meta, 
                                      annotation_colors = annotation_colors))

# plot PCAs --------------------------------------------------------------------
data_matrix_noNA <- cor_data$data_matrix
data_matrix_noNA[is.na(data_matrix_noNA)] <- 0 # replace NAs with 0

pca_matrix <- t(data_matrix_noNA)
sample_pca <- prcomp(pca_matrix)
(p_PCA <- autoplot(sample_pca, data = meta, x = 1, y = 2, size = 5, color = "Sex", shape = "Induction")+ theme_bw()+ 
    ggtitle("Principal Component Analysis")+ 
    scale_color_brewer(palette = "Set2")+ 
    geom_text(aes(label = IDs), vjust = -1, size = 5)+
    theme(text = element_text(size = 20)))
factoextra::fviz_contrib(sample_pca, choice = "var", 
                         axes = 1, top = 15, sort.val = c("desc"))

# merge data with metadata and normalize to original protein concentration ----------------------------------------------------- 
data_merged <- merge(data_sum, meta)
# normalize data to the sample concentration
data_norm <- data_merged %>% mutate(mean_norm = mean/sample_concentration)
data_norm$mean_norm_scaled <- scale(data_norm$mean_norm)[,1]
median(data_norm$mean_norm_scaled, na.rm = TRUE)

# plot cytokine levels in a heatmap for each Experiment ID separately ----------
#data_norm$ExperimentID <- factor(data_norm$ExperimentID, levels = c("BH06", "BH15", "BH20", "BH05"))

ggplot(data_norm, aes(x = IDs, y = Analyte, fill = mean_norm_scaled)) + #plot all experiments together
  geom_tile() +
  facet_wrap(MouseLine ~ factor(Induction, levels = c("EtOH", "TAM")),
             nrow = 4, scales = "free_x") +
  scale_fill_gradient2(low = "#4575B4", high = "#D73027", midpoint = 0, limits = range(data_norm$mean_norm_scaled, na.rm = TRUE)) +
  labs(title = "All Experiments") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

exp_list <- unique(data_norm$ExperimentID) #plot individual experiments

plot_list <- list()

for (exp in exp_list) { 
  df <- subset(data_norm, ExperimentID == exp)
  
  p <- ggplot(df, aes(x = IDs, y = Analyte, fill = mean_norm_scaled)) +
    geom_tile() +
    facet_grid(. ~ factor(Induction, levels = c("ETOH", "TAM")), space = "free", scales = "free") +
    scale_fill_gradient2(low="#4575B4", high="#D73027", midpoint = 0,
                         limits = range(data_norm$mean_norm_scaled, na.rm = TRUE)) +
    labs(title = paste("Experiment:", exp)) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
  
  # Store the plot in the list, named by experiment ID
  plot_list[[exp]] <- p
}

for (p in plot_list) print(p)

# Statistics --------------------------------------------------------------------
shapiro_EtOH <- shapiro.test(data_norm$mean_norm[data_norm$Induction=="EtOH"])$p.value
shapiro_TAM  <- shapiro.test(data_norm$mean_norm[data_norm$Induction=="TAM"])$p.value

if(shapiro_EtOH > 0.05 & shapiro_TAM > 0.05){
  chosen_test <- "t.test"
} else {
  chosen_test <- "wilcox.test"
}

print(chosen_test)

stats <- data_norm %>%
  group_by( Analyte) %>%
  summarise(
    test_result = list(
      if(chosen_test == "t.test"){
        t.test(mean_norm ~ Induction, data = cur_data())
      } else {
        wilcox.test(mean_norm ~ Induction, data = cur_data())
      }
    ),
    mean_EtOH = mean(mean_norm[Induction=="EtOH"], na.rm = TRUE),
    mean_TAM  = mean(mean_norm[Induction=="TAM"], na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(tidy = map(test_result, broom::tidy)) %>%
  unnest(tidy) %>%
  # Keep only the columns you want
  select( Analyte, mean_EtOH, mean_TAM, statistic, p.value)

# 2. Adjust p-values **per experiment** (across analytes)
stats <- stats %>%
  mutate(p_adj = p.adjust(p.value, method = "BH")) %>%
  ungroup()

stats <- stats %>%
  mutate(
    FC = mean_TAM / mean_EtOH,              
    log2FC = log2(mean_TAM / mean_EtOH)   
  )

ggplot(stats, aes(x = FC, y  = -log10(p.value))) + 
  geom_hex(alpha = 0.8)

ggplot(stats, aes(x = -log10(p.value))) + 
  geom_histogram()

# merge stats with data_norm ---------------------------------------------------
data_res <- merge(stats, data_norm)

# plot results -----------------------------------------------------------------
(p_foldchange <- ggplot(data_res, aes(x = MouseLine, y = fct_reorder(Analyte, log2FC), fill = log2FC)) +
  geom_tile(color = "white") +
  geom_text(aes(label = round(p_adj, 3)), color = "black", size = 3) +
  scale_fill_gradient2(low = "#4575B4", mid = "white", high = "#D73027", midpoint = 0) +
  theme_minimal() +
  labs(
    x = "Experiment",
    y = "Analyte",
    fill = "log2 Fold Change",
    title = "Analyte log2 Fold Changes (TAM / ETOH)"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text.y = element_text(size = 8)
  )
)
# plot individual analytes of interest ------------------------------------------
analytes <-unique(data_res$Analyte)
for(a in analytes){
  analytes_of_interest <- a
  #

  df_subset <- data_res %>%
   filter(Analyte %in% analytes_of_interest)

  df_stats <- data_norm %>%
   filter(Analyte %in% analytes_of_interest) %>%
    group_by( Analyte, Induction) %>%
    summarise(
      mean_val = mean(mean_norm, na.rm = TRUE),
      sd_val   = sd(mean_norm, na.rm = TRUE),
     n        = sum(!is.na(mean_norm)),
      se_val   = sd_val / sqrt(n),
      .groups = "drop"
   )

  df_error <- df_subset %>%
    group_by( Analyte, Induction) %>%
   summarise(
     mean_val = mean(mean_norm, na.rm = TRUE),
     sd_val   = sd(mean_norm, na.rm = TRUE),
      n        = sum(!is.na(mean_norm)),
      se_val   = sd_val / sqrt(n),
     .groups = "drop"
  )

  (p_individual <- ggplot(df_subset, aes(x = Induction, y = mean_norm)) +
    geom_point(size = 3, position = position_dodge(width = 0.5)) + 
    geom_errorbar(
     data = df_error,
     aes(x = Induction, y = mean_val, ymin = mean_val - sd_val, ymax = mean_val + sd_val),
     width = 0.2,
     position = position_dodge(width = 0.5)
   ) +
    #facet_grid(Analyte ~ MouseLine, scales = "free_y") +  # rows = analytes, consistent y-axis across experiments
    theme_bw() +
    labs(
     y = "Normalized Mean ± SD",
     title = "Individual Analytes by Induction and Experiment"
   ) +
   theme(axis.text.x = element_text(angle = 45, hjust = 1))
  )
  plot(p_individual)
}
# save generated plots -----------------------------------------------------------
plot_list <- lst(p_corHeatmap_ALL, p_corHeatmap_LIVER, p_foldchange, p_individual, p_PCA)
plot_names <- names(plot_list)
output_path <- "04_BH21.4_plots/"


for (i in 1:length(plot_list)) {
  file_name <- paste0(output_path, plot_names[i],".pdf")
  pdf(file_name, width = 12, height = 10)
  print(plot_list[[i]])
  dev.off()
}

for (i in 1:length(plot_list)) {
  file_name <- paste0(output_path, plot_names[i],".png")
  png(file_name, units = "in", height = 10, width = 12, res = 300)
  print(plot_list[[i]])
  dev.off()
}

for (i in 1:length(plot_list)) {
  file_name <- paste0(output_path, plot_names[i],".svg")
  svg(file_name, height = 10, width = 12)
  print(plot_list[[i]])
  dev.off()
}

