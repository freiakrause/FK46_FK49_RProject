
# libraries --------------------------------------------------------------------
library(tidyverse)
library(pheatmap) 
library(svDialogs)
# library(rstudioapi)
library(ggfortify)
library(factoextra)
library(broom)
library(purrr)
library('corrr')
library(ggcorrplot)
library("FactoMineR")
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
# define functions and annotation colors ---------------------------------------
annotation_colors <- list( #adjust this to your individual metadata
  Treatment = c(
    "TAM" = "darkred",
    "ctrl"  = "#666666"
  ),
  Sex = c(
    "male"   = "#3366FF",
    "female" = "#FF99FF"
  ),
  Batch = c(
    "1" = "lightblue",
    "2" = "darkblue"
  )
)

ExpId="FK49"

if (ExpId=="FK49") {
  setwd("C:/Users/b1084855/OneDrive - Universität Salzburg/Freigegebene Dokumente - AG_Tumorimmunologie/Data/Freia Krause/01_Experiments/FK49_CD-HFD_13wks/FK49_Analysis")

}  else if (ExpId == "FK46"){
  setwd("C:/Users/b1084855/OneDrive - Universität Salzburg/Freigegebene Dokumente - AG_Tumorimmunologie/Data/Freia Krause/01_Experiments/FK46_iALmice_high Fat diet 52 weeks 7d after injection/Analysis")
 } else{
  print("Let me set a folder Path and define ExpId and Folderpath")}

cytokines <-c("IL23","IL1a" ,"IFNy" ,"TNFa" , "MCP1" , "IL12p70" ,"IL1ß", "IL10" ,"IL6" ,"IL27" ,"IL17A" , "IFNß" ,"GMCSF") 
cytokine_pattern <- paste(cytokines, collapse = "|")

load(paste0("01_RawData/",ExpId,"_Data_prepared.Rda"))
data<-data%>%select(Animal, DOW, KILL.DATE,Sex,Treatment,BATCH,-TV,matches(cytokine_pattern),Tumor.no.yes,Ascites.no.yes,wks_diet)%>%
  filter(DOW== KILL.DATE,!Animal =="202",!Animal =="164") #202 has only NA/censored #164 is givin PCA1 maybe bc it had very low conc in liver lysate

# Columns with cytokine name are normalized to inpout protein concentration by formula: Cytokine[pg/mL]/Input[mg/mL] = pg cytokine /mg input
# Colum censoring gives info if detected values was in detection range FALSE or below/above detection range TRUE
# Column direction tell if it was below or above, if direction containes 2 values, duplicates of same samples had different values for censoring

# plot correlation heatmaps ----------------------------------------------------
#this function will take your data to calculate a pearson correlation between all samples
annotation <- data %>% filter(!Animal =="202",!Animal =="164")%>%
  select(Treatment, Sex, Animal,BATCH) %>%  #batch_lysis Injection,
  mutate(across(everything(), as.factor))   %>%arrange("T")%>%
  as.data.frame()

  d<-data%>% 
    arrange("Treatment")%>%
    select(-DOW,-KILL.DATE,-Sex,-BATCH,-Treatment,-ends_with("_censored"),
           -ends_with("_direction"),all_of(cytokines),-Ascites.no.yes,-Tumor.no.yes)
  colnames(d)
  d_transposed<-t(d)
  data_matrix <- as.matrix(d_transposed[, -1])  # rownames are preserved
  colnames(d_transposed) <- as.character(d_transposed[1, ])
  data_matrix <- as.matrix(d_transposed[-1,])  # colnames are preserved
  
  storage.mode(data_matrix) <- "numeric"
  
  cor_matrix <- cor(data_matrix,
                    use = "pairwise.complete.obs",
                    method = "spearman")
  
  rownames(annotation) <- colnames(cor_matrix)
  annotation$Animal <- NULL
  
  # ---- Plot heatmap ----
  heatmap_plot <- pheatmap(cluster_cols =T,
                           cluster_rows = T,
                                     cor_matrix,
                                     annotation_col = annotation,
                                     annotation_colors = annotation_colors)




# plot PCAs --------------------------------------------------------------------
  colSums(is.na(d))
  numerical_d<-d[,2:14]
  scaled_d <-scale(numerical_d)
  data.pca <- princomp(scaled_d)
  summary(data.pca)
  data.pca$loadings[, 1:2]
  fviz_eig(data.pca, addlabels = TRUE)
  fviz_pca_var(data.pca, col.var = "black")
  fviz_cos2(data.pca, choice = "var", axes = 1:2)
  fviz_pca_var(data.pca, col.var = "cos2",
               gradient.cols = c("black", "orange", "green"),
               repel = TRUE)
  fviz_contrib(data.pca, choice = "var", 
                            axes = 1, top = 15, sort.val = c("desc"))
  autoplot(data.pca, data = data, x = 1, y = 2, size = 5, color = "Treatment", shape = "Sex")+
   theme_bw()+ 
         ggtitle("Principal Component Analysis")+ 
         scale_color_brewer(palette = "Set2")+ 
         geom_text(aes(label = Animal), vjust = -1, size = 3)+
         theme(text = element_text(size = 20))



t_data_matrix<-t(data_matrix)
  heatmap_plot <- pheatmap(data_matrix,
                           scale= "row",
                           cluster_cols =T,
                           cluster_rows = T,
                           annotation_col = annotation,
                           annotation_colors = annotation_colors)
  
  data_sum <- data %>%
    group_by(Treatment, Sex) %>%
    summarise(
      across(
        all_of(cytokines),
        list(
          mean = ~mean(.x, na.rm = TRUE),
          sd   = ~sd(.x, na.rm = TRUE)
        ),
        .names = "{.col}_{.fn}"
      ),
      .groups = "drop"
    )
  
  # Set rownames using group info and keep only mean columns
  data_sum_matrix <- data_sum %>%
    unite("Group", Treatment, Sex, sep = "_") %>%   # e.g., Ctrl_male
    column_to_rownames("Group") %>%         # make it rownames
    select(ends_with("_mean")) %>%          # only mean values
    as.matrix()%>%t()
  
  pheatmap( data_sum_matrix,
    cluster_rows = T,
    cluster_cols = TRUE,
    scale = "row",   # optional: z-score each cytokine
    color = colorRampPalette(c("blue", "white", "red"))(100) )

  data <- data %>%
    mutate(  Sex_num = case_when(
        Sex %in% c("male", "M") ~ 1,
        Sex %in% c("female", "F") ~ 0,
        TRUE ~ NA_real_),
      Treatment_num = case_when(
        Treatment == "TAM" ~ 1,
        Treatment == "ctrl" ~ 0,
        TRUE ~ NA_real_)  )  
  cor_all <- data %>%
    select(all_of(cytokines),"Sex_num", "Treatment_num","wks_diet") %>%
    cor(use = "pairwise.complete.obs", method = "spearman")
  pheatmap(
    cor_all,
    color=colorRampPalette(c("blue", "grey90", "grey95", "red"))(100),
    breaks = seq(-1, 1, length.out = 101),  # 100 intervals  clustering_distance_rows = "correlation",
    clustering_distance_cols = "correlation",
    clustering_distance_rows = "correlation",
    main = "Correlation Cytokines and animal data long-term CD-HFD",
    display_numbers = TRUE,
    number_format = "%.1f"
  )  
  
  if (ExpId=="FK49") {
  Legendplex_Inflammation <- list(
    list(value="IL6",    y_title="IL6 [pg/mg]",    lowlimit=0.883),
    list(value="IL1a",   y_title="IL1a [pg/mg]",   lowlimit=0.076),
    list(value="IL1ß",    y_title="IL1ß [pg/mg]",    lowlimit=0.704),
    list(value="TNFa",   y_title="TNFa [pg/mg]",   lowlimit=1.935),
    list(value="MCP1",   y_title="MCP1 [pg/mg]",   lowlimit=1.77),
    list(value="IFNy",   y_title="IFNy [pg/mg]",   lowlimit=0.223),
    list(value="IL17A",  y_title="IL17A [pg/mg]",  lowlimit=0.381),
    list(value="IL12p70",y_title="IL12p70 [pg/mg]",lowlimit=0.599),
    list(value="IL10",   y_title="IL10 [pg/mg]",   lowlimit=5.345),
    list(value="IL23",   y_title="IL23 [pg/mg]",   lowlimit=9.61),
    list(value="IL27",   y_title="IL27 [pg/mg]",   lowlimit=11.82),
    list(value="IFNß",   y_title="IFNß [pg/mg]",   lowlimit=4.76),
    list(value="GMCSF",  y_title="GMCSF [pg/mg]",  lowlimit=0.296)
  )
  }else if(ExpId =="FK46"){
    Legendplex_Inflammation <- list(
      list(value="IL10",    y_title="IL10 [pg/mg]",    lowlimit=0.618),
      list(value="IL23",    y_title="IL23 [pg/mg]",    lowlimit=0.270),
      list(value="IL1a",    y_title="IL1a [pg/mg]",    lowlimit=0.058),
      list(value="IFNy",    y_title="IFNy [pg/mg]",    lowlimit=0.137),
      list(value="TNFa",    y_title="TNFa [pg/mg]",    lowlimit=0.135),
      list(value="MCP1",    y_title="MCP1 [pg/mg]",    lowlimit=0.235),
      list(value="IL12p70", y_title="IL12p70 [pg/mg]", lowlimit=0.197),
      list(value="IL6",     y_title="IL6 [pg/mg]",     lowlimit=0.070),
      list(value="IL1ß",    y_title="IL1ß [pg/mg]",    lowlimit=0.084),
      list(value="IL27",    y_title="IL27 [pg/mg]",    lowlimit=0.203),
      list(value="IL17A",   y_title="IL17A [pg/mg]",   lowlimit=0.074),
      list(value="IFNß",    y_title="IFNß [pg/mg]",    lowlimit=1.348),
      list(value="GMCSF",   y_title="GMCSF [pg/mg]",   lowlimit=0.061)
    )
    } else {print("Don't know which legenplex parameters you want")}
  
# Statistics --------------------------------------------------------------------
  load(paste0("01_RawData/",ExpId,"_Data_prepared.Rda"))
  data<-data%>%select(Animal, DOW, KILL.DATE,Sex,Treatment,BATCH,-TV,matches(cytokine_pattern),Tumor.no.yes,Ascites.no.yes,wks_diet)%>%
    filter(DOW== KILL.DATE,!Animal =="202",!Animal =="164") #202 has only NA/censored #164 is givin PCA1 maybe bc it had very low conc in liver lysate
  
#   data<-data%>%rename(Treatment=T) 
  do_Exigo <- function(inputdata, value, batch = "2", sex = "both", y_title, path_images,
                       normal_range = NULL, lowlimit = NULL, hilimit = NULL,
                       p_value_override = NULL) {
    # --- Data manipulation ---
    d <- inputdata %>%
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
    
    x_ctrl <- d[[value]][d$Treatment == "ctrl"]
    x_tam  <- d[[value]][d$Treatment == "TAM"]
    
    # Control group
    if (length(x_ctrl) > 3) {
      s1 <- shapiro.test(x_ctrl)
    } else {
      s1 <- list(p.value = 0.0)   # or any value you choose
    }
    
    # TAM group
    if (length(x_tam) > 3) {
      s2 <- shapiro.test(x_tam)
    } else {
      s2 <- list(p.value = 0.0)
    }
    
    
    # Helper function to check sufficient uncensored distinct values per group
    check_data_sufficient <- function(data, value_col, cens_col, group_col) {
      data %>%filter(!.data[[cens_col]]) %>%  # keep only uncensored
        group_by(.data[[group_col]]) %>%
        summarise(n_distinct_values = n_distinct(.data[[value_col]]), n_non_missing = sum(!is.na(.data[[value_col]]))) %>%
        pull(n_distinct_values) -> distinct_counts
      
      all(distinct_counts >= 2)}
    
    # Censored group comparison using cen2means, with safety check
    has_censored <- any(cens_data$cens_logical, na.rm = TRUE)
    
    # --- Censoring & statistical test ---
    
    if (has_censored) {
      if (all(cens_data$cens_logical, na.rm = TRUE)) {
        message("All values censored — skipping statistical test, plot only")
        p_value <- NA
        cor_title <- "All values are censored — no correlation computed"
        cor_p <- correlation <- cor_labels <- NULL
      } else {
        message("Censored data detected — using cen2means if possible")
        if (check_data_sufficient(cens_data, "value_numeric", "cens_logical", "Treatment")) {
          cen_result <- with(cens_data, cen2means(value_numeric, cens_logical, group = Treatment))
          p_value <- cen_result$pval
          message("Generated P value: ", p_value)
        } else {
          message("Insufficient uncensored values — skipping censored-data test")
          p_value <- NA
        }
        cor_title <- cor_p <- correlation <- cor_labels <- NULL
      }
    } else {
      message("No censored values — using t-test or Wilcoxon based on normality")
      test_result <- if (s1$p.value < 0.05 | s2$p.value < 0.05) {
        message("Data not normal — Wilcoxon test")
        wilcox.test(value_numeric ~ Treatment, data = cens_data)
      } else {
        message("Data normal — t-test")
        t.test(value_numeric ~ Treatment, data = cens_data)
      }
      p_value <- test_result$p.value
      cor_title <- cor_p <- correlation <- cor_labels <- NULL
    }
    if (!is.null(p_value_override)) {
      p_value <- p_value_override
    }
    # p_value <- test_result$p.value
    print(p_value)
    
    p_value_label <- ifelse(is.na(p_value),"No test only censored data",paste("p =", format.pval(p_value, digits = 3)))
    
    # Effect size
    cohen_d_result <- cohen.d(d[[value]] ~ d$Treatment)
    
    
    
    y_max <- max(d[[value]], na.rm = TRUE)
    y_min <- min(d[[value]], na.rm = TRUE)
    y_pos <- y_max + 0.15 * (y_max - y_min)
    # Plot 1: Treatment comparison
    p1 <- ggplot(d, aes(x = Treatment, y = .data[[value]])) +
      stat_summary(fun = mean, geom = "bar", aes(fill = Treatment),alpha = 0.5, width = 0.75, color = "black") +
      stat_summary(fun.data = mean_sdl, fun.args = list(mult = 1), geom = "errorbar", width = 0.2, color = "black") +
      scale_fill_manual(name = "Treatment", values = c("#4D4D4DBF", "#8B0000BF")) +
      guides(fill = guide_legend(order = 1, nrow = 1, byrow = TRUE)) +
      ggnewscale::new_scale_fill() +
      geom_point(aes(shape = event_status,  color = censor_status_combined), fill = "lightgrey",alpha=0.8,
                 position = position_jitter(width = 0.15), size = 5.3, stroke = 1.8) +
      scale_shape_manual(name = "Status", values = c(21, 22, 24, 25, 26)) +
      scale_color_manual(name = "Censoring", values = c("Below LOD" = "red",
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
             color = guide_legend(title = "Censoring", order = 4, nrow = 1, byrow = TRUE))#+
    #fill = guide_legend(title = "Sex", order = 3, nrow = 1, byrow = TRUE,
    #override.aes = list(shape = 21, color = NA, size = 5, fill = c("pink","lightblue"))
    # annotate("text", x = 1.5, y = y_pos,label = p_value_label, size = 6, color = "black", fontface = "italic")
    
    if (!is.null(normal_range)) {p1 <- p1 + geom_hline(yintercept = normal_range, linetype = "dashed",color = "darkgrey", linewidth = 0.7)}
    if (!is.null(lowlimit)) { p1 <- p1 + geom_hline(yintercept = lowlimit[1], linetype = "dashed",color = "red", linewidth = 0.7)}
    if (!is.null(hilimit)) {p1 <- p1 + geom_hline(yintercept = hilimit[1], linetype = "dashed", color = "darkblue", linewidth = 0.7)}
    
    # --- Save plots ---
    fname_val <- gsub("[^[:alnum:]_]", "_", value)
    filename1 <- paste0(ExpId,"_", fname_val, "_Batch", batch, "_", sex, ".png")
    # filename3 <- sub(".png$", ".svg", filename1)
    
    ggsave(filename = filename1, plot = p1, path = path_images, width = 4, height = 11, dpi = 300)
    # ggsave(filename = filename3, plot = p1, path = path_images, width = 4, height = 11)
    return(list(
      p_value_raw = p_value,
      plot_raw = p1
    ))
  }
  
run_all_LP_onecall <- function(inputdata, param_list, batch="2", sex="both",
                                 path_images, p_adjust_method="fdr") {

    n <- length(param_list)
    
    raw_pvals <- numeric(n)
    raw_plots <- vector("list", n)
    
    
    # PASS 1 — run do_LP once per parameter
    # (collect raw p-values + bare plots)
    
    message("STEP 1: running statistics & collecting raw plots")
    
    for (i in seq_along(param_list)) {
      p <- param_list[[i]]
      
      if (!(p$value %in% names(inputdata))) {
        message("Parameter not in data: ", p$value)
        raw_pvals[i]  <- NA
        raw_plots[[i]] <- NULL
        next
      }
      
      
      res <- do_Exigo(
        inputdata = inputdata,
        value     = p$value,
        batch     = batch,
        sex       = sex,
        y_title   = p$y_title,
        path_images = path_images,
        normal_range = p$normal_range %||% NULL,
        lowlimit     = p$lowlimit %||% NULL,
        hilimit      = p$hilimit %||% NULL
      )
      
      raw_pvals[i]  <- res$p_value_raw
      raw_plots[[i]] <- res$plot_raw
    }
    
    
    # PASS 2 — adjust p-values
    
    message("STEP 2: adjusting p-values using method: ", p_adjust_method)
    
    adj_pvals <- p.adjust(raw_pvals, method = p_adjust_method)
    # HEATMAP DATA (means + fold change)
    # ============================================================
    
    heat_df <- tibble(
      parameter = sapply(param_list, `[[`, "value"),
      p_adj = adj_pvals
    ) %>%
      rowwise() %>%
      mutate(
        ctrl_mean = mean(inputdata %>%
                           filter(Treatment == "ctrl") %>%
                           filter(case_when(sex == "both" ~ Sex %in% c("male","female"),
                                            TRUE ~ Sex == sex)) %>%
                           pull(parameter), na.rm = TRUE),
        tam_mean  = mean(inputdata %>%
                           filter(Treatment == "TAM") %>%
                           filter(case_when(sex== "both" ~ Sex %in% c("male","female"),
                                            TRUE ~ Sex == sex)) %>%
                           pull(parameter), na.rm = TRUE),
        log2FC = log2(tam_mean / ctrl_mean)
      ) %>%
      ungroup()
    param_order <- sapply(param_list, `[[`, "value")
    heat_df <- heat_df %>%
      mutate(parameter = factor(parameter, levels = rev(param_order)))
    heat_long <- heat_df %>%
      pivot_longer(cols = c(ctrl_mean, tam_mean),
                   names_to = "Treatment",
                   values_to = "Mean") %>%
      mutate(Treatment = recode(Treatment,
                                ctrl_mean = "ctrl",
                                tam_mean  = "TAM"))
    heat_long <- heat_long %>%
      mutate(parameter = factor(parameter, levels = rev(param_order)))
    p_heat_mean <- ggplot(heat_long,
                          aes(x = Treatment,
                              y = parameter,
                              fill = p_adj)) +
      geom_tile(color = "white", linewidth = 0.6) +
      # geom_text(aes(label = format.pval(p_adj, digits = 2)),
                # size = 3.5) +
     geom_text(aes(label = round(Mean,2)),
      size = 3.5) +
      scale_fill_gradientn(
        colors = c( "#fdae61", "white","blue"),
        values = scales::rescale(c(0, 0.01, 0.05, 1)),
        limits = c(0, 1),
        name = "adj p-value"
      )+
      theme_minimal(base_size = 13) +
      theme(axis.title = element_blank(),
            axis.text = element_text(face = "bold"),
            panel.grid = element_blank()) +
      labs(title = paste0("Legendplex – Mean values (", sex, ")"))
    
    ggsave(
      paste0(ExpId, "_Heatmap_Mean_", sex, ".png"),
      p_heat_mean,
      path = path_images,
      width = 6,
      height = 0.4 * nrow(heat_df) + 2,
      dpi = 300
    )
    
    p_heat_fc <- ggplot(heat_df,
                        aes(x = "TAM / ctrl",
                            y = parameter,
                            fill = p_adj)) +
      geom_tile(color = "white", linewidth = 0.6) +
      geom_text(aes(label = round(p_adj,3)),
                size = 3.5) +
      scale_fill_gradientn(
        colors = c( "#fdae61", "white","blue"),
        values = scales::rescale(c(0, 0.01, 0.05, 1)),
        limits = c(0, 1),
        name = "adj p-value"
      )+
      theme_minimal(base_size = 13) +
      theme(axis.title = element_blank(),
            axis.text = element_text(face = "bold"),
            panel.grid = element_blank()) +
      labs(title = paste0("Legendplex – Fold change (", sex, ")"))
    
    ggsave(
      paste0(ExpId, "_Heatmap_FoldChange_", sex, ".png"),
      p_heat_fc,
      path = path_images,
      width = 4,
      height = 0.4 * nrow(heat_df) + 2,
      dpi = 300
    )
    # PASS 3 — annotate plots & save
    
    message("STEP 3: writing plots with adjusted p-values")
    
    final_plots <- vector("list", n)
    
    for (i in seq_along(param_list)) {
      p <- param_list[[i]]
      
      pval_label <- paste0("adj p = ", format.pval(adj_pvals[i], digits = 3))
      
      # Position: ~15% above max
      y_max <- max(inputdata[[p$value]], na.rm = TRUE)
      y_min <- min(inputdata[[p$value]], na.rm = TRUE)
      y_pos <- y_max + 0.15*(y_max-y_min)
      
      # Add p-value annotation
      annotated_plot <- raw_plots[[i]] +
        annotate("text", x = 1.5, y = y_pos,
                 label = pval_label, size = 6,
                 fontface = "italic", color = "black")
      
      final_plots[[i]] <- annotated_plot
      
      # Save
      fname_val <- gsub("[^[:alnum:]_]", "_", p$value)
      file_name <- paste0(ExpId,"_", fname_val, "_Batch", batch, "_", sex, ".png")
      
      ggsave(file_name, plot = annotated_plot, path = path_images, width=4, height=11, dpi=300)
      
    }
    StatsOutput<-tibble(parameter = sapply(param_list, `[[`, "value"),
                        p_raw = raw_pvals,
                        p_adj = adj_pvals)
    write.csv2( StatsOutput, file = file.path(paste0(path_images, "/Exigo_StatsOutput_", sex,".csv")),
                row.names = FALSE,  na = "",  fileEncoding = "UTF-8")
    # --------------------------
    # Return everything as a table
    # --------------------------
    return(tibble(
      parameter = sapply(param_list, `[[`, "value"),
      p_raw = raw_pvals,
      p_adj = adj_pvals,
      plot = final_plots))
  }
  run_all_LP_onecall(inputdata = data, param_list = Legendplex_Inflammation,
                     batch = "ALL",sex = "male",  
                     path_images = "02_GeneratedData/Legenplex")
  run_all_LP_onecall(inputdata = data, param_list = Legendplex_Inflammation,
                     batch = "ALL",sex = "female",  
                     path_images = "02_GeneratedData/Legenplex")
  run_all_LP_onecall(inputdata = data, param_list = Legendplex_Inflammation,
                     batch = "ALL",sex = "both",  
                     path_images = "02_GeneratedData/Legenplex")
  