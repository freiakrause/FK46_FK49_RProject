rm(list=ls())
gc()

library(tidyverse)
library(pheatmap)
source("FK49_Definitions.R")
ExpId = "FK46"

# Read Raw Inputdata after general Data manipulation ------------------------------------------------------
load(file.path(PATHS$general_data[[paste0(ExpId,"_output")]], "01_RawData", paste0(ExpId,"_Data_prepared.Rda")))

data<-data%>%
  select(Animal, Sex, BATCH, Treatment,Score,wks_diet,wks_dead,Death)%>%
  mutate(wks_dead=round(wks_dead,1))
# Base animal metadata (only one row per animal)
animal_info <- data %>%
  select(Animal, Sex, BATCH, Treatment, wks_dead,Death) %>%
  distinct()

# All week × animal combinations
all_combos <- expand.grid(
  Animal = unique(data$Animal),
  wks_diet = sort(unique(data$wks_diet)))

# Join week × animal grid with animal info
data_full <- all_combos %>%
  left_join(animal_info, by = "Animal") %>%   # add stable info
  left_join(data %>% select(Animal, wks_diet, Score), 
          by = c("Animal", "wks_diet"))    # add Score when available
  

# 3. Assign Score_display for missing rows
data_full <- data_full %>%
  mutate(Score_display = case_when( is.na(Score) & wks_diet > wks_dead & Death == "1" ~ "\u2020",
                                    is.na(Score) & wks_diet > wks_dead & Death == "0" ~ "c",
                                    is.na(Score)& wks_diet < wks_dead ~ "-",
                                    TRUE ~ as.character(Score)),
         Score = case_when(Score_display == "\u2020"  ~ -2,
                           Score_display == "c"  ~ -3,
                           Score_display == "-"    ~ -1,
                           TRUE ~ Score  ))
  
  


# filtered <- inputdata %>%
#   filter(case_when(
#     sex == "female" ~ Sex == "female",
#     sex == "male" ~ Sex == "male",
#     sex == "both" ~ TRUE))
  common_timepoints <- data_full %>%       # Find common timepoints across both batches
    filter(BATCH %in% c(1, 2)) %>%
    group_by(wks_diet, BATCH) %>%
    summarise(n = n(), .groups = "drop") %>%
    group_by(wks_diet) %>%
    summarise(n_batches = n_distinct(BATCH)) %>%
    filter(n_batches == 2) %>%
    pull(wks_diet)
  
  data_full <- data_full %>%
    filter(BATCH %in% c(1, 2)) %>%
    filter(wks_diet!=24.3)
  
  data_full <- data_full %>%
    group_by(Animal, wks_diet,Treatment,Sex,Score_display,wks_dead,BATCH) %>%
    summarise(Score = first(Score), .groups = "drop") %>%# take first if duplicates
    arrange(Treatment,wks_dead,Sex) 
  
  # Numeric matrix for coloring
  mat <- data_full %>%
    select(Animal,wks_diet, Score) %>%
     arrange(wks_diet)%>%
    pivot_wider(names_from = wks_diet, values_from = Score) %>%
    column_to_rownames("Animal")%>%
    as.matrix()
  
  # Keep annotations
  ann <- data_full %>%
    select(Animal, Treatment,  Sex) %>%#,BATCH
    distinct() %>%
    column_to_rownames("Animal")
  ann <- ann[rownames(mat), ]  # make sure order matches

  
  # Transpose for pheatmap (animals as columns)
 # mat <- t(mat)
  
  
  # Create display matrix
  display_mat <- data_full %>%
    select(Animal, wks_diet, Score_display) %>%
    group_by(Animal,wks_diet)%>%
    summarise(Score_display=first(Score_display), .groups = "drop")%>%
    pivot_wider(names_from = wks_diet, values_from = Score_display) %>%
    column_to_rownames("Animal") %>%
    as.matrix()#%>%
    #t()  
  display_mat <- display_mat[rownames(mat), colnames(mat)]
  

  heat_color<-c( "grey70","grey80","#F0F0F0","white",colorRampPalette(c("white", "orange", "red"))(20))
                
  breaks <-c(-3,-2,-1, seq(0,max(mat,na.rm=TRUE),length.out = 21))   
  sex_colors <- Sex_colors  # from Definitions
  colnames(ann) <- c("T", "S")#, "B"
  ann$T <- as.factor(ann$T)
  ann$S <- as.factor(ann$S)
  #ann$B <- as.factor(ann$B) 
  treatment_colors <- Treatment_colors[c("Ctrl","TAM")]  # from Definitions
  #batch_colors<- c("1"= "#D9D9D9","2"= "#BFD8FF")
  annotation_colors=list(T = treatment_colors, S = sex_colors)#,B=batch_colors
  annotation_names_col <- c(T = "Treatment", S = "Sex") #, B = "BATCH"

  h <- pheatmap(mat,
                cluster_rows = FALSE,        
                cluster_cols = F,
                fontsize        = 8,
                fontsize_main   = 8,
                fontsize_row    = 6,
                fontsize_col    = 6,
                fontsize_number = 5,
                annotation_row  = ann,
                annotation_names_col =T,
                annotation_legend_side = "right", 
                color = heat_color,
                #gaps_col = 27,
                #gaps_row = 14,
                annotation_colors = annotation_colors,
                main = "Scoring Heatmap iAL on CD-HFD",
                labels_col = paste0(round(as.numeric(as.character(colnames(mat))),0)),
                labels_row = rownames(mat),
                border_color = "black",
                cellwidth = 10,
                cellheight = 10,
                angle_col = 0,
                display_numbers = display_mat,
                number_color = "grey30",
                # number_format = "%.0f",
                legend_breaks = c(0,5,10,15,20),
                legend_labels = c("0","5","10","15","20"))
  
ggsave(plot= h, filename= paste0(ExpId, "_Scores_HeatMap.png"), path = file.path(PATHS$general_data[[paste0(ExpId,"_output")]], "02_GeneratedData"), height = 20,width = 25, dpi = 300 )
h
