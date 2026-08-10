rm(list=ls())
gc()

library(dplyr)
library(ggplot2)
library(tibble)
library(pheatmap)

source("FK49_Definitions.R")
# Load Data and Statistics Data
# rename so that data columns and statistics columns that belong together have the same name
Animal_Level <- readRDS(file.path(PATHS$TEM$output,"CleanData/Animal_Dataframe.rds")) %>%
  select( -n_Images,-Mean_Perimeter,-Mean_Freq_StandardMitochondrium,-contains("SD"),-contains("nonHep"),-contains("Heps_per_Image"),
    -contains("Mito_per_Hep"),-contains("Hep_Area")) %>%
  rename_with(~ sub("^Mean_", "", .x), starts_with("Mean_")) %>%
  rename_with(~ sub("Mito_Length", "Length", .x), starts_with("Mito_Length")) %>%
  rename_with(~ sub("Mito_Width", "Width", .x), starts_with("Mito_Width")) %>%
  rename_with(~ sub("Mito_Area", "Area", .x), starts_with("Mito_Area")) %>%
  mutate(Treatment = factor(Treatment, levels = c("EtOH", "TAM")),
         Sex      = factor(Sex, levels = c("female", "male")))

levels(Animal_Level$Treatment) <- c("Ctrl","TAM")

Animal_values <- Animal_Level%>%
  select(-Animal,-Treatment,  -Sex,-Freq_Donut,-Freq_LostCristae,-starts_with("SD"))%>%
  colnames()
colnames(Animal_Level)
Statistics<-read.csv2(file=file.path(PATHS$TEM$output,"Statistics/01_Total_Stats.csv"))
Statistics$Variable

# Main Plots for Results ------
for(value in Animal_values){
  
#x generate position for p value label according to max y in plot  
  y_pos <- max(Animal_Level[[value]], na.rm = TRUE)
  y_pos <- y_pos + 0.03 * (y_pos - min(Animal_Level[[value]], na.rm = TRUE))
  adj_p <- Statistics %>%
    filter(Variable == value) %>%
    pull(adj_Treatment_p)
  adj_p <- ifelse(adj_p>=0.08,"ns",signif(adj_p, 3))
# Generate y label according to used variable  
  y_label <- case_when(
    value == "Mito_pro_Area"    ~ "Mitochondria/µm²",
    value == "Freq_damaged"     ~ "Damaged Mitochondria per Cell [%]",
    value == "Length"           ~ "Length [µm]",
    value == "Width"            ~ "Width [µm]",
    value == "Area"             ~ "Area [µm²]",
    value == "Circularity"      ~ "Circularity",
    TRUE                        ~ value
  )
#plot
  plot <- ggplot(Animal_Level, aes( x=Treatment, y=.data[[value]])) +
    stat_summary(fun.data = mean_sdl,fun.args = list(mult = 1), geom = "errorbar",width = 0.1,color = "black") + # mean+SD
    stat_summary(fun = mean, geom = "crossbar",linewidth= 0.2, width= 0.2) +                                    #crossbar for mean
    geom_jitter(aes(fill=Treatment,shape = Sex), width=0.1, size=4, alpha=0.5, color="black") +
    scale_fill_manual(name= "Treatment",values = Treatment_colors)+
    scale_shape_manual(name = "Sex", values = c(22, 24))+
    labs(y=y_label, x="") +
    #for freq i want to show them as perc because easier to read 
    scale_y_continuous(  labels = if (value %in% c("Freq_Donut","Freq_LostCristae", "Freq_damaged","Freq_StandardMitochondrium")) {
                                      function(x) paste0(x * 100, "%")
                                  } else {
                                    waiver()},
      expand = expansion(mult = c(0.05, 0.15)))+
    theme_classic() +
    theme(legend.position = "bottom",
          legend.box = "vertical",
          legend.spacing.y = unit(0, "mm"),
          legend.spacing.x = unit(0.5, "mm"),
          legend.box.spacing = unit(1, "mm"),
          legend.key.size = unit(2, "mm"),
          legend.key.width = unit(2, "mm"),
          legend.key.height = unit(2, "mm"),
          legend.key = element_blank(),
          legend.title = element_text(size = 8, face = "bold"),
          legend.text = element_text(size = 7),
          axis.line = element_line(color = "black", linewidth = 0.5),
          axis.ticks = element_line(color = "black", linewidth = 0.5),
          axis.title = element_text(size = 12, face = "bold"),
          axis.title.x = element_blank(),
          axis.text = element_text(size = 11, face = "bold"),
          strip.text = element_text(size = 10, hjust = 0.5) )+
    guides(shape = guide_legend(override.aes = list(  fill = "grey80",colour = "black", size = 4,alpha = 1), ncol=1),
           fill = guide_legend(override.aes = list(shape = 22,fill = Treatment_colors[c("Ctrl","TAM")],colour = "black", alpha = 0.5), ncol=1 ))+
    annotate("text",x = 1.5, y = y_pos, 
              label = paste0("adj.p = ", adj_p),
               size = 3, vjust = -0.5)
  print(plot)
  
  ggsave(plot=plot, filename = paste0(value,".png"),
         path= file.path(PATHS$TEM$output,"/Plots/"),
         dpi= 300, width=2, height=5.5 )
}
rm(Animal_values,Statistics)

#Fun Plots for Me -----
## Mito PCA -----
Mito_Level <- readRDS(file.path(PATHS$TEM$output,"CleanData/Mito_Dataframe.rds")) %>%
  filter(!Cell_Type== "nonHep")

PCA_data <- Mito_Level %>%select(Animal,Sex,Treatment, Cell_ID,Image_ID,Mito_ID, Area, Circularity, Width, Length) %>%
  filter( complete.cases(Area, Circularity,Width,Length))

PCA <- prcomp(PCA_data%>%select(Area,Circularity,Width,Length),scale. = TRUE)
PCA_scores <- as.data.frame(PCA$x) %>%
  bind_cols(PCA_data %>%select(Animal, Sex, Treatment,Image_ID,Cell_ID, Mito_ID))

Mito_PCA <- ggplot(PCA_scores, aes(x = PC1, y = PC2, color = Animal)) +
  geom_point(alpha = 0.3, size = 1) +
  theme_classic() +
  labs(x = paste0("PC1 (",round(summary(PCA)$importance[2, 1] * 100, 1),"%)"),
    y = paste0( "PC2 (",round(summary(PCA)$importance[2, 2] * 100, 1),"%)") )

round(PCA$rotation[, 1:2], 3)
ggsave(Mito_PCA,file= file.path(PATHS$TEM$output,"Plots/PCA_Mito.png"),dpi= 300,
       width = 5, height= 4.5,bg = 'white')
rm(Mito_PCA)

## Cell PCA -----
Cell_Level <- readRDS(file.path(PATHS$TEM$output,"CleanData/Cell_Dataframe.rds")) %>%
  filter(Cell_Type != "nonHep")%>%
  select(-contains("SD"),-Cell_Area)

PCA_data <- Cell_Level %>%
  select(Animal,Sex,Treatment,n_Mito, Mito_pro_Area, Mean_Mito_Area,
    Mean_Length,Mean_Width, Mean_Circularity, relMitoArea,Freq_damaged) %>%
  filter(complete.cases( n_Mito,  Mito_pro_Area,
      Mean_Mito_Area, Mean_Length,Mean_Width, Mean_Circularity,
      relMitoArea,Freq_damaged ))

PCA <- prcomp(PCA_data %>%select(Mito_pro_Area, Mean_Length,Mean_Width, Mean_Circularity, Mean_Mito_Area,Freq_damaged),  scale. = TRUE)
PCA_scores <- as.data.frame(PCA$x) %>%bind_cols(PCA_data %>%select(Animal, Sex, Treatment))

Cell_PCA <- ggplot(PCA_scores, aes(x = PC1, y = PC2, color = Animal)) +
  geom_point(alpha = 0.8, size = 1) +
  theme_classic() +
  labs(x = paste0("PC1 (",round(summary(PCA)$importance[2, 1] * 100, 1),"%)"),
      y = paste0( "PC2 (",round(summary(PCA)$importance[2, 2] * 100, 1),"%)") )
round(PCA$rotation[, 1:2], 3)
ggsave(Cell_PCA,file= file.path(PATHS$TEM$output,"Plots/PCA_Cell.png"),dpi= 300,
       width = 5, height= 4.5,bg = 'white')
rm(Cell_PCA)

## Image PCA -----
Image_Level <- readRDS(file.path(PATHS$TEM$output,"CleanData/Image_Dataframe.rds")) %>%
  select(-contains("SD"),-contains("nonHep_"),-Descriptor,-n_Cells,-n_Heps,-n_nonHeps,-Mean_Hep_Area,
         -Mean_Freq_LostCristae,-Mean_Freq_Donut,-Mean_Freq_StandardMitochondrium,-Mean_Perimeter)

PCA_data <- Image_Level %>%
  filter(complete.cases(Mean_Mito_Area,Mean_Mito_Length,Mean_Mito_Width, Mean_Circularity, 
                        Mean_Mito_pro_Area ,Mean_relMitoArea,Mean_Freq_damaged))

PCA <- prcomp(PCA_data %>%select(Mean_Mito_Area,Mean_Mito_Length,Mean_Mito_Width, Mean_Circularity, 
                                 Mean_Mito_pro_Area ,Mean_relMitoArea,Mean_Freq_damaged),  scale. = TRUE)
PCA_scores <- as.data.frame(PCA$x) %>%bind_cols(PCA_data %>%select(Animal, Sex, Treatment))

Image_PCA <- ggplot(PCA_scores, aes(x = PC1, y = PC2, color = Animal)) +
  geom_point(alpha = 0.8, size = 1) +
  theme_classic() +
  labs(x = paste0("PC1 (",round(summary(PCA)$importance[2, 1] * 100, 1),"%)"),
       y = paste0( "PC2 (",round(summary(PCA)$importance[2, 2] * 100, 1),"%)") )
round(PCA$rotation[, 1:2], 3)
ggsave(Image_PCA,file= file.path(PATHS$TEM$output,"Plots/PCA_Image.png"),dpi= 300,
       width = 5, height= 4.5,bg = 'white')
rm(Image_PCA)
## Animal PCA -----
PCA_data <- Animal_Level %>%
  select(Animal, Sex, Treatment, Length, Width, Circularity,Mito_pro_Area,Area,Freq_damaged) %>%
  filter(complete.cases(Length, Width, Circularity,Area, Mito_pro_Area,Freq_damaged))

PCA <- prcomp(PCA_data %>%select(Length, Width, Circularity,Area, Mito_pro_Area,Freq_damaged), scale. = TRUE)
PCA_scores <- as.data.frame(PCA$x) %>%bind_cols(PCA_data %>%select(Animal, Sex, Treatment))

Animal_PCA<-ggplot(PCA_scores, aes(x = PC1, y = PC2, color = Treatment)) +
  geom_point(alpha = 0.5, size = 4) +
  theme_classic() +
  scale_color_manual(values=Treatment_colors[c("Ctrl","TAM")])+
  labs(x = paste0("PC1 (",round(summary(PCA)$importance[2, 1] * 100, 1),"%)"),
      y = paste0( "PC2 (",round(summary(PCA)$importance[2, 2] * 100, 1),"%)") )
round(PCA$rotation[, 1:2], 3)

ggsave(Animal_PCA,file= file.path(PATHS$TEM$output,"Plots/PCA_Animal.png"),dpi= 300,
       width = 5, height= 4.5,bg = 'white')
rm(Animal_PCA)
## Animal Heatmap -----

Heatmap_data <- Animal_Level %>%
  select(Animal, Sex, Treatment, Length, Width, Circularity,Area, Mito_pro_Area,Freq_damaged) %>%
  column_to_rownames("Animal")

parameter_labels <- c(Mito_pro_Area    = "Mitochondria/µm²",
  Freq_Donut       = "Donut-shaped mitochondria [%]",
  Freq_LostCristae = "Mitochondria with lost cristae [%]",
  Freq_damaged = "Damaged Mitochondria per Cell [%]",
  Length           = "Length [µm]",
  Width            = "Width [µm]",
  Area             = "Area [µm²]",
  Perimeter        = "Perimeter [µm]",
  Circularity      = "Circularity")


Heatmap_matrix <- Heatmap_data %>%select(Length, Width, Circularity, Area, Mito_pro_Area,Freq_damaged)
Heatmap_matrix<-t(Heatmap_matrix)
rownames(Heatmap_matrix) <- parameter_labels[rownames(Heatmap_matrix)]
Heatmap_Annos_col<-Heatmap_data%>%select(Sex,Treatment)
Heatmap_colors <-list(Treatment = Treatment_colors[c("Ctrl","TAM")],
                      Sex         = Sex_colors)
Mito_Heatmap<-pheatmap(Heatmap_matrix, 
          cluster_rows = TRUE,
          cluster_cols = TRUE,
          scale = "row",
          annotation_col = Heatmap_Annos_col,
          annotation_colors=Heatmap_colors,
          show_rownames = TRUE, 
          show_colnames = FALSE)

ggsave(Mito_Heatmap,file= file.path(PATHS$TEM$output,"Plots/Heatmap.png"),dpi= 300,
       width = 7, height=6,bg = 'white')
