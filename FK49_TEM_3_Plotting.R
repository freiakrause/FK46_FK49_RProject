rm(list=ls())
gc()
#### bei nächsten mal hier weiter macehn.
#ANimal Level columns sind noch nicht vollständig und müssen beim preprocessing kontrolliert werden. 
#Was will ich für parameter plotten
#Dann muss geschuat werden, ob statisti und ANimal colmns die gleichen namen haben 
#damit der p value in den richtigen plot kommt
library(dplyr)
library(ggplot2)
library(tibble)
library(pheatmap)

source("FK49_Definitions.R")

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

# # color_map <- switch(plot_comp,
#                     Sex = Sex_colors,
#                     Treatment = Treatment_colors,
#                     T_D_S = T_D_S_colors,
#                     T_S = T_S_colors,
#                     Diet = Diet_colors,
#                     Timepoint = c("red","blue"),
#                     Time_Treat = c("red","blue","purple","royalblue"),
#                     stop("Unknown comparison for color mapping."))

# Main Plots for Results ------
for(value in Animal_values){
  
  y_pos <- max(Animal_Level[[value]], na.rm = TRUE)
  y_pos <- y_pos + 0.03 * (y_pos - min(Animal_Level[[value]], na.rm = TRUE))
  adj_p <- Statistics %>%
    filter(Variable == value) %>%
    pull(adj_Treatment_p)
  
  y_label <- case_when(
    value == "Mito_pro_Area"    ~ "Mitochondria/µm²",
    value == "Freq_damaged"     ~ "Visually damaged Mitochondria [%]",
    value == "Length"           ~ "Length [µm]",
    value == "Width"            ~ "Width [µm]",
    value == "Area"             ~ "Area [µm²]",
    value == "Circularity"      ~ "Circularity",
    TRUE                        ~ value
  )
  
  
  
  
  plot <- ggplot(Animal_Level, aes( x=Treatment, y=.data[[value]])) +
    stat_summary(fun.data = mean_sdl,fun.args = list(mult = 1), geom = "errorbar",width = 0.1,color = "black") +
    geom_jitter(aes(fill=Treatment,shape = Sex), width=0.1, size=4, alpha=0.5, color="black") +
    scale_fill_manual(name= "Treatment",values = Treatment_colors)+
    scale_shape_manual(name = "Sex", values = c(22, 24))+
    labs(y=y_label, x="") +
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
           fill = guide_legend(override.aes = list(shape = 22,fill = c("#4D4D4DBF", "#8B0000BF"),colour = "black", alpha = 0.8), ncol=1 ))+
    annotate("text",x = 1.5, y = y_pos, 
              label = paste0("adj.p = ", signif(adj_p, 3)),
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

PCA_data <- Mito_Level %>%select(Animal,Sex,Treatment, Mito_ID, Area, Circularity, Width, Length) %>%
  filter( complete.cases(Area, Circularity,Width,Length))

PCA <- prcomp(PCA_data%>%select(Area,Circularity,Width,Length),scale. = TRUE)
PCA_scores <- as.data.frame(PCA$x) %>%
  bind_cols(PCA_data %>%select(Animal, Sex, Treatment, Mito_ID))

ggplot(PCA_scores, aes(x = PC1, y = PC2, color = Treatment)) +
  geom_point(alpha = 0.3, size = 1) +
  theme_classic() +
  labs(x = paste0("PC1 (",round(summary(PCA)$importance[2, 1] * 100, 1),"%)"),
    y = paste0( "PC2 (",round(summary(PCA)$importance[2, 2] * 100, 1),"%)") )

round(PCA$rotation[, 1:2], 3)


## Cell PCA -----
Cell_Level <- readRDS(file.path(PATHS$TEM$output,"CleanData/Cell_Dataframe.rds")) %>%
  filter(Cell_Type != "nonHep")%>%
  select(-contains("SD"),-Image_ID,-Cell_ID,-Cell_Area)

PCA_data <- Cell_Level %>%
  select(Animal,Sex,Treatment,n_Mito, Mito_pro_Area, Mean_Mito_Area,
    Mean_Length,Mean_Width, Mean_Circularity, relMitoArea,Freq_Donut,
    Freq_LostCristae) %>%
  filter(complete.cases( n_Mito,  Mito_pro_Area,
      Mean_Mito_Area, Mean_Length,Mean_Width, Mean_Circularity,
      relMitoArea,Freq_Donut,  Freq_LostCristae ))

PCA <- prcomp(PCA_data %>%select( Mito_pro_Area, Mean_Length,Mean_Width, Mean_Circularity, Mean_Mito_Area),  scale. = TRUE)
PCA_scores <- as.data.frame(PCA$x) %>%bind_cols(PCA_data %>%select(Animal, Sex, Treatment))

ggplot(PCA_scores, aes(x = PC1, y = PC2, color = Treatment)) +
  geom_point(alpha = 0.3, size = 1) +
  theme_classic() +
  labs(x = paste0("PC1 (",round(summary(PCA)$importance[2, 1] * 100, 1),"%)"),
      y = paste0( "PC2 (",round(summary(PCA)$importance[2, 2] * 100, 1),"%)") )
round(PCA$rotation[, 1:2], 3)


## Animal PCA -----
PCA_data <- Animal_Level %>%
  select(Animal, Sex, Treatment, Length, Width, Circularity,Mito_pro_Area,Area,Freq_Donut, Freq_LostCristae) %>%
  filter(complete.cases(Length, Width, Circularity, Mito_pro_Area,Area,Freq_Donut,Freq_LostCristae))

PCA <- prcomp(PCA_data %>%select(Length, Width, Circularity,Area, Mito_pro_Area,Freq_Donut, Freq_LostCristae), scale. = TRUE)

PCA_scores <- as.data.frame(PCA$x) %>%bind_cols(PCA_data %>%select(Animal, Sex, Treatment))

ggplot(PCA_scores, aes(x = PC1, y = PC2, color = Treatment)) +
  geom_point(alpha = 0.3, size = 1) +
  theme_classic() +
  labs(x = paste0("PC1 (",round(summary(PCA)$importance[2, 1] * 100, 1),"%)"),
      y = paste0( "PC2 (",round(summary(PCA)$importance[2, 2] * 100, 1),"%)") )
round(PCA$rotation[, 1:2], 3)



## Animal Heatmap -----

Heatmap_data <- Animal_Level %>%
  select(Animal, Sex, Treatment, Length, Width, Circularity,Area, 
         Mito_pro_Area,Freq_damaged) %>%
  column_to_rownames("Animal")

parameter_labels <- c(Mito_pro_Area    = "Mitochondria/µm²",
  Freq_Donut       = "Donut-shaped mitochondria [%]",
  Freq_LostCristae = "Mitochondria with lost cristae [%]",
  Freq_damaged = "Mitochondria visually damaged [%]",
  Length           = "Length [µm]",
  Width            = "Width [µm]",
  Area             = "Area [µm²]",
  Perimeter        = "Perimeter [µm]",
  Circularity      = "Circularity")


Heatmap_matrix <- Heatmap_data %>%select(Length, Width, Circularity,Area, Mito_pro_Area,Freq_damaged)
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
          show_colnames = TRUE)

ggsave(Mito_Heatmap,file= file.path(PATHS$TEM$output,"Plots/Heatmap.png"),dpi= 300,
       width = 10, height= 9,bg = 'white')
