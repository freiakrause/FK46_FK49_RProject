rm(list=ls())
gc()
#### bei nächsten mal hier weiter macehn.
#ANimal Level columns sind noch nicht vollständig und müssen beim preprocessing kontrolliert werden. 
#Was will ich für parameter plotten
#Dann muss geschuat werden, ob statisti und ANimal colmns die gleichen namen haben 
#damit der p value in den richtigen plot kommt
library(dplyr)
library(ggplot2)
source("FK49_Definitions.R")

Animal_Level <- readRDS(file.path(PATHS$TEM$output,"CleanData/Animal_Dataframe.rds")) %>%
  select( -contains("nonHep"),
    -n_Images,-contains("SD"),
    -contains("Heps_per_Image"),
    -contains("Mito_per_Hep"),
    -contains("Hep_Area"),
    -contains("Mito_Area")) %>%
  rename_with(~ sub("^Mean_", "", .x), starts_with("Mean_")) %>%
  rename_with(~ sub("Mito_Length", "Length", .x), starts_with("Mito_Length")) %>%
  rename_with(~ sub("Mito_Width", "Width", .x), starts_with("Mito_Width")) %>%
  mutate( Treatment = factor(Treatment, levels = c("EtOH", "TAM")),
           Sex = factor(Sex, levels = c("female", "male")))

Animal_values <- Animal_Level%>%
  select(-Animal,-Treatment,  -Sex,-starts_with("SD"))%>%
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


for(value in Animal_values){
  
  y_pos <- max(Animal_Level[[value]], na.rm = TRUE)
  y_pos <- y_pos + 0.03 * (y_pos - min(Animal_Level[[value]], na.rm = TRUE))
  
  plot <- ggplot(Animal_Level, aes(x=Treatment, y=.data[[value]])) +
  stat_summary(fun.data = mean_sdl,fun.args = list(mult = 1), geom = "errorbar",width = 0.1,color = "black") +
    geom_jitter(aes(fill=Treatment,shape = Sex), width=0.1, size=4, alpha=0.5, color="black") +
    scale_fill_manual(name= "Treatment",values = c("#4D4D4DBF", "#8B0000BF"))+
    scale_shape_manual(name = "Sex", values = c(22, 24))+
    labs(y=value, x="") +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))) +
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
  geom_text(data=Animal_Level,aes(x = 1.5, y = y_pos, 
                label = paste0("adj.p = ", signif(Statistics$adj_Treatment_p[match(value, Animal_Level$value)], 3))),
            inherit.aes = FALSE, size = 3, vjust = -0.5)
  
  
  print(plot)
  
  ggsave(plot=plot, filename = paste0(value,".png"),
         path= file.path(PATHS$TEM$output,"/Plots/"),
         dpi= 300, width=2, height=5.5 )
}
