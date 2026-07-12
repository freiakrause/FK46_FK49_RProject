gc()
rm(list = ls())
home <- normalizePath("~") # bc Windows does not start inuser dir but in user/documents dir
parent <- dirname(home)
setwd(paste0(parent,"/OneDrive - Universität Salzburg/AG_Tumorimmunologie - Freia Krause/01_Experiments/FK46_FK49_RProject"))

library(dplyr)
library(tidyr)
library(ggplot2)
library('corrr')
library(ggcorrplot)
library("FactoMineR")
library(factoextra)
library(ggrepel)

source("FK49_Definitions.R")
home <- normalizePath("~") # bc Windows does not start inuser dir but in user/documents dir
parent <- dirname(home) # gives users dir works for me since i aways have the data from work in the onedrive at the same location and only local home dir changes
ExpId = "FK49"

#BA contains all animal but not 199 since it was massive outlier. technical problem, quality problem or Cholemia
#BA filtered does not contain 199, 186 and 187 since 186 and 187 were also a little bnit ouliters but not sure if this shows Cholemia in them. So I  want to analys both datasets and disucc with PI/Peers which to include.
BA         <-readRDS(file= paste0(parent,"/OneDrive - Universität Salzburg/AG_Tumorimmunologie - Freia Krause/01_Experiments/FK49_CD-HFD_13wks/FK49_Analysis/01_RawData/FK49_BA_preprocessed.rds"))
BA_filtered<-readRDS(file= paste0(parent,"/OneDrive - Universität Salzburg/AG_Tumorimmunologie - Freia Krause/01_Experiments/FK49_CD-HFD_13wks/FK49_Analysis/01_RawData/FK49_BAfiltered_preprocessed.rds"))

x<-readRDS(paste0(BApwd,"/BA_a/BA_a_filtered.rds"))
LMresults<-read.csv2(paste0(BApwd,"/BA_a/Stats/BA_a_filtered_LM_results.csv"),sep=";")
plot_BA <- function(data,
                    comparison = "Treatment",
                    top_n = 8,
                    data_type = c("log", "raw", "norm", "scaled"),
                    use_lm = FALSE,
                    XLIM=c(-9,9),
                    YLIM=c(0,7.5),
                    NAME= "",
                    method = NULL,
                    FOLDER = NULL
){
  # Load libraries
  library(dplyr)
  library(ggplot2)
  library(ggrepel)
  library(tidyr)
  library(pheatmap)
 
  PATH = BApwd
  # Determine which comparison to plot (first if multiple)
  plot_comp <- comparison[1]
  
  # Color mapping for plotting
  color_map <- switch(plot_comp,
                      Sex = Sex_colors,
                      Treatment = Treatment_colors,
                      T_D_S = T_D_S_colors,
                      T_S = T_S_colors,
                      Diet = Diet_colors,
                      Timepoint = c("red","blue"),
                      Time_Treat = c("red","blue","purple","royalblue"),
                      stop("Unknown comparison for color mapping."))
  
  
  

  
  top <- data %>% 
    group_by(Metabolite)%>%
    summarise(adj.p.value = first(adj.p.value))%>%
    arrange(adj.p.value)%>%
    slice_head(n=top_n)
  
  data<-data%>%filter(Metabolite %in% top$Metabolite)
  y_limits <- data %>%
    group_by(Metabolite) %>%
    summarise(
      y_min = min(value, na.rm = TRUE),
      y_max = max(value, na.rm = TRUE),
      adj.p.value = first(adj.p.value),
      adj_p_Trajectory = first(adj_p_Trajectory),
      .groups = "drop")
  
  if (all(c(-1, 11) %in% data$Timepoint)){
    p_spaghetti <- ggplot(data, aes(x = Timepoint, y = value, group = Animal, color = .data[[plot_comp]],)) +
      geom_line(alpha = 0.4, size= 0.1) +
      geom_point(size = 1.5,alpha = 0.5) +
      scale_color_manual(values = color_map) +
      facet_wrap(~Metabolite, scales = "free_y", ncol = 5,
                 labeller = label_wrap_gen(width = 15)) +
      labs(   y = "Concentration in Serum [µmol/L]",
              x = "Timepoint",title = paste0("Top ", top_n, " Metabolites (paired changes)") ) +
      theme_classic() +
      theme(strip.text = element_text(size = 10, hjust = 0.5))+
      stat_summary(aes(group = Treatment),fun = mean,geom = "line", size = 1.2, alpha = 0.9)+
    geom_text(data = y_limits,aes(x = 1.5,  y = y_max + 0.02 * (y_max - y_min), 
                                  label = paste0("adj.p = ", signif(adj_p_Trajectory, 3))),
              inherit.aes = FALSE,  size = 3,  vjust = 0)
    
    ggsave(plot= p_spaghetti, 
           filename= paste0(NAME,"_spaghetti",top_n,".png"),
           path =paste0(PATH,"/",FOLDER,"/" ),
           width= 18, height = 10, dpi = 300 )
  }
  
  p_bar_alltp <- ggplot(data,aes(x = Timepoint, y = value, fill = Treatment)) +
    stat_summary(fun = mean, geom = "bar", position = position_dodge(width = 0.8),width = 0.7) +
    stat_summary(fun.data = mean_sdl,
                 fun.args = list(mult = 1),  
                 geom = "errorbar", 
                 position = position_dodge(width = 0.8),  
                 width = 0.2,
                 color = "black") +
    geom_point(position = position_dodge(width = 0.8), shape = 21, color = "black", size = 1.5, alpha = 0.5) +
    scale_fill_manual(values = color_map) +
    facet_wrap(~Metabolite, scales = "free_y", ncol = 11, labeller = label_wrap_gen(width = 15)) +
    labs(  y = "Concentration in Serum [µmol/L]",x = "Timepoint",  title = paste0("Top ", top_n, " Metabolites by ", plot_comp)) +
    theme_classic() +
    theme(strip.text = element_text(size = 10, hjust = 0.5)) +
    geom_text(data = y_limits,aes(x = 1.5,  y = y_max + 0.05 * (y_max - y_min), label = paste0("adj.p = ", signif(adj.p.value, 3))),
              inherit.aes = FALSE,  size = 3,  vjust = 0)
  
  print(p_bar_alltp)
  
  ggsave( plot = p_bar_alltp, filename = paste0(NAME, "_barstps", top_n, ".png"),
          width = 18, height = 6, dpi = 300,path = paste0(PATH, "/", FOLDER, "/"))
  
  data<-data%>%filter(Timepoint==11)
  y_limits <- data %>%
    group_by(Metabolite) %>%
    summarise(y_min = min(value, na.rm = TRUE), y_max = max(value, na.rm = TRUE))
  
  # ---- Violin plot ----
  p_violin <- ggplot(data, aes(x=.data[[plot_comp]], y=value, fill=.data[[plot_comp]])) +
    geom_violin(trim=FALSE, alpha=0.5) +
    geom_boxplot(width=0.1, outlier.shape=NA) +
    geom_jitter(width=0.1, size=1.5, alpha=0.7) +
    scale_fill_manual(values=color_map) +
    facet_wrap(~Metabolite, scales="free_y", ncol = 5, labeller = label_wrap_gen(width = 15)) +
    labs(y="Concentration in Serum [µmol/L]", x=plot_comp, title=paste0("Top ", top_n, " Metabolites by ", plot_comp)) +
    theme_classic() +
    theme(strip.text = element_text(size = 10, hjust = 0.5)) +
    geom_blank(data = y_limits, aes(y = y_min, yend = y_max, group = Metabolite), inherit.aes = FALSE) +
    geom_text(data = y_limits,  aes(x = 1.5, y = y_max + 0.05 * (y_max - y_min), label = paste0("adj.p = ", signif(data$adj.p.value[match(Metabolite, data$Metabolite)], 3))),
              inherit.aes = FALSE, size = 3, vjust = -0.5)
  
  print(p_violin)
  ggsave(plot= p_violin, filename= paste0(NAME,"_violins",top_n,".png"), width= 18, height = 9, dpi = 300, path =paste0(PATH,"/",FOLDER,"/" ))
  
  
  p_bar <- ggplot(data ,aes(x=.data[[plot_comp]], y=value, fill=.data[[plot_comp]])) +
    stat_summary(fun = mean, geom = "bar", position = position_dodge()) +
    stat_summary(fun.data = mean_sdl,fun.args = list(mult=1), geom = "errorbar",width= 0.2, color= "black") +
    geom_jitter(aes(fill=.data[[plot_comp]]),color="black", shape=21,width=0.1, size=1.5, alpha=0.5) +
    scale_fill_manual(values=color_map) +
    facet_wrap(~Metabolite, scales="free_y", ncol = 11, labeller = label_wrap_gen(width = 15)) +
    labs(y="Concentration in Serum [µmol/L]",
         x= plot_comp, 
         title=paste0("Top ", top_n, " Bile Acids by ", plot_comp)) +
    theme_classic() +
    theme(strip.text = element_text(size = 10, hjust = 0.5)) +
    geom_blank(data = y_limits, aes(y = y_min, yend = y_max + 0.10 * (y_max - y_min), group = Metabolite), inherit.aes = FALSE) +
    geom_text(data = y_limits, 
              aes(x = 1.5, y = y_max + 0.03 * (y_max - y_min), 
                  label = paste0("adj.p = ", signif(data$adj.p.value[match(Metabolite, data$Metabolite)], 3))),
              inherit.aes = FALSE, size = 3, vjust = -0.5)
  print(p_bar)
  ggsave(plot= p_bar, filename= paste0(NAME,"_bars",top_n,".png"), width= 18, height = 6, dpi = 300, path =paste0(PATH,"/",FOLDER,"/" ))
  
}


plot_BA (data = x, 
         FOLDER = "BA_a",  NAME ="BA_a_filtered", comparison = c("Treatment"),   
         top_n = 9, data_type = "log", use_lm = TRUE, XLIM= c(-10,10), YLIM= c(0,5), 
         method = "targeted")

plot_BA (data = x, 
         FOLDER = "BA_a",  NAME ="BA_a_filtered", comparison = c("Treatment"),   
         top_n = 21, data_type = "log", use_lm = TRUE, XLIM= c(-10,10), YLIM= c(0,5), 
         method = "targeted")

# ---- Volcano plot ----
p_volcano <- ggplot(LMresults,aes(x = log2FC, y = negLog10FDR)) +
  geom_point(aes(fill = direction), alpha = 0.5, size = 3,stroke = 0.5,
             position=position_jitter(width= 0.08), shape=21,color= "black") +
  scale_fill_manual(values = c("blue","grey60", "firebrick")) +
  geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed", color = "grey80") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey80") +
  labs( title = paste0("Volcano plot - Treatment"),
    x = expression(paste("FC [", log[2], "]")),       
    y = expression(paste("-log"[10], "(adj.p.value)")))+
  theme_classic()+  
  theme(panel.grid= element_line(color ="grey90", linewidth = 0.1))+
  geom_text_repel(data = LMresults %>% filter(significant == TRUE), 
                  aes(label = Metabolite),
                  size = 4,
                  max.overlaps = 20) +
  coord_cartesian(xlim = c(-4, 4), ylim = c(0, 3))
  
print(p_volcano)
ggsave(plot= p_volcano, filename= paste0("BA_a_volcano.png"), width= 6, height = 9, dpi = 300, 
       path =paste0(BApwd,"/Ba_a/" ))

# ---- Heatmap ----
annotation_colors <- list(
  Sex = Sex_colors,
  Treatment = Treatment_colors,
  Diet = Diet_colors,
  Timepoint =c("-1"="red","11"="blue")
)
meta_cols <- c("Animal","Sample","Sex","Treatment","Diet","ExpID","T_D_S","T_D" ,"T_S","Time_Treat","Timepoint")
BA_heat<-BA_filtered%>%filter(Timepoint== 11)%>%select(-GDCA,-GCDCA)
matrix_heat<-BA_heat%>%select(-all_of(meta_cols))
matrix_heat <-as.matrix(matrix_heat)
rownames(matrix_heat)<-BA_heat$Sample
ann <- data.frame(
  Sex = BA_heat$Sex,
  Treatment = BA_heat$Treatment,
  Diet = BA_heat$Diet,
  Timepoint = BA_heat$Timepoint)
rownames(ann) <- BA_heat$Sample
ann$Treatment <- factor(ann$Treatment, levels = c("Ctrl", "TAM"))

ord <- order(ann$Treatment)
ann <- ann[ord, ]
matrix_heat<-t(matrix_heat)
matrix_heat <- matrix_heat[, ord]
heatmap_height <- nrow(matrix_heat)/8 + 3

# clustered
p_heat <- pheatmap::pheatmap(matrix_heat,
                             scale = "row",
                             cluster_rows = TRUE, 
                             cluster_cols = TRUE,
                             annotation_col  = ann,
                             annotation_colors = annotation_colors)
ggsave(plot= p_heat, filename= paste0("Ba_a_heatmap_clustered.png"),limitsize = FALSE, 
       width= 10, height = heatmap_height, dpi = 300,bg = "white", path =paste0(BApwd,"/Ba_a/" ))

# non-clustered
p_heat <- pheatmap::pheatmap(matrix_heat,
                             scale = "row",
                             cluster_rows = TRUE, 
                             cluster_cols = FALSE,
                             annotation_col  = ann,
                             annotation_colors = annotation_colors)
ggsave(plot= p_heat, filename= paste0("Ba_a_heatmap.png"),limitsize = FALSE, 
       width= 10, height = heatmap_height, dpi = 300,bg = "white", path =paste0(BApwd,"/Ba_a/" ))

# Stacke Barplot for relative abundance of Bile Acids
a<-x%>%
  filter(Timepoint == 11)%>%
   mutate(    Origin= factor(case_when(Metabolite %in% BA_primary ~ "1°",
                                        Metabolite %in% BA_secondary ~ "2°",
                                        TRUE ~ NA),
                             levels = c("1°","2°")),
              Conjugation= factor(case_when(Metabolite %in% BA_uncon ~ "unconjugated",
                                     Metabolite %in% BA_con ~ "conjugated",
                                     TRUE ~ NA),
                                  levels = c("unconjugated","conjugated")),
              OriConjugation= factor(case_when(Metabolite %in% BA_primary_uncon ~ "1°uncon",
                                        Metabolite %in% BA_primary_con ~ "1°con",
                                        Metabolite %in% BA_secondary_uncon ~ "2°uncon",
                                        Metabolite %in% BA_secondary_con ~ "2°con",
                                        TRUE ~ NA),
                                     levels = c("1°uncon","1°con","2°uncon","2°con")))




plot_stack <- function(data, group_var, value_var, fill_var, ylab, position_type = "fill",NAME,WIDTH,HEIGHT) {
  plot<-ggplot(data, aes(x = {{group_var}}, y = {{value_var}}, fill = {{fill_var}})) +
    geom_bar(stat = "identity", position = position_type) +
    labs(y = ylab, x = "") +
    theme_classic()
  
  # if(position_type == "stack"){
  #   plot<-plot+geom_text(aes(label = round({{value_var}}, 2)),position = "jitter", vjust = 1, size = 3)
  # }
  type<-case_when(position_type == "fill" ~ "rel",
                  position_type == "stack"~ "abs",
                  TRUE ~ NA)
  ggsave(plot, file =paste0(BApwd,"/BA_a/Ba_a_stacked_",NAME,"_",type,".png"), dpi=300,width=WIDTH, height=HEIGHT )
}
ALL <- a %>%
  group_by(Treatment,Metabolite) %>%
  summarise(sum = sum(value), .groups = "drop")
plot_stack(ALL, Treatment, sum, Metabolite, "Relative Abundance [%]", "fill",NAME="all",WIDTH= 10, HEIGHT=10)
plot_stack(ALL, Treatment, sum, Metabolite, "Concentration [µmol/L]", "stack",NAME="all",WIDTH= 10, HEIGHT=10)

Origin <- a %>%
  group_by(Treatment, Origin) %>%
  summarise(sumOrigin = sum(value), .groups = "drop")

plot_stack(Origin, Treatment, sumOrigin, Origin, "Relative Abundance [%]", "fill",NAME="origin",WIDTH= 5, HEIGHT=10)
plot_stack(Origin, Treatment, sumOrigin, Origin, "Concentration [µmol/L]", "stack",NAME="origin",WIDTH= 5, HEIGHT=10)

Conjugation <- a %>%
  group_by(Treatment, Conjugation) %>%
  summarise(sumConj = sum(value), .groups = "drop")

plot_stack(Conjugation, Treatment, sumConj, Conjugation, "Relative Abundance [%]", "fill",NAME="conj",WIDTH= 5, HEIGHT=10)
plot_stack(Conjugation, Treatment, sumConj, Conjugation, "Concentration [µmol/L]", "stack",NAME="conj",WIDTH= 5, HEIGHT=10)

OriginConjugation <- a %>%
  group_by(Treatment, OriConjugation) %>%
  summarise(sumOC = sum(value), .groups = "drop")

plot_stack(OriginConjugation, Treatment, sumOC, OriConjugation, "Relative Abundance [%]", "fill",NAME="oc",WIDTH= 5, HEIGHT=10)
plot_stack(OriginConjugation, Treatment, sumOC, OriConjugation, "Concentration [µmol/L]", "stack",NAME="oc",WIDTH= 5, HEIGHT=10)