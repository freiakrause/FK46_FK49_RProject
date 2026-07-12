
rm(list = ls())
gc()
home <- normalizePath("~") # bc Windows does not start inuser dir but in user/documents dir
parent <- dirname(home)
setwd(paste0(parent,"/OneDrive - Universität Salzburg/AG_Tumorimmunologie - Dokumente/Data/Freia Krause/01_Experiments/FK46_FK49_RProject"))

library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)
library(patchwork)
library(waffle)
library(pheatmap)
library(lubridate)
library(tidyverse)
library(ggstatsplot)
library(grid)
#library(superb)
#library(rstatix)

# Read Raw Inputdata and general Data manipulation ------------------------------------------------------
ExpId="FK46"

if (ExpId=="FK49") {
  setwd(paste0(parent,"/OneDrive - Universität Salzburg/AG_Tumorimmunologie - Dokumente/Data/Freia Krause/01_Experiments/FK49_CD-HFD_13wks/FK49_Analysis"))
  }  else if (ExpId == "FK46"){
      setwd(paste0(parent,"/OneDrive - Universität Salzburg/AG_Tumorimmunologie - Dokumente/Data/Freia Krause/01_Experiments/FK46_iALmice_high Fat diet 52 weeks 7d after injection/Analysis"))
      } else{
        print("Let me set a folder Path and define ExpId and Folderpath")}
load(paste0("01_RawData/",ExpId,"_Data_prepared.Rda"))
##### Try to do NASH display

# Waffel Plots -----

tiles_per_block <-10  # Set fixed block size for ctrl per NASH category

## Data Summary for Waffels -----
nash_summary <- data %>%
  #filter(Sex == "male")%>%
  #filter(Ascites.no.yes==1)%>%
  filter(BATCH %in% c(1, 2)) %>%
  filter(complete.cases(NASH_S2)) %>%
  mutate(NASH_Category = case_when(
    is.na(NASH_S2) ~ NA_character_,
    NASH_S == "0" ~ "No NASH",
    NASH_S != "0" & NASH_B == "0" ~ "Only Steatosis",
    NASH_S != "0" & NASH_B != "0" & NASH_I <= "0.5" ~ "Borderline NASH",
    NASH_S != "0" & NASH_B != "0" & NASH_I > "0" ~ "Definite NASH",
    TRUE ~ NA_character_)) %>%
  #filter(!is.na(NASH_Category)) %>%
  mutate(NASH_Category = factor(NASH_Category,levels = c("Definite NASH","Borderline NASH", "Only Steatosis","No NASH"))) %>%
  group_by(NASH_Category, Treatment) %>%
  summarise(n = n(),.groups = "drop") %>%
  complete(NASH_Category,
           Treatment = c("ctrl", "TAM"),
           fill = list(n = 0))

tiles_per_block <-max(nash_summary$n)  # Set fixed block size for ctrl per NASH category to make plot align to the left for ctrl and tam part.
## Generate Blocks for Waffle -----
ctrl_blocks <- nash_summary %>%
  filter(Treatment == "ctrl") %>%
  mutate(fill_type = "ctrl",value = n)

dummy_Ctrl_blocks<- ctrl_blocks %>%
  mutate( value = pmax(tiles_per_block - value, 0),  # dummy blocks add up onto ctrl blocks to make ctrl and tam in all categroies use left / right side
    fill_type = "dummy" )

tam_blocks <- nash_summary %>%
  filter(Treatment == "TAM") %>%
  mutate( fill_type = "TAM",value = n)
# dummy_TAM_blocks<- tam_blocks %>%
#   mutate( value = pmax(tiles_per_block - value, 0),  # dummy blocks add up onto ctrl blocks to make ctrl and tam in all categroies use left / right side
#           fill_type = "dummyT" )
#Combine all blocks
waffle_data <- bind_rows(ctrl_blocks, dummy_Ctrl_blocks, tam_blocks) %>% #,dummy_TAM_blocks
  select(NASH_Category, fill_type, value) %>%
  filter(value > 0) %>%
  uncount(value)%>%mutate(fill_type=factor(fill_type, levels=c("ctrl","TAM","dummy")))  # one row per tile #,"dummyT"

## Fishers exact test ----
NASH_Fishers<-nash_summary%>%pivot_wider(names_from = Treatment, values_from = n)
NASH_test <- as.matrix(NASH_Fishers[, -1]) #fishers needs to have first colmn remove 
test<-fisher.test(NASH_test)
test$p.value

## Waffle Plot -----
nash_p <- ggplot(waffle_data, aes(fill = fill_type, values = 1)) +
  geom_waffle(aes(colour = fill_type), size = 0.5, n_rows = 1, height = 1.6, width = 0.9,
              radius = unit(1, "pt"), alpha = 0.5) +
  scale_colour_manual(values = c("black", "black" ,"white","white"),guide="none") +
  scale_fill_manual(values = c("ctrl" = "#4D4D4DBF", "TAM" =  "#8B0000BF", "dummy" = "white"), #,"dummyT" = "white"
                    name = "Treatment", labels = c("ctrl"="Control","TAM"="TAM","dummy"="")) + #,"dummyT"=""
  coord_equal(expand = TRUE) +
  facet_grid(NASH_Category ~ ., switch = "y", space = "free_y",  labeller = label_wrap_gen(3)) +
  theme(panel.spacing = unit(0.5 ,"lines"),
    plot.tag.position = c(0, -0.2),
    plot.tag = element_text(hjust = -0.2, size = 10),
    strip.text = element_text(size = 13, face = "bold"),
    strip.text.y = element_text(angle = 90),
    legend.position = "top",
    strip.background = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.background = element_blank(),
    panel.background = element_rect(fill = "white", color = "black", linewidth = 1) ) +
  labs(tag = paste0("Fisher's Exact Test: p ", round(test$p.value, 5),"\nTreatment does influence NASH severity"))

ggsave(filename = paste0(ExpId,"_NASHScore_both.png"), plot = nash_p, path = "02_GeneratedData/NASH_Score", width = 10, height =12 ,dpi = 300)

#Show
nash_p

rm(ctrl_blocks,dummy_blocks,tam_blocks, waffle_data,nash_p,tiles_per_block)

df_long <- NASH_Fishers %>%
  pivot_longer( cols = c(ctrl, TAM), names_to = "Treatment",   values_to = "count" ) %>%
  uncount(count)%>%
  mutate(NASH_Category=factor(NASH_Category,levels=c("No NASH","Only Steatosis","Borderline NASH","Definite NASH")) )              # replicate rows according to count

## Barplot ----
p1<-ggbarstats( df_long, NASH_Category, Treatment,results.subtitle = FALSE, subtitle = paste0("Fisher's exact test", ", p-value = ",round(test$p.value, 5)))+
  scale_fill_manual(values = c("No NASH" ="#A8E6B1", "Only Steatosis" =  "#FFDAB9", "Borderline NASH" = "#FDBA74","Definite NASH" = "#FF9999"))

ggsave(filename = "FK49_NASHScore_bar_both.png", plot = p1, path = "02_GeneratedData/NASH_Score", width = 5, height =10 ,dpi = 300)

# Preapre Dataframe for BoxPlots and HeatMap -----
d<-data%>%filter(BATCH %in% c(1, 2)) %>%
  filter(complete.cases(NASH_S2))%>%
  select(Animal,Treatment,NASH_S,NASH_B,NASH_I,NASH_SAF)%>%
  mutate(across(starts_with("NASH_"), as.character))%>%
  mutate(across(starts_with("NASH_"), as.numeric))%>%
  arrange(Treatment,desc(NASH_SAF),desc(NASH_I),desc(NASH_B),desc(NASH_S))

# Boxplots of indiviual NASH Parameters -----
## Prepare object for loops for boxplots and p value -----
wilcox_results <- list()
plots <- list()
pvalues <- c()   # numeric vector
vars <- c("NASH_S","NASH_B","NASH_I","NASH_SAF")

for (i in vars) {
  a <- wilcox.test(  as.formula(paste(i,"~Treatment")),  data = d,  exact = FALSE,  correct = FALSE,   conf.int = FALSE) # ist mann whitney u
  wilcox_results[[i]] <- a
  pvalues[i] <- a$p.value
}

padjusted <- p.adjust(pvalues, method = "BH")

for (i in vars) {
   plots[[i]] <- ggplot(d, aes(x = Treatment, y = .data[[i]], fill = Treatment)) +
    geom_boxplot() +
    scale_fill_manual(values = c("#4D4D4DBF","#8B0000BF")) +
    annotate("text",  x = 1.5,y = max(d[[i]], na.rm = TRUE),  size = 5,
              label = paste0( "Mann-Whitney U:\n","adj p ",signif(padjusted[i], 3))) +
    ggtitle(i)
  
   ggsave( filename = paste0(ExpId,"_NASHScore_box_", i, ".png"),plot = plots[[i]],  path = "02_GeneratedData/NASH_Score",
    width = 5,  height = 10,  dpi = 300 )
}

# Do HeatMap of NASH Scores -----
## Preapre Matrix and Data for HeatMap -----

  mat <- as.matrix(d[, c("NASH_I","NASH_B","NASH_S","NASH_SAF")])#,
  rownames(mat) <- d$Animal
  ann <- data.frame(Treatment = d$Treatment)
  rownames(ann) <- d$Animal
  mat<-t(mat) #rotates transposes the matrix
  row_annot <- data.frame(adj.p = round(padjusted[c("NASH_I","NASH_B","NASH_S",  "NASH_SAF")],5))
  rownames(row_annot) <- c("NASH_I", "NASH_B","NASH_S", "NASH_SAF")
  
## Heatmap  -----
  h <- pheatmap(mat,
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    fontsize        = 10,
    fontsize_main   = 7,
    fontsize_row    = 9,
    fontsize_col    = 9,
    fontsize_number = 8,
    annotation_col  = ann,
    annotation_row  = row_annot,       
    color = colorRampPalette(c("white", "orange", "red"))(50),
    annotation_colors = list(Treatment = c("ctrl" = "#4D4D4DBF", "TAM" = "#8B0000BF")),
    main = "NASH Scoring Heatmap iAL after 13wks CD-HFD",
    labels_col = "Animals",
    labels_row = c("Inflammation","Ballooning","Steatosis",  "Total Score"),
    border_color = "black",
    cellwidth = 10,
    cellheight = 10,
    angle_col = 0,
    gaps_row = 3,
    display_numbers = FALSE,
    number_format = "%.3f",
    legend_breaks = c(0,2,4),
    legend_labels = c("0","2 ", "4"))
  
ggsave(filename = paste0(ExpId,"_NASHScore_HeatMap.png"), plot = h, path = "02_GeneratedData/NASH_Score", width = 10, height =5 ,dpi = 300)
#Show plot
h
