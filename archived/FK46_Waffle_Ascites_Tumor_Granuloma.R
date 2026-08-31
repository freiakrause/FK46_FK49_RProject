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
# Waffle Plots of Tumor, Ascites and Granuloma Incidence ------------------
# setwd("C:/Users/b1084855/OneDrive - Universität Salzburg/Freigegebene Dokumente - AG_Tumorimmunologie/Data/Freia Krause/01_Experiments/FK49_CD-HFD_13wks/FK49_Analysis")
# load("01_RawData/FK49_Data_prepared.Rda")
do_waffle <- function(data, variable, sex_to_use=c("both"),params, order=NULL, path = "02_GeneratedData/Ascites_Tumor",Abb=NULL,nice_name) {
  
  if(sex_to_use == "both"){
    sex_to_filter = c("female","male")
  }else{sex_to_filter = sex_to_use}
  
  variable <- rlang::as_string(rlang::ensym(variable))
  yes_lab <- paste0(Abb)
  no_lab  <- paste0("no",Abb)
  
  
  d <- data %>%
    select(Animal,Sex, Treatment,all_of(variable)) %>%
    filter(Sex %in% sex_to_filter)%>%
    filter(complete.cases(.data[[variable]])) %>%
    group_by(Animal,Treatment, .data[[variable]]) %>%
    summarise(n = n()) %>%
    ungroup() %>%
    group_by(Treatment) %>%
    summarise(  !!yes_lab := sum(.data[[variable]] == 1),
                !!no_lab  := sum(.data[[variable]] == 0),    .groups = "drop") %>%
    pivot_longer(!Treatment, names_to = "Category", values_to = "Count")
  
  if(sex_to_use != "both"){
    tiles_per_block_ctrl = 5} else{tiles_per_block_ctrl= max(d$Count)}
  
  # Build waffle blocks
  ctrl_blocks <- d %>%
    filter(Treatment == "ctrl") %>%
    mutate(value = Count, fill_type = "ctrl")
  
  dummy_ctrl <- ctrl_blocks %>%
    mutate(value = pmax(tiles_per_block_ctrl - value, 0),   fill_type = "dummy")
  
  tam_blocks <- d %>%
    filter(Treatment == "TAM") %>%
    mutate(value = Count, fill_type = "TAM")
  
  if(sex_to_use != "both"){
    tiles_per_block_tam = 8
  }else{tiles_per_block_tam = 14}
  
  dummy_tam <- tam_blocks %>%
      mutate(value = pmax(tiles_per_block_tam - value, 0),   fill_type = "dummy")
  category_order <-order
  
  waffle_data <- bind_rows(ctrl_blocks, dummy_ctrl, tam_blocks,dummy_tam) %>%
    mutate(Category=factor(Category, levels = category_order))%>%
    uncount(value) %>%
    mutate(fill_type = factor(fill_type, levels = c("ctrl", "TAM", "dummy")))
  
  waffle_data$Category<-droplevels(waffle_data$Category)
  
  # Fisher test table
  fishers <- d %>%
    filter(Category %in% params) %>%
    select(Category, Treatment, Count) %>%
    pivot_wider(names_from = Treatment, values_from = Count)
  
  fish_matrix <- as.matrix(fishers[, c("ctrl", "TAM")])
  fisher_res <- fisher.test(fish_matrix)
  tag_part1 <-paste0("Fisher's Exact Test:\n")
  tag_part2 <- paste0("\nOdds Ratio = ", round(fisher_res$estimate,4))
  tag_part3 <- paste0("\n95% CI [",round(fisher_res$conf.int[1],2),",",round(fisher_res$conf.int[2],2),"]")
  tag_part4 <-if (fisher_res$p.value < 0.05){
    paste0("\n    p = ", signif(fisher_res$p.value, 3), "\nTreatment does influence occurence of ", nice_name, ".")
  } else  {paste0("\n       p = ", signif(fisher_res$p.value, 3), " \nTreatment does not influence occurence of ", nice_name, ".")} 
  
  # Waffle plot
  plot <- ggplot(waffle_data, aes(fill = fill_type, values = 1)) +
    geom_waffle(aes(color = fill_type), n_rows = 1, size = 0.5,
                height = 1.6, width = 0.9,  radius = unit(1, "pt"), alpha = 0.5) +
    scale_colour_manual(values = c("black", "black", "white"), guide = "none") +
    scale_fill_manual(values = c("ctrl" = "#4D4D4DBF", "TAM" = "#8B0000BF", "dummy" = "white"),
                      labels = c("ctrl"="Control","TAM"="TAM","dummy"=""),
                      name="Treatment") +
    coord_equal() +
    facet_grid(Category ~ ., switch = "y",   space = "free_y", labeller = label_wrap_gen(3)) +
    theme(plot.tag.position = c(0, -0.2),
      plot.tag = element_text(hjust = -0.5, size = 9),
      panel.spacing = unit(0.5, "lines"),
      strip.text = element_text(size = 13, face = "bold"),
      strip.text.y = element_text(angle = 90),
      legend.position = "top",
      strip.background = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      panel.grid = element_blank(),
      panel.background = element_rect(fill = "white", color = "black", linewidth = 1), 
      plot.margin = margin(t = 0,  # Top margin
                           r = 0,  # Right margin
                           b = 0,  # Bottom margin
                           l = 0))+ # Left margin  
    labs( tag = paste0(tag_part1,tag_part2,tag_part3,tag_part4))
                       
  # Save plot
 # filename1 <- paste0("FK46_", nice_name, "_",sex_to_use,".svg")
  filename2 <- paste0("FK46_", nice_name, "_",sex_to_use,".png")
  
  if(sex_to_use != "both"){
    w = 10
    h = 8
  } else{
      w = 10 
      h = 5.5}
  
 # ggsave(filename = filename1, plot = plot, path = path,  width = w, height = h, dpi = 300)
  ggsave(filename = filename2, plot = plot, path = path,  width = w, height = h, dpi = 300)
  
  return(list(fishers_table = fish_matrix,
    fisher_p = fisher_res$p.value,
    fisher =fisher_res,
    plot = plot))
}
do_waffle(data = data,variable = "Ascites.no.yes",  Abb  = "A", sex_to_use="both",params = c("A", "noA"), order = c( "noA","A"), nice_name= "Ascites ")
do_waffle(data = data,variable = "Tumor.no.yes",    Abb  = "T", sex_to_use="both", params = c("T", "noT"), order = c( "noT","T"), nice_name= "Tumor ")
do_waffle(data = data,variable = "Granuloma",Abb  = "G", sex_to_use="both", params = c("G", "noG"), order = c("noG","G" ), nice_name= "Granuloma")

do_waffle(data = data,variable = "Ascites.no.yes",  Abb  = "A", sex_to_use="male",params = c("A", "noA"), order = c( "noA","A"), nice_name= "Ascites ")
do_waffle(data = data,variable = "Tumor.no.yes",    Abb  = "T", sex_to_use="male", params = c("T", "noT"), order = c( "noT","T"), nice_name= "Tumor ")
do_waffle(data = data,variable = "Granuloma",Abb  = "G", sex_to_use="male", params = c("G", "noG"), order = c("noG","G" ), nice_name= "Granuloma")

do_waffle(data = data,variable = "Ascites.no.yes",  Abb  = "A", sex_to_use="female",params = c("A", "noA"), order = c( "noA","A"), nice_name= "Ascites ")
do_waffle(data = data,variable = "Tumor.no.yes",    Abb  = "T", sex_to_use="female", params = c("T", "noT"), order = c( "noT","T"), nice_name= "Tumor ")
do_waffle(data = data,variable = "Granuloma",Abb  = "G", sex_to_use="female", params = c("G", "noG"), order = c("noG","G" ), nice_name= "Granuloma")