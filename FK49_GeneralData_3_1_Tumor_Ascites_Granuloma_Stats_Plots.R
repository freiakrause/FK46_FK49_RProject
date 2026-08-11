rm(list=ls())
gc()
library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)
library(rlang)
source("FK49_Definitions.R")

# hier morgen weitermachen, code schick und kommentare. und statistik in methods. dann fine
output_pwd <-file.path(PATHS$TAG$output)
load(file.path(PATHS$TAG$input,"FK49_Data_prepared.Rda"))
exigo <- c("ALB", "TP", "GLOB","A.G", "TB", "GGT", "AST", "ALT", "ALP", "AMY","Crea","UA","BUN","GLU","TC","TG")

do_waffle <- function(data, variable, sex_to_use=c("both"),params, order=NULL, path = output_pwd,Abb=NULL,nice_name) {
  
  if(sex_to_use == "both"){
    sex_to_filter = c("female","male")
  }else{sex_to_filter = sex_to_use}
  
  variable <- rlang::as_string(rlang::ensym(variable))
  yes_lab <- paste0(Abb)
  no_lab  <- paste0("no",Abb)
  
  d <- data %>%
    select(Animal,Sex, Treatment,all_of(variable)) %>%
    filter(Sex %in% sex_to_filter)%>%
    filter(complete.cases(.)) %>%
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
    filter(Treatment == "Ctrl") %>%
    mutate(value = Count, fill_type = "Ctrl")
  
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
    mutate(fill_type = factor(fill_type, levels = c("Ctrl", "TAM", "dummy")))
  
  waffle_data$Category<-droplevels(waffle_data$Category)
  
  # Create tile positions without waffle package
  
  # waffle_data <- waffle_data %>%
  #   group_by(Treatment, Category) %>%
  #   mutate(tile_id = row_number(),
  #          x = ((tile_id - 1) %% if_else(Treatment == "Ctrl",
  #                                        tiles_per_block_ctrl,
  #                                        tiles_per_block_tam)) + 1,
  #          y = ceiling(tile_id / if_else(Treatment == "Ctrl",
  #                                        tiles_per_block_ctrl,
  #                                        tiles_per_block_tam))) %>%
  #   ungroup() %>%
  #   mutate(
  #     x = if_else(
  #       Treatment == "TAM",
  #       x + tiles_per_block_ctrl,
  #       x
  #     )
  #   )
 waffle_data <- waffle_data %>%
    group_by(Treatment, Category) %>%
    mutate(tile_id = row_number(),
           x = ((tile_id - 1) %% 7) + 1,
           y = ceiling(tile_id / 7)) %>%
    ungroup() %>%
    mutate(
      x = if_else(
        Treatment == "TAM",
        x + 7,
        x
      )
    )
  # Fisher test table
 fishers <- d %>%
   filter(Category %in% params) %>%
   select(Category, Treatment, Count) %>%
   pivot_wider(names_from = Treatment, values_from = Count)
 
 fish_matrix <- as.matrix(fishers[, c("Ctrl", "TAM")])
 fisher_res <- fisher.test(fish_matrix)
 interpretation <- if (fisher_res$p.value < 0.05) {
                   paste0("Treatment has an effect on ", nice_name, " outcome.")
                   } else {
                  paste0("Treatment has no evidence of an effect on ", nice_name, " outcome.")
                       }
 
 fisher_results <- tibble(
   Parameter = c(
     "Variable",
     "Sex",
     "Category 1",
     "Category 2",
     "Ctrl Category 1",
     "Ctrl Category 2",
     "TAM Category 1",
     "TAM Category 2",
     "Odds Ratio",
     "95% CI lower",
     "95% CI upper",
     "p-value",
     "Alternative",
     "Method",
     "Interpretation"
   ),
   Value = c(
     nice_name,
     sex_to_use,
     params[1],
     params[2],
     fish_matrix[1, "Ctrl"],
     fish_matrix[2, "Ctrl"],
     fish_matrix[1, "TAM"],
     fish_matrix[2, "TAM"],
     unname(fisher_res$estimate),
     fisher_res$conf.int[1],
     fisher_res$conf.int[2],
     fisher_res$p.value,
     fisher_res$alternative,
     fisher_res$method,
     interpretation
   )
 )
 
 write.csv2(
   fisher_results,
   file.path(path, paste0("FK46_", nice_name, "_", sex_to_use, "_Fisher.csv")),
   row.names = FALSE
 )
 
 write.csv2(
   fisher_results,
   file.path(path, paste0("FK46_", nice_name, "_", sex_to_use, "_Fisher.csv")),
   row.names = FALSE
 )
  # Waffle plot
  
  plot <- ggplot(waffle_data, aes(x = x, y = y, fill = fill_type)) +
    geom_tile(color = "white", linewidth = 0.8, width = 1, height = 1) +
    scale_fill_manual(values = c(Treatment_colors[c("Ctrl","TAM")], "dummy" = "white"),
                      labels = c("Ctrl"="Control","TAM"="TAM","dummy"=""),
                      name="Treatment") +
    scale_x_continuous(breaks = c(3.5, 10.5),labels = c("Ctrl", "TAM")) +
    scale_y_continuous(breaks = NULL, labels = NULL ) +
    coord_equal() +
    facet_grid(Category ~ ., switch = "y", space = "free_y", labeller = label_wrap_gen(3)) +
    theme(plot.tag.position = c(0, -0.4),
          plot.tag = element_text(hjust = -0.5, size = 9),
          panel.spacing = unit(0.5, "lines"),
          strip.text = element_text(size = 13, face = "bold"),
          strip.text.y = element_text(angle = 90),
          legend.position = "top",
          strip.background = element_blank(),
          axis.text.y = element_blank(),
          axis.text.x = element_text(size = 13, face = "bold"),
          axis.ticks = element_blank(),
          axis.title.x =element_blank(),
          axis.title.y = element_blank(),
          panel.grid = element_blank(),
          panel.background = element_rect(fill = "white", color = "black", linewidth = .8),
          plot.margin = margin(t = 0,  # Top margin
                               r = 0,  # Right margin
                               b = 0,  # Bottom margin
                               l = 0)) # Left margin
 
  
  # Save plot
  
  # filename1 <- paste0("FK46_", nice_name, "_",sex_to_use,".svg")
  
  filename2 <- paste0("FK46_", nice_name, "_",sex_to_use,".png")
  
  if(sex_to_use != "both"){
    w = 5
    h = 2
  } else{
    w = 6
    h = 3}
  
  # ggsave(filename = filename1, plot = plot, path = path,  width = w, height = h, dpi = 300)
  
  ggsave(filename = filename2, plot = plot, path = path,  width = w, height = h, dpi = 300)
  
  return(list(fishers_table = fish_matrix,
              fisher_p = fisher_res$p.value,
              fisher =fisher_res,
              plot = plot))
}



do_waffle(data = data,variable = "Tumor.no.yes",    
           Abb  = "T", sex_to_use="both", 
           params = c("T", "noT"), order = c( "T","noT"), nice_name= "Tumor ")



do_waffle(data = data,variable = "Ascites.no.yes",    
          Abb  = "A", sex_to_use="both", 
          params = c("A", "noA"), order = c( "A","noA"), nice_name= "Ascites ")

do_waffle(data = data,variable = "Granuloma.no.yes",    
          Abb  = "G", sex_to_use="both", 
          params = c("G", "noG"), order = c( "noG","G"), nice_name= "Granuloma ")

