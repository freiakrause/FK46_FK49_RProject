
rm(list = ls())
gc()
library(tidyverse)
library(ggstatsplot)
library(pheatmap)
library(rstatix)
library(grid)
library(scales)
source("FK49_Definitions.R")


# Read Raw Inputdata and general Data manipulation ------------------------------------------------------
ExpId="FK49"

if (ExpId=="FK49") {
  
  load(file.path(PATHS$MASH$FK49_input,"FK49_Data_prepared.Rda"))
  output_pwd <-file.path(PATHS$MASH$FK49_output)
  }  else if (ExpId == "FK46"){

    load(file.path(PATHS$MASH$FK46_input,"FK46_Data_prepared.Rda"))
    output_pwd <-file.path(PATHS$MASH$FK46_output)
      } else{
        print("Let me set a folder Path and define ExpId and Folderpath")}
data<-data %>% rename_with( ~ str_replace(.x, "^NASH_", "MASH_"), starts_with("NASH_"))# it's not called NASH anymore but MASH

# Waffel Plots -----

## Data Summary for Waffels -----
MASH_summary <- data %>%select(Animal,Sex, Treatment,starts_with("MASH_"),-MASH_N1,-MASH_N2)%>%
  filter(!is.na(MASH_S)) %>%
  mutate( across(starts_with("MASH_"),~ as.numeric(as.character(.)))) %>% 
  mutate(MASH_Category = case_when(
    is.na(MASH_S) ~ NA_character_,
    MASH_S == 0 ~ "No MASH",
    MASH_S != 0 & MASH_B == 0 & MASH_I == 0 ~ "Only Steatosis",
    MASH_S != 0 & MASH_B != 0 & MASH_I <= 0.5 ~ "Borderline MASH",
    MASH_S != 0 & MASH_B != 0 & MASH_I > 0.5 ~ "Definite MASH",
    TRUE ~ NA_character_)) %>%
   mutate(MASH_Category = factor(MASH_Category,levels = c("Definite MASH","Borderline MASH", "Only Steatosis","No MASH"))) %>%
  group_by(MASH_Category, Treatment) %>%
  summarise(n = n(),.groups = "drop") %>%
  complete(MASH_Category, Treatment = c("Ctrl", "TAM"), fill = list(n = 0))



# Fishers exact test ----
MASH_Fishers<-MASH_summary%>%pivot_wider(names_from = Treatment, values_from = n)
MASH_test <- as.matrix(MASH_Fishers[, -1]) #fishers needs to have first colmn remove 
test<-fisher.test(MASH_test)
test$p.value
tiles_per_block <-max(MASH_summary$n)  # Set fixed block size for Ctrl per MASH category to make plot align to the left for Ctrl and tam part.



# Generate Blocks for Tile plot for Categories -----
Ctrl_blocks <- MASH_summary %>%
  filter(Treatment == "Ctrl") %>%
  mutate(fill_type = "Ctrl",value = n)

dummy_Ctrl_blocks<- Ctrl_blocks %>%
  mutate( value = pmax(tiles_per_block - value, 0),  # dummy blocks add up onto Ctrl blocks to make Ctrl and tam in all categroies use left / right side
    fill_type = "dummy" )

tam_blocks <- MASH_summary %>%
  filter(Treatment == "TAM") %>%
  mutate( fill_type = "TAM",value = n)

dummy_tam_blocks<- tam_blocks %>%
   mutate( value = pmax(tiles_per_block - value, 0),  # dummy blocks add up onto TAM blocks to make Ctrl and tam in all categroies use left / right side
           fill_type = "dummy" )

# combine all blocks -----
waffle_data <- bind_rows(Ctrl_blocks, dummy_Ctrl_blocks,tam_blocks,dummy_tam_blocks) %>% 
  mutate(MASH_Category=factor(MASH_Category))%>%
  uncount(value) %>%# generates individual rows for all future blocks/tiles
  mutate(fill_type = factor(fill_type, levels = c("Ctrl", "TAM", "dummy")))%>%  
  group_by(Treatment, MASH_Category) %>%
  mutate(tile_id = row_number(),
         x = ((tile_id - 1) %% max(MASH_summary$n) ) + 1, #%% gibt ganz zahligen rest der division, 2%%9 2 da 2 nicht durch 7 teilbar.Bis 7 kommt immer die zahl raus, die getilet werden soll
         y = ceiling(tile_id / max(MASH_summary$n) )) %>% # if there are more than 9 it will get y biggger 1 so be in the next row
  ungroup() %>%
  mutate( x = if_else(Treatment == "TAM",x + max(MASH_summary$n) , x ))


## Plot the tile Plot -----
MASH_p <- ggplot(waffle_data, aes(fill = fill_type, values = 1)) +
  geom_tile(aes(x=x,y=y),color = "white", linewidth = 0.8, width = 1, height =1) +
  scale_colour_manual(values = c("black", "black" ,"white"),guide="none") +
  scale_fill_manual(values = c(alpha(Treatment_colors[c("Ctrl","TAM")],0.4), "dummy" = "white"),
                    labels = c("Ctrl"="Ctrl","TAM"="TAM","dummy"=""),
                     name= "Treatment")+
  scale_x_continuous(breaks = c(4.5, 13.5),labels = c("Ctrl", "TAM")) +
 # coord_equal(expand = TRUE) +
  facet_grid(MASH_Category ~ ., switch = "y", space = "free_y",  labeller = label_wrap_gen(8)) +
  theme(panel.spacing = unit(0.5 ,"lines"),
        plot.tag.position = c(0, -0.2),
        plot.tag = element_text(hjust = -0.2, size = 9),
        strip.text = element_text(size = 9, face = "bold"),
        strip.text.y = element_text(angle = 90),
        legend.position = "top",
        strip.background = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.text.x = element_text(size = 13, face = "bold"),
        axis.ticks = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.background = element_blank(),
        panel.background = element_rect(fill = "white", color = "black", linewidth = 1) ) #+
#  labs(tag = paste0("Fisher's Exact Test: p ", round(test$p.value, 5),"\nTreatment does influence MASH severity"))

ggsave(filename = paste0(ExpId,"_MASHScore_both.png"), plot = MASH_p, path = output_pwd, 
       width = 8, height =6,dpi = 300,bg="white")


rm(Ctrl_blocks,dummy_tam_blocks,tam_blocks, waffle_data,MASH_p,tiles_per_block)
# Do Barplot as different viszualisation -----
df_long <- MASH_Fishers %>%
  pivot_longer( cols = c(Ctrl, TAM), names_to = "Treatment",   values_to = "count" ) %>%
  uncount(count)%>%
  mutate(MASH_Category=factor(MASH_Category,levels=c("No MASH","Only Steatosis","Borderline MASH","Definite MASH")) )              # replicate rows according to count

## Barplot ----
p1<-ggbarstats(df_long, MASH_Category, Treatment,results.subtitle = FALSE, subtitle = paste0("Fisher's exact test", ", p-value = ",round(test$p.value, 5)))+
  scale_fill_manual(values = c("No MASH" ="#A8E6B1", "Only Steatosis" =  "#FFDAB9", "Borderline MASH" = "#FDBA74","Definite MASH" = "#FF9999"))
ggsave(filename = "FK49_MASHScore_bar_both.png", plot = p1, path =output_pwd, width = 5, height =10 ,dpi = 300)



# Preapre Dataframe for BoxPlots and HeatMap of MASH Scores -----
d<-data%>%filter(BATCH %in% c(1, 2)) %>%
  filter(complete.cases(MASH_S2))%>%
  select(Animal,Sex,Treatment,MASH_S,MASH_B,MASH_I,MASH_SAF)%>%
  mutate(across(starts_with("MASH_"), as.character))%>%
  mutate(across(starts_with("MASH_"), as.numeric))%>%
  arrange(Treatment,desc(MASH_SAF),desc(MASH_I),desc(MASH_B),desc(MASH_S))

# Boxplots of indiviual MASH Parameters -----
## Prepare object for loops for boxplots and p value -----
wilcox_results <- list()
plots <- list()
pvalues <- c()   # numeric vector
method <- c() 
data.name <- c() 
alternative <- c() 
vars <- c("MASH_S","MASH_B","MASH_I","MASH_SAF")
Summary_W_stats<- c()
# Paired Wilcoxon Rank Sum test on MASH Scores -----
for (i in vars) {
  Summary_W_stats<-d%>%group_by(Treatment)%>%get_summary_stats()
  a <- wilcox.test(as.formula(paste0(i,"~ Treatment")),  data = d,exact = FALSE,  correct = FALSE,   
                   conf.int = FALSE) # ist mann whitney u
  #ist Mann-Whitney-U test da piared=FALSE
  wilcox_results[[i]] <- a
  pvalues[i] <- a$p.value
  method[i] <- a$method
  data.name[i] <- a$data.name
  alternative[i] <- a$alternative
  
}

padjusted <- p.adjust(pvalues, method = "fdr")
wilcox_table <- data.frame(
  Variable = vars,
  p_value = pvalues,
  p_adjusted = padjusted,
  data.name =data.name,
  alternative = alternative,
  method= method)
rownames( wilcox_table ) <- NULL
MASH_Category_stats <- bind_rows( MASH_Fishers%>%mutate(method= "Summary_Category"), as.data.frame(test[1:4]))%>%
  select(method,colnames(MASH_Fishers),colnames(as.data.frame(test[1:4])))

MASH_Score_stats <- bind_rows( Summary_W_stats%>%mutate(method= "Summary_Scores"),wilcox_table)%>%
  select(method,colnames(Summary_W_stats),colnames(wilcox_table))

write.csv2(MASH_Category_stats,paste0(output_pwd,"/Statistics/",ExpId,"_MASH_Category_Stats.csv"))
write.csv2(MASH_Score_stats,paste0(output_pwd,"/Statistics/",ExpId,"_MASH_Score_Stats.csv"))


# Boxplots of MASH Scores -----
for (i in vars) {
   plots[[i]] <- ggplot(d, aes(x = Treatment, y = .data[[i]], fill = Treatment)) +
    geom_boxplot() +
    scale_fill_manual(values = c(alpha(Treatment_colors[c("Ctrl","TAM")],0.4)) )+
    annotate("text",  x = 1.5,y = max(d[[i]], na.rm = TRUE),  size = 5,
              label = paste0( "Mann-Whitney U:\n","adj p ",signif(padjusted[i], 3))) +
    ggtitle(i)
    ggsave( filename = paste0(ExpId,"_MASHScore_box_", i, ".png"),plot = plots[[i]],  path = output_pwd,
    width = 5,  height = 10,  dpi = 300 )
}

# Do HeatMap of MASH Scores -----
## Prepare Matrix and Data for HeatMap -----

  mat <- as.matrix(d[, c("MASH_I","MASH_B","MASH_S","MASH_SAF")])#,
  rownames(mat) <- d$Animal
  ann <- data.frame(Treatment = d$Treatment)
  rownames(ann) <- d$Animal
  mat<-t(mat) #rotates transposes the matrix
  row_annot <- data.frame(adj.p = round(padjusted[c("MASH_I","MASH_B","MASH_S",  "MASH_SAF")],5))
  rownames(row_annot) <- c("MASH_I", "MASH_B","MASH_S", "MASH_SAF")
  
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
    annotation_colors = list(Treatment = Treatment_colors[c("Ctrl","TAM")]),
    main = "MASH Scoring Heatmap iAL after 13wks CD-HFD",
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
  
ggsave(filename = paste0(ExpId,"_MASHScore_HeatMap.png"),
       plot = h, path = output_pwd, width = 10, height =5 ,dpi = 300,bg="white")

