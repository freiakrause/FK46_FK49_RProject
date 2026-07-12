gc()
library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)
source("FK49_Definitions.R")
data <- read.csv2(paste0(TEM_in_pwd, "/FK49_TEM_test.csv")) %>%
  mutate(Classification = gsub("Ring or Donut Mito","Donut",Classification),
    Area = as.numeric(gsub(",", ".", Area.µm.2)),
    Perimeter = as.numeric(gsub(",", ".", Perimeter.µm)),
    "Length.µm" = as.numeric(gsub(",", ".", Length.µm)))%>%
  select(Image, Classification, Area, Perimeter, Length.µm)

head(data)
df<-data%>%
  group_by(Image)%>% 
  mutate(Image_ID=cur_group_id())%>%ungroup()
# Generate Pedigree ----
# to allocate mesurments to specific mitochondria in Specific cells in specific images 
# Image numbering
df <- df %>%
  mutate(Image_ID = as.numeric(factor(Image)))

# Initialize columns
df$Cell_ID <- NA
df$Mito_ID <- NA

current_image <- 0
current_cell <- 0
current_mito <- 0

for(i in 1:nrow(df)) {
  
  # New image?
  if(df$Image[i] != current_image){
    current_image <- df$Image[i]
    current_cell <- 0
    current_mito <- 0
  }
  
  # New hepatocyte
  if(df$Classification[i] == "Hepatocyte"){
    current_cell <- current_cell + 1
    current_mito <- 0
    df$Cell_ID[i] <- current_cell
  }
  
  # Everything inside a hepatocyte gets same Cell_ID
  if(current_cell > 0){
    df$Cell_ID[i] <- current_cell
  }
  
  # New mitochondrion
  if(df$Classification[i] == "Donut" |
        df$Classification[i] == "Mitochondrium") {
    current_mito <- current_mito + 1
    df$Mito_ID[i] <- current_mito
  }
  
  # Width and Length belong to current mitochondrion
  if(df$Classification[i] %in% c("Width","Length")){
    df$Mito_ID[i] <- current_mito
  }
}

# Create pedigree IDs
df <- df %>%
  mutate(
    Pedigree_ID = ifelse(
      !is.na(Mito_ID),
      paste(Image_ID, Cell_ID, Mito_ID, sep="_"),
      ifelse(!is.na(Cell_ID),
             paste(Image_ID, Cell_ID, sep="_"),
             as.character(Image_ID))
    )
  )
## Get Analysis Order Id ----
#just to check if the order of analysis does influence/correclates with measuremtns . 
#if i learn overtime and changed

analysis_order <- df %>%
  distinct(Image) %>%
  mutate(Animal = as.numeric(str_extract(Image, "^[0-9]+"))) %>%
  arrange(Animal, Image) %>%       # image names increase within animal
  group_by(Animal) %>%
  mutate(Image_rank = row_number()) %>%
  ungroup() %>%
  arrange(Image_rank, Animal) %>%  # reproduces analysis order
  mutate(Analysis_ID = as.factor(row_number()))

df <- df %>%
  left_join(analysis_order %>%
            select(Image, Analysis_ID,Animal),
            by = "Image")%>%
  select(-"Image")


str(df)
# Look at  problematic mitochondria and improve them in qupath! No mito shold be problem!
wrong_mitos <- df %>%
  filter(Classification %in% c("Width", "Length")) %>%
  group_by(Image_ID, Pedigree_ID) %>%
  summarize(n = n(), .groups = "drop") %>%
  filter(n != 2)
## Check Lengths and Widths ----
# Length should be defined as longer than width. If this was accidentally wrongly assigned during measurements, we want to change it here
# Extract width and length for each mitochondrion
# Continue with only valid mitochondria
dims <- df %>%
  filter(Classification %in% c("Width", "Length")) %>%
  anti_join(wrong_mitos, by = c("Image_ID", "Pedigree_ID")) %>%
  select(Pedigree_ID, Classification, Length.µm) %>%
  pivot_wider(names_from = Classification,
              values_from = Length.µm) %>%
  mutate(Width_new = pmin(Width, Length),
         Length_new = pmax(Width, Length))

# Join corrected values back
df <- df %>%
  left_join(dims %>% select(Pedigree_ID, Width_new, Length_new),
    by = "Pedigree_ID")

# Replace values in Width rows
df$Length.µm[df$Classification == "Width"] <-
  df$Width_new[df$Classification == "Width"]

# Replace values in Length rows
df$Length.µm[df$Classification == "Length"] <-
  df$Length_new[df$Classification == "Length"]

# Remove temporary columns
df <- df %>%
  select(-Width_new, -Length_new)

rm(dims)
df<-df%>%mutate(Circularity = 4* pi *Area/(Perimeter*Perimeter))
# Reorder Columns to my liking
df<-df[,c("Animal","Image_ID","Cell_ID","Mito_ID","Pedigree_ID","Classification","Area","Perimeter","Length.µm","Circularity","Analysis_ID")]
df<-df%>%mutate(Image_ID = as.factor(Image_ID))
df_long <- df %>%
  pivot_longer( cols = c(Area, Perimeter, Length.µm, Circularity),
    names_to = "Variable",
    values_to = "Value",
    values_drop_na = TRUE )
#x<-df%>%group_by(Image_ID)%>%summarize(fraction= n(Mitochondrium)/(n(Mitochondrium)+n(Donut)))
n_mito_image <- df %>%
  filter(Classification %in% c("Mitochondrium", "Donut")) %>%
  group_by(Image_ID) %>%
  summarise(n_mito = n())

n_mito_cell <- df %>%
  filter(Classification %in% c("Mitochondrium", "Donut")) %>%
  group_by(Image_ID, Cell_ID) %>%
  summarise(n_mito = n())

mito_area_image <- df %>%
  filter(Classification %in% c("Mitochondrium", "Donut")) %>%
  group_by(Image_ID) %>%
  summarise(total_mito_area = sum(Area, na.rm = TRUE))

mito_area_cell <- df %>%
  filter(Classification %in% c("Mitochondrium", "Donut")) %>%
  group_by(Image_ID, Cell_ID) %>%
  summarise(total_mito_area = sum(Area, na.rm = TRUE))

hepatocyte_area <- df %>%
  filter(Classification == "Hepatocyte") %>%
  select(Image_ID, Cell_ID, cell_area = Area)

mito_density_cell <- n_mito_cell %>%
  left_join(hepatocyte_area,
            by = c("Image_ID", "Cell_ID")) %>%
  mutate(mito_per_area = n_mito / cell_area)

mito_fraction_cell <- mito_area_cell %>%
  left_join(hepatocyte_area,
            by = c("Image_ID", "Cell_ID")) %>%
  mutate(mito_area_fraction = total_mito_area / cell_area)



df %>%
  filter(Classification == "Length") %>%
  ggplot(aes(x = Analysis_ID, y = Length.µm)) +
  geom_jitter(width = 0.2, size = 1, alpha = 0.6) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  labs(x = "Analysis_ID", y = "Mitochondrial length (µm)")
  
  df %>%
    filter(Classification == "Width") %>%
    ggplot(aes(x = Image_ID, y = Length.µm)) +
    geom_jitter(width = 0.2, size = 1, alpha = 0.6) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  labs(x = "Image_ID", y = "Mitochondrial width (µm)")
  
  df %>%
    filter(Classification %in% c("Mitochondrium", "Donut")) %>%
    ggplot(aes(x = Image_ID, y = Area)) +
    geom_jitter(width = 0.2, size = 1, alpha = 0.6) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
    labs(x = "Image_ID", y = "Mitochondrial Area (µm2)")
  
  df %>%
    filter(Classification %in% c("Mitochondrium", "Donut")) %>%
    ggplot(aes(x = Image_ID, y = Perimeter)) +
    geom_jitter(width = 0.2, size = 1, alpha = 0.6) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
    labs(x = "Image_ID", y = "Mitochondrial Perimeter (µm)")
  
 df %>%
    filter(Classification %in% c("Mitochondrium", "Donut")) %>%
    ggplot(aes(x = Image_ID, y = Perimeter)) +
    geom_jitter(width = 0.2, size = 1, alpha = 0.6) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
    labs(x = "Image_ID", y = "Mitochondrial Perimeter (µm)")
  
  mitos_circle<-df %>%
    filter(Classification %in% c("Mitochondrium", "Donut")) %>%
    ggplot(aes(x = Image_ID, y = Circularity)) +
    geom_jitter(width = 0.2, size = 1, alpha = 0.6) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
    labs(x = "Image_ID", y = "Mitochondrial Circularity")
  ggsave(mitos_circle,path = TEM_out_pwd,file="Mito_Circularity.png",width = 8, height=4)
mito_density_cell %>%
    ggplot(aes(x = as.character(Image_ID), y = mito_per_area)) +
    geom_jitter(width = 0.1, size = 1, alpha = 0.6) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
    labs(x = "Image_ID", y = "Mitochondria/area/cell")

mito_fraction_cell %>%
  ggplot(aes(x = as.character(Image_ID), y = mito_area_fraction*100)) +
  geom_jitter(width = 0.2, size = 1, alpha = 0.6) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  labs(x = "Image_ID", y = "Mitochondria/area/cell")
 