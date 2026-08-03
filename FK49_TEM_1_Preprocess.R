rm(list=ls())
gc()
library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)
library(gridExtra)
library(grid)
source("FK49_Definitions.R")
###### Muss noch treatment, Animal real und Sex überall mitziehen 03.08.26
data <- read.csv2(file.path(PATHS$TEM$input, "/FK49_TEM_test.csv"),sep=",") %>%
  filter(!Classification== "Granules",
          !  Classification == "FUN",
           !Classification == "Unsure",
            !Classification == "Nucleus") %>%
  mutate(Classification = gsub("Ring or Donut Mito","Donut",Classification),
         Classification = gsub("LostCristMito","LostCristae",Classification),
         Classification = ifelse(Classification == "",Name,Classification),
         Classification = gsub("imae_descriptor","image_descriptor",Classification),
         Classification = gsub("image_dscriptor","image_descriptor",Classification),
         Classification = gsub("Lymphocyte\\?","nonHepatocyte",Classification),
      Area = as.numeric(gsub(",", ".", Area.µm.2)),
    Perimeter = as.numeric(gsub(",", ".", Perimeter.µm)),
    "Length.µm" = as.numeric(gsub(",", ".", Length.µm)),
    Descriptor = Num.points)%>%
  select(Image, Classification, Area, Perimeter, Length.µm,Descriptor)

animal_metadata<-read.csv2(file.path(PATHS$TEM$input, "/TEM_Animal_Metadata.csv"),sep=";")

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
  if(df$Classification[i] == "Hepatocyte" |
     df$Classification[i] == "nonHepatocyte" ){
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
     df$Classification[i] == "LostCristae" |
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
df <- df%>%mutate(Cell_ID = ifelse(Classification== "image_descriptor", NA,Cell_ID))

## Get Analysis Order Id ----
#just to check if the order of analysis does influence/correlates with measurements . 
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
# Look at  problematic mitochondria and improve them in qupath! No mito should be problem!
# problematic ==length or width measurement missing
#QC check input data for completeness and correctness ----
  ## Do all Mitos have 2  measurments -----
   #in pedigree below them (lenght and width)
wrong_mitos <- df %>%
  filter(Classification %in% c("Width", "Length")) %>%
  group_by(Image_ID, Cell_ID,Mito_ID,Pedigree_ID) %>%
  summarize(n = n(), .groups = "drop") %>%
  filter(n != 2)
### Checked and corrected this QuPath on 03.08.26 and is fine now.
 ## Do I have 10 animals and 10 images per Animal -----
df <- df %>%
  left_join(animal_metadata %>%
      select( Animal_blind, Animal_real = Animal,Treatment,Sex),
    by = c("Animal" = "Animal_blind") )
total_numbers <- df %>%
  summarise(Total_Animals   = n_distinct(Animal),
    Total_Images    = n_distinct(Image_ID),
    Total_TAM       = n_distinct(Animal[Treatment == "TAM"]),
    Total_EtOH      = n_distinct(Animal[Treatment == "EtOH"]),
    Total_Female    = n_distinct(Animal[Sex == "Female"]),
    Total_Male      = n_distinct(Animal[Sex == "Male"]),
    Total_Cells     = sum(Classification %in% c("Hepatocyte", "nonHepatocyte")),
    Hepatocytes     = sum(Classification == "Hepatocyte"),
    nonHepatocytes  = sum(Classification == "nonHepatocyte"),
    Total_AllMitos  = sum(Classification %in% c("Mitochondrium","Donut","LostCristae")),
    StandardMitos   = sum(Classification == "Mitochondrium"),
    Donuts          = sum(Classification == "Donut"),
    LostCristae     = sum(Classification == "LostCristae"),
    Mito_TAM        = sum(Classification %in% c("Mitochondrium","Donut","LostCristae") &
                            Treatment == "TAM"),
    Mito_EtOH       = sum(Classification %in% c("Mitochondrium","Donut","LostCristae") &
                            Treatment == "EtOH"),
    Mito_Female     = sum(Classification %in% c("Mitochondrium","Donut","LostCristae") &
                            Sex == "Female"),
    Mito_Male       = sum(Classification %in% c("Mitochondrium","Donut","LostCristae") &
                            Sex == "Male")) %>%
  pivot_longer( cols = everything(),names_to = "Metric",values_to = "Count" )

write.csv(total_numbers,file=file.path(PATHS$TEM$output, "Background_TotalNumbers.csv"))
png(file.path(PATHS$TEM$output,"Background_TotalNumbers.png"),  width = 800,  height = 1200,  res = 300)
grid.table(total_numbers)
dev.off()

## AAnimal summary
animal_summary <- df %>%
  group_by(Animal, Treatment, Sex) %>%
  summarise(# Sampling depth
    Images = n_distinct(Image_ID),
    Total_Cells = sum(Classification %in% c("Hepatocyte", "nonHepatocyte")),
    Hepatocytes = sum(Classification == "Hepatocyte"),
    nonHepatocytes = sum(Classification == "nonHepatocyte"),
    # Cell area
    Hep_Area = sum(   Area[Classification == "Hepatocyte"],   na.rm = TRUE),
    Mean_Hep_Area = mean(  Area[Classification == "Hepatocyte"],   na.rm = TRUE),
    # Mito numbers
    AllMitos = sum(  Classification %in% c("Mitochondrium","Donut","LostCristae")),
    Mito_per_HepA = round(AllMitos / Hep_Area,4  ),
    Mean_Mitos_per_Image = AllMitos / Images,
    Mean_Mitos_per_Hepatocyte = AllMitos / Hepatocytes,
    # Mito morphology
    StandardMitos = sum(Classification == "Mitochondrium"),
    Donuts = sum(Classification == "Donut"),
    LostCristae = sum(Classification == "LostCristae"),
    Freq_StandardMito = StandardMitos / AllMitos,
    Freq_Donut = Donuts / AllMitos,
    Freq_LostCristae = LostCristae / AllMitos,
    Freq_Damaged = (Donuts + LostCristae) / AllMitos,
    .groups = "drop"
  )
write.csv(animal_summary,file=file.path(PATHS$TEM$output, "Background_AnimalNumbers.csv"))
png(file.path(PATHS$TEM$output,"Background_AnimalNumbers.png"),  width = 3000,  height = 1200,  res = 300)
grid.table(animal_summary)
dev.off()

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
df<-df[,c("Animal","Image_ID","Cell_ID","Mito_ID","Pedigree_ID","Classification","Area","Perimeter","Length.µm","Circularity","Descriptor","Analysis_ID")]
df<-df%>%mutate(Image_ID = as.factor(Image_ID))
df_long <- df %>%
  pivot_longer( cols = c(Area, Perimeter, Length.µm, Circularity),
    names_to = "Variable",
    values_to = "Value",
    values_drop_na = TRUE )
str(df)
head(df)
unique(df$Classification)




##  Dataframe für Mitochondrien erstellen, säubern und speichern
Mito_Dataframe <- df %>%
  filter(!Classification %in% c("Hepatocyte", "nonHepatocyte", "image_descriptor")) %>%
  mutate( Shape = if_else(Classification %in% c("Donut", "Mitochondrium", "LostCristae"),
                          Classification, NA_character_  )) %>%
  group_by(Pedigree_ID) %>%
  fill(Shape, Area, Perimeter, Circularity, .direction = "down") %>%
  pivot_wider(id_cols = c(Animal, Image_ID, Cell_ID, Mito_ID, Pedigree_ID,Shape, Area, Perimeter, Circularity, Analysis_ID),
    names_from = Classification,values_from = Length.µm) %>%
  select(-Donut,-Mitochondrium,-LostCristae)%>%
  ungroup()                           
saveRDS(Mito_Dataframe,file.path(PATHS$TEM$output,"Mito_Dataframe.rds"))
write.csv(Mito_Dataframe,file.path(PATHS$TEM$output,"Mito_Dataframe.csv"))
## Dataframe für Zellen erstellen säubern und speichern

Cell_Metadata <- df %>%
  filter(Classification %in% c("Hepatocyte", "nonHepatocyte")) %>%
  select( Animal,  Image_ID,Cell_ID,Cell_Type = Classification, Cell_Area = Area) %>%
  distinct()

Cell_Dataframe <- Mito_Dataframe %>%
  left_join(Cell_Metadata,by = c("Animal", "Image_ID", "Cell_ID") )%>%
  group_by(Animal, Image_ID,Cell_ID) %>%
  summarise(n_Mito = n(),
            Cell_Area= mean(Cell_Area),
            Cell_Type = first(Cell_Type),
    Mito_pro_Area = n_Mito/Cell_Area,
    Mean_Mito_Area = mean(Area, na.rm = TRUE),
    SD_Mito_Area = sd(Area, na.rm = TRUE),
    Mean_Perimeter = mean(Perimeter, na.rm = TRUE),
    SD_Perimeter = sd(Perimeter, na.rm = TRUE),
    Mean_Length = mean(Length, na.rm = TRUE),
    SD_Length = sd(Length, na.rm = TRUE),
    Mean_Width = mean(Width, na.rm = TRUE),
    SD_Width = sd(Width, na.rm = TRUE),
    Mean_Circularity = mean(Circularity, na.rm = TRUE),
    sd_Circularity = sd(Circularity, na.rm = TRUE),
    Total_Mito_Area = sum(Area, na.rm = TRUE),
    relMitoArea =sum(Area, na.rm = TRUE)/Cell_Area,
    Freq_Donut = sum(Shape == "Donut")/sum(Shape %in% c("Donut","LostCristae","Mitochondrium")),
    Freq_StandardMitochondrium = sum(Shape == "Mitochondrium")/sum(Shape %in% c("Donut","LostCristae","Mitochondrium")),
    Freq_LostCristae = sum(Shape == "LostCristae")/sum(Shape %in% c("Donut","LostCristae","Mitochondrium")),
    Freq_damaged=sum(Shape %in% c("LostCristae","Donut")) / sum(Shape %in% c("Donut","LostCristae","Mitochondrium")),
    .groups = "drop"
  )
saveRDS(Cell_Dataframe,file.path(PATHS$TEM$output,"Cell_Dataframe.rds"))
write.csv(Cell_Dataframe,file.path(PATHS$TEM$output,"Cell_Dataframe.csv"))
## Image Dataframe erstellen und speichern ---- 

Image_Metadata <- df %>%
  filter(Classification == "image_descriptor") %>%
  select( Animal,  Image_ID,Descriptor) %>%
  distinct()

Image_Dataframe <- Cell_Dataframe %>%
  group_by(Animal, Image_ID) %>%
  summarise(n_Cells = n(),
    Mean_Cell_Area = mean(Cell_Area, na.rm = TRUE),
    SD_Cell_Area = sd(Cell_Area, na.rm = TRUE),
    Mean_n_Mito = mean(n_Mito, na.rm = TRUE),
    SD_n_Mito = sd(n_Mito, na.rm = TRUE),
    Mean_Mito_Area = mean(Mean_Mito_Area, na.rm = TRUE),
    SD_Mito_Area = sd(Mean_Mito_Area, na.rm = TRUE),
    Mean_Mito_Length = mean(Mean_Length, na.rm = TRUE),
    SD_Mito_Length = sd(Mean_Length, na.rm = TRUE),
    Mean_Mito_Width = mean(Mean_Width, na.rm = TRUE),
    SD_Mito_Width = sd(Mean_Width, na.rm = TRUE),
    Mean_Circularity = mean(Mean_Circularity, na.rm = TRUE),
    SD_Circularity = sd(Mean_Circularity, na.rm = TRUE),
    Mean_relMitoArea = mean(relMitoArea, na.rm = TRUE),
    Mean_Freq_Donut = mean(Freq_Donut, na.rm = TRUE),
    Mean_Freq_StandardMitochondrium = mean(Freq_StandardMitochondrium, na.rm = TRUE),
    Mean_Freq_LostCristae = mean(Freq_LostCristae, na.rm = TRUE),
    Mean_Freq_damaged = mean(Freq_damaged, na.rm = TRUE),
    .groups = "drop"
  )
Image_Dataframe <- Image_Dataframe %>%left_join(Image_Metadata, by = c("Animal", "Image_ID"))
Animal_Dataframe <- Image_Dataframe %>%
  group_by(Animal) %>%
  summarise(
    n_Images = n(),
    Mean_Cells_per_Image = mean(n_Cells, na.rm = TRUE),
    SD_Cells_per_Image = sd(n_Cells, na.rm = TRUE),
    Mean_Cell_Area = mean(Mean_Cell_Area, na.rm = TRUE),
    SD_Cell_Area = sd(Mean_Cell_Area, na.rm = TRUE),
    Mean_Mito_per_Cell = mean(Mean_n_Mito, na.rm = TRUE),
    SD_Mito_per_Cell = sd(Mean_n_Mito, na.rm = TRUE),
    Mean_Mito_Area = mean(Mean_Mito_Area, na.rm = TRUE),
    SD_Mito_Area = sd(Mean_Mito_Area, na.rm = TRUE),
    Mean_Mito_Length = mean(Mean_Mito_Length, na.rm = TRUE),
    SD_Mito_Length = sd(Mean_Mito_Length, na.rm = TRUE),
    Mean_Mito_Width = mean(Mean_Mito_Width, na.rm = TRUE),
    SD_Mito_Width = sd(Mean_Mito_Width, na.rm = TRUE),
    Mean_Circularity = mean(Mean_Circularity, na.rm = TRUE),
    SD_Circularity = sd(Mean_Circularity, na.rm = TRUE),
    Mean_relMitoArea = mean(Mean_relMitoArea, na.rm = TRUE),
    Mean_Freq_Donut = mean(Mean_Freq_Donut, na.rm = TRUE),
    Mean_Freq_StandardMitochondrium = mean(Mean_Freq_StandardMitochondrium, na.rm = TRUE),
    Mean_Freq_LostCristae = mean(Mean_Freq_LostCristae, na.rm = TRUE),
    Mean_Freq_damaged = mean(Mean_Freq_damaged, na.rm = TRUE),
    .groups = "drop"
  )


# 
# Animal_Metadata <- data.frame(
#   Animal = c(1,2,3),
#   Group = c("Ctrl","Ctrl","Treatment"),
#   Sex = c("M","F","M")
# )
# 
# Animal_Dataframe <- Animal_Dataframe %>%left_join(Animal_Metadata,by = "Animal" )
saveRDS(Animal_Dataframe,file.path(PATHS$TEM$output,"Animal_Dataframe.rds"))
write.csv(Animal_Dataframe,file.path(PATHS$TEM$output,"Animal_Dataframe.csv"))



# #x<-df%>%group_by(Image_ID)%>%summarize(fraction= n(Mitochondrium)/(n(Mitochondrium)+n(Donut)))
# n_mito_image <- df %>%
#   filter(Classification %in% c("Mitochondrium", "Donut")) %>%
#   group_by(Image_ID) %>%
#   summarise(n_mito = n())
# 
# n_mito_cell <- df %>%
#   filter(Classification %in% c("Mitochondrium", "Donut")) %>%
#   group_by(Image_ID, Cell_ID) %>%
#   summarise(n_mito = n())
# 
# mito_area_image <- df %>%
#   filter(Classification %in% c("Mitochondrium", "Donut")) %>%
#   group_by(Image_ID) %>%
#   summarise(total_mito_area = sum(Area, na.rm = TRUE))
# 
# mito_area_cell <- df %>%
#   filter(Classification %in% c("Mitochondrium", "Donut")) %>%
#   group_by(Image_ID, Cell_ID) %>%
#   summarise(total_mito_area = sum(Area, na.rm = TRUE))
# 
# hepatocyte_area <- df %>%
#   filter(Classification == "Hepatocyte") %>%
#   select(Image_ID, Cell_ID, cell_area = Area)
# 
# mito_density_cell <- n_mito_cell %>%
#   left_join(hepatocyte_area,
#             by = c("Image_ID", "Cell_ID")) %>%
#   mutate(mito_per_area = n_mito / cell_area)
# 
# mito_fraction_cell <- mito_area_cell %>%
#   left_join(hepatocyte_area,
#             by = c("Image_ID", "Cell_ID")) %>%
#   mutate(mito_area_fraction = total_mito_area / cell_area)
# 
# 
# 
# df %>%
#   filter(Classification == "Length") %>%
#   ggplot(aes(x = Analysis_ID, y = Length.µm)) +
#   geom_jitter(width = 0.2, size = 1, alpha = 0.6) +
#   theme_bw() +
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
#   labs(x = "Analysis_ID", y = "Mitochondrial length (µm)")
#   
#   df %>%
#     filter(Classification == "Width") %>%
#     ggplot(aes(x = Image_ID, y = Length.µm)) +
#     geom_jitter(width = 0.2, size = 1, alpha = 0.6) +
#     theme_bw() +
#     theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
#   labs(x = "Image_ID", y = "Mitochondrial width (µm)")
#   
#   df %>%
#     filter(Classification %in% c("Mitochondrium", "Donut")) %>%
#     ggplot(aes(x = Image_ID, y = Area)) +
#     geom_jitter(width = 0.2, size = 1, alpha = 0.6) +
#     theme_bw() +
#     theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
#     labs(x = "Image_ID", y = "Mitochondrial Area (µm2)")
#   
#   df %>%
#     filter(Classification %in% c("Mitochondrium", "Donut")) %>%
#     ggplot(aes(x = Image_ID, y = Perimeter)) +
#     geom_jitter(width = 0.2, size = 1, alpha = 0.6) +
#     theme_bw() +
#     theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
#     labs(x = "Image_ID", y = "Mitochondrial Perimeter (µm)")
#   
#  df %>%
#     filter(Classification %in% c("Mitochondrium", "Donut")) %>%
#     ggplot(aes(x = Image_ID, y = Perimeter)) +
#     geom_jitter(width = 0.2, size = 1, alpha = 0.6) +
#     theme_bw() +
#     theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
#     labs(x = "Image_ID", y = "Mitochondrial Perimeter (µm)")
#   
#   mitos_circle<-df %>%
#     filter(Classification %in% c("Mitochondrium", "Donut")) %>%
#     ggplot(aes(x = Image_ID, y = Circularity)) +
#     geom_jitter(width = 0.2, size = 1, alpha = 0.6) +
#     theme_bw() +
#     theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
#     labs(x = "Image_ID", y = "Mitochondrial Circularity")
#   ggsave(mitos_circle,path = TEM_out_pwd,file="Mito_Circularity.png",width = 8, height=4)
# mito_density_cell %>%
#     ggplot(aes(x = as.character(Image_ID), y = mito_per_area)) +
#     geom_jitter(width = 0.1, size = 1, alpha = 0.6) +
#     theme_bw() +
#     theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
#     labs(x = "Image_ID", y = "Mitochondria/area/cell")
# 
# mito_fraction_cell %>%
#   ggplot(aes(x = as.character(Image_ID), y = mito_area_fraction*100)) +
#   geom_jitter(width = 0.2, size = 1, alpha = 0.6) +
#   theme_bw() +
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
#   labs(x = "Image_ID", y = "Mitochondria/area/cell")
#  