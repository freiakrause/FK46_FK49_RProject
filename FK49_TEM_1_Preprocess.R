rm(list=ls())
gc()

library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)
library(gridExtra)
library(grid)
source("FK49_Definitions.R")

# Load Data and Clean up colums --------------------------------
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

# Load Animal Metadata ----------------------------------------------------

animal_metadata<-read.csv2(file.path(PATHS$TEM$input, "/TEM_Animal_Metadata.csv"),sep=";")


# Generate Pedigree -------------------------------------------------------
# to allocate measurements to specific mitochondria in Specific cells in specific images 
df<-data%>% mutate(Image_ID = as.numeric(factor(Image)))

# Initialize columns
df$Cell_ID <- NA
df$Mito_ID <- NA

# Set start to 0
current_image <- 0
current_cell <- 0
current_mito <- 0
#Over all rows sum up Image ID, Cell ID, Mito ID
for(i in 1:nrow(df)) {
  
  # New image
  if(df$Image[i] != current_image){
    current_image <- df$Image[i]
    current_cell <- 0
    current_mito <- 0
  }
  
  # New cell
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

## Get Analysis Order Id -------------------------------------------------------------
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

#QC check input data for completeness and correctness ----------------------
  ## Do all Mitos have 2  measurments -----
  #in pedigree below them (lenght and width)
  #Look at  problematic mitochondria and improve them in qupath! No mito should be problem!
   # problematic ==length or width measurement missing

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
    by = c("Animal" = "Animal_blind") )%>%
  mutate(Animal_blind = Animal, 
         Animal = Animal_real)%>%
  select(-Animal_blind,-Animal_real)

total_numbers <- df %>%
  summarise(Total_Animals   = n_distinct(Animal),
    Total_Images    = n_distinct(Image_ID),
    Total_TAM       = n_distinct(Animal[Treatment == "TAM"]),
    Total_EtOH      = n_distinct(Animal[Treatment == "EtOH"]),
    Total_Female    = n_distinct(Animal[Sex == "female"]),
    Total_Male      = n_distinct(Animal[Sex == "male"]),
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
                            Sex == "female"),
    Mito_Male       = sum(Classification %in% c("Mitochondrium","Donut","LostCristae") &
                            Sex == "male")) %>%
  pivot_longer( cols = everything(),names_to = "Metric",values_to = "Count" )

write.csv(total_numbers,file=file.path(PATHS$TEM$output, "Background/1_TotalNumbers.csv"))
png(file.path(PATHS$TEM$output,"Background/1_TotalNumbers.png"),  
    width = 800,  height = 1500,  res = 300)
grid.table(total_numbers)
dev.off()

## AAnimal summary
animal_summary <- df %>%
  group_by(Animal, Treatment, Sex) %>%
  summarise(# Sampling depth
    Images = n_distinct(Image_ID),
    Total_Cells = sum(Classification %in% c("Hepatocyte", "nonHepatocyte")),
    Hepatocytes = sum(Classification == "Hepatocyte"),
    nonHeps = sum(Classification == "nonHepatocyte"),
    # Cell area
    Hep_Area = sum(Area[Classification == "Hepatocyte"],   na.rm = TRUE),
    Mean_Hep_Area = mean(Area[Classification == "Hepatocyte"],   na.rm = TRUE),
    # Mito numbers
    AllMitos = sum(Classification %in% c("Mitochondrium","Donut","LostCristae")),
    Mito_per_HepA = round(AllMitos / Hep_Area,4  ),
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
write.csv(animal_summary,file=file.path(PATHS$TEM$output, "Background/1_AnimalNumbers.csv"))
png(file.path(PATHS$TEM$output,"Background/1_AnimalNumbers.png"),  
    width = 6000,  height = 1200,  res = 300)
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
df<-df[,c("Animal","Sex","Treatment","Image_ID","Cell_ID","Mito_ID","Pedigree_ID","Classification","Area","Perimeter","Length.µm","Circularity","Descriptor","Analysis_ID")]
df<-df%>%mutate(Image_ID = as.factor(Image_ID))
df_long <- df %>%
  pivot_longer( cols = c(Area, Perimeter, Length.µm, Circularity),
    names_to = "Variable",
    values_to = "Value",
    values_drop_na = TRUE )
str(df)
head(df)
unique(df$Classification)



# Sum up data in Dataframes of different levels ----------------------------------------
##  Dataframe Mitos ------------------
Mito_Dataframe <- df %>%
  filter(!Classification %in% c("Hepatocyte", "nonHepatocyte", "image_descriptor")) %>%
  mutate( Shape = if_else(Classification %in% c("Donut", "Mitochondrium", "LostCristae"),
                          Classification, NA_character_  )) %>%
  group_by(Pedigree_ID) %>%
  fill(Shape, Area, Perimeter, Circularity, .direction = "down") %>%
  pivot_wider(id_cols = c(Animal, Sex,Treatment,Image_ID, Cell_ID, Mito_ID, Pedigree_ID,Shape, Area, Perimeter, Circularity, Analysis_ID),
    names_from = Classification,values_from = Length.µm) %>%
  select(-Donut,-Mitochondrium,-LostCristae)%>%
  ungroup()                           
saveRDS(Mito_Dataframe,file.path(PATHS$TEM$output,"CleanData/Mito_Dataframe.rds"))
write.csv(Mito_Dataframe,file.path(PATHS$TEM$output,"CleanData/Mito_Dataframe.csv"))


##  Dataframe Cells ------------------
Cell_Metadata <- df %>%
  filter(Classification %in% c("Hepatocyte", "nonHepatocyte", "image_descriptor")) %>%
  select( Animal,Image_ID,Cell_ID,Cell_Type = Classification, Cell_Area = Area) %>%
  distinct()

Cell_Dataframe <- Mito_Dataframe %>%
  left_join(Cell_Metadata,by = c("Animal", "Image_ID", "Cell_ID") )%>%
  group_by(Animal, Image_ID,Cell_ID) %>%
  summarise(Sex=first(Sex),
            Treatment=first(Treatment),
            n_Mito = n(),
            Cell_Type = first(Cell_Type),
            Cell_Area= mean(Cell_Area),
            
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
    SD_Circularity = sd(Circularity, na.rm = TRUE),
    Total_Mito_Area = sum(Area, na.rm = TRUE),
    relMitoArea =sum(Area, na.rm = TRUE)/Cell_Area,
    Freq_Donut = sum(Shape == "Donut")/sum(Shape %in% c("Donut","LostCristae","Mitochondrium")),
    Freq_StandardMitochondrium = sum(Shape == "Mitochondrium")/sum(Shape %in% c("Donut","LostCristae","Mitochondrium")),
    Freq_LostCristae = sum(Shape == "LostCristae")/sum(Shape %in% c("Donut","LostCristae","Mitochondrium")),
    Freq_damaged=sum(Shape %in% c("LostCristae","Donut")) / sum(Shape %in% c("Donut","LostCristae","Mitochondrium")),
    .groups = "drop"
  )
saveRDS(Cell_Dataframe,file.path(PATHS$TEM$output,"CleanData/Cell_Dataframe.rds"))
write.csv(Cell_Dataframe,file.path(PATHS$TEM$output,"CleanData/Cell_Dataframe.csv"))


##  Dataframe Iamges ------------------
Image_Metadata <- df %>%
  filter(Classification == "image_descriptor") %>%
  select( Animal, Image_ID,Descriptor) %>%
  distinct()

Image_Dataframe <- Cell_Dataframe %>%
  group_by(Animal, Image_ID)%>%
  summarize(Sex=first(Sex),
    Treatment=first(Treatment),
    n_Cells = n(),
    n_Heps= sum (Cell_Type== "Hepatocyte"),
    n_nonHeps= sum (Cell_Type== "nonHepatocyte"),
    
    Mean_Hep_Area = mean(Cell_Area[Cell_Type== "Hepatocyte"], na.rm = TRUE),
    SD_Hep_Area = sd(Cell_Area[Cell_Type== "Hepatocyte"], na.rm = TRUE),
    Mean_n_Mito = mean(n_Mito[Cell_Type== "Hepatocyte"], na.rm = TRUE),
    SD_n_Mito = sd(n_Mito[Cell_Type== "Hepatocyte"], na.rm = TRUE),
    Mean_Mito_Area = mean(Mean_Mito_Area[Cell_Type== "Hepatocyte"], na.rm = TRUE),
    SD_Mito_Area = sd(Mean_Mito_Area[Cell_Type== "Hepatocyte"], na.rm = TRUE),
    Mean_Mito_Length = mean(Mean_Length[Cell_Type== "Hepatocyte"], na.rm = TRUE),
    SD_Mito_Length = sd(Mean_Length[Cell_Type== "Hepatocyte"], na.rm = TRUE),
    Mean_Mito_Width = mean(Mean_Width[Cell_Type== "Hepatocyte"], na.rm = TRUE),
    SD_Mito_Width = sd(Mean_Width[Cell_Type== "Hepatocyte"], na.rm = TRUE),
    Mean_Circularity = mean(Mean_Circularity[Cell_Type== "Hepatocyte"], na.rm = TRUE),
    SD_Circularity = sd(Mean_Circularity[Cell_Type== "Hepatocyte"], na.rm = TRUE),
    Mean_relMitoArea = mean(relMitoArea[Cell_Type== "Hepatocyte"], na.rm = TRUE),
    SD_relMitoArea = sd(relMitoArea[Cell_Type== "Hepatocyte"], na.rm = TRUE),
    Mean_Freq_Donut = mean(Freq_Donut[Cell_Type== "Hepatocyte"], na.rm = TRUE),
    SD_Freq_Donut = sd(Freq_Donut[Cell_Type== "Hepatocyte"], na.rm = TRUE),
    Mean_Freq_StandardMitochondrium = mean(Freq_StandardMitochondrium[Cell_Type== "Hepatocyte"], na.rm = TRUE),
    SD_Freq_StandardMitochondrium = sd(Freq_StandardMitochondrium[Cell_Type== "Hepatocyte"], na.rm = TRUE),
    Mean_Freq_LostCristae = mean(Freq_LostCristae[Cell_Type== "Hepatocyte"], na.rm = TRUE),
    SD_Freq_LostCristae = sd(Freq_LostCristae[Cell_Type== "Hepatocyte"], na.rm = TRUE),
    Mean_Freq_damaged = mean(Freq_damaged[Cell_Type== "Hepatocyte"], na.rm = TRUE),
    SD_Freq_damaged = sd(Freq_damaged[Cell_Type== "Hepatocyte"], na.rm = TRUE),
    
    Mean_nonHep_Area = mean(Cell_Area[Cell_Type== "nonHepatocyte"], na.rm = TRUE),
    SD_nonHep_Area = sd(Cell_Area[Cell_Type== "nonHepatocyte"], na.rm = TRUE),
    Mean_nonHep_n_Mito = mean(n_Mito[Cell_Type== "nonHepatocyte"], na.rm = TRUE),
    SD_nonHep_n_Mito = sd(n_Mito[Cell_Type== "nonHepatocyte"], na.rm = TRUE),
    Mean_nonHep_Mito_Area = mean(Mean_Mito_Area[Cell_Type== "nonHepatocyte"], na.rm = TRUE),
    SD_nonHep_Mito_Area = sd(Mean_Mito_Area[Cell_Type== "nonHepatocyte"], na.rm = TRUE),
    Mean_nonHep_Mito_Length = mean(Mean_Length[Cell_Type== "nonHepatocyte"], na.rm = TRUE),
    SD_nonHep_Mito_Length = sd(Mean_Length[Cell_Type== "nonHepatocyte"], na.rm = TRUE),
    Mean_nonHep_Mito_Width = mean(Mean_Width[Cell_Type== "nonHepatocyte"], na.rm = TRUE),
    SD_nonHep_Mito_Width = sd(Mean_Width[Cell_Type== "nonHepatocyte"], na.rm = TRUE),
    Mean_nonHep_Circularity = mean(Mean_Circularity[Cell_Type== "nonHepatocyte"], na.rm = TRUE),
    SD_nonHep_Circularity = sd(Mean_Circularity[Cell_Type== "nonHepatocyte"], na.rm = TRUE),
    Mean_nonHep_relMitoArea = mean(relMitoArea[Cell_Type== "nonHepatocyte"], na.rm = TRUE),
    SD_nonHep_relMitoArea = sd(relMitoArea[Cell_Type== "nonHepatocyte"], na.rm = TRUE),
    Mean_nonHep_Freq_Donut = mean(Freq_Donut[Cell_Type== "nonHepatocyte"], na.rm = TRUE),
    SD_nonHep_Freq_Donut = sd(Freq_Donut[Cell_Type== "nonHepatocyte"], na.rm = TRUE),
    Mean_nonHep_Freq_StandardMitochondrium = mean(Freq_StandardMitochondrium[Cell_Type== "nonHepatocyte"], na.rm = TRUE),
    SD_nonHep_Freq_StandardMitochondrium = sd(Freq_StandardMitochondrium[Cell_Type== "nonHepatocyte"], na.rm = TRUE),
    Mean_nonHep_Freq_LostCristae = mean(Freq_LostCristae[Cell_Type== "nonHepatocyte"], na.rm = TRUE),
    SD_nonHep_Freq_LostCristae = sd(Freq_LostCristae[Cell_Type== "nonHepatocyte"], na.rm = TRUE),
    Mean_nonHep_Freq_damaged = mean(Freq_damaged[Cell_Type== "nonHepatocyte"], na.rm = TRUE),
    SD_nonHep_Freq_damaged = sd(Freq_damaged[Cell_Type== "nonHepatocyte"], na.rm = TRUE),
    .groups = "drop"
  )
Image_Dataframe <- Image_Dataframe %>%
  left_join(Image_Metadata, by = c("Animal", "Image_ID"))

saveRDS(Image_Dataframe,file.path(PATHS$TEM$output,"CleanData/Image_Dataframe.rds"))
write.csv(Image_Dataframe,file.path(PATHS$TEM$output,"CleanData/Image_Dataframe.csv"))


##  Dataframe Animals------------------
Animal_Dataframe <- Image_Dataframe %>%
  group_by(Animal) %>%
  summarise(Sex=first(Sex),
            Treatment=first(Treatment),
            n_Images = n(),
            
            Mean_Heps_per_Image = mean(n_Heps, na.rm = TRUE),
            SD_Heps_per_Image = sd(n_Heps, na.rm = TRUE),
            Mean_Hep_Area = mean(Mean_Hep_Area, na.rm = TRUE),
            SD_Hep_Area = sd(Mean_Hep_Area, na.rm = TRUE),
            Mean_Mito_per_Hep = mean(Mean_n_Mito, na.rm = TRUE),    
            SD_Mito_per_Hep = sd(Mean_n_Mito, na.rm = TRUE),
            Mean_Mito_Area = mean(Mean_Mito_Area, na.rm = TRUE),
            SD_Mito_Area = sd(Mean_Mito_Area, na.rm = TRUE),
            Mean_Mito_Length = mean(Mean_Mito_Length, na.rm = TRUE),
            SD_Mito_Length = sd(Mean_Mito_Length, na.rm = TRUE),
            Mean_Mito_Width = mean(Mean_Mito_Width, na.rm = TRUE),
            SD_Mito_Width = sd(Mean_Mito_Width, na.rm = TRUE),
            Mean_Circularity = mean(Mean_Circularity, na.rm = TRUE),
            SD_Circularity = sd(Mean_Circularity, na.rm = TRUE),
            Mean_relMitoArea = mean(Mean_relMitoArea, na.rm = TRUE),
            SD_relMitoArea = sd(Mean_relMitoArea, na.rm = TRUE),
            Mean_Freq_Donut = mean(Mean_Freq_Donut, na.rm = TRUE),
            SD_Freq_Donut = sd(Mean_Freq_Donut, na.rm = TRUE),
            Mean_Freq_StandardMitochondrium = mean(Mean_Freq_StandardMitochondrium, na.rm = TRUE),
            SD_Freq_StandardMitochondrium = sd(Mean_Freq_StandardMitochondrium, na.rm = TRUE),
            Mean_Freq_LostCristae = mean(Mean_Freq_LostCristae, na.rm = TRUE),
            SD_Freq_LostCristae = sd(Mean_Freq_LostCristae, na.rm = TRUE),
            Mean_Freq_damaged = mean(Mean_Freq_damaged, na.rm = TRUE),
            SD_Freq_damaged = sd(Mean_Freq_damaged, na.rm = TRUE),
            
            Mean_nonHeps_per_Image = mean(n_nonHeps, na.rm = TRUE),
            SD_nonHeps_per_Image = sd(n_nonHeps, na.rm = TRUE),
            Mean_nonHep_Area = mean(Mean_nonHep_Area, na.rm = TRUE),
            SD_nonHep_Area = sd(Mean_nonHep_Area, na.rm = TRUE),
            Mean_Mito_per_nonHep = mean(Mean_nonHep_n_Mito, na.rm = TRUE),
            SD_Mito_per_nonHep = sd(Mean_nonHep_n_Mito, na.rm = TRUE),
            Mean_nonHep_Mito_Area = mean(Mean_nonHep_Mito_Area, na.rm = TRUE),
            SD_nonHep_Mito_Area = sd(Mean_nonHep_Mito_Area, na.rm = TRUE),
            Mean_nonHep_Mito_Length = mean(Mean_nonHep_Mito_Length, na.rm = TRUE),
            SD_nonHep_Mito_Length = sd(Mean_nonHep_Mito_Length, na.rm = TRUE),
            Mean_nonHep_Mito_Width = mean(Mean_nonHep_Mito_Width, na.rm = TRUE),
            SD_nonHep_Mito_Width = sd(Mean_nonHep_Mito_Width, na.rm = TRUE),
            Mean_nonHep_Circularity = mean(Mean_nonHep_Circularity, na.rm = TRUE),
            SD_nonHep_Circularity = sd(Mean_nonHep_Circularity, na.rm = TRUE),
            Mean_nonHep_relMitoArea = mean(Mean_nonHep_relMitoArea, na.rm = TRUE),
            SD_nonHep_relMitoArea = sd(Mean_nonHep_relMitoArea, na.rm = TRUE),
            Mean_nonHep_Freq_Donut = mean(Mean_nonHep_Freq_Donut, na.rm = TRUE),
            SD_nonHep_Freq_Donut = sd(Mean_nonHep_Freq_Donut, na.rm = TRUE),
            Mean_nonHep_Freq_StandardMitochondrium = mean(Mean_nonHep_Freq_StandardMitochondrium, na.rm = TRUE),
            SD_nonHep_Freq_StandardMitochondrium = sd(Mean_nonHep_Freq_StandardMitochondrium, na.rm = TRUE),
            Mean_nonHep_Freq_LostCristae = mean(Mean_nonHep_Freq_LostCristae, na.rm = TRUE),
            SD_nonHep_Freq_LostCristae = sd(Mean_nonHep_Freq_LostCristae, na.rm = TRUE),
            Mean_nonHep_Freq_damaged = mean(Mean_nonHep_Freq_damaged, na.rm = TRUE),
            SD_nonHep_Freq_damaged = sd(Mean_nonHep_Freq_damaged, na.rm = TRUE),
            .groups = "drop"
  )

saveRDS(Animal_Dataframe,file.path(PATHS$TEM$output,"CleanData/Animal_Dataframe.rds"))
write.csv(Animal_Dataframe,file.path(PATHS$TEM$output,"CleanData/Animal_Dataframe.csv"))


