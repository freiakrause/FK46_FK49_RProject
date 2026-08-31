rm(list = ls())
gc()
library(tidyverse)
source("FK49_Definitions.R")

ExpID <- "FK49"

if (ExpID == "FK49") {
  output_pwd <- file.path(PATHS$exigo$FK49_output)
  Exigo_cols <-PARAMETERS$EXIGO$FK49_Exigo_cols
  load(file.path( PATHS$exigo$FK49_input, "FK49_Data_prepared.Rda"))
   d1 <- data %>%select( Animal, Sex,Treatment,BATCH,-TV,starts_with(Exigo_cols)) %>%
    filter(!is.na(ALB)) %>% mutate( Animal = as.character(Animal),BATCH = as.character(BATCH) )
  
  
  load(file.path( PATHS$BH_baseline$input,"BH15_Data_prepared.Rda"))
  baseline_data <- data %>% select( Animal, Sex,Treatment,BATCH,-TV,starts_with(Exigo_cols)) %>%
    filter(!is.na(ALB)) %>%
    mutate( Animal = as.character(Animal),
            BATCH = as.character(BATCH),
            Treatment = gsub("ctrl", "Ctrl", Treatment) )
  
  exp_and_BaseLine <- bind_rows(d1, baseline_data)
  
    # Save prepared data -----
  save(d1, baseline_data,exp_and_BaseLine, file = file.path(PATHS$exigo$FK49_input,  "FK49_Exigo_prepared.Rda"))
  
} else if (ExpID == "FK46") {
  output_pwd <- file.path(PATHS$exigo$FK46_output)
  Exigo_cols <-PARAMETERS$EXIGO$FK46_Exigo_cols
  load(file.path(PATHS$exigo$FK46_input,"FK46_Data_prepared.Rda"))
  d1 <- data %>%select( Animal, Sex,Treatment,BATCH,-TV,starts_with(Exigo_cols)) %>%
    filter(!is.na(ALB)) %>% mutate( Animal = as.character(Animal),BATCH = as.character(BATCH) )
  
  
  load(file.path( PATHS$BH_baseline$input,"BH15_Data_prepared.Rda"))
  baseline_data <- data %>% select( Animal, Sex,Treatment,BATCH,-TV,starts_with(Exigo_cols)) %>%
    filter(!is.na(ALB)) %>%
    mutate( Animal = as.character(Animal),
            BATCH = as.character(BATCH),
            Treatment = gsub("ctrl", "Ctrl", Treatment) )
  
  exp_and_BaseLine <- bind_rows(d1, baseline_data)
  
  # Save prepared data -----
  save(d1, baseline_data,exp_and_BaseLine, file = file.path(PATHS$exigo$FK46_input,  "FK46_Exigo_prepared.Rda"))
  
} else {stop("Unknown ExpID.")
}

load(file.path( PATHS$exigo$FK49_input, "FK49_Data_prepared.Rda"))

