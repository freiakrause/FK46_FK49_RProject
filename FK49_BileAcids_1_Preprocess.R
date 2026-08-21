
rm(list = ls())
gc()
library(dplyr)
library(tidyr)
library(ggplot2)
library('corrr')
library(ggcorrplot)
library("FactoMineR")
library(factoextra)
library(ggrepel)
source("FK49_Definitions.R")
ExpId = "FK49"

  process_BA <- function(df, meta_cols) {
    df %>%
      mutate(Sample = as.factor(Sample),
        T_D_S = paste0(Treatment, "_", Diet, "_", Sex),
        T_D   = paste0(Treatment, "_", Diet),
        T_S   = paste0(Treatment, "_", Sex),
        Treatment = case_when(Treatment == "EtOH" ~ "Ctrl", TRUE ~ Treatment),
        across(where(is.character), ~ gsub(",", ".", .))) %>%
      pivot_longer( all_of(BA_to_test), names_to = "parameter", values_to = "value_raw" ) %>%
      mutate(
        censored = str_detect(value_raw, "[<>]"),
        direction = case_when(
                    str_detect(value_raw, "^>") ~ ">",
                    str_detect(value_raw, "^<") ~ "<",
                    TRUE ~ NA_character_  ),
        numeric_value = as.numeric(str_remove_all(value_raw, "[^0-9\\.]") )
      ) %>%
      group_by(parameter) %>%
      mutate( numeric_value = if_else(
                              value_raw == "<LOD",min(numeric_value[!censored & !is.na(numeric_value)], na.rm = TRUE) / 2,
                              numeric_value )
      ) %>%
      ungroup() %>%
      pivot_wider(names_from = parameter, values_from = c(value_raw,censored, direction,numeric_value ),names_glue = "{parameter}_{.value}" ) %>%
     mutate(
        Treatment = factor(Treatment, levels = c("Ctrl", "TAM")),
        Sex = factor(Sex, levels = c("female", "male")),
        Diet = factor(Diet, levels = c("ND", "CDHFD13")),
        Timepoint = factor(Timepoint, levels = c("-1", "11")),
        Time_Treat = paste0(Timepoint, "_", Treatment)
      ) %>%

      rename_with(~ str_replace(.x, "_numeric_value$", ""), ends_with("_numeric_value")) %>%
      rename_with(~ str_replace(.x, "_value_raw$", "_raw"), ends_with("_value_raw")) %>%
      mutate(
        SUM = rowSums(select(., all_of(BA_to_test)), na.rm = TRUE),
        SUMprim = rowSums(select(., all_of(BA_primary)), na.rm = TRUE),
        SUMsec = rowSums(select(., all_of(BA_secondary)), na.rm = TRUE)
      )
  }


if (ExpId=="FK49") {
  meta_cols <-PARAMETERS$BA$meta_cols
  BA_to_test<-PARAMETERS$BA$BA_sort
  BA_secondary<-PARAMETERS$BA$BA_secondary
  BA_primary<-PARAMETERS$BA$BA_primary
  output_pwd <- PATHS$BA$output
  BA <- read.csv(PATHS$BA$input,  sep=";", stringsAsFactors=FALSE, check.names=FALSE) %>%
    dplyr::select(-`Sample Code`, -`Sample No`) %>%
    process_BA(meta_cols)
  
  saveRDS(BA, file= paste0(dirname(dirname(output_pwd)),"/01_RawData/FK49_BA_raw.rds"))
 

  # Define the output folder names
  #output_folders <- c("CDHFD_a", "CDHFD_m", "CDHFD_f", "ND_f","NDvsCDHFD")
  
  
  # Create the output folders under both targeted and untargeted paths
  #create_output_folders(targeted_pwd, output_folders)

}  else if (ExpId == "FK46"){
  print("You dont have data for this Experiment")
} else{
  print("You dont have data for this Experiment")
}


