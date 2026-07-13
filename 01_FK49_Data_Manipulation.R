library(tidyr)
library(dplyr)
library(stringr)
library(lubridate)
library(ggplot2 )

# Read Raw Inputdata and general Data manipulation ------------------------------------------------------
ExpId="FK49"

  
data_manipulation_FK49 <- function(d = df,ExpId=NULL) {
  #set working directory. If there is a new (unknown to me) experiment, path needs to be included here
  if (ExpId=="FK49") {
    setwd("C:/Users/b1084855/OneDrive - Universität Salzburg/Freigegebene Dokumente - AG_Tumorimmunologie/Data/Freia Krause/01_Experiments/FK49_CD-HFD_13wks/FK49_Analysis")
  }  else if (ExpId == "FK46"){
    setwd("C:/Users/b1084855/OneDrive - Universität Salzburg/Freigegebene Dokumente - AG_Tumorimmunologie/Data/Freia Krause/01_Experiments/FK46_iALmice_high Fat diet 52 weeks 7d after injection/Analysis")
  } else{
    print("Let me set a folder Path and define ExpId and Folderpath")}
  
  # ------ Startweight and Endpoint Weight ------
  startweight <- d %>% filter(DOW == START.Diet) %>% mutate(Startweight = Weight) %>% select(Animal, Startweight)
  EP_weight <- d %>%filter(DOW == KILL.DATE) %>%mutate(EP_weight = Weight) %>%  select(Animal, EP_weight)
  exigo <- c("ALB", "TP", "GLOB","A.G", "TB", "GGT", "AST", "ALT", "ALP", "AMY","Crea","UA","BUN","GLU","TC","TG")
  
  # ------ Main Data Manipulation ------
  data <- d %>%
    left_join(startweight, by = "Animal") %>%
    left_join(EP_weight, by = "Animal") %>%
    mutate(across(c(Weight, Startweight, EP_weight), ~ gsub("x", NA_character_, .))) %>%
    mutate(Sex = as.factor(Sex)) %>%
    mutate(across(c(KILL.DATE, START.I, START.Diet, DOB, DOW),~ as.Date(.x, format = "%d.%m.%Y"))) %>%
    mutate( wks_dead       = round((as.numeric(KILL.DATE - START.Diet) / 7), 2),
            wks_exp_total  = round((as.numeric(KILL.DATE - START.I) / 7), 1),
            age_total      = round((as.numeric(KILL.DATE - DOB) / 7), 1),
            wks_diet       = round((as.numeric(DOW - START.Diet) / 7), 2),
            days_diet      = round((as.numeric(DOW - START.Diet)), 0),
            age_start      = round((as.numeric(START.Diet - DOB) / 7), 1)) %>%
    mutate(across(c(Liver, Spleen, Weight, LNc, LNld, LNm, LivFACS, Fat, Food, Water, Score, Startweight, EP_weight),
                  ~ as.numeric(gsub(",", ".", .)))) %>%
    mutate(rel.weight    = (Weight / Startweight) * 100,
           rel_EP_weight = (EP_weight / Startweight) * 100,
           Spleen        = Spleen * 1000,
           Liver_rel     = Liver / Weight * 100,
           Spleen_rel    = ((Spleen / 1000) / Weight) * 100,
           Fat_rel       = Fat / Weight * 100,
           NASH_S1 = as.numeric(NASH_S1),
           NASH_S2 = as.numeric(NASH_S2),
           NASH_B1 = as.numeric(NASH_B1),
           NASH_B2 = as.numeric(NASH_B2),
           NASH_I1 = as.numeric(NASH_I1),
           NASH_I2 = as.numeric(NASH_I2)) %>%
    rename( BATCH = Batch ) %>%
    mutate(Treatment = factor(Treatment,  levels = c("EtOH", "TAM"),  labels = c("ctrl", "TAM"))) %>%
    mutate(across(all_of(exigo), as.character)) %>%
    mutate(DFactor = (EXIGOSample + EXIGOBuffer) / EXIGOSample) %>%
    rowwise() %>%
    mutate( NASH_S = median(c_across(c(NASH_S1, NASH_S2)), na.rm = TRUE),
            NASH_B = median(c_across(c(NASH_B1, NASH_B2)), na.rm = TRUE),
            NASH_I = median(c_across(c(NASH_I1, NASH_I2)), na.rm = TRUE),
            NASH_SAF = if_all(c(NASH_S, NASH_B, NASH_I), is.na) %>%
              if_else(., NA_real_, sum(c_across(c(NASH_S, NASH_B, NASH_I)), na.rm = TRUE))) %>%
    ungroup() %>%
    mutate(across(c(NASH_S, NASH_B, NASH_I),  ~ factor(., 
                                                       levels = c(0, 0.25, 0.5, 0.75, 1, 1.25, 1.5, 1.75, 2, 2.25, 2.5, 2.75, 3), 
                                                       ordered = TRUE)))%>%
    # ------ Pivoting EXIGO values ------
  pivot_longer(all_of(exigo),  names_to = "parameter", values_to = "value") %>%
    mutate( censored = str_detect(value, "[<>]"),
            direction = case_when(
              str_detect(value, "^>") ~ ">",
              str_detect(value, "^<") ~ "<",
              TRUE ~ NA_character_),
            numeric_value = as.numeric(str_remove_all(value, "[^0-9\\.]")),
            numeric_value_diluted = numeric_value * DFactor) %>%
    pivot_wider(names_from  = parameter,
                values_from = c(value, censored, direction, numeric_value, numeric_value_diluted),
                names_glue  = "{parameter}_{.value}") %>%
    rename_with(~ str_replace(., "_numeric_value_diluted$", ""))
    
  
  animal_count <- data %>%  group_by(DOW, Cage) %>%  summarise(n_animals = n(), .groups = "drop")
  data <- data %>%left_join(animal_count, by = c("DOW", "Cage"))
  

  
  rm(animal_count,EP_weight,startweight)
  ## Prepare Food Measurements ------------------------------------------------------
  df_food <- data %>%
    arrange(Cage, DOW) %>%
    group_by(Cage, DOW) %>%
    summarise(Food = first(Food), n_animals = n(),.groups = "drop") %>%
    arrange(Cage, DOW) %>%
    group_by(Cage) %>%
    mutate( DayDiff = as.numeric(DOW - lag(DOW, default = first(DOW))),
            NewPeriod = if_else(is.na(DayDiff) | DayDiff > 2, 1, 0),
            Block = cumsum(NewPeriod),
            Food_consumed_cage = lag(Food) - Food,
            Food_consumed = Food_consumed_cage / n_animals,
            Block = as.factor(Block),
            DOW_days = cumsum(DayDiff)) %>%
    ungroup() %>%
    mutate(Food_consumed = case_when(        # easy fix to get NA/0 on first day per block: 
      DOW_days == 7  ~ NA,                     # otherwise it compares with last weighing day of block before. i dont want that
      DOW_days == 63 ~ NA,
      TRUE ~ Food_consumed))
  
  ## Prepare Water Measurements ------------------------------------------------------
  df_water <- data %>%
    arrange(Cage, DOW) %>%
    group_by(Cage, DOW) %>%
    summarise(Water = first(Water), n_animals = n(),.groups = "drop") %>%
    arrange(Cage, DOW) %>%
    group_by(Cage) %>%
    mutate( DayDiff = as.numeric(DOW - lag(DOW, default = first(DOW))),
            NewPeriod = if_else(is.na(DayDiff) | DayDiff > 2, 1, 0),
            Block = cumsum(NewPeriod),
            Water_consumed_cage = lag(Water) - Water,
            Water_consumed = Water_consumed_cage / n_animals,
            Block = as.factor(Block),
            DOW_days = cumsum(DayDiff)) %>%
    mutate(Water_consumed = case_when(        # easy fix to get NA/0 on first day per block: 
      DOW_days == 7  ~ NA,                     # otherwise it compares with last weighing day of block before. i dont want that
      DOW_days == 63 ~ NA,
      TRUE ~ Water_consumed))
  
  data <- data %>%
    left_join(df_water %>% select(Cage, DOW, Water_consumed), by = c("Cage", "DOW"))%>%
    left_join(df_food %>% select(Cage, DOW, Food_consumed, Block), by = c("Cage", "DOW"))
  ## Check Food and Water Loss in empty cages (manually look at it)------------------------------------------------------
  
  data<-data%>%
    mutate(Diet = case_when(
      Animal == "EC2" ~ "Normal",
      Animal == "EC3" ~ "Normal",
      Animal == "EC1" ~ "CD-HFD",
      !Animal %in% c("EC1", "EC2","EC3") & Block == "0" ~ "Normal",       #  Block 0 was before CD-HFD
      !Animal %in% c("EC1", "EC2","EC3") & Block != "0"~ "CD-HFD",        # After Block 0 CD-HFD started for animals
      TRUE ~ "WTF" ))
  
  empty_cage_loss_food <- data %>%
    filter(Animal == "EC1" | Animal == "EC2"| Animal == "EC3") %>% 
    filter(!is.na(Food_consumed)) %>%
    mutate(Diet = as.character(Diet))%>%
    group_by( Diet,Block,days_diet) %>%
    summarise(mean_empty_loss = mean(Food_consumed, na.rm = TRUE), .groups = "drop")
  
  plot<-ggplot(filter(empty_cage_loss_food), aes(x=Diet, y = mean_empty_loss)) +
    geom_boxplot(position = position_dodge(width = 0.8)) +
    theme_bw()
  ggsave(filename = "FK49_Food_Empty_Cages.png", plot = plot, path = "02_GeneratedData/FoodIntake/Background", width = 4.5, height = 9, dpi = 300)
  
  empty_cage_loss_water <- data %>%
    filter(Animal == "EC1" | Animal == "EC2"| Animal == "EC3") %>% 
    filter(!is.na(Water_consumed)) %>%
    group_by( Block,days_diet) %>%
    summarise(mean_empty_loss = mean(Water_consumed, na.rm = TRUE), .groups = "drop")
  water_loss<-empty_cage_loss_water%>%summarise(mean_empty_loss = mean(mean_empty_loss, na.rm = TRUE))
  
  plot<-ggplot(filter(empty_cage_loss_water), aes( x=Block, y = mean_empty_loss)) +
    geom_boxplot(position = position_dodge(width = 0.8)) +
    theme_bw()
  ggsave(filename = "FK49_Water_Empty_Cages.png", plot = plot, path = "02_GeneratedData/FoodIntake/Background", width = 4.5, height = 9, dpi = 300)
  
  
  #Decided not to subtract food and water values of empty cages from animal cages since it is quite variable,
  #also i get negative values for animal water conspution and that cannot be correct
  #I am unsure if to subtract the daily values or the mean values I figure all animal cages should have the same mean loss due to handling
  rm(empty_cage_loss_water,empty_cage_loss_food,water_loss,plot)
  
  ## Add Food and Water Consumption to Dataframe ------------------------------------------------------
  rm(df_water,df_food)
  
  
  data_sum_for_check1<-data%>%
    group_by(wks_diet,BATCH)%>%
    summarise(weight = mean(Weight, na.rm = TRUE), 
              n = n(),                                # Are there the expected numbers of animals per week/day
              sd = sd(Weight, na.rm = TRUE))         # Yes to the main weeks there are 10 and 18 animals for Batch 1 and 2. 
  #But we also see some weighing dates that only happend in Batch 1or 2
  
  # week -0.86 and -0.71 was forgotten to weight in Batch 2 since sheet was hiding in animal room :'D
  # week 1, 1.14, 1.29, 1.57 was measured in Batch 2 but not 1 because animal caretakes got confused and thought they had to weigh food and water and mice everyday every week
  # week 4.29 in Batch 1 should have been 4.43 but animals were weight 1 day before due to Brandmeldeanalgen Test in Uni and closure of uni. It's 1 day difference
  # week 8, 8.14,8.29,8.57  was measured in Batch 1 but not 2 because animal caretaker got confused and thought they had to weigh food and water and mice everyday every week
  # week 12.29  12.43 13.14 and 13.29 becuase date of prep was different for Batch 1 an dBatch 2. batch 2 was prepr a bit later in their course of diet.
  
  ## Mutate individual time points ------------------------------------------------------
  # which would mess up automated mutation to summaries the different batches. 
  # 1 day real life difference is ignored here
  
  data<-data%>%mutate(wks_diet= case_when(
    wks_diet== 4.29 ~ 4.43,
    wks_diet== 12.29 ~ 12.43,
    TRUE ~wks_diet),
    days_diet=case_when(
      days_diet == 30 ~31,
      days_diet == 86~87,
      TRUE ~days_diet)
  )
  
  ## Checkpoint 2  for wks_diet if manual mutation worked ------------------------------------------------------
  data_sum_for_check2<-data%>%
    group_by(wks_diet,days_diet,BATCH)%>%
    summarise(weight = mean(Weight, na.rm = TRUE),
              n= n(), sd = sd(Weight, na.rm = TRUE))
  
  ## Automatically mutate individual timepoints ------------------------------------------------------
  # was commented out here becuase mutating manually individual points worked out fine
  #Since I sometimes measured BATCH 1 and BATCH 2 not exactly on the day after 
  #the same time but sometimes 1 or 2 days before or after, 
  #i want to group some timepoints to the mean of the timepoint
  # group_close_timepoints <- function(inputdata, tolerance = 0.1) {
  #   unique_times <- sort(unique(inputdata$wks_diet))
  #   groups <- list()
  #   
  #     while (length(unique_times) > 0) {
  #     ref <- unique_times[1]
  #     close_vals <- unique_times[abs(unique_times - ref) <= tolerance]
  #     groups[[length(groups) + 1]] <- close_vals
  #     unique_times <- setdiff(unique_times, close_vals)}
  #   
  #     replacements <- lapply(groups, function(g) {
  #     rep(mean(g), length(g))}) %>% unlist()
  #   
  #   value_map <- data.frame(
  #     original = unlist(groups),
  #     new = replacements)
  #   
  #    inputdata <- inputdata %>%
  #     left_join(value_map, by = c("wks_diet" = "original")) %>%
  #     mutate(wks_diet = ifelse(is.na(new), wks_diet, new)) %>%
  #     select(-new)%>%
  #     mutate(wks_diet = round(wks_diet,digits=1))
  #   return(inputdata)
  # }
  # 
  # # Save in data2 at this moment, to not mess up data since i might do troubleshooting
  # data2<-group_close_timepoints(data)
  # 
  # ## Checkpoint 3  for wks_diet if automated mutation worked ------------------------------------------------------
  # data_sum_for_check3<-data2%>%
  #                     group_by(wks_diet,BATCH)%>%
  #                     summarise(weight = mean(Weight, na.rm = TRUE),
  #                      n = n(), sd = sd(Weight, na.rm = TRUE))
  #
  # rm( data_sum_for_check3,group_close_timepoints,data2)
  
  rm(data_sum_for_check1,data_sum_for_check2,d)
  Legendplex_data <-NULL
  if (ExpId=="FK49") {
    pw<-"C:/Users/b1084855/OneDrive - Universität Salzburg/Freigegebene Dokumente - AG_Tumorimmunologie/Data/Freia Krause/01_Experiments/FK49_CD-HFD_13wks/FK49_Legendplex/"
    Legendplex_data <-read.csv(paste0(pw,"FK49_Legenplex_clean.csv") )%>%
      mutate(Animal = as.character(Animal))
  }  else if (ExpId == "FK46"){
    pw<-"C:/Users/b1084855/OneDrive - Universität Salzburg/Freigegebene Dokumente - AG_Tumorimmunologie/Data/Freia Krause/01_Experiments/FK46_iALmice_high Fat diet 52 weeks 7d after injection/FK46_Legendplex/"
    Legendplex_data <-read.csv("FK46_Legendplex_clean.csv")%>%
      mutate(Animal = as.character(Animal))
  } else{
    print("I dont know if there is clean legenplex data")}  
  
  if(is.null(Legendplex_data)){
    data <- data
  }  else{
    LD <- data %>%
      filter(DOW == KILL.DATE)%>%
      select(Animal,DOW)%>%
      distinct()%>%
      left_join(Legendplex_data,by = "Animal")
    data<-data %>% left_join(LD, by = c("Animal","DOW") )    }
  gc()
  save(data,file = paste0("01_RawData/",ExpId,"_Data_prepared.Rda"))
  return(data=data)
}
df <- read.csv("C:/Users/b1084855/OneDrive - Universität Salzburg/Freigegebene Dokumente - AG_Tumorimmunologie/Data/Freia Krause/01_Experiments/FK49_CD-HFD_13wks/FK49_CD-HFD_13wks_Data.csv", sep = ";") %>%
  filter(!is.na(Animal), trimws(Animal) != "") 


#df <- read.csv("C:/Users/b1084855/OneDrive - Universität Salzburg/Freigegebene Dokumente - AG_Tumorimmunologie/Data/Freia Krause/01_Experiments/FK49_CD-HFD_13wks/FK49_Data_Ctrl_animals.csv", sep = ";") %>%
#  filter(!is.na(Animal), trimws(Animal) != "")
d<-data_manipulation_FK49(d=df,ExpId="FK49")
gc()
