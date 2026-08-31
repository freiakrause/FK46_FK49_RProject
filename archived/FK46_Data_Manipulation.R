library(dplyr)
library(tidyr)
library(stringr)
library(lubridate)
ExpId="FK46"
df<-read.csv("C:/Users/b1084855/OneDrive - Universität Salzburg/Freigegebene Dokumente - AG_Tumorimmunologie/Data/Freia Krause/01_Experiments/FK46_iALmice_high Fat diet 52 weeks 7d after injection/Raw Data/FK46_Data_for_R.csv",sep=";")
setwd("C:/Users/b1084855/OneDrive - Universität Salzburg/Freigegebene Dokumente - AG_Tumorimmunologie/Data/Freia Krause/01_Experiments/FK46_iALmice_high Fat diet 52 weeks 7d after injection/Analysis")
startweight<-df%>%filter(DOW ==START.Diet)%>%mutate(Startweight=Weight)%>%select(Animal,Startweight)
EP_weight<-df%>%filter(DOW ==KILL.DATE)%>%mutate(EP_weight=Weight)%>%select(Animal,EP_weight)

exigo <- c("ALB", "TP", "GLOB", "TB", "AST", "ALT", "GGT", "ALP", "TBA", "TC","A.G")  #Parametes LiverPanel Exigo
data <- df %>%
  merge(startweight, by = "Animal") %>%
  merge(EP_weight, by = "Animal") %>%
  mutate(Treatment = factor(Treatment, levels = c("EtOH","TAM"), labels = c("ctrl", "TAM")))%>%
  mutate(Sex = as.factor(Sex)) %>%
  mutate(Animal = as.character(Animal)) %>%
  mutate(across(c(KILL.DATE, START.I, START.Diet, START, DOB, DOW), ~ as.Date(., "%d.%m.%Y"))) %>%
  mutate(wks_dead = round((as.numeric(KILL.DATE - START.Diet) / 7), digits = 2),
         wks_exp_total = round((as.numeric(KILL.DATE - START.I) / 7), digits = 1),
         age_total = round((as.numeric(KILL.DATE - DOB) / 7), digits = 1),
         wks_diet = round((as.numeric(DOW - START.Diet) / 7), digits = 2),
         age_start = round((as.numeric(START - DOB) / 7), digits = 1),
         Weight = as.numeric(gsub(",", ".", gsub("x", "NA", Weight))),
         Startweight = as.numeric(gsub(",", ".", gsub("x", "NA", Startweight))),
         EP_weight = as.numeric(gsub(",", ".", gsub("x", "NA", EP_weight))),
         rel.weight = (Weight / Startweight) * 100,
         rel_EP_weight=(EP_weight/Startweight)* 100,
         Score = as.numeric(Score),
         Treatment =gsub("EtOH", "Ctrl",Treatment),
         NASH_S1 = as.numeric(NASH_S1),
         NASH_S2 = as.numeric(NASH_S2),
         NASH_B1 = as.numeric(NASH_B1),
         NASH_B2 = as.numeric(NASH_B2),
         NASH_I1 = as.numeric(NASH_I1),
         NASH_I2 = as.numeric(NASH_I2),)%>%
  mutate(across(c(Liver,Spleen,Weight), ~ as.numeric(gsub(",", ".", .)))) %>%
  mutate(Spleen = as.numeric(Spleen) * 1000) %>%
  mutate(Liver_rel = Liver / Weight * 100,
         Spleen_rel = ((Spleen / 1000) / Weight) * 100)%>%
  mutate(across(all_of(exigo), as.character)) %>%
  mutate(DFactor = (USEPANEL + USEBUFFER) / USEPANEL) %>%
  pivot_longer(all_of(exigo), names_to = "parameter", values_to = "value") %>%
  mutate(censored = str_detect(value, "[<>]"),
         direction = case_when(str_detect(value, "^>") ~ ">",  str_detect(value, "^<") ~ "<",TRUE ~ NA_character_),
         numeric_value = as.numeric(str_remove_all(value, "[^0-9\\.]")),
         numeric_value_diluted = numeric_value * DFactor)  %>%         
  pivot_wider(names_from  = parameter,
              values_from = c(value, censored, direction, numeric_value, numeric_value_diluted),
              names_glue  = "{parameter}_{.value}") %>%
  rename_with(~ str_replace(., "_numeric_value_diluted$", "")) %>%
  rowwise() %>%
  mutate( NASH_S = median(c_across(c(NASH_S1, NASH_S2)), na.rm = TRUE),
          NASH_B = median(c_across(c(NASH_B1, NASH_B2)), na.rm = TRUE),
          NASH_I = median(c_across(c(NASH_I1, NASH_I2)), na.rm = TRUE),
          NASH_SAF = if_all(c(NASH_S, NASH_B, NASH_I), is.na) %>%
            if_else(., NA_real_, sum(c_across(c(NASH_S, NASH_B, NASH_I)), na.rm = TRUE))) %>%
  ungroup() %>%
  mutate(across(c(NASH_S, NASH_B, NASH_I),  ~ factor(., 
                                                     levels = c(0, 0.25, 0.5, 0.75, 1, 1.25, 1.5, 1.75, 2, 2.25, 2.5, 2.75, 3), 
                                                     ordered = TRUE)))


#Checkpoint 1  for wks_diet
data_sum_for_check1<-data%>%group_by(wks_diet,BATCH)%>%
  summarise(weight = mean(Weight, na.rm = TRUE), n = n(), sd = sd(Weight, na.rm = TRUE))

#Mutate individual timepoints, which would mess up automated mutation
data<-data%>%mutate(wks_diet= case_when(wks_diet== 7.57 ~ 7.835,
                                        wks_diet== 8.14 ~ 7.835,
                                        wks_diet== 25.43 ~ 25.715,
                                        wks_diet== 26.00 ~ 25.715,
                                        TRUE ~wks_diet))

#Checkpoint 2  for wks_diet if manual mutation worked
data_sum_for_check2<-data%>%group_by(wks_diet,BATCH)%>%
  summarise(weight = mean(Weight, na.rm = TRUE),n = n(), sd = sd(Weight, na.rm = TRUE))

#Automatically mutate individual timepoints
#Since I sometimes measured BATCH 1 and BATCH 2 not exactly on the day after 
#the same time but sometimes 1 or 2 days before or after, 
#i want to group some timepoints to the mean of the timepoint
group_close_timepoints <- function(inputdata, tolerance = 0.4) {
  unique_times <- sort(unique(inputdata$wks_diet))
  groups <- list()
  
  # Create groups of close values
  while (length(unique_times) > 0) {
    ref <- unique_times[1]
    close_vals <- unique_times[abs(unique_times - ref) <= tolerance]
    groups[[length(groups) + 1]] <- close_vals
    unique_times <- setdiff(unique_times, close_vals)}
  
  # Create a lookup table for replacements
  replacements <- lapply(groups, function(g) {
    rep(mean(g), length(g))}) %>% unlist()
  
  value_map <- data.frame(
    original = unlist(groups),
    new = replacements)
  
  # Join with original data
  inputdata <- inputdata %>%
    left_join(value_map, by = c("wks_diet" = "original")) %>%
    mutate(wks_diet = ifelse(is.na(new), wks_diet, new)) %>%
    select(-new)%>%
    mutate(wks_diet = round(wks_diet,digits=1))
  return(inputdata)
}

# Save in data2 at this moment, to not mess up data since i might do troubleshooting
data2<-group_close_timepoints(data)

#Checkpoint 3  for wks_diet if automated mutation worked
data_sum_for_check3<-data2%>%group_by(wks_diet,BATCH)%>%
  summarise(weight = mean(Weight, na.rm = TRUE), n = n(), sd = sd(Weight, na.rm = TRUE))

#Since everything worked as expected, remove data2, save only data
data<-data2
rm(data2,data_sum_for_check1,data_sum_for_check2,data_sum_for_check3,df)
gc()
Legendplex_data <-NULL
if (ExpId=="FK49") {
  pw<-"C:/Users/b1084855/OneDrive - Universität Salzburg/Freigegebene Dokumente - AG_Tumorimmunologie/Data/Freia Krause/01_Experiments/FK49_CD-HFD_13wks/FK49_Legendplex/"
  Legendplex_data <-read.csv(paste0(pw,"FK49_Legenplex_clean.csv") )%>%
    mutate(Animal = as.character(Animal))
}  else if (ExpId == "FK46"){
  #pw<-"C:/Users/b1084855/OneDrive - Universität Salzburg/Freigegebene Dokumente - AG_Tumorimmunologie/Data/Freia Krause/01_Experiments/FK46_iALmice_high Fat diet 52 weeks 7d after injection/FK46_Legendplex/"
  pw<-"C:/Users/b1084855/OneDrive - Universität Salzburg/Freigegebene Dokumente - AG_Tumorimmunologie/Data/Freia Krause/01_Experiments/FK49_CD-HFD_13wks/FK49_Legendplex/"
  Legendplex_data <-read.csv(paste0(pw,"FK46_Legenplex_clean.csv") )%>%
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
rm(list = ls())
gc()
