rm=(list=ls())
gc()
library(tidyverse)
library(NADA2)
library(emmeans)
source("FK49_Definitions.R")
ExpID= "FK49"   # Decide if you want to load data from FK46 or FK49
 
if(ExpID == "FK49"){
   load(file = file.path(PATHS$exigo$FK49_input,  "FK49_Exigo_prepared.Rda"))
   param_list=PARAMETERS$EXIGO$FK49_Exigo_Comprehensive_Panel
  
   }else if(ExpID == "FK46") {
     load(file = file.path(PATHS$exigo$FK46_input,  "FK46_Exigo_prepared.Rda"))
     param_list=  param_list=PARAMETERS$EXIGO$FK46_Exigo_Liver_Panel
   }else{print("Give me an exisiting Experiment ID to load the correct data from the correct path.")
     }

analyze_exigo_parameter <- function(inputdata,value,batch = "ALL",reference_batch = NULL) {
 #in preprocessing i generated colms for censoring and direction 
 ## to not write a list for their names it is easier to give the names to the function this way 
  censored_col <- paste0(value, "_censored")
  direction_col <- paste0(value, "_direction")
  d <- inputdata %>%filter(!is.na(.data[[value]]))
  
  if (!identical(batch, "ALL")) {
    d <- d %>%filter(BATCH %in% batch)
  }
  
  #Censoring information -----
  d <- d %>%
    mutate( value_numeric = as.numeric(.data[[value]]),
            censored = as.character(.data[[censored_col]]),
            direction = as.character(.data[[direction_col]]),
            cens_logical = censored == "TRUE" )
  
  #summary/basic stats
  n_total <- nrow(d)
  n_ctrl <- sum(d$Treatment == "Ctrl",na.rm = TRUE)
  n_tam <- sum(  d$Treatment == "TAM",  na.rm = TRUE)
  n_female <- sum(d$Sex == "female",na.rm = TRUE)
  n_male <- sum(d$Sex == "male",na.rm = TRUE)
  n_censored <- sum(d$cens_logical, na.rm = TRUE)
  n_below_lod <- sum( d$cens_logical & d$direction == "<",na.rm = TRUE)
  n_above_uloq <- sum(d$cens_logical & d$direction == ">",na.rm = TRUE)
  
  # set up results table
  result <- tibble(
    parameter = value,
    method = NA_character_,
    n_total = n_total,
    n_ctrl = n_ctrl,
    n_tam = n_tam,
    n_female = n_female,
    n_male = n_male,
    n_censored = n_censored,
    n_below_LOD = n_below_lod,
    n_above_ULOQ = n_above_uloq,
    mean_ctrl = NA_real_,
    sd_ctrl = NA_real_,
    mean_tam = NA_real_,
    sd_tam = NA_real_,
    p_treatment = NA_real_,
    p_adj = NA_real_,
    p_sex = NA_real_,
    p_interaction = NA_real_,
    p_raw_female = NA_real_,
    p_raw_male = NA_real_,
    treatment_effect_female = NA_real_,
    CI_low_female = NA_real_,
    CI_high_female = NA_real_,
    treatment_effect_male = NA_real_,
    CI_low_male= NA_real_,
    CI_high_male = NA_real_
  )
  #Statistics for censored data -----
  if (n_censored > 0) {
    uncensored <- d %>%filter(!cens_logical) #the single values that got FALSE in cens logical are not censored
    sufficient <- nrow(uncensored) >= 4 & # if more or equal to 4 values are uncensored test can be run
      all(uncensored %>% group_by(Treatment) %>%
            summarise(n_dist = n_distinct(value_numeric), .groups = "drop") %>%
            pull(n_dist) >= 2)                 # 2 values per treatment
    
    if (sufficient) {
      cen_result <- suppressWarnings(
        with(d,cen2means( value_numeric,cens_logical, group = Treatment)))
      result$method <- "cen2means"
      result$p_treatment <- cen_result$pval
    } else {
      result$method <- "cen2means_not_performed"
      }
    
  } else {
    # Statistics for uncensored data -----
    model <- lm(value_numeric ~ Treatment * Sex, data = d)
    # ANOVA
    anova_result <- anova(model)
    p_treatment <- anova_result["Treatment", "Pr(>F)"]
    p_sex <- anova_result["Sex", "Pr(>F)"]
    p_interaction <- anova_result["Treatment:Sex", "Pr(>F)"]
    # Treatment effect separately within each Sex
    emm <- emmeans(model, ~ Treatment | Sex)
    contrast_result <- contrast(emm,method = "revpairwise")
    contrast_summary <- summary(contrast_result,infer = TRUE)
    result$method <- "linear_model_Treatment_x_Sex"
    result$p_treatment <- anova_result["Treatment", "Pr(>F)"]
    result$p_sex <- anova_result["Sex", "Pr(>F)"]
    result$p_interaction <- anova_result["Treatment:Sex", "Pr(>F)"]  
    result$p_raw_female <- contrast_summary$p.value[contrast_summary$Sex == "female"]
    result$p_raw_male <- contrast_summary$p.value[contrast_summary$Sex == "male"]
    result$treatment_effect_female <- contrast_summary$estimate[contrast_summary$Sex == "female"]
    result$CI_low_female <- contrast_summary$lower.CL[contrast_summary$Sex == "female"]
    result$CI_high_female <- contrast_summary$upper.CL[ contrast_summary$Sex == "female"]
    result$treatment_effect_male <- contrast_summary$estimate[contrast_summary$Sex == "male"]
    result$CI_low_male <- contrast_summary$lower.CL[contrast_summary$Sex == "male"]
    result$CI_high_male <- contrast_summary$upper.CL[ contrast_summary$Sex == "male"]
    result$mean_ctrl <- mean( d$value_numeric[d$Treatment == "Ctrl"],  na.rm = TRUE  )
    result$sd_ctrl <- sd( d$value_numeric[d$Treatment == "Ctrl"],  na.rm = TRUE )
    result$mean_tam <- mean(d$value_numeric[d$Treatment == "TAM"],na.rm = TRUE)
    result$sd_tam <- sd(d$value_numeric[d$Treatment == "TAM"],na.rm = TRUE
    )
  }
  result
}

if (ExpID=="FK49") {
}  else if (ExpID == "FK46"){
} else{
    print("I don't know which Exigo you used.")}


results_list <- lapply(param_list,  function(x) {analyze_exigo_parameter(inputdata = d1, value = x$value, batch = "ALL" ) }) # perform the function over the paramter list and put it in list
StatsOutput <- bind_rows(results_list) # make tibble out of the list 
StatsOutput <- StatsOutput %>%mutate(p_adj = p.adjust(p_treatment,method = "fdr")) #adjust pvalues
write.csv2( StatsOutput,file = file.path(output_pwd, "FK49_Exigo_Statistics.csv"),row.names = FALSE)




results_list <- lapply(param_list,  function(x) {analyze_exigo_parameter(inputdata = baseline_data, value = x$value, batch = "ALL" ) }) # perform the function over the paramter list and put it in list
StatsOutput <- bind_rows(results_list) # make tibble out of the list 
StatsOutput <- StatsOutput %>%mutate(p_adj = p.adjust(p_treatment,method = "fdr")) #adjust pvalues
write.csv2( StatsOutput,file = file.path(output_pwd, "BH15_Exigo_Statistics.csv"),row.names = FALSE)

rm=(list=ls())
gc()