rm(list=ls())
gc()
library(tidyverse)
library(NADA2)
library(emmeans)
source("FK49_Definitions.R")
ExpID= "FK49"   # Decide if you want to load data from FK46 or FK49


if(ExpID == "FK49"){
 data<- readRDS (file = file.path(dirname(dirname(PATHS$legendplex$FK49_output)),  "01_RawData/FK49_Legendplex_clean.Rds"))
  cytokine_list<-PARAMETERS$Legendplex$cytokine_list  
  output_pwd = file.path(PATHS$legendplex$FK49_output)
}else if(ExpID == "FK46") {
  # load(file = file.path(PATHS$exigo$FK46_input,  "FK46_Exigo_prepared.Rda"))
  # param_list=  param_list=PARAMETERS$EXIGO$FK46_Exigo_Liver_Panel
  # output_pwd = file.path(PATHS$exigo$FK46_output)
  }else{print("Give me an exisiting Experiment ID to load the correct data from the correct path.")
}

analyze_legendplex_parameter <- function(inputdata,value,batch = "ALL",reference_batch = NULL) {
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
      if (sufficient) {
        
        d <- d %>% mutate(Treatment = factor(Treatment, levels = c("Ctrl", "TAM")),
                          Sex = factor(Sex, levels = c("female", "male")))
        
        cen_result <- capture.output(suppressWarnings(with(d,cen2way(ifelse(cens_logical, value_numeric * 2, value_numeric),
                                                                     cens_logical, Treatment, Sex, LOG = TRUE, interact = TRUE))))
        p_lines <- cen_result[grep("Treatment|Sex|interaction", cen_result)]
        p_values <- as.numeric(sub(".*\\s([0-9]+\\.[0-9]+)$", "\\1", p_lines))
        
        result$method <- "cen2way"
        result$p_treatment <- p_values[3]
        
        #  Effect size: Geometric Mean Ratio (TAM / Ctrl)-----
        y1 <- ifelse(d$cens_logical, d$value_numeric * 2, d$value_numeric)
        y2 <- d$cens_logical
        e <- ifelse(d$Treatment == "Ctrl", 1, -1)
        s <- ifelse(d$Sex == "female", 1, -1)
        int <- e * s
        
        lnvar <- log(y1)
        fconst <- max(lnvar, na.rm = TRUE)
        flip.log <- fconst + 1 - lnvar
        detect <- !y2
        
        logCensData <- survival::Surv(flip.log, detect, type = "right")
        cen_model <- survival::survreg(logCensData ~ e + s + int, dist = "gaussian")
        
        # Reverse the flipping used by cen2way()
        beta <- -coef(cen_model)
        beta[1] <- fconst + 1 + beta[1]
        
        # log geometric means, marginal over Sex
        b0 <- beta["(Intercept)"]
        bT <- beta["e"]
        bS <- beta["s"]
        bI <- beta["int"]
        
        log_GM_ctrl <- mean(c(b0 + bT + bS + bI, b0 + bT - bS - bI))
        log_GM_tam <- mean(c(b0 - bT + bS - bI, b0 - bT - bS + bI))
        
        # geometric means
        GM_ctrl <- exp(log_GM_ctrl)
        GM_tam <- exp(log_GM_tam)
        
        log_GMR <- -2 * beta["e"]
        GMR <- exp(log_GMR)
        # 95% CI for log(GMR)
        se_b_tr <- sqrt(vcov(cen_model)["e", "e"])
        se_log_GMR <- 2 * se_b_tr
        
        # Back-transform CI to GMR scale
        CI_low <- exp(log_GMR - 1.96 * se_log_GMR)
        CI_high <- exp(log_GMR + 1.96 * se_log_GMR)
        
        # Store results
        result$effect_CI_low <- CI_low
        result$effect_CI_high <- CI_high
        
        # store results
        result$mean_ctrl <- GM_ctrl
        result$mean_tam <- GM_tam
        
        result$effect_size <- GMR
        result$effect_size_type <- "GMR"
      }
      
    } else {
      result$method <- "cen2way_not_performed"
    }
    
  } else {
    # Statistics for uncensored data -----
    model <- lm(value_numeric ~ Treatment * Sex, data = d)
    # ANOVA
    anova_result <- anova(model)
    #ANOVA effect size die model struktur berücksichtigt
    
    p_treatment <- anova_result["Treatment", "Pr(>F)"]
    p_sex <- anova_result["Sex", "Pr(>F)"]
    p_interaction <- anova_result["Treatment:Sex", "Pr(>F)"]
    # Treatment effect separately within each Sex
    emm <- emmeans(model, ~ Treatment | Sex) #gesamt model
    contrast_result <- contrast(emm,method = "revpairwise")
    contrast_summary <- summary(contrast_result, infer = TRUE)
    
    emm_treatment <- emmeans(model, ~ Treatment) #gemittelter treatment effect um effekt size zeigen zu können
    contrast_treatment <- contrast( emm_treatment,  method = "revpairwise")
    summary_treatment <- summary(  contrast_treatment,infer = TRUE  )
    eff <- eff_size(emm_treatment,sigma = sigma(model),  edf = df.residual(model) )
    eff_summary <- summary(eff)
    result$effect_size_calc <-"emmeans(model ~ Treatment) eff_size() "
    result$effect_size_type <-"standardized model effect"
    
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
    result$sd_tam <- sd(d$value_numeric[d$Treatment == "TAM"],na.rm = TRUE)
    result$effect_size <- -eff_summary$effect.size #with . so that x-1 so that i can sho it next to GMR that shoes TAM/ctrl effect sze. here it was calculcate Ctrl-TAM
    result$effect_CI_low <- -eff_summary$lower.CL
    result$effect_CI_high <- -eff_summary$upper.CL
    
  }
  result
}


results_list <- lapply(cytokine_list,  function(x) {analyze_legendplex_parameter(inputdata = data, value = x$value, batch = "ALL" ) }) # perform the function over the paramter list and put it in list
StatsOutput <- bind_rows(results_list) # make tibble out of the list 
StatsOutput <- StatsOutput %>%mutate(p_adj = p.adjust(p_treatment,method = "fdr")) #adjust pvalues
write.csv2(StatsOutput,file = file.path(output_pwd, paste0(ExpID,"_Legendplex_Statistics.csv")),row.names = FALSE)

rm=(list=ls())
gc()
