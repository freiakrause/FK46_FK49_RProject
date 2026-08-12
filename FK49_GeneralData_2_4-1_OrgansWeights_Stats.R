rm(list=ls())
gc()
library(knitr)
library(car)
library(dplyr)
library(broom)  # for tidy() function
library(tidyr)
library(dplyr)
library(survival)
library(lubridate)
library(tidyverse)
library(ggsurvfit)
library(survminer)
library(patchwork)
library(superb)
library(ggbreak)
library(tibble)
library(rstatix)
library(lmerTest)
library(emmeans)
library(grid)
library(ggnewscale)
library(NADA2)
library(effsize)
source("FK49_Definitions.R")
# Read Raw Inputdata after general Data manipulation ------------------------------------------------------
output_pwd <- file.path(PATHS$Organs$output)

load(file.path(PATHS$Organs$input, "FK49_Data_prepared.Rda"))

#Stat Analysis Organ Weight-----
analyze_organ_weight <- function(inputdata, value) {
  
  d <- inputdata %>%filter(complete.cases(.data[[value]])) %>%
    mutate( Sex = factor(Sex),Treatment = factor(Treatment),BATCH = factor(BATCH))
  d[[value]] <- as.numeric(d[[value]])
  
## Linear model -----
  formula <- as.formula(paste(value, "~ Treatment * Sex + BATCH")) # Is there T effect, is there sex efect, isthere tretment and sex interaction, is there batch effect
  model <- lm(formula, data = d)
  anova_table <- anova(model)
  p_treatment  <- anova_table["Treatment", "Pr(>F)"]
  p_sex        <- anova_table["Sex", "Pr(>F)"]
  p_interaction <- anova_table["Treatment:Sex", "Pr(>F)"]
  p_batch      <- anova_table["BATCH", "Pr(>F)"]
  
  if (!is.na(p_interaction) && p_interaction < 0.05) {emm <- emmeans(model, ~ Treatment | Sex) #diff in treatment effect by sex, therefor mmeans per sex sex
   } else {  emm <- emmeans(model, ~ Treatment) # no diff in treatment effect by sex, therefor mmeans averaged over sex
   }
  contrasts <- contrast(emm,method = "pairwise",adjust = "none")%>%as.data.frame()
  #
  list(
    model = model,
    anova = anova_table,
    contrasts = contrasts,
    p_treatment = p_treatment,
    p_sex = p_sex,
    p_interaction = p_interaction,
    p_batch = p_batch
  )
}



absolute_organ_variables <- c("Liver","Spleen", "Fat")
rel_organ_variables <- c("Liver_rel","Spleen_rel","Fat_rel")

absolute_results <- lapply( absolute_organ_variables, function(x) analyze_organ_weight(data, x))
rel_results <- lapply(rel_organ_variables,  function(x) analyze_organ_weight(data, x))
names(absolute_results) <- absolute_organ_variables
names(rel_results) <- rel_organ_variables

# make absoulte results into df and add interpretation -----
absolute_stats <- map2_dfr( absolute_results, names(absolute_results),
                            ~ tibble(Variable = .y,
                                     #model = .x$model,
                                     anova = .x$anova_table,
                                     contrasts = .x$contrasts,
                                     p_Treatment = .x$p_treatment,
                                     p_Sex = .x$p_sex,
                                     p_Interaction = .x$p_interaction,
                                     p_Batch = .x$p_batch ))%>%
  mutate(p_Treatment_adj = p.adjust( p_Treatment,method = "fdr"),
         p_Sex_adj = p.adjust( p_Sex,method = "fdr"),
         p_Interaction_adj = p.adjust(p_Interaction,method = "fdr"),
         p_Batch_adj = p.adjust(p_Batch,method = "fdr"),
         Type = "Absolute")%>%
  mutate(
    Interpretation_Treatment = case_when(
    is.na(p_Treatment_adj) ~ NA_character_, p_Treatment_adj < 0.05 ~ "Significant",TRUE ~ "Not significant"),
    Interpretation_Sex = case_when(
      is.na(p_Sex_adj) ~ NA_character_,p_Sex_adj < 0.05 ~ "Significant", TRUE ~ "Not significant"),
    Interpretation_Interaction = case_when(
      is.na(p_Interaction_adj) ~ NA_character_,p_Interaction_adj < 0.05 ~ "Significant",TRUE ~ "Not significant" ),
    Interpretation_Batch = case_when(
      is.na(p_Batch_adj) ~ NA_character_,p_Batch_adj < 0.05 ~ "Significant",TRUE ~ "Not significant")
         ) %>%
  rowwise()%>%
  mutate(Overall_Interpretation = paste0(
      if (is.na(p_Treatment_adj)) {"Treatment effect could not be determined"
      } else if (p_Treatment_adj < 0.05) {"Organ weight differs between treatments"
      } else {"Organ weight does not differ between treatments"
      },
      ", ",
      if (is.na(p_Sex_adj)) {"sex effect could not be determined"
      } else if (p_Sex_adj < 0.05) {"differs between sexes"
      } else {"does not differ between sexes"
      },
      ", and ",
      if (is.na(p_Interaction_adj)) {"the treatment-by-sex interaction could not be determined"
      } else if (p_Interaction_adj < 0.05) {"the treatment effect differs between sexes"
      } else {"the treatment effect does not differ between sexes"
      },
      ". ",
      if (is.na(p_Batch_adj)) {"The batch effect could not be determined."
      } else if (p_Batch_adj < 0.05) {"A significant batch effect was observed."
      } else {"No significant batch effect was observed."
      }
    ) ) %>%
  ungroup()

rel_stats <- map2_dfr(rel_results,names(rel_results),
  ~ tibble(
    Variable = .y,
    #model = .x$model,
    anova = .x$anova_table,
    contrasts = .x$contrasts,
    p_Treatment = .x$p_treatment,
    p_Sex = .x$p_sex,
    p_Interaction = .x$p_interaction,
    p_Batch = .x$p_batch))%>%
  mutate(
    p_Treatment_adj = p.adjust(p_Treatment,method = "fdr"),
    p_Sex_adj = p.adjust(p_Sex,method = "fdr"),
    p_Interaction_adj = p.adjust(p_Interaction,method = "fdr"),
    p_Batch_adj = p.adjust(p_Batch,method = "fdr"),
    Type = "Relative" )%>%
  mutate( Interpretation_Treatment = case_when(is.na(p_Treatment_adj) ~ NA_character_,
                                               p_Treatment_adj < 0.05 ~ "Significant",
                                               TRUE ~ "Not significant"),
          
          Interpretation_Sex = case_when(is.na(p_Sex_adj) ~ NA_character_,
                                         p_Sex_adj < 0.05 ~ "Significant",
                                         TRUE ~ "Not significant"),
          Interpretation_Interaction = case_when(is.na(p_Interaction_adj) ~ NA_character_,
                                                 p_Interaction_adj < 0.05 ~ "Significant",
                                                 TRUE ~ "Not significant" ),
          
          Interpretation_Batch = case_when(is.na(p_Batch_adj) ~ NA_character_,
                                           p_Batch_adj < 0.05 ~ "Significant",
                                           TRUE ~ "Not significant")
  ) %>%
  
  rowwise() %>%
  mutate(Overall_Interpretation = paste0(
      if (is.na(p_Treatment_adj)) {"Organ weight could not be evaluated for treatment"
      } else if (p_Treatment_adj < 0.05) {"Organ weight differs between treatments"
      } else {"Organ weight does not differ between treatments"
      },
      ", ",
      
      if (is.na(p_Sex_adj)) {"sex effect could not be determined"
      } else if (p_Sex_adj < 0.05) {"differs between sexes"
      } else {"does not differ between sexes"
      },
      ", and ",
      if (is.na(p_Interaction_adj)) {"the treatment-by-sex interaction could not be determined"
      } else if (p_Interaction_adj < 0.05) {"the treatment effect differs between sexes"
      } else {"the treatment effect does not differ between sexes"
      },
      ". ",
      if (is.na(p_Batch_adj)) {"The batch effect could not be determined."
      } else if (p_Batch_adj < 0.05) {"A significant batch effect was observed."
      } else {"No significant batch effect was observed."
      }
    ) ) %>%
  ungroup()


organ_stats <- bind_rows(absolute_stats,rel_stats)
organ_stats <- organ_stats %>%
  select(Variable,Type,p_Treatment,p_Treatment_adj,Interpretation_Treatment,
        p_Sex,p_Sex_adj,Interpretation_Sex, p_Interaction,p_Interaction_adj,Interpretation_Interaction,
        p_Batch,p_Batch_adj,Interpretation_Batch, Overall_Interpretation,contrasts)


write.csv2( organ_stats,file = file.path(output_pwd,"Statistics/FK49_Organ_Weight_Statistics.csv"),row.names = FALSE)
saveRDS(organ_stats,file = file.path(output_pwd,"Statistics/FK49_Organ_Weight_Statistics.rds"))


           