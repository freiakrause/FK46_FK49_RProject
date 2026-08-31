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
library(waffle)
library(rstatix)
library(lmerTest)
library(emmeans)
library(grid)
library(ggnewscale)
library(NADA2)
library(effsize)
# Read Raw Inputdata after general Data manipulation ------------------------------------------------------
setwd("C:/Users/b1084855/OneDrive - Universität Salzburg/Freigegebene Dokumente - AG_Tumorimmunologie/Data/Freia Krause/01_Experiments/FK46_iALmice_high Fat diet 52 weeks 7d after injection/Analysis/")
load("01_RawData/FK46_Data_prepared.Rda")
#Plot Survival Function --------------------------------------------------------
do_surv_plot <- function(inputdata, sex = "both"){
  d <-inputdata %>%
    select(Animal, KILL.DATE,wks_diet,START.Diet, Sex, Treatment, Death, BATCH) %>%
    filter(case_when(
      sex == "female" ~ Sex == "female",
      sex == "male" ~ Sex == "male",
      sex == "both" ~ TRUE))
  
  surv_df <- d %>%
    group_by(Animal) %>%
    summarise(
      Sex = first(Sex),
      Treatment = first(Treatment),
      time = as.numeric(difftime(first(KILL.DATE), first(START.Diet), units = "weeks")),
      Death = first(Death))
  
  
  
  fit <- survfit(Surv(time, Death) ~ Treatment, data = surv_df)      # Create survival obectwith Surv() and Fit the survival model
  survdiff(Surv(time, Death) ~ Treatment, data = surv_df)
  
  # Plot using ggsurvfit from survminer (not ggplot2)
  survival <-ggsurvplot(fit, 
                        data = surv_df, 
                        type = "survival", # Kaplan-Meier survival curve,
                        surv.median.line = "none",
                        risk.table = FALSE,# Include risk table
                        tables.height = 0.1,
                        tables.theme = theme_cleantable(),
                        pval = TRUE, # Show p-value for the log-rank test
                        pval.coord= c(15,0.35),
                        pval.size = 4,
                        pval.method = TRUE,
                        pval.method.size =4,
                        pval.method.coord =c(15,0.4) ,
                        conf.int = TRUE,  # Include confidence intervals
                        conf.int.alpha = 0.1,
                        palette = c("grey30", "palevioletred4"),  # Customize colors for Treatment
                        legend.title = "Treatment",  # Legend title
                        legend.labs = c("Ctrl", "TAM"), # Custom labels for the legend
                        censor.shape=124,
                        cumevents = T,
                        cumcensor = T,
                        fontsize = 4,
                        ncensor.plot = FALSE,
                        font.main = c(12, "bold", "black"),
                        font.x = c(12, "bold", "black"),
                        font.y = c(12, "bold", "black"),
                        font.tickslab = c(10, "plain", "black"))
  
  survival$plot<-survival$plot+
    scale_x_continuous(name = "Time on CD-HFD [wks]", 
                       limits = c(0,41),  
                       breaks = seq(0,40,4),  
                       minor_breaks = seq(0, 41, by = 1))+
    guides(x = guide_axis(cap = "upper",minor.ticks = TRUE), y = guide_axis(cap = "upper"))+
    ggtitle(paste0("Survival Curve of ",sex," iAL mice on CD-HFD"))
  
  print(survival)
  
  ggsave(filename = paste0("FK46_Survival_",sex,".png"),plot = survival$plot, path = "02_GeneratedData", 
         width = 8, height= 5, dpi=300)
  rm(fit,survival,d)
  gc()
}

do_surv_plot(data,"female")
do_surv_plot(data,"male")
do_surv_plot(data,"both")
