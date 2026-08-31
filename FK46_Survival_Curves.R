library(tidyverse)
library(survminer)
library(survival)

source("FK49_Definitions.R")
ExpId <- "FK46"

load(file.path(PATHS$general_data[[paste0(ExpId, "_output")]], "01_RawData",
               paste0(ExpId, "_Data_prepared.Rda")))

do_surv_plot <- function(inputdata, sex = "both") {
  
  d <- inputdata %>%
    select(Animal, KILL.DATE, START.Diet, Sex, Treatment, Death) %>%
    filter(case_when(
      sex == "female" ~ Sex == "female",
      sex == "male" ~ Sex == "male",
      sex == "both" ~ TRUE
    ))
  
  surv_df <- d %>%
    group_by(Animal) %>%
    summarise(
      Sex = first(Sex),
      Treatment = first(Treatment),
      time = as.numeric(difftime(first(KILL.DATE), first(START.Diet), units = "weeks")),
      Death = first(Death),
      .groups = "drop"
    )
  
  # Kaplan-Meier
  fit <- survfit(Surv(time, Death) ~ Treatment, data = surv_df)
  
  # Log-rank test
  logrank <- survdiff(Surv(time, Death) ~ Treatment, data = surv_df)
  logrank_p <- 1 - pchisq(logrank$chisq, length(logrank$n) - 1)
  
  # Cox model
  cox <- coxph(Surv(time, Death) ~ Treatment, data = surv_df)
  cs <- summary(cox)
  HR <- exp(coef(cox))
  HR_CI <- exp(confint(cox))
  HR_p <- cs$coefficients[, "Pr(>|z|)"]
  
  # PH assumption
  ph_test <- cox.zph(cox)
  
  # Median survival
  median_survival <- surv_median(fit)
  
  # Survival probabilities
  s <- summary(fit, times = c(20, 30, 40), extend = TRUE)
  survival_at_times <- data.frame(
    Treatment = s$strata, Time_weeks = s$time, Survival = s$surv,
    CI_lower = s$lower, CI_upper = s$upper
  )
  
  # Statistical summary
  tab <- table(surv_df$Treatment, surv_df$Death)
  med <- setNames(median_survival$median, gsub("Treatment=", "", median_survival$strata))
  
  statistics <- data.frame(
    Sex = sex,
    n_Ctrl = sum(surv_df$Treatment == "Ctrl"),
    n_TAM = sum(surv_df$Treatment == "TAM"),
    Deaths_Ctrl = ifelse("Ctrl" %in% rownames(tab) && "1" %in% colnames(tab), tab["Ctrl","1"], 0),
    Deaths_TAM = ifelse("TAM" %in% rownames(tab) && "1" %in% colnames(tab), tab["TAM","1"], 0),
    Censored_Ctrl = ifelse("Ctrl" %in% rownames(tab) && "0" %in% colnames(tab), tab["Ctrl","0"], 0),
    Censored_TAM = ifelse("TAM" %in% rownames(tab) && "0" %in% colnames(tab), tab["TAM","0"], 0),
    Logrank_chisq = logrank$chisq,
    Logrank_df = length(logrank$n) - 1,
    Logrank_p = logrank_p,
    HR_TAM_vs_Ctrl = unname(HR),
    HR_CI_lower = unname(HR_CI[1]),
    HR_CI_upper = unname(HR_CI[2]),
    Cox_p = unname(HR_p),
    PH_test_p = ph_test$table[1, "p"],
    Median_survival_Ctrl = ifelse("Ctrl" %in% names(med), med["Ctrl"], NA),
    Median_survival_TAM = ifelse("TAM" %in% names(med), med["TAM"], NA)
  )
  
  # Save statistics
  write.csv2(
    statistics,
    file.path(PATHS$general_data[[paste0(ExpId, "_output")]],
              "02_GeneratedData",
              paste0(ExpId, "_Survival_", sex, "_Statistics.csv")),
    row.names = FALSE
  )
  
  # Plot
  survival <- ggsurvplot(
    fit, data = surv_df, type = "survival",
    surv.median.line = "none", risk.table = FALSE,
    pval = TRUE, pval.coord = c(15, 0.35), pval.size = 4,
    pval.method = TRUE, pval.method.size = 4,
    pval.method.coord = c(15, 0.4),
    conf.int = TRUE, conf.int.alpha = 0.1,
    palette = unname(Treatment_colors[c("Ctrl", "TAM")]),
    legend.title = "Treatment", legend.labs = c("Ctrl", "TAM"),
    censor.shape = 124, cumevents = TRUE, cumcensor = TRUE,
    fontsize = 4, ncensor.plot = FALSE,
    font.main = c(12, "bold", "black"),
    font.x = c(12, "bold", "black"),
    font.y = c(12, "bold", "black"),
    font.tickslab = c(10, "plain", "black")
  )
  
  survival$plot <- survival$plot +
    scale_x_continuous(
      name = "Time on CD-HFD [wks]", limits = c(0, 41),
      breaks = seq(0, 40, 4), minor_breaks = seq(0, 41, 1)
    ) +
    guides(x = guide_axis(cap = "upper", minor.ticks = TRUE),
           y = guide_axis(cap = "upper")) +
    ggtitle(paste0("Survival Curve of ", sex, " iAL mice on CD-HFD"))
  
  print(survival)
  
  ggsave(
    filename = paste0(ExpId, "_Survival_", sex, ".png"),
    plot = survival$plot,
    path = file.path(PATHS$general_data[[paste0(ExpId, "_output")]],
                     "02_GeneratedData"),
    width = 8, height = 5, dpi = 300
  )
  
  return(list(
    data = surv_df, KM = fit, logrank = logrank, logrank_p = logrank_p,
    Cox = cox, HR = HR, HR_CI = HR_CI, HR_p = HR_p,
    PH_test = ph_test, median_survival = median_survival,
    survival_at_times = survival_at_times, statistics = statistics
  ))
}

stats_female <- do_surv_plot(data, "female")
stats_male <- do_surv_plot(data, "male")
stats_both <- do_surv_plot(data, "both")