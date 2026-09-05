###############################################################################
# FK49_Proteomics_2_Statistics.R
#
# Statistical analysis of FK49 proteomics data using limma.
###############################################################################

rm(list = ls())
gc()

library(readxl)
library(dplyr)
library(tidyr)
library(tibble)
library(limma)

source("FK49_Definitions.R")

proteom_input_pwd <- PATHS$proteomics$input
proteom_output_pwd <- PATHS$proteomics$output
stats_pwd <- file.path(proteom_output_pwd, "Statistics")


# Data ------------------------------------------------------------------------

Proteins <- read_excel( file.path(proteom_input_pwd, "FK49_TAM_EtOH_STATISTICAL_results.xlsx"), sheet = 1) %>%
  mutate(  Protein.Names = gsub("_MOUSE", "", Protein.Names)) %>%
  dplyr::rename(
    Name = Protein.Names,
    logFC = logFC...6,
    pValue = pValue...7,
    adj_pvalue = BH_corrected_pvalue...8,
    Subset_1 = Subset...9,
    Direction = Regulation...10,
    logFC_2 = logFC...33,
    pValue_2 = pValue...34,
    BH_corrected_pvalue_2 = BH_corrected_pvalue...35,
    Subset_2 = Subset...36,
    Regulation_2 = Regulation...37) %>%
    dplyr::select( -Regulation_2, -logFC_2, -pValue_2, -BH_corrected_pvalue_2,
    -Subset_2, -pValue, -adj_pvalue, -Direction, -Subset_1)
colnames(Proteins) <- gsub("EtOH", "Ctrl", colnames(Proteins))

meta <- read_excel( file.path(proteom_input_pwd, "FK49_FK46_BH_SampleList_Proteomics_Shubam.xlsx")) %>%
  filter( ExpID == "FK49",  Batch == "2") %>%
  dplyr::select( Animal,Sex, Treatment,PotentialShubhamNames )%>%
  mutate( Treatment = gsub("EtOH", "Ctrl", Treatment),
    PotentialShubhamNames = gsub("EtOH", "Ctrl", PotentialShubhamNames) )

protein_matrix <- Proteins %>%
  dplyr::select(
    F_TAM_1, F_TAM_2,
    M_Ctrl_1, M_Ctrl_2, M_Ctrl_3, M_Ctrl_4, M_Ctrl_5,
    F_TAM_3, F_TAM_4, F_TAM_5,
    M_TAM_1, M_TAM_2, M_TAM_3,
    F_Ctrl_1, F_Ctrl_2, F_Ctrl_3, F_Ctrl_4,
    M_TAM_4) %>%
  dplyr::mutate(across(everything(), as.numeric)) %>%
  as.matrix()

rownames(protein_matrix) <- Proteins$Name

meta <- meta[ match(colnames(protein_matrix),  meta$PotentialShubhamNames),]

meta$Sex <- factor(meta$Sex)
meta$Treatment <- factor(meta$Treatment)


# limma model -----------------------------------------------------------------
design <- model.matrix( ~ Treatment * Sex,  data = meta)
fit <- lmFit( protein_matrix, design)
fit <- eBayes(fit)
coef_names <- colnames(design)
treatment_coef <- which(coef_names == "TreatmentTAM")
interaction_coef <- which(coef_names == "TreatmentTAM:Sexmale")


# Contrasts -------------------------------------------------------------------

C_overall <- rep(0, ncol(design))
C_overall[treatment_coef] <- 1
C_overall[interaction_coef] <- 0.5

C_interaction <- rep(0, ncol(design))
C_interaction[interaction_coef] <- 1

C_female <- rep(0, ncol(design))
C_female[treatment_coef] <- 1

C_male <- rep(0, ncol(design))
C_male[treatment_coef] <- 1
C_male[interaction_coef] <- 1

fit_contrasts <- contrasts.fit(
  fit,
  contrasts = cbind(
    Treatment_overall = C_overall,
    Treatment_x_Sex = C_interaction,
    Treatment_female = C_female,
    Treatment_male = C_male
  )
)

fit_contrasts <- eBayes(fit_contrasts)


# Results ---------------------------------------------------------------------

results_Treatment <- topTable(fit_contrasts,coef = "Treatment_overall",  number = Inf,adjust.method = "BH") %>%
  rownames_to_column("Name") %>%
 dplyr::select(Name, logFC, P.Value, adj.P.Val) %>%
  dplyr::rename( logFC_Treatment = logFC,
    pValue_Treatment = P.Value,
    adj_pvalue_Treatment = adj.P.Val)

results_Interaction <- topTable(fit_contrasts,coef = "Treatment_x_Sex", number = Inf, adjust.method = "BH") %>%
  rownames_to_column("Name") %>%
  dplyr::select(Name, logFC, P.Value, adj.P.Val) %>%
  dplyr::rename( logFC_Treatment_x_Sex = logFC,
    pValue_Treatment_x_Sex = P.Value,
    adj_pvalue_Treatment_x_Sex = adj.P.Val)

results_Female <- topTable( fit_contrasts,coef = "Treatment_female",number = Inf, adjust.method = "BH") %>%
  rownames_to_column("Name") %>%
  dplyr::select(Name, logFC, P.Value, adj.P.Val) %>%
  dplyr::rename( logFC_Treatment_female = logFC,
    pValue_Treatment_female = P.Value,
    adj_pvalue_Treatment_female = adj.P.Val )

results_Male <- topTable( fit_contrasts,coef = "Treatment_male",number = Inf, adjust.method = "BH") %>%
  rownames_to_column("Name") %>%
  dplyr::select(Name, logFC, P.Value, adj.P.Val) %>%
  dplyr::rename(logFC_Treatment_male = logFC,pValue_Treatment_male = P.Value,adj_pvalue_Treatment_male = adj.P.Val )


# Save individual statistical results -----------------------------------------

write.csv2( results_Treatment,file.path(proteom_output_pwd, "Statistics/01_LIMMA_Treatment_overall.csv"), row.names = FALSE)
write.csv2( results_Interaction,file.path(proteom_output_pwd, "Statistics/01_LIMMA_Treatment_x_Sex.csv"), row.names = FALSE)
write.csv2( results_Female, file.path(proteom_output_pwd, "Statistics/01_LIMMA_Treatment_female.csv"),row.names = FALSE)
write.csv2( results_Male, file.path(proteom_output_pwd, "Statistics/01_LIMMA__Treatment_male.csv"),row.names = FALSE)


# Combined results ------------------------------------------------------------

Proteins_withstats <- Proteins %>%
  left_join(results_Treatment, by = "Name") %>%
  left_join(results_Interaction, by = "Name") %>%
  left_join(results_Female, by = "Name") %>%
  left_join(results_Male, by = "Name") %>%
  mutate( Direction = case_when(adj_pvalue_Treatment < 0.05 & logFC_Treatment > 1 ~ "Up",
                                adj_pvalue_Treatment < 0.05 & logFC_Treatment < -1 ~ "Down",
                                TRUE ~ "NS" ),
          Direction_Male = case_when(adj_pvalue_Treatment_male < 0.05 &logFC_Treatment_male > 1 ~ "Up",
                                     adj_pvalue_Treatment_male < 0.05 &logFC_Treatment_male < -1 ~ "Down",
                                     TRUE ~ "NS"),
          Direction_Female = case_when(adj_pvalue_Treatment_female < 0.05 &logFC_Treatment_female > 1 ~ "Up",
                                 adj_pvalue_Treatment_female < 0.05 &logFC_Treatment_female < -1 ~ "Down",
                                 TRUE ~ "NS"),
          Direction_Interaction = case_when(adj_pvalue_Treatment_x_Sex < 0.05 &abs(logFC_Treatment_x_Sex) > 1 ~ "Different",
                                            TRUE ~ "NS")
  )


# Save combined results and matrix --------------------------------------------

write.csv2( Proteins_withstats, file.path(proteom_output_pwd, "Statistics/02_LIMMA_combined_stats.csv"),  row.names = FALSE)
saveRDS(Proteins_withstats, file.path(proteom_output_pwd, "Statistics/02_LIMMA_combined_stats.rds"))

saveRDS(protein_matrix, file.path(proteom_output_pwd, "Data/01_protein_matrix.rds"))
saveRDS(meta, file.path(proteom_output_pwd, "Data/01_metadata.rds") )

### Calculate Correlation between Animals ------------------------------------
### Calculate and save correlations ------------------------------------------

protein_matrix_centered <- t(
  scale(t(protein_matrix), center = TRUE, scale = FALSE)
)

animal_cor_spearman_centered <- cor(
  protein_matrix_centered,
  method = "spearman",
  use = "pairwise.complete.obs"
)

animal_cor_spearman <- cor(
  protein_matrix,
  method = "spearman",
  use = "pairwise.complete.obs"
)

animal_cor_pearson <- cor(
  protein_matrix,
  method = "pearson",
  use = "pairwise.complete.obs"
)

animal_cor_pearson_centered <- cor(
  protein_matrix_centered,
  method = "pearson",
  use = "pairwise.complete.obs"
)


### Significant proteins -----------------------------------------------------

Proteins_for_correlation <- Proteins_withstats %>%
  filter(
    adj_pvalue_Treatment < 0.05,
    abs(logFC_Treatment) > 1
  ) %>%
  dplyr::select(
    Name,
    F_TAM_1, F_TAM_2, M_Ctrl_1, M_Ctrl_2,
    M_Ctrl_3, M_Ctrl_4, M_Ctrl_5,
    F_TAM_3, F_TAM_4, F_TAM_5,
    M_TAM_1, M_TAM_2, M_TAM_3,
    F_Ctrl_1, F_Ctrl_2, F_Ctrl_3,
    F_Ctrl_4, M_TAM_4
  )%>%as.data.frame()

rownames(Proteins_for_correlation) <- Proteins_for_correlation$Name

Proteins_for_correlation <- Proteins_for_correlation %>%
  dplyr::select(-Name) %>%
  as.matrix() %>%
  t()

proteins_for_correlation_centered <- t(scale(Proteins_for_correlation,center = TRUE,scale = FALSE ))

protein_cor_spearman <- cor(Proteins_for_correlation,  method = "spearman",  use = "pairwise.complete.obs")

animal_sig_cor_spearman <- cor(  t(Proteins_for_correlation),  method = "spearman",  use = "pairwise.complete.obs")

animal_sig_cor_spearman_centered <- cor(  proteins_for_correlation_centered,  method = "spearman",  use = "pairwise.complete.obs")

animal_sig_cor_pearson <- cor(  t(Proteins_for_correlation),  method = "pearson",  use = "pairwise.complete.obs")

animal_sig_cor_pearson_centered <- cor( proteins_for_correlation_centered,method = "pearson", use = "pairwise.complete.obs")


### Save correlation matrices -----------------------------------------------

saveRDS(
  list(
    animal_cor_spearman = animal_cor_spearman,
    animal_cor_spearman_centered = animal_cor_spearman_centered,
    animal_cor_pearson = animal_cor_pearson,
    animal_cor_pearson_centered = animal_cor_pearson_centered,
    protein_cor_spearman = protein_cor_spearman,
    animal_sig_cor_spearman = animal_sig_cor_spearman,
    animal_sig_cor_spearman_centered = animal_sig_cor_spearman_centered,
    animal_sig_cor_pearson = animal_sig_cor_pearson,
    animal_sig_cor_pearson_centered = animal_sig_cor_pearson_centered
  ),
  file.path(proteom_output_pwd, "Statistics", "03_Correlation_matrices.rds")
)

write.csv2( animal_cor_spearman,file.path(proteom_output_pwd, "Statistics", "03_Correlation_Animal_Spearman.csv"))
write.csv2( animal_cor_spearman_centered, file.path(proteom_output_pwd, "Statistics", "03_Correlation_Animal_Spearman_centered.csv"))
write.csv2( animal_cor_pearson, file.path(proteom_output_pwd, "Statistics", "03_Correlation_Animal_Pearson.csv"))
write.csv2(  animal_cor_pearson_centered,file.path(proteom_output_pwd, "Statistics", "03_Correlation_Animal_Pearson_centered.csv"))
write.csv2( protein_cor_spearman, file.path(proteom_output_pwd, "Statistics", "03_Correlation_Protein_Spearman.csv"))
write.csv2(  animal_sig_cor_spearman,  file.path(proteom_output_pwd, "Statistics", "03_Correlation_Animal_SignificantProteins_Spearman.csv"))
write.csv2( animal_sig_cor_spearman_centered, file.path(proteom_output_pwd, "Statistics", "03_Correlation_Animal_SignificantProteins_Spearman_centered.csv"))
write.csv2( animal_sig_cor_pearson,  file.path(proteom_output_pwd, "Statistics", "03_Correlation_Animal_SignificantProteins_Pearson.csv"))
write.csv2(  animal_sig_cor_pearson_centered,  file.path(proteom_output_pwd, "Statistics", "03_Correlation_Animal_SignificantProteins_Pearson_centered.csv"))