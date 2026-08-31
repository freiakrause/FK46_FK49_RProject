###############################################################################
# FK49_Microbiome_2_Statistics.R
#
# Statistical analysis of FK49 16S microbiome data.
# Loads preprocessed phyloseq objects from Script 1.
#
# Analyses:
#   2a. Alpha diversity — LMM with Sex as covariate/interaction
#   2b. Beta diversity — PERMANOVA (adonis2) + betadisper
#   2c. Differential abundance — ANCOM-BC2 (longitudinal, compositional)
#   2d. Differential abundance — ALDEx2 (cross-sectional complement)
#   2e. Consensus DA table (ANCOM-BC2 + ALDEx2)
#   2f. Targeted Lactobacillus analysis (hypothesis-driven)
#
# Sex strategy: Sex included as covariate in primary models. Sex-stratified
# post-hoc analyses run ONLY if a Sex interaction is significant.
#
# Inputs:  ps.rds, ps_genus.rds, ps_clr.rds, metadata.rds (from Script 1)
# Outputs: CSV tables in PATHS$microbiome$output_stats
###############################################################################

rm(list = ls())
gc()

library(tidyverse)
library(microbiome)
library(phyloseq)
library(lme4)
library(lmerTest)
library(emmeans)
library(car)
library(vegan)
library(ALDEx2)

# ANCOMBC requires R >= 4.5.0; load if available, otherwise skip gracefully
ancombc_available <- requireNamespace("ANCOMBC", quietly = TRUE)
if (ancombc_available) library(ANCOMBC)

source("FK49_Definitions.R")

# ============================================================
# CONFIGURATION
# ============================================================
mb_params <- PARAMETERS$microbiome
output_dir <- PATHS$microbiome$output
stats_dir  <- PATHS$microbiome$output_stats

# Ensure output folder exists
if (!dir.exists(stats_dir)) dir.create(stats_dir, recursive = TRUE)

# Load preprocessed data
ps         <- readRDS(file.path(output_dir, "ps.rds"))
ps_genus   <- readRDS(file.path(output_dir, "ps_genus.rds"))
ps_genus_rel <- readRDS(file.path(output_dir, "ps_genus_rel.rds"))
ps_clr     <- readRDS(file.path(output_dir, "ps_clr.rds"))
metadata   <- readRDS(file.path(output_dir, "metadata.rds"))

fdr_thresh <- mb_params$fdr_threshold
cat("=== FK49 Microbiome Statistics ===\n")
cat("Samples:", nsamples(ps), "| Taxa:", ntaxa(ps), "| Genera:", ntaxa(ps_genus), "\n\n")

# Helper: significance stars
sig_stars <- function(p) {
  ifelse(is.na(p), "NA",
         ifelse(p < 0.001, "***",
                ifelse(p < 0.01, "**",
                       ifelse(p < 0.05, "*", "ns"))))
}

# ============================================================
# 2a. ALPHA DIVERSITY (LMM)
# ============================================================
cat("--- 2a. Alpha Diversity ---\n")

alpha_div <- estimate_richness(ps, measures = mb_params$alpha_measures) %>%
  rownames_to_column(var = "SampleID") %>%
  mutate(SampleID = gsub("^X", "", SampleID)) %>%
  left_join(metadata, by = "SampleID")

alpha_results <- list()
sex_interaction_significant <- FALSE

for (a in mb_params$alpha_measures) {
  cat("  Fitting", a, "...\n")
  
  formula_full <- as.formula(paste0(a, " ~ Treatment * Sex * Feces + (1|Animal)"))
  formula_reduced <- as.formula(paste0(a, " ~ Treatment * Feces + Sex + (1|Animal)"))
  
  model <- tryCatch(
    lmer(formula_full, data = alpha_div),
    error = function(e) {
      cat("    Full model failed, using reduced model (Sex as covariate)\n")
      lmer(formula_reduced, data = alpha_div)
    }
  )
  
  model_type <- ifelse(as.character(formula(model))[3] ==
                         "Treatment * Sex * Feces + (1 | Animal)", "full_3way", "reduced_additive")
  
  anova_tab <- Anova(model, type = 3)
  anova_tab_df <- as.data.frame(anova_tab) %>%
    rownames_to_column("Term") %>%
    mutate(Measure = a, ModelType = model_type)
  
  sex_terms <- anova_tab_df %>%
    filter(grepl("Sex", Term) & grepl(":", Term))
  if (any(sex_terms$`Pr(>Chisq)` < fdr_thresh, na.rm = TRUE)) {
    sex_interaction_significant <- TRUE
  }
  
  emm <- emmeans(model, ~ Treatment | Feces)
  pwc <- contrast(emm, method = "pairwise", adjust = "BH")
  pwc_df <- as.data.frame(pwc) %>%
    mutate(Measure = a, Contrast = "Ctrl_vs_TAM") %>%
    mutate(significance = sig_stars(p.value))
  
  alpha_results[[paste0(a, "_anova")]] <- anova_tab_df
  alpha_results[[paste0(a, "_pwc")]] <- pwc_df
  
  png(file.path(stats_dir, paste0("alpha_diagnostics_", a, ".png")),
      width = 800, height = 400)
  par(mfrow = c(1, 2))
  plot(model, main = paste(a, "residuals vs fitted"))
  qqnorm(resid(model), main = paste(a, "QQ"))
  qqline(resid(model))
  par(mfrow = c(1, 1))
  dev.off()
}

alpha_anova_all <- bind_rows(alpha_results[grepl("_anova", names(alpha_results))])
alpha_pwc_all   <- bind_rows(alpha_results[grepl("_pwc", names(alpha_results))])

write.csv(alpha_anova_all, file.path(stats_dir, "alpha_diversity_anova.csv"), row.names = FALSE)
write.csv(alpha_pwc_all,   file.path(stats_dir, "alpha_diversity_posthoc.csv"), row.names = FALSE)

cat("  Sex interaction significant:", sex_interaction_significant, "\n")

if (sex_interaction_significant) {
  cat("  Running sex-stratified post-hoc alpha diversity...\n")
  alpha_stratified <- list()
  
  for (a in mb_params$alpha_measures) {
    for (s in c("male", "female")) {
      data_sub <- alpha_div %>% filter(Sex == s)
      model_sub <- lmer(as.formula(paste0(a, " ~ Treatment * Feces + (1|Animal)")),
                        data = data_sub)
      emm_sub <- emmeans(model_sub, ~ Treatment | Feces)
      pwc_sub <- as.data.frame(contrast(emm_sub, method = "pairwise", adjust = "BH"))
      pwc_sub$Measure <- a
      pwc_sub$Sex <- s
      pwc_sub$significance <- sig_stars(pwc_sub$p.value)
      alpha_stratified[[paste(a, s, sep = "_")]] <- pwc_sub
    }
  }
  alpha_stratified_all <- bind_rows(alpha_stratified)
  write.csv(alpha_stratified_all, file.path(stats_dir, "alpha_diversity_stratified.csv"),
            row.names = FALSE)
}

# ============================================================
# 2b. BETA DIVERSITY (PERMANOVA + BETADISPER)
# ============================================================
cat("\n--- 2b. Beta Diversity ---\n")

ps_rel_mat <- as(otu_table(ps_genus_rel), "matrix")
if (!taxa_are_rows(ps_genus_rel)) ps_rel_mat <- t(ps_rel_mat)
bc_dist <- vegdist(t(ps_rel_mat), method = "bray")

clr_mat <- as(otu_table(ps_clr), "matrix")
if (!taxa_are_rows(ps_clr)) clr_mat <- t(clr_mat)
aitchison_dist <- dist(t(clr_mat), method = "euclidean")

meta_df <- as(sample_data(ps_genus_rel), "data.frame")
meta_df$Animal <- as.factor(meta_df$Animal)

set.seed(42)
permanova_bc <- adonis2(bc_dist ~ Treatment * Feces + Sex,
                        data = meta_df, strata = meta_df$Animal,
                        permutations = 999)
cat("  Bray-Curtis PERMANOVA R2:", round(permanova_bc$R2, 3), "\n")

permanova_bc_df <- data.frame(
  Term = rownames(permanova_bc),
  Df = permanova_bc$Df,
  SumSqs = permanova_bc$SumOfSqs,
  R2 = permanova_bc$R2,
  F = permanova_bc$F,
  p = permanova_bc$`Pr(>F)`,
  Distance = "Bray-Curtis"
)

set.seed(42)
permanova_ait <- adonis2(aitchison_dist ~ Treatment * Feces + Sex,
                         data = meta_df, strata = meta_df$Animal,
                         permutations = 999)

permanova_ait_df <- data.frame(
  Term = rownames(permanova_ait),
  Df = permanova_ait$Df,
  SumSqs = permanova_ait$SumOfSqs,
  R2 = permanova_ait$R2,
  F = permanova_ait$F,
  p = permanova_ait$`Pr(>F)`,
  Distance = "Aitchison"
)

permanova_all <- bind_rows(permanova_bc_df, permanova_ait_df)
write.csv(permanova_all, file.path(stats_dir, "beta_diversity_permanova.csv"), row.names = FALSE)

disp_bc <- betadisper(bc_dist, meta_df$Treatment)
disp_test <- permutest(disp_bc, permutations = 999)
betadisper_df <- data.frame(
  Test = "Betadisper (Bray-Curtis, Treatment)",
  F = disp_test$tab$F[1],
  p = disp_test$tab$`Pr(>F)`[1]
)
write.csv(betadisper_df, file.path(stats_dir, "betadisper_results.csv"), row.names = FALSE)

cat("  Betadisper p =", round(disp_test$tab$`Pr(>F)`[1], 4), "\n")

if (sex_interaction_significant) {
  cat("  Running sex-stratified PERMANOVA...\n")
  perm_stratified <- list()
  for (s in c("male", "female")) {
    idx <- meta_df$Sex == s
    bc_sub <- as.dist(as.matrix(bc_dist)[idx, idx])
    meta_sub <- meta_df[idx, ]
    set.seed(42)
    perm_sub <- adonis2(bc_sub ~ Treatment * Feces, data = meta_sub,
                        strata = meta_sub$Animal, permutations = 999)
    perm_df <- data.frame(
      Term = rownames(perm_sub), Df = perm_sub$Df, SumSqs = perm_sub$SumOfSqs,
      R2 = perm_sub$R2, F = perm_sub$F, p = perm_sub$`Pr(>F)`,
      Distance = "Bray-Curtis", Sex = s
    )
    perm_stratified[[s]] <- perm_df
  }
  write.csv(bind_rows(perm_stratified),
            file.path(stats_dir, "beta_diversity_permanova_stratified.csv"), row.names = FALSE)
}

# ============================================================
# 2c. DIFFERENTIAL ABUNDANCE — ANCOM-BC2
# ============================================================

cat("\n--- 2c. Differential Abundance: ANCOM-BC2 ---\n")

ancom_treatment <- NULL
ancom_global <- NULL
ancom_pairwise <- NULL

if (ancombc_available) {
  
  meta_ancom <- as(sample_data(ps_genus), "data.frame")
  meta_ancom$Animal <- factor(meta_ancom$Animal)
  meta_ancom$Treatment <- factor(meta_ancom$Treatment, levels = c("Ctrl", "TAM"))
  meta_ancom$Feces <- factor(meta_ancom$Feces, levels = mb_params$timepoints)
  meta_ancom$Sex <- factor(meta_ancom$Sex, levels = c("female", "male"))
  
  set.seed(42)
  
  ancom_res <- ancombc2(
    data = ps_genus,
    fix_formula = "Treatment * Feces + Sex",
    rand_formula = "(1 | Animal)",
   group = "Treatment",
    p_adj_method = mb_params$fdr_method,
    prv_cut = mb_params$prevalence_threshold,
    lib_cut = mb_params$min_reads,
    global = TRUE,
    pairwise = TRUE
  )
  
  # --- Model coefficients ---
  res_ancom <- ancom_res$res
  
  if (!is.null(res_ancom)) {
    ancom_treatment <- res_ancom %>%
      dplyr::select(
        OTU = taxon,
        # F1: Treatment effect at reference timepoint F1
        log2FC_Treatment_F1 = lfc_TreatmentTAM,
        SE_Treatment_F1 = se_TreatmentTAM,
        W_Treatment_F1 = W_TreatmentTAM,
        p_Treatment_F1 = p_TreatmentTAM,
        q_Treatment_F1 = q_TreatmentTAM,
        # F2: Treatment × Feces interaction
        log2FC_Treatment_F2_interaction = `lfc_TreatmentTAM:FecesF2`,
        SE_Treatment_F2_interaction = `se_TreatmentTAM:FecesF2`,
        W_Treatment_F2_interaction = `W_TreatmentTAM:FecesF2`,
        p_Treatment_F2_interaction = `p_TreatmentTAM:FecesF2`,
        q_Treatment_F2_interaction = `q_TreatmentTAM:FecesF2`,
        # F3: Treatment × Feces interaction
        log2FC_Treatment_F3_interaction = `lfc_TreatmentTAM:FecesF3`,
        SE_Treatment_F3_interaction = `se_TreatmentTAM:FecesF3`,
        W_Treatment_F3_interaction = `W_TreatmentTAM:FecesF3`,
        p_Treatment_F3_interaction = `p_TreatmentTAM:FecesF3`,
        q_Treatment_F3_interaction = `q_TreatmentTAM:FecesF3`,
        # F4: Treatment × Feces interaction
        log2FC_Treatment_F4_interaction = `lfc_TreatmentTAM:FecesF4`,
        SE_Treatment_F4_interaction = `se_TreatmentTAM:FecesF4`,
        W_Treatment_F4_interaction = `W_TreatmentTAM:FecesF4`,
        p_Treatment_F4_interaction = `p_TreatmentTAM:FecesF4`,
        q_Treatment_F4_interaction = `q_TreatmentTAM:FecesF4`
      ) %>%
      left_join(
        as.data.frame(tax_table(ps_genus)) %>%
          rownames_to_column("OTU") %>%
          dplyr::select(OTU, Genus, Family, Order, Phylum),
        by = "OTU"
      ) %>%
      mutate(
        interpretation_F1 = ifelse(
          q_Treatment_F1 < fdr_thresh,
          "Significant treatment effect at timepoint 1",
          "No significant treatment effect at timepoint 1"
        ),
        interpretation_F2 = ifelse(
          q_Treatment_F2_interaction < fdr_thresh,
          "Significant change in treatment effect at timepoint 2 vs timepoint 1",
          "No significant change in treatment effect at timepoint 2 vs timepoint 1"
        ),
        interpretation_F3 = ifelse(
          q_Treatment_F3_interaction < fdr_thresh,
          "Significant change in treatment effect at timepoint 3 vs timepoint 1",
          "No significant change in treatment effect at timepoint 3 vs timepoint 1"
        ),
        interpretation_F4 = ifelse(
          q_Treatment_F4_interaction < fdr_thresh,
          "Significant change in treatment effect at timepoint 4 vs timepoint 1",
          "No significant change in treatment effect at timepoint 4 vs timepoint 1"
        )
      )
    
    write.csv(
      ancom_treatment,
      file.path(stats_dir, "ancombc2_model_coefficients.csv"),
      row.names = FALSE
    )
  }
  
  # --- Overall Treatment effect ---
  if (!is.null(ancom_res$res_global)) {
    ancom_global <- ancom_res$res_global %>%
      dplyr::select(
        OTU = taxon,
        W_Global_Treatment = W,
        p_Global_Treatment = p_val,
        q_Global_Treatment = q_val,
        diff_abn_Global_Treatment = diff_abn,
        passed_ss_Global_Treatment = passed_ss,
        diff_robust_abn_Global_Treatment = diff_robust_abn) %>%
      mutate(interpretation = ifelse(q_Global_Treatment < fdr_thresh,
          "Significant overall treatment effect",
          "No significant overall treatment effect"))
    
    write.csv(ancom_global, file.path(stats_dir, "ancombc2_global_treatment.csv"),row.names = FALSE  )
  }
  
  # --- Pairwise Treatment comparisons ---
  if (!is.null(ancom_res$res_pair)) {
    ancom_pairwise <- ancom_res$res_pair
    write.csv(  ancom_pairwise, file.path(stats_dir, "ancombc2_pairwise.csv"),  row.names = FALSE  )
    cat("\nPairwise columns:\n")
    print(colnames(ancom_pairwise))
  }
  
  saveRDS(ancom_res,file.path(stats_dir, "ancombc2_full_results.rds"))
  
} else {
  cat("  ANCOMBC package not available (requires R >= 4.5.0).\n")
  cat("  Skipping ANCOM-BC2. ALDEx2 will be used as primary DA method.\n")
  cat("  Install ANCOMBC on your machine with: BiocManager::install('ANCOMBC')\n")
}

# --- Sex-stratified ANCOM-BC2 (only if Sex interaction significant AND ANCOMBC available) ---
if (sex_interaction_significant && ancombc_available) {
  cat("  Running sex-stratified ANCOM-BC2...\n")
  for (s in c("male", "female")) {
    ps_sub <- subset_samples(ps_genus, Sex == s)
    ps_sub <- prune_taxa(taxa_sums(ps_sub) > 0, ps_sub)
    meta_sub <- as(sample_data(ps_sub), "data.frame")
    meta_sub$Animal <- as.factor(meta_sub$Animal)
    meta_sub$Treatment <- factor(meta_sub$Treatment, levels = c("Ctrl", "TAM"))
    meta_sub$Feces <- factor(meta_sub$Feces, levels = mb_params$timepoints)
    
    res_sub <- tryCatch({
      ancombc2(data = ps_sub, fix_formula = "Treatment * Feces",
               rand_formula = "~ 1 | Animal", group = "Treatment",
               p_adj_method = mb_params$fdr_method,
               prv_cut = mb_params$prevalence_threshold,
               lib_cut = mb_params$min_reads)
    }, error = function(e) {
      cat("    ANCOM-BC2 failed for", s, ":", conditionMessage(e), "\n")
      NULL
    })
    
    if (!is.null(res_sub) && !is.null(res_sub$res)) {
      ancom_sub <- res_sub$res %>%
        dplyr::select(OTU = id, Genus = tax_genus, Family = tax_family,
                      log2FC = diff_Treatment_TAM_vs_Ctrl,
                      SE = se_Treatment_TAM_vs_Ctrl,
                      W = W_Treatment_TAM_vs_Ctrl,
                      p = p_Treatment_TAM_vs_Ctrl,
                      q = q_Treatment_TAM_vs_Ctrl) %>%
        mutate(Sex = s, significance = sig_stars(q))
      write.csv(ancom_sub,
                file.path(stats_dir, paste0("ancombc2_results_", s, ".csv")),
                row.names = FALSE)
    }
  }
}

# ============================================================
# 2d. DIFFERENTIAL ABUNDANCE — ALDEx2 (cross-sectional complement)
# ============================================================
cat("\n--- 2d. Differential Abundance: ALDEx2 ---\n")

run_aldex <- function(ps_sub, label) {
  counts <- as(otu_table(ps_sub), "matrix")
  if (taxa_are_rows(ps_sub)) counts <- t(counts)
  counts <- t(counts)
  
  conditions <- as.character(sample_data(ps_sub)$Treatment)
  
  cat("    ALDEx2:", label, "- n =", ncol(counts), "samples\n")
  
  aldex_out <- aldex(counts, conditions, mc.samples = 128, test = "t", effect = TRUE,
                     verbose = FALSE)
  
  aldex_df <- aldex_out %>%
    rownames_to_column("OTU") %>%
    mutate(
      label = label,
      q_welch = p.adjust(we.ep, method = mb_params$fdr_method),
      q_wilcox = p.adjust(wi.ep, method = mb_params$fdr_method),
      significance_welch = sig_stars(q_welch)
    )
  
  tax_info <- as.data.frame(tax_table(ps_sub)) %>%
    rownames_to_column("OTU") %>%
    dplyr::select(OTU, Genus, Family, Order, Phylum)
  aldex_df <- aldex_df %>% left_join(tax_info, by = "OTU")
  
  return(aldex_df)
}

aldex_all <- list()

for (f in mb_params$timepoints) {
  ps_sub <- subset_samples(ps_genus, Feces == f)
  ps_sub <- prune_taxa(taxa_sums(ps_sub) > 0, ps_sub)
  aldex_all[[paste0("per_timepoint_", f)]] <- run_aldex(ps_sub, paste0("Feces_", f))
}

ps_cdhfd <- subset_samples(ps_genus, Feces %in% c("F3", "F4"))
ps_cdhfd <- prune_taxa(taxa_sums(ps_cdhfd) > 0, ps_cdhfd)
aldex_all[["cdhfd_pooled"]] <- run_aldex(ps_cdhfd, "CDHFD_pooled_F3F4")

ps_male_cdhfd <- subset_samples(ps_genus, Sex == "male" & Feces %in% c("F3", "F4"))
ps_male_cdhfd <- prune_taxa(taxa_sums(ps_male_cdhfd) > 0, ps_male_cdhfd)
aldex_all[["male_cdhfd"]] <- run_aldex(ps_male_cdhfd, "Male_CDHFD_pooled")

aldex_combined <- bind_rows(aldex_all)
write.csv(aldex_combined, file.path(stats_dir, "aldex2_results.csv"), row.names = FALSE)

aldex_sig <- aldex_combined %>% filter(q_welch < fdr_thresh)
write.csv(aldex_sig, file.path(stats_dir, "aldex2_significant.csv"), row.names = FALSE)
cat("  ALDEx2 significant taxa (Welch q <", fdr_thresh, "):", nrow(aldex_sig), "\n")

# ============================================================
# 2e. CONSENSUS DA TABLE
# ============================================================
cat("\n--- 2e. Consensus DA Table ---\n")

if (!is.null(ancom_treatment)) {
  ancom_for_consensus <- ancom_treatment %>%
    dplyr::select(Genus, Phylum, Family, Order,
                  ancom_log2FC = log2FC_Treatment_F1,
                  ancom_q = q_Treatment_F1,
                  ancom_W = W_Treatment_F1)
  
  aldex_cdhfd <- aldex_all[["cdhfd_pooled"]] %>%
    dplyr::select(Genus, aldex_effect = effect,
                  aldex_q_welch = q_welch,
                  aldex_q_wilcox = q_wilcox)
  
  consensus <- ancom_for_consensus %>%
    full_join(aldex_cdhfd, by = "Genus") %>%
    mutate(
      ancom_sig = ancom_q < fdr_thresh,
      aldex_sig = aldex_q_welch < fdr_thresh,
      consensus = case_when(
        ancom_sig & aldex_sig ~ "consensus_significant",
        ancom_sig | aldex_sig ~ "suggestive",
        TRUE ~ "ns"
      )
    ) %>%
    arrange(consensus, ancom_q)
} else {
  consensus <- aldex_all[["cdhfd_pooled"]] %>%
    dplyr::select(Genus, Phylum, Family, Order,
                  aldex_effect = effect,
                  aldex_q_welch = q_welch,
                  aldex_q_wilcox = q_wilcox) %>%
    mutate(
      aldex_sig = aldex_q_welch < fdr_thresh,
      consensus = ifelse(aldex_sig, "aldex_significant", "ns")
    ) %>%
    arrange(consensus, aldex_q_welch)
}

write.csv(consensus, file.path(stats_dir, "da_consensus.csv"), row.names = FALSE)
cat("  Consensus significant:", sum(consensus$consensus == "consensus_significant"), "\n")
cat("  Suggestive:", sum(consensus$consensus == "suggestive"), "\n")

# ============================================================
# 2f. TARGETED LACTOBACILLUS ANALYSIS (hypothesis-driven)
# ============================================================
cat("\n--- 2f. Targeted Lactobacillus Analysis ---\n")

ps_lacto <- subset_taxa(ps_genus, Genus == "Lactobacillus")
if (ntaxa(ps_lacto) > 0) {
  lacto_counts <- as(otu_table(ps_lacto), "matrix")
  if (!taxa_are_rows(ps_lacto)) lacto_counts <- t(lacto_counts)
  lacto_total <- colSums(lacto_counts)
  
  lacto_clr <- log(lacto_total + 0.5) - mean(log(lacto_total + 0.5))
  
  lacto_df <- data.frame(
    SampleID = names(lacto_total),
    Lacto_counts = lacto_total,
    Lacto_CLR = lacto_clr
  ) %>%
    left_join(metadata, by = "SampleID")
  
  formula_full_lacto <- as.formula("Lacto_CLR ~ Treatment * Sex * Feces + (1|Animal)")
  formula_reduced_lacto <- as.formula("Lacto_CLR ~ Treatment * Feces + Sex + (1|Animal)")
  
  model_lacto <- tryCatch(
    lmer(formula_full_lacto, data = lacto_df),
    error = function(e) {
      cat("    Lactobacillus full 3-way model failed, using reduced model\n")
      lmer(formula_reduced_lacto, data = lacto_df)
    }
  )
  
  anova_lacto <- Anova(model_lacto, type = 3)
  anova_lacto_df <- as.data.frame(anova_lacto) %>%
    rownames_to_column("Term")
  
  emm_lacto <- emmeans(model_lacto, ~ Treatment | Feces)
  pwc_lacto <- as.data.frame(contrast(emm_lacto, method = "pairwise", adjust = "BH"))
  pwc_lacto$significance <- sig_stars(pwc_lacto$p.value)
  pwc_lacto$Contrast <- "Ctrl_vs_TAM"
  
  pwc_lacto_sex_list <- list()
  for (s in c("male", "female")) {
    lacto_sub <- lacto_df %>% filter(Sex == s)
    model_lacto_sub <- lmer(Lacto_CLR ~ Treatment * Feces + (1|Animal),
                            data = lacto_sub)
    emm_sub <- emmeans(model_lacto_sub, ~ Treatment | Feces)
    pwc_sub <- as.data.frame(contrast(emm_sub, method = "pairwise", adjust = "BH"))
    pwc_sub$Sex <- s
    pwc_sub$significance <- sig_stars(pwc_sub$p.value)
    pwc_lacto_sex_list[[s]] <- pwc_sub
  }
  pwc_lacto_male <- bind_rows(pwc_lacto_sex_list)
  
  emm_summary <- as.data.frame(emm_lacto)
  
  lacto_results <- list(
    anova = anova_lacto_df,
    posthoc_pooled = pwc_lacto,
    posthoc_by_sex = pwc_lacto_male,
    emmeans = emm_summary,
    data = lacto_df
  )
  
  write.csv(anova_lacto_df, file.path(stats_dir, "lactobacillus_anova.csv"), row.names = FALSE)
  write.csv(pwc_lacto, file.path(stats_dir, "lactobacillus_posthoc_pooled.csv"), row.names = FALSE)
  write.csv(pwc_lacto_male, file.path(stats_dir, "lactobacillus_posthoc_by_sex.csv"), row.names = FALSE)
  
  cat("  Lactobacillus ANOVA:\n")
  print(anova_lacto_df)
  cat("\n  Post-hoc (pooled, Ctrl vs TAM per Feces):\n")
  print(pwc_lacto %>% dplyr::select(Feces, estimate, p.value, significance))
  cat("\n  Post-hoc (male, Ctrl vs TAM per Feces):\n")
  print(pwc_lacto_male %>% filter(Sex == "male") %>%
          dplyr::select(Feces, estimate, p.value, significance))
  
} else {
  cat("  WARNING: No Lactobacillus taxa found at genus level\n")
}

# ============================================================
# SAVE SEX INTERACTION FLAG FOR PLOT SCRIPT
# ============================================================
sex_flag <- data.frame(
  sex_interaction_significant = sex_interaction_significant,
  note = ifelse(sex_interaction_significant,
                "Sex-faceted plots will be generated",
                "No Sex interaction significant; combined-sex plots only")
)
write.csv(sex_flag, file.path(stats_dir, "sex_interaction_flag.csv"), row.names = FALSE)

cat("\n=== Statistics complete ===\n")
cat("Results saved to:", stats_dir, "\n")