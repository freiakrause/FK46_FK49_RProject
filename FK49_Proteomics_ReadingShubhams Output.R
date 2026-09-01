#FK49 Read Proteomics what we got from Shubham - files with p values, FCs and some transformed protein abundance data. Dont know the exact transormation
rm(list=ls())
gc()
library(readxl)
library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)
library(ggrepel)
library(janitor)
library(tibble)
library(limma)
library(pheatmap)
source("FK49_Definitions.R")

proteom_input_pwd<-PATHS$proteomics$input
proteom_output_pwd<-PATHS$proteomics$output
Proteins<-read_excel(file.path(proteom_input_pwd,"/FK49_TAM_EtOH_STATISTICAL_results.xlsx"), sheet = 1)%>%
  mutate(Protein.Names= (gsub("_MOUSE", "", Protein.Names)))%>%
  rename(Name = Protein.Names,
              logFC = logFC...6,
              pValue = pValue...7,
              adj_pvalue = BH_corrected_pvalue...8,
              Subset_1 = Subset...9,
              Direction = Regulation...10,
              logFC_2 = logFC...33,
              pValue_2 = pValue...34,
              BH_corrected_pvalue_2 = BH_corrected_pvalue...35,
              Subset_2 = Subset...36,
              Regulation_2 = Regulation...37)%>%
  select(-Regulation_2,-logFC_2,-pValue_2,-BH_corrected_pvalue_2,-Subset_2,-pValue,-adj_pvalue,-Direction,-Subset_1,-Subset_2)
meta<-read_excel(file.path(proteom_input_pwd,"FK49_FK46_BH_SampleList_Proteomics_Shubam.xlsx"))%>%
                   filter(ExpID=="FK49", Batch == "2")%>%select(Animal,Sex,Treatment,PotentialShubhamNames)
protein_matrix <- Proteins %>%
  select(  F_TAM_1, F_TAM_2,  M_EtOH_1, M_EtOH_2, M_EtOH_3, M_EtOH_4, M_EtOH_5, 
           F_TAM_3, F_TAM_4, F_TAM_5,  M_TAM_1, M_TAM_2, M_TAM_3,  F_EtOH_1, 
           F_EtOH_2, F_EtOH_3, F_EtOH_4,  M_TAM_4 ) %>%
  as.matrix()
rownames(protein_matrix) <- Proteins$Name
meta <- meta[match(colnames(protein_matrix),
                   meta$PotentialShubhamNames), ]
meta$Sex <- factor(meta$Sex)
meta$Treatment <- factor(meta$Treatment)
design <- model.matrix(~ Treatment * Sex, data = meta)
fit <- lmFit(protein_matrix, design)
fit <- eBayes(fit)

coef_names <- colnames(design)
treatment_coef <- which(coef_names == "TreatmentTAM")
interaction_coef <- which(coef_names == "TreatmentTAM:Sexmale")

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

fit_contrasts <- contrasts.fit( fit,contrasts = cbind(  Treatment_overall = C_overall,
                                                        Treatment_x_Sex = C_interaction,
                                                        Treatment_female = C_female,
                                                        Treatment_male = C_male  ))

fit_contrasts <- eBayes(fit_contrasts)
results_Treatment <- topTable(fit_contrasts,coef = "Treatment_overall", number = Inf, adjust.method = "BH") %>%
  rownames_to_column("Name") %>%
  select(Name, logFC, P.Value, adj.P.Val) %>%
  rename( logFC_Treatment = logFC,
    pValue_Treatment = P.Value,
    adj_pvalue_Treatment = adj.P.Val)

results_Interaction <- topTable(fit_contrasts,coef = "Treatment_x_Sex", number = Inf, adjust.method = "BH")%>%
  rownames_to_column("Name") %>%
  select(Name, logFC, P.Value, adj.P.Val) %>%
  rename(  logFC_Treatment_x_Sex = logFC,
    pValue_Treatment_x_Sex = P.Value,
    adj_pvalue_Treatment_x_Sex = adj.P.Val)

results_Female <- topTable( fit_contrasts,  coef = "Treatment_female",number = Inf, adjust.method = "BH")%>%
  rownames_to_column("Name") %>%
  select(Name, logFC, P.Value, adj.P.Val) %>%
  rename(  logFC_Treatment_female = logFC,
    pValue_Treatment_female = P.Value,
    adj_pvalue_Treatment_female = adj.P.Val  )

results_Male <- topTable( fit_contrasts,coef = "Treatment_male",number = Inf, adjust.method = "BH")%>%
  rownames_to_column("Name") %>%
  select(Name, logFC, P.Value, adj.P.Val) %>%
  rename(  logFC_Treatment_male = logFC,
    pValue_Treatment_male = P.Value,
    adj_pvalue_Treatment_male = adj.P.Val
  )
write.csv2(results_Treatment, file.path(proteom_output_pwd, "results_Treatment_overall.csv"), row.names = FALSE)
write.csv2(results_Interaction, file.path(proteom_output_pwd, "results_Treatment_x_Sex.csv"), row.names = FALSE)
write.csv2(results_Female, file.path(proteom_output_pwd, "results_Treatment_female.csv"), row.names = FALSE)
write.csv2(results_Male, file.path(proteom_output_pwd, "results_Treatment_male.csv"), row.names = FALSE)

Proteins <- Proteins %>%
  left_join(results_Treatment, by = "Name") %>%
  left_join(results_Interaction, by = "Name") %>%
  left_join(results_Female, by = "Name") %>%
  left_join(results_Male, by = "Name") %>%
  mutate( Direction = case_when(
      adj_pvalue_Treatment < 0.05 & logFC_Treatment > 1 ~ "Up",
      adj_pvalue_Treatment < 0.05 & logFC_Treatment < -1 ~ "Down",
      TRUE ~ "NS" ),
    Direction_Male = case_when( adj_pvalue_Treatment_male < 0.05 & logFC_Treatment_male > 1 ~ "Up",
      adj_pvalue_Treatment_male < 0.05 & logFC_Treatment_male < -1 ~ "Down",TRUE ~ "NS" ),
    Direction_Female = case_when(
      adj_pvalue_Treatment_female < 0.05 & logFC_Treatment_female > 1 ~ "Up",
      adj_pvalue_Treatment_female < 0.05 & logFC_Treatment_female < -1 ~ "Down",
      TRUE ~ "NS"  ),
    Direction_Interaction = case_when(
      adj_pvalue_Treatment_x_Sex < 0.05 & abs(logFC_Treatment_x_Sex) > 1 ~ "Different",
      TRUE ~ "NS"  )
  )
# Volcano -----
p_volcano <- ggplot(Proteins,aes(x = logFC_Treatment, y = -log10(pValue_Treatment))) +
  geom_point(aes(fill = Direction), alpha = 0.5, size = 3,stroke = 0.5,
              shape=21,color= "black") +
  scale_fill_manual(values = c("blue","grey60","firebrick")) +
  geom_vline(xintercept = c(-1,1), linetype = "dashed", color = "grey80") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey80") +
  labs( title = paste0("Volcano plot - Treatment"),
        x = expression(paste("FC [", log[2], "]")),       
        y = expression(paste("-log"[10], "(adj.p.value)")))+
  theme_classic()+  
  theme(panel.grid= element_line(color ="grey90", linewidth = 0.1))+
  geom_text_repel(data = Proteins %>% filter(Direction %in% c("Up","Down")), 
                  aes(label = Name),
                  size = 2.5,
                  max.overlaps = 25) +
  coord_cartesian(xlim = c(-4.5, 4.5), ylim = c(0, 15))

p_volcano
ggsave(plot= p_volcano, filename= paste0("Prots_a_volcano.png"), width= 9, height = 9, dpi = 300, 
       path =proteom_output_pwd)

# Heatmaps -----
## Heatmap up -----
values_for_heatmap_up<-as.data.frame(Proteins%>% filter(adj_pvalue_Treatment<0.05&logFC_Treatment>1)%>%
                                    select(Name,Genes,contains(c("EtOH","TAM"))) )
rownames(values_for_heatmap_up)<-values_for_heatmap_up$Name
values_for_heatmap_up$Name <-NULL
ann_col <- data.frame(Sample = colnames(values_for_heatmap_up %>% select(-Genes))) %>%
  separate( Sample, into = c("Sex", "Treatment", "Replicate"),sep = "_") %>%
  select(-Replicate) %>%
  mutate(Sex = case_when(Sex == "F" ~ "female", Sex == "M" ~ "male" ),
         Treatment = case_when(Treatment == "EtOH" ~ "Ctrl", TRUE ~ Treatment ),
         Sex = factor( Sex, levels = c("female", "male")),
         Treatment = factor(Treatment, levels = c("Ctrl", "TAM"))
  )
values_for_heatmap_up$Genes <-NULL
rownames(ann_col) <- colnames(values_for_heatmap_up)


ph<-pheatmap(values_for_heatmap_up,
             cluster_rows = T,
             cluster_cols = F,
             show_colnames = F,
             scale="row",
            annotation_col = ann_col,
             annotation_colors = list(Sex=Sex_colors,Treatment= Treatment_colors[c("Ctrl","TAM")]))
ggsave(ph, file= "02_Heat_sig_up.png",
       path=proteom_output_pwd, width = 12,
       height=1+nrow(values_for_heatmap_up)/5,dpi=300, bg="white",limitsize = FALSE)

## Heatmapt down -----
values_for_heatmap_down<-as.data.frame(Proteins%>% filter(adj_pvalue_Treatment<0.05 & logFC_Treatment< -1)%>%
                                       select(Name,Genes,contains(c("EtOH","TAM"))) )
rownames(values_for_heatmap_down)<-values_for_heatmap_down$Name
values_for_heatmap_down$Name <-NULL
ann_col <- data.frame(Sample = colnames(values_for_heatmap_down %>% select(-Genes))) %>%
  separate( Sample, into = c("Sex", "Treatment", "Replicate"),sep = "_") %>%
  select(-Replicate) %>%
  mutate(Sex = case_when(Sex == "F" ~ "female", Sex == "M" ~ "male" ),
         Treatment = case_when(Treatment == "EtOH" ~ "Ctrl", TRUE ~ Treatment ),
         Sex = factor( Sex, levels = c("female", "male")),
         Treatment = factor(Treatment, levels = c("Ctrl", "TAM")))

values_for_heatmap_down$Genes <-NULL
rownames(ann_col) <- colnames(values_for_heatmap_down)


ph<-pheatmap(
  values_for_heatmap_down,
  cluster_rows = T,
  cluster_cols =F,
  scale = "row",
  show_colnames = F,
    annotation_col = ann_col,
  annotation_colors = list(
    Sex = Sex_colors,
    Treatment = Treatment_colors[c("Ctrl", "TAM")]
  )
)

ggsave(ph, file= "02_Heat_sig_down.png",path=proteom_output_pwd,width = 12,
       height=1+nrow(values_for_heatmap_down)/5,dpi=300, bg="white",limitsize = FALSE)

## HeatmaptTop50 up and Top50 DOwn inone -----
values_for_heatmap_top100 <- as.data.frame(
  Proteins %>%
    filter(adj_pvalue_Treatment < 0.05) %>%
    arrange(logFC_Treatment) %>%
    slice_head(n = 50) %>%
    bind_rows(
      Proteins %>%
        filter(adj_pvalue_Treatment < 0.05) %>%
        arrange(desc(logFC_Treatment)) %>%
        slice_head(n = 50)
    ) %>%
    select(Name, Genes, contains(c("EtOH", "TAM")))
)

rownames(values_for_heatmap_top100) <- values_for_heatmap_top100$Name
values_for_heatmap_top100$Name <- NULL

ann_col <- data.frame(
  Sample = colnames(values_for_heatmap_top100 %>% select(-Genes))
) %>%
  separate(Sample, into = c("Sex", "Treatment", "Replicate"), sep = "_") %>%
  select(-Replicate) %>%
  mutate(
    Sex = case_when(Sex == "F" ~ "female", Sex == "M" ~ "male"),
    Treatment = case_when(Treatment == "EtOH" ~ "Ctrl", TRUE ~ Treatment),
    Sex = factor(Sex, levels = c("female", "male")),
    Treatment = factor(Treatment, levels = c("Ctrl", "TAM"))
  )

values_for_heatmap_top100$Genes <- NULL
rownames(ann_col) <- colnames(values_for_heatmap_top100)

ph <- pheatmap(
  values_for_heatmap_top100,
  cluster_rows = TRUE,
  cluster_cols = T,
  scale = "row",
  show_colnames = FALSE,
  annotation_col = ann_col,
  annotation_colors = list(
    Sex = Sex_colors,
    Treatment = Treatment_colors[c("Ctrl", "TAM")]
  )
)

ggsave( ph, file = "02_Heat_top50up_top50down.png",
  path = proteom_output_pwd,width = 12, height = 2 + nrow(values_for_heatmap_top100) / 8,
  dpi = 300,bg = "white", limitsize = FALSE)
# PCA ----
## PCA decriptive/exploraory -all proteins -----

pca <- prcomp(
  t(protein_matrix),
  center = TRUE,
  scale. = TRUE
)

# PCA scores
pca_df <- as.data.frame(pca$x) %>%
  tibble::rownames_to_column("Sample") %>%
  tidyr::separate(Sample,  into = c("Sex", "Treatment", "Replicate"), sep = "_" ) %>%
  mutate(  Sex = case_when(Sex == "F" ~ "female",Sex == "M" ~ "male"  ),
    Treatment = case_when(Treatment == "EtOH" ~ "Ctrl",
      TRUE ~ Treatment),Sex = factor(Sex, levels = c("female", "male")),
    Treatment = factor(Treatment, levels = c("Ctrl", "TAM")))

# PCA Plot
pca_plot <- ggplot( pca_df, aes(x = PC1, y = PC2, shape = Sex, color = Treatment, fill = Treatment)) +
  geom_point(size = 4,alpha = 0.5) +
  stat_ellipse(   aes(x = PC1, y = PC2, group = Treatment, fill = Treatment),
    type = "norm",level = 0.95,  geom = "polygon",  alpha = 0.05,
    linewidth = 0) +
  stat_ellipse(  aes(x = PC1, y = PC2, group = Treatment, fill = Treatment), type = "norm",level = 0.95,
    alpha = 0.2,linewidth = 0.8) +
  scale_color_manual( values = Treatment_colors[c("Ctrl", "TAM")] ) +
  scale_fill_manual(values = Treatment_colors[c("Ctrl", "TAM")] ) +
  scale_shape_manual( values = Sex_shape) +
  theme_classic() +
  labs( x = paste0("PC1 (", round(100 * summary(pca)$importance[2, 1], 1),"%)"),
    y = paste0( "PC2 (",round(100 * summary(pca)$importance[2, 2], 1),"%)") )
pca_plot
ggsave( pca_plot, file = "03_PCA_Proteins_all.png",
        path = proteom_output_pwd,width = 5, height = 5,
        dpi = 300,bg = "white", limitsize = FALSE)
## PCA significant proteins-----
pca_proteins <- Proteins %>%
  filter( adj_pvalue_Treatment < 0.05,
    abs(logFC_Treatment) > 1) %>% pull(Name)

protein_matrix_sig <- protein_matrix[pca_proteins, ]
pca_sig <- prcomp(
  t(protein_matrix_sig),
  center = TRUE,
  scale. = TRUE)

# PCA scores
pca_sig_df <- as.data.frame(pca_sig$x) %>%
  tibble::rownames_to_column("Sample") %>%
  tidyr::separate(Sample,  into = c("Sex", "Treatment", "Replicate"), sep = "_" ) %>%
  mutate(  Sex = case_when(Sex == "F" ~ "female",Sex == "M" ~ "male"  ),
           Treatment = case_when(Treatment == "EtOH" ~ "Ctrl",
                                 TRUE ~ Treatment),Sex = factor(Sex, levels = c("female", "male")),
           Treatment = factor(Treatment, levels = c("Ctrl", "TAM")))

# PCA Plot
pca_sig_plot <- ggplot( pca_sig_df, aes(x = PC1, y = PC2, shape = Sex, color = Treatment, fill = Treatment)) +
  geom_point(size = 4,alpha = 0.5) +
  stat_ellipse(  aes(x = PC1, y = PC2, group = Treatment, fill = Treatment),
                 type = "norm",level = 0.95,  geom = "polygon",  alpha = 0.05,
                 linewidth = 0) +
  stat_ellipse(  aes(x = PC1, y = PC2, group = Treatment, fill = Treatment), type = "norm",level = 0.95,
                alpha = 0.2,linewidth = 0.8) +
  scale_color_manual( values = Treatment_colors[c("Ctrl", "TAM")] ) +
  scale_fill_manual(values = Treatment_colors[c("Ctrl", "TAM")] ) +
  scale_shape_manual( values = Sex_shape) +
  theme_classic() +
  labs( x = paste0("PC1 (", round(100 * summary(pca_sig)$importance[2, 1], 1),"%)"),
        y = paste0( "PC2 (",round(100 * summary(pca_sig)$importance[2, 2], 1),"%)") )
pca_sig_plot
ggsave( pca_sig_plot, file = "03_PCA_Proteins_sig.png",
        path = proteom_output_pwd,width = 5, height = 5,
        dpi = 300,bg = "white", limitsize = FALSE)

 # Exploratory -----
biological_oxidations=c("Nnmt","Cyp4a10","Gsta","Cyp2a4","Gstt2",
                        "Ugt1a9","Cyp2c29","Cyp3a11","Gstt1","Fmo2",
                        "Ugp2","Mgst3","Sult1c2","Sult1b1","Adh4")
metabolism <- c("Adh4","Abcc3","Acot1","Ttr","Alb","C3","Ca3","Fabp2",
                "Cox7a1", "Apoa4","Fabp5","Hsd3b5","Acsl6","Lpgat1","Tymp",
                "Slco1a4","Acot2","Mvk","Fga","Orm2","Orm1","Gas6","Itih3",
                "Fgb","Fgg")

phaseII_conjugation_of_compounds <- c("Nnmt","Gsta","Ugp2","Ugt1a9","Gstt2","Gstt1","Mgst3","Sult1c2","Sult1b1")
metabolism_of_lipids <- c("Fabp2","Cox7a1","Apoa4","Fabp5","Hsd3b5","Acsl6",
                          "Lpgat1", "Acot2","Mvk")
platelet_degranulation <- c("Fga","Fgb","Fgg","Orm1","Orm2","Gas6","Itih3")
response_to_elevated_platelet_cytosolic_Ca2 <- c( "Fga","Fgb","Fgg","Gas6","Itih3")
drug_ADME <- c( "Cyp4a10","Cyp2a4","Ugt1a9","Cyp2c29","Cyp3a11","Fmo2", 
                "Mgst3","Sult1c2","Adh4")
ciprofloxacin_ADME <- c("Fga","Fgb","Fgg","Gas6")
GRB2_SOS_integrin_MAPK <- c("Itih3","Gas6","Fga")
p130Cas_integrin_MAPK <- c( "Itih3","Gas6","Fga")


values_for_heatmap<-as.data.frame(Proteins%>%
  select(Name,Genes, contains(c("Etoh","Tam")))%>%
  filter(Genes %in% c(biological_oxidations,metabolism_of_lipids,
                      metabolism,phaseII_conjugation_of_compounds)))

values_for_heatmap<-as.data.frame(Proteins%>%
                                    select(Name,Genes, contains(c("Etoh","Tam")))%>%
                                    filter(Genes %in% c(bile_acids_production)))
rownames(values_for_heatmap)<-values_for_heatmap$Name
values_for_heatmap<-values_for_heatmap

  
  

  ann <- data.frame(Sample = colnames(values_for_heatmap)) %>%
    separate(Sample, into = c("Sex", "Treatment", "Replicate"), sep = "_") %>%
    select(-Replicate)%>%mutate(
      Sex = case_when(Sex == "F" ~ "female",  Sex == "M" ~ "male")  )
  
  annotation_row <- data.frame(
    Metabolism = values_for_heatmap$Genes %in% metabolism,
    bile_acids_production=values_for_heatmap$Genes %in% bile_acids_production,
    Biological_oxidations = values_for_heatmap$Genes %in% biological_oxidations,
    Metabolism_of_lipids = values_for_heatmap$Genes %in% metabolism_of_lipids,
    PhaseII_conjugation = values_for_heatmap$Genes %in% phaseII_conjugation_of_compounds,
    row.names = values_for_heatmap$Genes
  )
  annotation_row[] <- lapply(annotation_row, factor, levels = c(FALSE, TRUE))
  
  annotation_colors <- list(
    Sex = Sex_colors,
    Treatment = Treatment_colors,
    Diet = Diet_colors )

  annotation_colors$Biological_oxidations <- c("FALSE" = "grey95", "TRUE" = "#8dd3c7")
  annotation_colors$Metabolism <- c("FALSE" = "grey95", "TRUE" = "#ffffb3")
  annotation_colors$Metabolism_of_lipids <- c("FALSE" = "grey95", "TRUE" = "#fb8072")
  annotation_colors$PhaseII_conjugation <- c("FALSE" = "grey90", "TRUE" = "#bebada")
  
  annotation_colors$PhaseII_conjugation <- c("FALSE" = "grey90", "TRUE" = "forestgreen")
  
  rownames(annotation_row) <- rownames(values_for_heatmap)
  rownames(ann) <- colnames(values_for_heatmap)
  heatmap_height <-nrow(heatmap_data)/7+3
  values_for_heatmap<-values_for_heatmap%>%select(-Name,-Genes)
  p_heat<- pheatmap::pheatmap(values_for_heatmap,
                              scale = "row",
                              cluster_rows = TRUE, 
                              cluster_cols = TRUE,
                              annotation_col  = ann,
                              annotation_row = annotation_row,
                              annotation_names_row = FALSE,
                              annotation_colors = annotation_colors,
                              show_colnames = FALSE,
                              treeheight_row = 0,  
                              treeheight_col = 5,   
                              main = "Top4 Enriched Processes")
  
  print(p_heat)
  ggsave(plot= p_heat, filename= "Heatmap_clustered.png",limitsize = FALSE, 
         width= 9, height = heatmap_height, dpi = 500,bg = "white",path=paste0(proteom_output_pwd))
#### HeatMap for Bile Acids
 Uptake<- c("Slc10a1","Slco1a1","Slco1b2","Slco2b1")
 Secretion<-c("Abcb11","Abcc2","Abcg5","Abcg8")
 Synthesis<-c("Cyp7a1","Cyp8b1","Cyp27a1","Cyp7b1","Hsd3b7","Akr1d1","Baat")
 TF<-c("Nr1h4","Nr0b2","Rxra","Hnf4a","Ppara","Ppard","Pparg")
 Conjugation <- c(
   "Slc27a5",  # bile acid activation
   "Baat",     # taurine/glycine conjugation
   "Sult2a1",  # sulfation
   "Sult2a2",
   "Ugt1a1",   # glucuronidation
   "Ugt1a6",
   "Ugt2b1"
 )

 values_for_heatmap<-as.data.frame(Proteins%>%
                                     select(Name,Genes, contains(c("Etoh","Tam")))%>%
                                     filter(Genes %in% c(Uptake,Secretion,Synthesis,TF,Conjugation)))
 rownames(values_for_heatmap)<-values_for_heatmap$Name
 row_category <- case_when(
   values_for_heatmap$Genes %in% Uptake ~ "Uptake",
   values_for_heatmap$Genes %in% Secretion ~ "Secretion",
   values_for_heatmap$Genes %in% Synthesis ~ "Synthesis",
   values_for_heatmap$Genes %in% TF ~ "TF",
   values_for_heatmap$Genes %in% Conjugation ~ "Conjugation"
   
 )
 
 row_category <- factor(
   row_category,
   levels = c("Synthesis",   
     "Secretion",
     "Conjugation",
     "Uptake",
     "TF"
   )
 )

 
 
 
 
 ann <- data.frame(Sample = colnames(values_for_heatmap)) %>%
   separate(Sample, into = c("Sex", "Treatment", "Replicate"), sep = "_") %>%
   select(-Replicate)%>%mutate(
     Sex = case_when(Sex == "F" ~ "female",  Sex == "M" ~ "male")  )
 
 annotation_row <- data.frame(
   Uptake=values_for_heatmap$Genes %in% Uptake,
   Secretion = values_for_heatmap$Genes %in% Secretion,
   Synthesis = values_for_heatmap$Genes %in% Synthesis,
   TF = values_for_heatmap$Genes %in% TF,
   Conjugation = values_for_heatmap$Genes %in% Conjugation,
   row.names = values_for_heatmap$Genes
 )
 annotation_row[] <- lapply(annotation_row, factor, levels = c(FALSE, TRUE))
 
 annotation_colors <- list(
   Sex = Sex_colors,
   Treatment = Treatment_colors,
   Diet = Diet_colors )
 
 annotation_colors$Uptake <- c("FALSE" = "grey95", "TRUE" = "#8dd3c7")
 annotation_colors$Secretion <- c("FALSE" = "grey95", "TRUE" = "#ffffb3")
 annotation_colors$Synthesis <- c("FALSE" = "grey90", "TRUE" = "#bebada")
 annotation_colors$TF <- c("FALSE" = "grey90", "TRUE" = "forestgreen")
 annotation_colors$Conjugation <- c("FALSE" = "grey90", "TRUE" = "brown3")
 
 rownames(annotation_row) <- rownames(values_for_heatmap)
 rownames(ann) <- colnames(values_for_heatmap)
 heatmap_height <-nrow(heatmap_data)/7+3
 
 ord <- order(row_category)
 
 values_for_heatmap <- values_for_heatmap[ord, ]
 annotation_row <- annotation_row[ord, ]
 values_for_heatmap<-values_for_heatmap%>%select(-Name,-Genes)
 gaps_row <- cumsum(table(row_category[ord]))[-length(table(row_category))]
 p_heat<- pheatmap::pheatmap(values_for_heatmap,
                             scale = "row",
                             cluster_rows = FALSE, 
                             cluster_cols = FALSE,
                             annotation_col  = ann,
                             annotation_row = annotation_row,
                             annotation_names_row = FALSE,
                             annotation_colors = annotation_colors,
                             show_colnames = FALSE,
                             treeheight_row = 0,  
                             treeheight_col = 5,  
                             gaps_row = gaps_row,
                             main = "Bile Acid Metabolism")
 
 print(p_heat)
 ggsave(plot= p_heat, filename= "Heatmap_BA_clustered.png",limitsize = FALSE, 
        width= 9, height = heatmap_height, dpi = 500,bg = "white",path=paste0(proteom_output_pwd))
 
 #### HeatMap for PPARs
 PPARa_FAoxidation<-c( "Ppara","Acox1",  "Acot1",  "Acot2",  "Acsl1",  "Acadm",
   "Acadl",   "Cpt1a",   "Cyp4a10",   "Cyp4a14",   "Ehhadh",   "Fabp1",
   "Fgf21",   "Hmgcs2",   "Pdk4" )
 
 PPARg_lipogenesis<-c( "Pparg","Cd36", "Fabp4",  "Plin2",  "Lpl",  "Scd1", "Fasn", "Adipoq" )
 #PPARdb_FAoxidation_Energy<-c("Ppard","Pparb", "Pdk4", "Cpt1a","Ucp2", "Angptl4","Acadl")
 PPAR_Coactivators<- c(  "Rxra",     # heterodimer partner
   "Ppargc1a", # PGC-1α
   "Ppargc1b", # PGC-1β
   "Nr1h4",    # FXR
   "Srebf1",
   "Hnf4a"
 )
 values_for_heatmap<-as.data.frame(Proteins%>%
                                     select(Name,Genes, contains(c("Etoh","Tam")))%>%
                                     filter(Genes %in% c(PPARa_FAoxidation,PPARg_lipogenesis,PPARdb_FAoxidation_Energy,PPAR_Coactivators)))
 rownames(values_for_heatmap)<-values_for_heatmap$Name
 row_category <- case_when(
   values_for_heatmap$Genes %in% PPARa_FAoxidation ~ "PPARa_FAoxidation",
   #values_for_heatmap$Genes %in% PPARdb_FAoxidation_Energy ~ "PPARdb_FAoxidation_Energy",
   values_for_heatmap$Genes %in% PPAR_Coactivators ~ "PPAR_Coactivators",
   values_for_heatmap$Genes %in% PPARg_lipogenesis ~ "PPARg_lipogenesis")
 
 row_category <- factor(
   row_category,
   levels = c("PPARa_FAoxidation","PPAR_Coactivators",   #,"PPARdb_FAoxidation_Energy"
              "PPARg_lipogenesis"   ) )
 
 
 ann <- data.frame(Sample = colnames(values_for_heatmap)) %>%
   separate(Sample, into = c("Sex", "Treatment", "Replicate"), sep = "_") %>%
   select(-Replicate)%>%mutate(
     Sex = case_when(Sex == "F" ~ "female",  Sex == "M" ~ "male")  )
 
 
 

 
 annotation_row <- data.frame(
   PPARa_FAoxidation = values_for_heatmap$Genes %in% PPARa_FAoxidation,
   PPARg_lipogenesis = values_for_heatmap$Genes %in% PPARg_lipogenesis,
   #PPARdb_FAoxidation_Energy = values_for_heatmap$Genes %in% PPARdb_FAoxidation_Energy,
   PPAR_Coactivators = values_for_heatmap$Genes %in% PPAR_Coactivators,
   row.names = values_for_heatmap$Genes
 )
 annotation_row[] <- lapply(annotation_row, factor, levels = c(FALSE, TRUE))
 
 annotation_colors <- list(
   Sex = Sex_colors,
   Treatment = Treatment_colors,
   Diet = Diet_colors )
 
 annotation_colors$PPARg_lipogenesis <- c("FALSE" = "grey95", "TRUE" = "#8dd3c7")
 annotation_colors$PPARa_FAoxidation <- c("FALSE" = "grey95", "TRUE" = "#ffffb3")
 annotation_colors$PPAR_Coactivators <- c("FALSE" = "grey95", "TRUE" = "orange2")
 #annotation_colors$PPARdb_FAoxidation_Energy <- c("FALSE" = "grey95", "TRUE" = "darkblue")
 
 rownames(annotation_row) <- rownames(values_for_heatmap)
 rownames(ann) <- colnames(values_for_heatmap)
 heatmap_height <-nrow(heatmap_data)/7+3
 
 ord <- order(row_category)
 
 values_for_heatmap <- values_for_heatmap[ord, ]
 annotation_row <- annotation_row[ord, ]
 values_for_heatmap<-values_for_heatmap%>%select(-Name,-Genes)
 gaps_row <- cumsum(table(row_category[ord]))[-length(table(row_category))]
 p_heat<- pheatmap::pheatmap(values_for_heatmap,
                             scale = "row",
                             cluster_rows = FALSE, 
                             cluster_cols = TRUE,
                             annotation_col  = ann,
                             annotation_row = annotation_row,
                             annotation_names_row = FALSE,
                             annotation_colors = annotation_colors,
                             show_colnames = FALSE,
                             treeheight_row = 0,  
                             treeheight_col = 5,  
                             gaps_row = gaps_row,
                             main = "PPARs")
 
 print(p_heat)
 ggsave(plot= p_heat, filename= "Heatmap_PPAR_clustered.png",limitsize = FALSE, 
        width= 9, height = heatmap_height, dpi = 500,bg = "white",path=paste0(proteom_output_pwd))
 ########
 Metabolic_Pathways <- list(
   
   "Purine metabolism" = c(
     "Ppat","Gart","Pfas","Paics","Adsl","Atic",
     "Impdh1","Impdh2","Gmps","Ampd1","Ampd2","Ampd3",
     "Nt5e","Ada","Xdh","Hprt","Aprt","Adk","Ak1","Ak2",
     "Nme1","Nme2","Pnp","Entpd1","Entpd2","Entpd5",
     "Gda","Adenosine kinase","Itpa","Enpp1"          ))
 #   ),
 #   
 #   "Pyrimidine metabolism" = c(
 #     "Cad","Dhodh","Umps","Tk1","Tk2","Tyms","Dtymk",
 #     "Cmpk1","Cmpk2","Nme1","Nme2","Uck1","Uck2",
 #     "Upp1","Upp2","Dck","Nt5c","Nt5c2","Cda","Dpyd",
 #     "Dctd","Rrm1","Rrm2","Rrm2b"
 #   ),
 #   
 #   "D-Amino acid metabolism" = c(
 #     "Dao","Ddo","Got1","Got2","Gpt","Gpt2","Ddost"
 #   ),
 #   
 #   "Arginine and proline metabolism" = c(
 #     "Arg1","Arg2","Ass1","Asl","Otc","Cps1",
 #     "Oat","Prodh","Pycr1","Pycr2","Pycrl","Aldh18a1",
 #     "Nos1","Nos2","Nos3","Sat1","Smox","Sms","Srm",
 #     "Amd1","Mat1a","Mat2a"
 #   ),
 #   
 #   "Glycine, serine and threonine metabolism" = c(
 #     "Phgdh","Psat1","Psph",
 #     "Shmt1","Shmt2",
 #     "Gldc","Amt","Gcsh",
 #     "Bhmt","Bhmt2","Gnmt","Sardh","Dmgdh",
 #     "Mthfd1","Mthfd2","Mthfd1l","Mthfr","Dhfr",
 #     "Thnsl1","Thnsl2"
 #   ),
 #   
 #   "Pyruvate metabolism" = c(
 #     "Pdha1","Pdhb","Dlat","Dld",
 #     "Pdk1","Pdk2","Pdk3","Pdk4",
 #     "Pdp1","Pdp2",
 #     "Pc","Pcx",
 #     "Ldha","Ldhb",
 #     "Me1","Me2","Me3",
 #     "Acly","Acss1","Acss2"
 #   ),
 #   
 #   "Butanoate metabolism" = c(
 #     "Acss1","Acss2","Oxct1","Bdh1","Bdh2",
 #     "Hmgcs2","Hmgcl","Hadha","Hadhb",
 #     "Echs1","Acat1","Acat2","Acadm","Acadl","Acadvl"
 #   ),
 #   
 #   "Starch and sucrose metabolism" = c(
 #     "Gaa","Pygm","Pygl","Pygb",
 #     "Gys1","Gys2",
 #     "Agl","Gbe1",
 #     "Hk1","Hk2","Hk3","Gpi1",
 #     "Pgm1","Pgm2","Ugp2"
 #   ),
 #   
 #   "Phenylalanine metabolism" = c(
 #     "Pah","Tat","Hpdl",
 #     "Got1","Got2",
 #     "Aldh1a1","Aldh3a1","Aldh3b1",
 #     "Maoa","Maob","Ddc","Comt"
 #   ),
 #   
 #   "Glyoxylate and dicarboxylate metabolism" = c(
 #     "Agxt","Grhpr","Hao1","Hao2",
 #     "Ldha","Mdh1","Mdh2",
 #     "Got1","Got2","Gpt","Gpt2"
 #   ),
 #   
 #   "One carbon pool by folate" = c(
 #     "Shmt1","Shmt2",
 #     "Mthfd1","Mthfd1l","Mthfd2",
 #     "Mthfr","Dhfr","Tyms","Atic",
 #     "Gart","Mtr","Mtrr","Fpgs"
 #   ),
 #   
 #   "Histidine metabolism" = c(
 #     "Hal","Uroc1","Amdhd1",
 #     "Ftcd","Hnmt","Abp1","Aoc1"
 #   ),
 #   
 #   "Pentose phosphate pathway" = c(
 #     "G6pdx","Pgls","Pgd",
 #     "Rpia","Rpe","Tkt","Taldo1",
 #     "Prps1","Prps2","Rbks"
 #   ),
 #   
 #   "Citrate cycle (TCA cycle)" = c(
 #     "Cs",
 #     "Aco1","Aco2",
 #     "Idh1","Idh2","Idh3a","Idh3b","Idh3g",
 #     "Ogdh","Dlst","Dld",
 #     "Suclg1","Sucla2","Suclg2",
 #     "Sdha","Sdhb","Sdhc","Sdhd",
 #     "Fh1",
 #     "Mdh1","Mdh2",
 #     "Pdha1","Pdhb","Dlat",
 #     "Pc","Me1","Me2","Me3"
 #   ),
 #   
 #   "Nitrogen metabolism" = c(
 #     "Glul","Gls","Gls2",
 #     "Glud1",
 #     "Car2","Car3","Ca2","Ca4",
 #     "Ass1","Asl","Cps1","Otc"
 #   ),
 #   
 #   "Alanine, aspartate and glutamate metabolism" = c(
 #     "Got1","Got2",
 #     "Gpt","Gpt2",
 #     "Glud1",
 #     "Asns",
 #     "Gls","Gls2",
 #     "Glul",
 #     "Aspa","Asl","Ass1","Ddo"
 #   ),
 #   
 #   "Valine, leucine and isoleucine biosynthesis" = c(
 #     "Bcat1","Bcat2",
 #     "Bckdha","Bckdhb","Dbt","Dld",
 #     "Acadm","Acad8","Ivd","Hmgcs2"
 #   ),
 #   
 #   "Phenylalanine, tyrosine and tryptophan metabolism" = c(
 #     "Pah","Tat","Tyr",
 #     "Ddc","Comt",
 #     "Maoa","Maob",
 #     "Tph1","Tph2","Ido1","Ido2","Tdo2","Kmo"
 #   ),
 #   
 #   "Arginine biosynthesis" = c(
 #     "Cps1","Nags","Otc",
 #     "Ass1","Asl",
 #     "Arg1","Arg2",
 #     "Nos1","Nos2","Nos3",
 #     "Oat","Aldh18a1"
 #   )
 #   
 # )

 Metabolic_Pathways <- unique(unlist(Metabolic_Pathways))
 values_for_heatmap <- as.data.frame( Proteins %>%
     select(Name, Genes, contains(c("Etoh", "Tam"))) %>%
     filter(Genes %in% Metabolic_Pathways))
 
 rownames(values_for_heatmap) <- values_for_heatmap$Name
 

 
 row_category <- sapply(values_for_heatmap$Genes, function(x){
   names(Metabolic_Pathways)[sapply(Metabolic_Pathways, function(y) x %in% y)][1]
 }
 )
 
 row_category <- factor(row_category,levels = names(Metabolic_Pathways))
 annotation_row <- data.frame( row.names = values_for_heatmap$Genes)
   
   for(i in names(Metabolic_Pathways)){
     annotation_row[[i]] <- factor( values_for_heatmap$Genes %in% Metabolic_Pathways[[i]],
                                     levels = c(FALSE, TRUE))
   }
 annotation_row<-
   rownames(annotation_row) <- rownames(values_for_heatmap)
   library(RColorBrewer)
   
   pathway_cols <- colorRampPalette(
     brewer.pal(12, "Set3")
   )(length(Metabolic_Pathways))
   
   annotation_colors <- list(
     Sex = Sex_colors,
     Treatment = Treatment_colors,
     Diet = Diet_colors
   )
   
   for(i in seq_along(names(Metabolic_Pathways))){
     annotation_colors[[names(Metabolic_Pathways)[i]]] <-
       c("FALSE" = "grey95",
         "TRUE" = pathway_cols[i])
   }
   ord <- order(row_category)
   
   values_for_heatmap <- values_for_heatmap[ord, ]
   annotation_row <- annotation_row[ord, ]
   
   values_for_heatmap <- values_for_heatmap %>%
     select(-Name, -Genes)
   ann <- data.frame(Sample = colnames(values_for_heatmap)) %>%
     separate(Sample, into = c("Sex", "Treatment", "Replicate"), sep = "_") %>%
     select(-Replicate)%>%mutate(
       Sex = case_when(Sex == "F" ~ "female",  Sex == "M" ~ "male")  )
   #gaps_row <- cumsum(table(row_category[ord]))[-length(table(row_category))]
 
   p_metabolism <- pheatmap::pheatmap(
     values_for_heatmap,
     scale = "row",
     cluster_rows = FALSE,
     cluster_cols = TRUE,
     annotation_col = ann,
     #annotation_row = annotation_row,
     annotation_names_row = FALSE,
     annotation_colors = annotation_colors,
     show_colnames = FALSE,
     treeheight_row = 0,
     treeheight_col = 5,
     #gaps_row = gaps_row,
     main = "Metabolic pathways"
   )
   
   print(p_metabolism)
   
   ggsave(
     plot = p_metabolism,
     filename = "Heatmap_Metabolic_pathways_clustered.png",
     limitsize = FALSE,
     width = 10,
     height = nrow(values_for_heatmap)/7 + 3,
     dpi = 500,
     bg = "white",
     path = proteom_output_pwd
   )