gc()
rm(list=ls())
library(tidyverse)
library(phyloseq)
library(tools)
library(lme4)      
library(lmerTest)   
library(car)        
library(emmeans)
library(gridExtra)
library(DESeq2)
library(pheatmap)
# Load Data and Prepare Phyloseq Object------

setwd("C:/Users/b1084855/OneDrive - Universität Salzburg/AG_Tumorimmunologie - Dokumente/Data/Freia Krause/01_Experiments/FK49_CD-HFD_13wks/FK49_Analysis")
raw_data_pw <- "C:/Users/b1084855/OneDrive - Universität Salzburg/AG_Tumorimmunologie - Dokumente/Data/Freia Krause/01_Experiments/FK49_CD-HFD_13wks/FK49_Microbiome/"

#Give Info were files are, read files and combine into dataframe
files <- list.files(paste0(raw_data_pw,"16S"), full.names = TRUE)

df_list <- lapply(files, function(f) {
  read.csv(f, sep = ";",header = TRUE,stringsAsFactors = FALSE,   check.names = FALSE  )})


df_list <- Map(function(df, fname) {
  df$SampleID <- file_path_sans_ext(basename(fname))
  df  },
  df_list, files)
all_df <- do.call(rbind, df_list)

colnames(all_df)[1] <- "Sequence" #dname first col sequence 
#sample id was defined from file name. but filnema contained date and extra info
#i only want to keep iAL(number)_Fnumber
all_df<- all_df %>% mutate(SampleID = gsub(".*iAL([0-9]+)-([A-Z0-9]+)-.*", "\\1_\\2", SampleID)) 
clean_df <- all_df %>%
  group_by(SampleID) %>%
  filter(n() > 20) %>%  # 3 Samples have below 20 Genus/OTUs
  ungroup()

# Define Count Matrix, Taconomy Matrix and Meta Data fophylosseq object
##Count Matrix
count_mat <- clean_df %>%
  select(Sequence, SampleID, abundance) %>%
  pivot_wider(names_from = SampleID,
              values_from = abundance,
              values_fill = 0)%>%
  column_to_rownames(var = "Sequence") %>%
  as.matrix()
## later on we us apply as.numeric onmatrix, this removes rownames(this are the seqs).
##I want to save the seqs and give them a shorter ID for handling here and fuse them back to the count matric after as.numeric 

sequences<-rownames(count_mat)
sequence_ids <- paste0("OTU", seq_along(sequences))
seq_info <- data.frame(
  seq_ids = sequence_ids,
  Sequence = sequences,
  stringsAsFactors = FALSE)
write.csv(seq_info,file= "02_GeneratedData/FK49_Microbiome/Seq_Info.csv" )
rm(seq_info)

count_mat<-  apply(count_mat, 2, as.numeric)
rownames(count_mat) <- sequence_ids

##Taxonomy Matrix
tax_mat_unique <- clean_df %>%
  select(Sequence, Kingdom:Species) %>%
  distinct(Sequence, .keep_all = TRUE)%>%
  column_to_rownames(var = "Sequence") %>%
    as.matrix()

#some seqs were used with different taxonomy data. I dont know about that. 
#So i kept the seqs with tax data it first appears with
tax_mat_unique <- tax_mat_unique[sequences, ]  # ensure the order matches sequences
rownames(tax_mat_unique) <- sequence_ids

##Meta Data
metadata <- read.csv(paste0(raw_data_pw,"/FK49_CD-HFD_13wks_Microbiome_Meta.csv"),sep = ";")%>%
  mutate(T_T_D_S= paste0(Treatment,"_",Feces,"_",Diet,"_",Sex))%>%
  mutate(T_T_D= paste0(Treatment,"_",Feces,"_",Diet))%>%
  mutate(T_T= paste0(Treatment,"_",Feces))%>%
  mutate(T_S= paste0(Treatment,"_",Sex))
rownames(metadata) <- metadata$SampleID
## Combine Counts, Tax and Meta in phyloseq object
ps <- phyloseq(
  otu_table(count_mat, taxa_are_rows = TRUE),
  tax_table(tax_mat_unique),
  sample_data(metadata) )

ps<- prune_taxa(taxa_sums(ps) > 0, ps) # taxa with remove 0 counts
# Do Quality Check or Something? ------

table(tax_table(ps)[, "Phylum"], exclude = NULL)

#for some measures we need absolute counts and for some we need relative counts /abundance ps and ps_rel
ps <- subset_taxa(ps, !is.na(Phylum) & !Phylum %in% c("", "uncharacterized"))
ps_rel <- transform_sample_counts(ps, function(x) x / sum(x))

rm(all_df,clean_df,count_mat,df_list,tax_mat_unique,files,sequence_ids,sequences)
#https://bioconductor.org/help/course-materials/2017/BioC2017/Day1/Workshops/Microbiome/MicrobiomeWorkflowII.html#using_phyloseq
#### Compute prevalence of each feature, store as data.frame ----
#Test something so that i can do quality check. 
  prevdf = apply(X = otu_table(ps),
                 MARGIN = ifelse(taxa_are_rows(ps), yes = 1, no = 2),
                 FUN = function(x){sum(x > 0)})
  # Add taxonomy and total read counts to this data.frame
  prevdf = data.frame(Prevalence = prevdf,
                      TotalAbundance = taxa_sums(ps),
                      tax_table(ps))
  plyr::ddply(prevdf, "Phylum", function(df1){cbind(mean(df1$Prevalence),sum(df1$Prevalence))})
  # Subset to the remaining phyla
  prevdf1 = subset(prevdf, Phylum %in% get_taxa_unique(ps, "Phylum"))
  ggplot(prevdf1, aes(TotalAbundance, Prevalence / nsamples(ps),color=Phylum)) +
   geom_hline(yintercept = 0.025
              , alpha = 0.5, linetype = 2) +  geom_point(size = 2, alpha = 0.7) +
   scale_x_log10() +  xlab("Total Abundance") + ylab("Prevalence [Frac. Samples]") +
   facet_wrap(~Phylum) + theme(legend.position="none")

  #Define prevalence threshold as 5% of total samples
  prevalenceThreshold = 0.025 * nsamples(ps)
  prevalenceThreshold
  keepTaxa = rownames(prevdf1)[(prevdf1$Prevalence >= prevalenceThreshold)]
  ps2 = prune_taxa(keepTaxa, ps)
  # How many genera would be present after filtering?
  length(get_taxa_unique(ps2, taxonomic.rank = "Genus"))
  ps3 = tax_glom(ps2, "Genus", NArm = TRUE)

  plot_abundance = function(physeq,title = "",
                            Facet = "Order", Color = "Phylum"){
    # Arbitrary subset, based on Phylum, for plotting
    p1f = subset_taxa(physeq, Phylum %in% c("Firmicutes"))
    mphyseq = psmelt(p1f)
    mphyseq <- subset(mphyseq, Abundance > 0)
    ggplot(data = mphyseq, mapping = aes_string(x = "Sex",y = "Abundance",
                                                  color = Color, fill = Color)) +
      geom_violin(fill = NA) +
      geom_point(size = 1, alpha = 0.3,
                position = position_jitter(width = 0.3)) +
     facet_wrap(facets = Facet) + scale_y_log10()+
     theme(legend.position="none")
  }
  # Transform to relative abundance. Save as new object.
  ps3ra = transform_sample_counts(ps3, function(x){x / sum(x)})
  plotBefore = plot_abundance(ps3,"")
  plotAfter = plot_abundance(ps3ra,"")
  # Combine each plot into one graphic.
  grid.arrange(nrow = 2,  plotBefore, plotAfter)
  psOrd = subset_taxa(ps3ra, Order == "Lactobacillales")
  plot_abundance(psOrd, Facet = "Genus", Color = NULL)
  qplot(sample_data(ps)$Week, geom = "histogram",binwidth=1) + xlab("Week")
  qplot(log10(rowSums(otu_table(ps))),binwidth=0.2) +
    xlab("Logged counts-per-sample")
rm(plotAfter,plotBefore,prevalenceThreshold,keepTaxa,psOrd,ps3ra,ps3,ps2,prevdf,prevdf1,plot_abundance)
gc()

#https://academic.oup.com/proteincell/article/14/10/713/7147618
#Diverstiy Analysis
##alpha: Richness(number) or distribution (evenness) of the sample, rarefaction curve for richness
##beta similarity or dissmilarity of two communities
###bray-curtis dissimilariy (size as ovrall abundance per sample and shape as abundance of each taxon of the communites)
#Difference Analysis
#Biomarker Identification
#Network Analysis
#Function predicton

# Alpha diverstiy has to be calulated on the raw counts, not normalized counts
a_diversity_measures <-c("Observed", "Chao1", "ACE", "Shannon", "Simpson", "InvSimpson")
alpha_diversity <- estimate_richness(ps,   measures = a_diversity_measures)

alpha_diversity <- alpha_diversity %>%
  rownames_to_column(var = "SampleID") %>%   #sampleID was shitty formatted here with X in fromnt. dont know why. formatted it to my style
  mutate(SampleID = gsub("^X", "", SampleID))
#Join with meta data since i dont know how to flexible tell ggplot otherwise which sample has which info. Might be in phylo object
alpha_with_meta <- alpha_diversity %>%
  left_join(metadata, by = "SampleID")

hist(alpha_with_meta$Shannon)# looks to me like roughly normal so no transformation
rm(alpha_diversity)
#Plot alpha diversity calculated with different methods -----
for(a in a_diversity_measures){
  
  #Generate Summary statistics for the plot Mean, SD 
  alpha_summary <- alpha_with_meta %>%
    group_by(Feces, Treatment) %>%
    summarise(mean = mean(.data[[a]], na.rm = TRUE), 
              sd = sd(.data[[a]], na.rm = TRUE),
              n=n() ,
              .groups = "drop") 
  # Fit the linear mixed-effects model 
  model_alpha <- lmer(formula = as.formula(paste0(a, " ~ Treatment * Feces + (1|Animal)")),
                    data = alpha_with_meta)
  
  anova_table_alpha <-anova(model_alpha, type = 3)
  
  anova_label_alpha <- paste0("ANOVA over linear mixed-effects model\n",
                        "Treatment: F = ", round(anova_table_alpha$F[1], 2), ", p = ", signif(anova_table_alpha$`Pr(>F)`[1], 3), "\n",
                        "Time: F = ", round(anova_table_alpha$F[2], 2), ", p < ", format.pval(anova_table_alpha$`Pr(>F)`[2], digits = 1), "\n",
                        "Interaction: F = ", round(anova_table_alpha$F[3], 2), ", p = ", signif(anova_table_alpha$`Pr(>F)`[3], 3))


# Calculate the estimated marginal means for both Treatment and wks_diet
  emm_alpha <- emmeans(model_alpha, ~ Treatment | Feces)

  plot(model_alpha)          # residuals vs fitted
  qqnorm(resid(model_alpha))
  qqline(resid(model_alpha))


# Perform pairwise contrasts between Treatment levels (ctrl vs TAM) at all levels of wks_diet
  pwc_alpha <- contrast(emm_alpha, method = "pairwise", adjust = "bonferroni")
  posthoc_label_alpha<- "Post Hoc: Pairwise with Bonferroni correction"
  test_label<- paste0(anova_label_alpha,"\n",posthoc_label_alpha)

  pwc_df_alpha <- as.data.frame(pwc_alpha)
  pwc_df_rounded_alpha <- pwc_df_alpha %>%
   mutate(rounded_p_value = ifelse(is.na(p.value), "NA", round(p.value, 3))) %>%
   mutate(significance = case_when(
     is.na(p.value) ~ "NA",                     # For NA p-values
     p.value < 0.001 ~ "***",                   # p < 0.001 is highly significant
     p.value >= 0.001 & p.value < 0.01 ~ "**",  # 0.001 ≤ p < 0.01 is significant
     p.value >= 0.01 & p.value < 0.05 ~ "*",    # 0.01 ≤ p < 0.05 is moderately significant
     p.value >= 0.05 ~ "NS",                    # p ≥ 0.05 is not significant
     TRUE ~ "NA"   ))  %>%                      # Default case
   select(Feces, rounded_p_value, significance)


##Plot and save all a diverstiy measures -----

    p<-ggplot() +
    geom_jitter(data=alpha_with_meta, 
                aes(x = Feces,  y = .data[[a]],   color = Treatment),
                    width = 0.1, size = 3, alpha = 0.8) +  
    geom_line(data = alpha_summary,
            aes(x = Feces, y = mean, group = Treatment,color= Treatment),
             size = 1 ) +
    geom_ribbon(data = alpha_summary,
      aes(x = Feces, fill = Treatment,group = Treatment,
          ymin = mean-sd,   
          ymax =mean+sd),  
      alpha = 0.1) +
   geom_text(data = pwc_df_rounded_alpha,
    aes(x = c(1,2,3,4,5),
      y = alpha_summary$mean[c(2, 3, 6, 8,10)] + 
          alpha_summary$sd[c(2, 3, 6, 8,10)]+
          0.02*alpha_summary$mean[c(2, 3, 6, 8,10)]  ,
      label = significance),  size = 2.5,  fontface = "italic") +
  labs(  title = "Diversity Over Time",  x = "Time (Feces)",  y = a)+
  scale_color_manual(values = c("EtOH" = "#4D4D4DBF", "TAM" = "#8B0000BF")) +
    scale_fill_manual(values = c("EtOH" = "#4D4D4DBF", "TAM" = "#8B0000BF")) +  
    theme_minimal() +
      theme(legend.position = "bottom",
            axis.line = element_line(color = "black", linewidth = 0.5),
            axis.ticks = element_line(color = "black", linewidth = 0.5),
            axis.title = element_text(size = 20, face = "bold"),
            axis.title.x = element_blank(),
            axis.text = element_text(size = 19, face = "bold"),
            plot.title = element_blank(),
            legend.title = element_text(size = 10),
            legend.text = element_text(size = 10),
            panel.grid = element_blank())+
      annotate("text", x = 1, y = min(alpha_summary$mean)-2*min(alpha_summary$sd),  # 5% above the bottom
               label = test_label,
               size = 2,
               hjust = 0,
               color = "black", 
               fontface = "italic")
      
      
  
  ggsave(plot=p,filename= paste0("02_GeneratedData/FK49_Microbiome/FK49_Microbiome_aDiversity_",a,".png"),
         width = 9, height = 6, dpi = 300)
  rm(p,alpha_summary,anova_label_alpha,anova_table_alpha,emm_alpha,model_alpha,pwc_df_alpha,pwc_alpha,pwc_df_rounded_alpha)  
}


rm(alpha_with_meta)
gc()
# Plot beta diversity -----


## Function to PLot PCoAs ----
plot_PCoA<-function(dataset,
                    shapedots=NULL,
                    colordots=NULL,
                    linetypeellipse=NULL,
                    groupellipse=NULL,
                    colorellipse=NULL,
                    limitscolor=c("ND","CDHFD","EtOH", "TAM"),
                    valuescolor= NULL,
                    limitsshape= c("female","male"),
                    valuesshape=c( female =8, male   = 17),
                    valueslinetype = c( female = "dashed",  male   = "solid"),
                    limitslinetype= c("female", "male")
                    ){
  message("Dataset received: ", dataset)

  ## Subsetting dataset in different interesting sets ----
  fullset <-switch(
    dataset,
   "ps_rel"              = ps_rel,
   "psrel_CDHFD"         = ps_rel,
    "psrel_ND"            = ps_rel,
    "psrel_male"          = ps_rel,
    "psrel_male_CDHFD"    = ps_rel,
    "psrel_female"        = ps_rel,
    "psrel_female_CDHFD"  = ps_rel,
    "ps"                 = ps,
     "ps_CDHFD"            = ps,
     "ps_ND"               = ps,
    "ps_male"             = ps,
     "ps_male_CDHFD"       = ps,
     "ps_female"           = ps,
     "ps_female_CDHFD"     = ps,
    stop("Dataset not recognized")
  )

  DATA <- switch(
    dataset,
    "ps_rel"              = ps_rel,
    "psrel_CDHFD"         = prune_samples(sample_data(fullset)$Diet_short  %in% c("CDHFD"),fullset),
    "psrel_ND"            = prune_samples(sample_data(fullset)$Diet_short  %in% c("ND")  ,fullset),
    "psrel_male"          = prune_samples(sample_data(fullset)$Sex         %in% c("male"),fullset),
    "psrel_male_CDHFD"    = prune_samples(sample_data(fullset)$Sex         %in% c("male") & 
                                                  sample_data(fullset)$Diet_short     %in% c("CDHFD"),fullset),
    "psrel_female"        = prune_samples(sample_data(fullset)$Sex         %in% c("female"),fullset),
    "psrel_female_CDHFD"  = prune_samples(sample_data(fullset)$Sex         %in% c("female") & 
                                                     sample_data(fullset)$Diet_short %in% c("CDHFD"),fullset),
    "ps"                  = ps,
    "ps_CDHFD"            = prune_samples(sample_data(fullset)$Diet_short  %in% c("CDHFD"),fullset),
    "ps_ND"               = prune_samples(sample_data(fullset)$Diet_short  %in% c("ND")  ,fullset),
    "ps_male"             = prune_samples(sample_data(fullset)$Sex         %in% c("male"),fullset),
    "ps_male_CDHFD"       = prune_samples(sample_data(fullset)$Sex         %in% c("male") & 
                                                     sample_data(fullset)$Diet_short  %in% c("CDHFD"),fullset),
    "ps_female"           = prune_samples(sample_data(fullset)$Sex         %in% c("female"),fullset),
    "ps_female_CDHFD"     = prune_samples(sample_data(fullset)$Sex         %in% c("female") & 
                                                     sample_data(fullset)$Diet_short  %in% c("CDHFD"), fullset),
     stop("Dataset not recognized")
     )

  message("Dataset after pruning: ", class(DATA))
  ordination <-ordinate(DATA, method="PCoA", distance="bray")
  eig <- ordination$values$Relative_eig
  pc1 <- round(eig[1] * 100, 2)
  pc2 <- round(eig[2] * 100, 2)

  PCoA_plot<- plot_ordination(DATA, ordination, color = colordots, shape = shapedots)+
    geom_point(size = 2.5) +
    stat_ellipse(aes_string(color= colorellipse,  linetype = linetypeellipse ,       group = groupellipse), 
                 type = "norm", linewidth = 0.8)+
    labs(title = "PCoA",  x = paste0("PCoA1 [",pc1,"%]"),  y = paste0("PCoA2 [",pc2,"%]"))+
    theme_minimal() +
    theme(legend.position = "bottom",
          axis.line = element_line(color = "black", linewidth = 0.5),
          axis.ticks = element_line(color = "black", linewidth = 0.5),
          axis.title = element_text(size = 16, face = "bold"),
          axis.text = element_text(size = 16, face = "bold"),
          plot.title = element_blank(),
          legend.title = element_text(size = 10),
          legend.text = element_text(size = 10),
          panel.grid = element_blank())+
    scale_color_manual( limits = limitscolor, values = valuescolor)+
    scale_shape_manual(values = valuesshape,  limits =limitsshape) +
    scale_linetype_manual(values =valueslinetype ,limits = limitslinetype)
  
  ggsave(plot=PCoA_plot,
       filename= paste0("02_GeneratedData/FK49_Microbiome/BDiversity_PCoA_",
                        dataset,
                        "_.png"),
       width = 9, height = 6, dpi = 300)
  return(PCoA_plot)
  rm(DATA,fullset, PCoA_plot,eig, pc1, pc2, ordination)
  gc()
}
## Plot PCoA of different data subsets -----
plot_PCoA(dataset="ps_rel", shapedots=NULL,  colordots="Treatment",colorellipse = "Diet_short",    groupellipse = "Diet_short",
          linetypeellipse= NULL,  limitscolor= c("ND","CDHFD","EtOH", "TAM"), valuescolor= c("darkorange3","darkviolet",adjustcolor(c("#4D4D4D", "#8B0000"), alpha.f = 1)),
          limitsshape= c("female","male"),valuesshape= c( female =8, male   = 17),valueslinetype = c( female = "dashed",  male   = "solid"),  limitslinetype= c("female", "male") )

plot_PCoA(dataset="psrel_CDHFD", shapedots="Sex",colordots="Treatment",colorellipse = "Treatment", groupellipse = "T_S",
  linetypeellipse= "Sex",  limitscolor = c("EtOH","TAM"),  valuescolor= c(adjustcolor(c("#4D4D4D", "#8B0000"), alpha.f = 0.7)),
   limitsshape=c("female","male") ,  valuesshape= c( female =8, male   = 17),  valueslinetype = c( female = "dashed",  male   = "solid"),  limitslinetype= c("female", "male"))  

plot_PCoA(  dataset="psrel_male_CDHFD",  shapedots="Feces",  colordots="Treatment",  colorellipse = "Treatment",  groupellipse = "Treatment",
  linetypeellipse= NULL,  limitscolor = c("EtOH","TAM"),  valuescolor= c(adjustcolor(c("#4D4D4D", "#8B0000"), alpha.f = 0.7)),
   limitsshape=c("F3","F4","F5") ,  valuesshape= c( "F3" =15, "F4" = 16,"F5"= 8),  valueslinetype = NULL,  limitslinetype= NULL)       

plot_PCoA(  dataset="psrel_female",  shapedots="Feces",  colordots="Treatment",  colorellipse = "Treatment",  groupellipse = "Treatment",
  linetypeellipse= NULL,  limitscolor = c("EtOH","TAM"),  valuescolor= c(adjustcolor(c("#4D4D4D", "#8B0000"), alpha.f = 0.7)),
  limitsshape=c("F3","F4","F5") ,  valuesshape= c( "F3" =15, "F4" = 16,"F5"= 8),  valueslinetype = NULL,  limitslinetype= NULL)   



# Barplots 



plot_Micro_Bars <-function (dataset=NULL, GROUP_MEAN= FALSE, subsetTAX = NULL, RANK = NULL, FILL = NULL, savename = NULL, WRAP = NULL, WIDTH = 6.5, HEIGHT= 9){
  message("Dataset received: ", dataset)
  
  ## Subsetting dataset in different interesting sets ----
  fullset <-switch(
    dataset,
    "ps_rel"              = ps_rel,
    "psrel_CDHFD"         = ps_rel,
    "psrel_ND"            = ps_rel,
    "psrel_male"          = ps_rel,
    "psrel_male_CDHFD"    = ps_rel,
    "psrel_female"        = ps_rel,
    "psrel_female_CDHFD"  = ps_rel,
    "ps"                 = ps,
    "ps_CDHFD"            = ps,
    "ps_ND"               = ps,
    "ps_male"             = ps,
    "ps_male_CDHFD"       = ps,
    "ps_female"           = ps,
    "ps_female_CDHFD"     = ps,
    stop("Dataset not recognized")
  )
  
  DATA <- switch(
    dataset,
    "ps_rel"              = ps_rel,
    "psrel_CDHFD"         = prune_samples(sample_data(fullset)$Diet_short  %in% c("CDHFD"),fullset),
    "psrel_ND"            = prune_samples(sample_data(fullset)$Diet_short  %in% c("ND")  ,fullset),
    "psrel_male"          = prune_samples(sample_data(fullset)$Sex         %in% c("male"),fullset),
    "psrel_male_CDHFD"    = prune_samples(sample_data(fullset)$Sex         %in% c("male") & 
                                            sample_data(fullset)$Diet_short     %in% c("CDHFD"),fullset),
    "psrel_female"        = prune_samples(sample_data(fullset)$Sex         %in% c("female"),fullset),
    "psrel_female_CDHFD"  = prune_samples(sample_data(fullset)$Sex         %in% c("female") & 
                                            sample_data(fullset)$Diet_short %in% c("CDHFD"),fullset),
    "ps"                  = ps,
    "ps_CDHFD"            = prune_samples(sample_data(fullset)$Diet_short  %in% c("CDHFD"),fullset),
    "ps_ND"               = prune_samples(sample_data(fullset)$Diet_short  %in% c("ND")  ,fullset),
    "ps_male"             = prune_samples(sample_data(fullset)$Sex         %in% c("male"),fullset),
    "ps_male_CDHFD"       = prune_samples(sample_data(fullset)$Sex         %in% c("male") & 
                                            sample_data(fullset)$Diet_short  %in% c("CDHFD"),fullset),
    "ps_female"           = prune_samples(sample_data(fullset)$Sex         %in% c("female"),fullset),
    "ps_female_CDHFD"     = prune_samples(sample_data(fullset)$Sex         %in% c("female") & 
                                            sample_data(fullset)$Diet_short  %in% c("CDHFD"), fullset),
    stop("Dataset not recognized")
  )
  
  if(!is.null(subsetTAX)){
    DATA <- switch (subsetTAX, 
                    "Firmicutes" =  subset_taxa(DATA, Phylum == "Firmicutes"),
                    "Lactobacillales" =  subset_taxa(DATA, Order == "Lactobacillales"),
                    "Lactobacillus" =  subset_taxa(DATA, Genus  == "Lactobacillus"),
                    "Clostridiales" =  subset_taxa(DATA, Order  == "Clostridiales")
                    
           )
  }
 
  message("Dataset after pruning: ", class(DATA))
  DATA = tax_glom(DATA, taxrank= RANK , NArm=FALSE)
  if(GROUP_MEAN == TRUE){
    desired_order <- c(
      "EtOH_F1", "TAM_F1",
      "EtOH_F2", "TAM_F2",
      "EtOH_F3", "TAM_F3",
      "EtOH_F4", "TAM_F4",
      "EtOH_F5", "TAM_F5")
    
  df <- psmelt(DATA)
  df_sum <- df %>%
    group_by(T_T, .data[[FILL]], Treatment,Feces) %>%
    summarise(
      mean_abund = mean(Abundance),
      sd_abund   = sd(Abundance),
      .groups = "drop"
    )
  df_sum$T_T <- factor(df_sum$T_T, levels = desired_order)
  
  # Facet order by max mean abundance
  facet_order <- df_sum %>%
    group_by(.data[[FILL]]) %>%
    summarise(max_abund = max(mean_abund, na.rm = TRUE)) %>%
    arrange(desc(max_abund)) %>%
    pull(.data[[FILL]])
  
  df_sum[[FILL]] <- factor(df_sum[[FILL]], levels = facet_order)
  df[[FILL]]     <- factor(df[[FILL]],     levels = facet_order)
  
  barplot <- ggplot(df_sum,  aes(x = Feces,   y = mean_abund,  fill = Treatment ) ) +
            geom_point(data=df,aes(y= Abundance,color = Treatment), position = position_jitter(width= 0.1),     size = 3)+
        geom_line(data = df_sum,  aes(x = Feces, y = mean_abund, group = Treatment,color= Treatment), size = 1 ) +
    geom_ribbon(data = df_sum,  aes(x = Feces, fill = Treatment, group = Treatment, ymin = mean_abund-sd_abund,   
                    ymax =mean_abund+sd_abund),alpha = 0.1) +
            #geom_col(  position = "identity",    color = "black",   linewidth = 0.2) +
    scale_fill_manual( values = c("TAM"="#8B0000BF","EtOH"="#4D4D4DBF"))+
    scale_color_manual( values = c("TAM"="#8B0000BF","EtOH"="#4D4D4DBF"))+
        facet_wrap(  vars(.data[[FILL]]), scales = "free_y") +
        theme_classic() +
    theme(legend.position = "right",
          axis.text.x = element_blank(),
          axis.text = element_text(size = 14, face = "bold"),
      axis.title = element_text(size = 16, face = "bold") ) +
      ylab("Relative abundance")
    
  
  } else {
     barplot<- plot_bar(DATA, fill = FILL)+
      facet_wrap(vars(.data[[WRAP]]), scales="free_x", nrow=1)+
        theme(legend.position = "right",
          axis.line = element_line(color = "black", linewidth = 0.5),
          axis.ticks = element_line(color = "black", linewidth = 0.5),
          axis.title = element_text(size = 16, face = "bold"),
          axis.text = element_text(size = 16, face = "bold"),
          plot.title = element_blank(),
           axis.text.x = element_blank(),
          legend.title = element_text(size = 10),
          legend.text = element_text(size = 10),
          panel.grid = element_blank())

 
  }
  ggsave(plot=barplot,
         filename= paste0("02_GeneratedData/FK49_Microbiome/Abundance_",
                          savename,
                          "_.png"),
         width = WIDTH, height = HEIGHT, dpi = 300,limitsize = FALSE)
  return(barplot)
}
# you would not publieh absolute count values but i want to look at them do tunderstand my data and data quality better.
plot_Micro_Bars(dataset= "ps", RANK= "Phylum", FILL= "Phylum", WRAP= "Feces", savename = "1ALL_P") #absolute
plot_Micro_Bars(dataset= "ps_rel", RANK= "Phylum", FILL= "Phylum", WRAP= "Feces", savename = "1ALL_Prel") #relative

plot_Micro_Bars(dataset= "ps_ND", RANK= "Phylum", FILL= "Phylum", WRAP= "T_T", savename = "2ND_P")
plot_Micro_Bars(dataset= "psrel_ND", RANK= "Phylum", FILL= "Phylum", WRAP= "T_T", savename = "2ND_Prel")

plot_Micro_Bars(dataset= "psrel_ND", RANK= "Order", FILL= "Order", WRAP= "Treatment", savename = "3ND_Orel")
plot_Micro_Bars(dataset= "ps_ND", RANK= "Order", FILL= "Order", WRAP= "Treatment", savename = "3ND_O")

plot_Micro_Bars(dataset= "psrel_CDHFD", RANK= "Phylum", FILL= "Phylum", WRAP= "T_T", savename = "4CDHFD_Prel")
plot_Micro_Bars(dataset= "ps_CDHFD", RANK= "Phylum", FILL= "Phylum", WRAP= "T_T", savename = "4CDHFD_P")

plot_Micro_Bars(dataset= "psrel_CDHFD", RANK= "Order", FILL= "Order", WRAP= "Treatment", savename = "5CDHFD_Orel")
plot_Micro_Bars(dataset= "ps_CDHFD", RANK= "Order", FILL= "Order", WRAP= "Treatment", savename = "5CDHFD_O")

plot_Micro_Bars(dataset= "psrel_male", RANK= "Phylum", FILL= "Phylum", WRAP= "Treatment", savename = "6male_Prel")
plot_Micro_Bars(dataset= "ps_male", RANK= "Phylum", FILL= "Phylum", WRAP= "Treatment", savename = "6male_P")

plot_Micro_Bars(dataset= "psrel_male", RANK= "Order", FILL= "Order", WRAP= "T_T", savename = "7male_Orel", WIDTH=8)
plot_Micro_Bars(dataset= "ps_male", RANK= "Order", FILL= "Order", WRAP= "T_T", savename = "7male_O", WIDTH = 8)

plot_Micro_Bars(dataset= "psrel_male", RANK= "Class", FILL= "Class", WRAP= "Treatment", savename = "8male_Crel")
plot_Micro_Bars(dataset= "ps_male", RANK= "Class", FILL= "Class", WRAP= "Treatment", savename = "8male_C")

plot_Micro_Bars(dataset= "psrel_male_CDHFD", RANK= "Phylum", FILL= "Phylum", WRAP= "Treatment", savename = "9CDHFD_male_Prel")
plot_Micro_Bars(dataset= "ps_male_CDHFD", RANK= "Phylum", FILL= "Phylum", WRAP= "Treatment", savename = "9CDHFD_male_P")

plot_Micro_Bars(dataset= "psrel_male_CDHFD", RANK= "Order", FILL= "Order", WRAP= "Treatment", savename = "10CDHFD_male_Orel")
plot_Micro_Bars(dataset= "ps_male_CDHFD", RANK= "Order", FILL= "Order", WRAP= "Treatment", savename = "10CDHFD_male_O")

plot_Micro_Bars(dataset= "psrel_male_CDHFD", RANK= "Class", FILL= "Class", WRAP= "Treatment", savename = "11CDHFD_male_Crel")
plot_Micro_Bars(dataset= "ps_male_CDHFD", RANK= "Class", FILL= "Class", WRAP= "Treatment", savename = "11CDHFD_male_C")

plot_Micro_Bars(dataset= "psrel_male_CDHFD",subsetTAX="Firmicutes", RANK= "Family", FILL= "Family", WRAP= "Treatment", savename = "12CDHFD_male_Firmi_rel")
plot_Micro_Bars(dataset= "ps_male_CDHFD",subsetTAX="Firmicutes" ,RANK= "Family", FILL= "Family", WRAP= "Treatment", savename = "12CDHFD_male_Firmi")


plot_Micro_Bars(dataset= "psrel_male_CDHFD",subsetTAX="Lactobacillales", RANK= "Family", FILL= "Family", WRAP= "Treatment", savename = "13CDHFD_male_Lacto_rel")
plot_Micro_Bars(dataset= "ps_male_CDHFD",subsetTAX="Lactobacillales" ,RANK= "Family", FILL= "Family", WRAP= "Treatment", savename = "13CDHFD_male_Lacto")
plot_Micro_Bars(dataset= "psrel_male_CDHFD",subsetTAX="Lactobacillales", RANK= "Genus", FILL= "Genus", WRAP= "Treatment", savename = "13_2CDHFD_male_Lacto_rel")
plot_Micro_Bars(dataset= "ps_male_CDHFD",subsetTAX="Lactobacillales" ,RANK= "Genus", FILL= "Genus", WRAP= "Treatment", savename = "13_2CDHFD_male_Lacto")


plot_Micro_Bars(dataset= "psrel_male",subsetTAX="Lactobacillales", RANK= "Genus", FILL= "Genus", WRAP= "T_T", savename = "14male_Lacto_rel",WIDTH=12)
plot_Micro_Bars(dataset= "ps_male",subsetTAX="Lactobacillales" ,RANK= "Genus", FILL= "Genus", WRAP= "T_T", savename = "14male_Lacto",WIDTH=12)


plot_Micro_Bars(dataset= "psrel_female_CDHFD",subsetTAX="Lactobacillales", RANK= "Genus", FILL= "Genus", WRAP= "Treatment", savename = "15_2CDHFD_female_Lacto_rel")
plot_Micro_Bars(dataset= "ps_female_CDHFD",subsetTAX="Lactobacillales" ,RANK= "Genus", FILL= "Genus", WRAP= "Treatment", savename = "15_2CDHFD_female_Lacto")


plot_Micro_Bars(dataset= "psrel_female",subsetTAX="Lactobacillales", RANK= "Genus", FILL= "Genus", WRAP= "T_T", savename = "16female_Lacto_rel",WIDTH=12)
plot_Micro_Bars(dataset= "ps_female",subsetTAX="Lactobacillales" ,RANK= "Genus", FILL= "Genus", WRAP= "T_T", savename = "16female_Lacto",WIDTH=12)

plot_Micro_Bars(dataset= "psrel_male",GROUP_MEAN = FALSE,subsetTAX="Lactobacillus", RANK= "Genus", FILL= "Genus", WRAP= "T_T", savename = "17male_Lactobacillus_rel",WIDTH=12)

plot_Micro_Bars(dataset= "psrel_male",GROUP_MEAN = TRUE, subsetTAX="Lactobacillus", RANK= "Genus", FILL= "Genus", WRAP= "T_T", savename = "18male_Lactobacillus_MEANrel",WIDTH=9)
plot_Micro_Bars(dataset= "ps_male",GROUP_MEAN = TRUE, subsetTAX="Lactobacillus", RANK= "Genus", FILL= "Genus", WRAP= "T_T", savename = "19male_Lactobacillus_MEAN",WIDTH=9)

plot_Micro_Bars(dataset= "psrel_female",GROUP_MEAN = TRUE, subsetTAX="Lactobacillus", RANK= "Genus", FILL= "Genus", WRAP= "T_T", savename = "20female_Lactobacillus_MEANrel",WIDTH=9)
plot_Micro_Bars(dataset= "ps_female",GROUP_MEAN = TRUE, subsetTAX="Lactobacillus", RANK= "Genus", FILL= "Genus", WRAP= "T_T", savename = "21female_Lactobacillus_MEAN",WIDTH=9)

plot_Micro_Bars(dataset= "ps_rel",GROUP_MEAN = TRUE, subsetTAX="Lactobacillus", RANK= "Genus", FILL= "Genus", WRAP= "T_T", savename = "22all_Lactobacillus_MEANrel",WIDTH=9)
plot_Micro_Bars(dataset= "ps",GROUP_MEAN = TRUE, subsetTAX="Lactobacillus", RANK= "Genus", FILL= "Genus", WRAP= "T_T", savename = "23all_Lactobacillus_MEAN",WIDTH=9)

plot_Micro_Bars(dataset= "psrel_male",GROUP_MEAN = TRUE, subsetTAX="Firmicutes", RANK= "Class", FILL= "Class", WRAP= "T_T", savename = "24male_Firmi_ClassMEANrel",WIDTH=25,HEIGHT=25)
plot_Micro_Bars(dataset= "psrel_male",GROUP_MEAN = TRUE, subsetTAX="Firmicutes", RANK= "Order", FILL= "Order", WRAP= "T_T", savename = "25male_Firmi_OrderMEANrel",WIDTH=50,HEIGHT=50)
plot_Micro_Bars(dataset= "psrel_male",GROUP_MEAN = TRUE, subsetTAX="Firmicutes", RANK= "Family", FILL= "Family", WRAP= "T_T", savename = "26male_Firmi_FamilyMEANrel",WIDTH=50,HEIGHT=50)
plot_Micro_Bars(dataset= "psrel_male",GROUP_MEAN = TRUE, subsetTAX="Firmicutes", RANK= "Genus", FILL= "Genus", WRAP= "T_T", savename = "27male_Lactobacillus_MEANrel",WIDTH=50,HEIGHT=50)

plot_Micro_Bars(dataset= "psrel_male",GROUP_MEAN = TRUE, subsetTAX="Clostridiales", RANK= "Genus", FILL= "Genus", WRAP= "T_T", savename = "28male_Clostridiales_GenusMEANrel",WIDTH=50,HEIGHT=50)
plot_Micro_Bars(dataset= "psrel_male",GROUP_MEAN = TRUE, subsetTAX="Ruminococcaceae", RANK= "Genus", FILL= "Genus", WRAP= "T_T", savename = "28male_Clostridiales_MEAN",WIDTH=25,HEIGHT=25)
plot_Micro_Bars(dataset= "ps_male",GROUP_MEAN = TRUE, subsetTAX="Firmicutes", RANK= "Class", FILL= "Class", WRAP= "T_T", savename = "29male_Firmi_ClassMEAN",WIDTH=50,HEIGHT=50)

plot_Micro_Bars(dataset= "psrel_male",GROUP_MEAN = TRUE, subsetTAX=NULL, RANK= "Genus", FILL= "Genus", WRAP= "T_T", savename = "30male_male_ALLGenusMEANrel",WIDTH=50,HEIGHT=50)
plot_Micro_Bars(dataset= "ps_male",GROUP_MEAN = TRUE, subsetTAX=NULL, RANK= "Genus", FILL= "Genus", WRAP= "T_T", savename = "31male_male_ALLGenusMEAN",WIDTH=50,HEIGHT=50)



# Barplots are only for descriptive thigns not for statistics. I need to create object that is statistically testabel
# Define taxonomic levels you want
tax_levels <- c("Phylum", "Class", "Order","Family", "Genus")

# Define the comparisons you want
res_list <- list(
  Sex = "Sex_male_vs_female",
  Treatment = "Treatment_TAM_vs_EtOH",
  Feces_F2_vs_F1 = "Feces_F2_vs_F1",
  Feces_F3_vs_F1 = "Feces_F3_vs_F1",
  Feces_F4_vs_F1 = "Feces_F4_vs_F1",
  Feces_F5_vs_F1 = "Feces_F5_vs_F1",
  Sex_Treatment = "Sexmale.TreatmentTAM"
)

# Initialize empty list to store all results
all_tax_res <- list()

# Loop over taxonomic levels
for(tax in tax_levels){
  
  # Agglomerate phyloseq object at this taxonomic level
  ps_tax <- tax_glom(prune_taxa(taxa_sums(ps) > 0, ps), taxrank = tax, NArm = TRUE)
  
  # Convert to DESeq2
  dds <- phyloseq_to_deseq2(ps_tax, ~ Animal + Sex * Treatment * Feces)
  
  # Run DESeq2
  dds <- DESeq(dds, test = "Wald", fitType = "parametric")
  
  # Loop over comparisons
  for(name in names(res_list)){
    res <- results(dds, name = res_list[[name]])
    res_df <- as.data.frame(res)
    
    # Add taxonomic info and level
    res_df$Taxon <- as.character(tax_table(ps_tax)[rownames(res_df), tax])
    res_df$TaxonomicLevel <- tax
    res_df$Comparison <- name
    
    # Store in master list
    all_tax_res[[paste(tax, name, sep = "_")]] <- res_df
  }
}

# Combine everything into one big data frame
combined_tax_res <- do.call(rbind, all_tax_res)
write.csv(combined_tax_res, "02_GeneratedData/FK49_Microbiome/DESeq2_all_tax_levels.csv", row.names = TRUE)



all_tax_res_male <- list()

for(tax in tax_levels){
  
  # Subset samples first
  ps_male_cdhfd <- subset_samples(ps, Sex=="male" & Feces %in% c("F3","F4","F5"))
  
  # Agglomerate
  ps_male_cdhfd <- tax_glom(prune_taxa(taxa_sums(ps_male_cdhfd) > 0, ps_male_cdhfd), taxrank = tax, NArm = TRUE)
  
  # DESeq2
  dds_male_cdhfd <- phyloseq_to_deseq2(ps_male_cdhfd, ~ Animal + Treatment + Feces)
  dds_male_cdhfd <- DESeq(dds_male_cdhfd)
  
  # Define only valid contrasts
  res_list_male <- list(
    Treatment = "Treatment_TAM_vs_EtOH",
    Feces_F4_vs_F3 = "Feces_F4_vs_F3",
    Feces_F5_vs_F3 = "Feces_F5_vs_F3"
  )
  
  for(name in names(res_list_male)){
    res <- results(dds_male_cdhfd, name = res_list_male[[name]])
    res_df <- as.data.frame(res)
    
    # Add taxonomic info
    res_df$Taxon <- as.character(tax_table(ps_male_cdhfd)[rownames(res_df), tax])
    res_df$TaxonomicLevel <- tax
    res_df$Comparison <- name
    
    all_tax_res_male[[paste(tax, name, sep = "_")]] <- res_df
  }
}

combined_tax_res_male <- do.call(rbind, all_tax_res_male)
write.csv(combined_tax_res_male, "02_GeneratedData/FK49_Microbiome/DESeq2_all_tax_males.csv", row.names = TRUE)


# DIY function
plot_abundance <- function(dataset = "ps_rel", # phyloseq data
                           subsetTAX = NULL,
                           title = "",
                           Facet = "Phylum", # taxa rank for facets
                           Category = "Feces", # categorical features for x axis
                           Color = "Phylum",
                           legend = "none",
                           nameadd= NULL, WIDTH = 9, HEIGHT=9) { 
  ## Subsetting dataset in different interesting sets ----
  fullset <-switch(
    dataset,
    "ps_rel"              = ps_rel,
    "psrel_CDHFD"         = ps_rel,
    "psrel_ND"            = ps_rel,
    "psrel_male"          = ps_rel,
    "psrel_male_CDHFD"    = ps_rel,
    "psrel_female"        = ps_rel,
    "psrel_female_CDHFD"  = ps_rel,
    "ps"                 = ps,
    "ps_CDHFD"            = ps,
    "ps_ND"               = ps,
    "ps_male"             = ps,
    "ps_male_CDHFD"       = ps,
    "ps_female"           = ps,
    "ps_female_CDHFD"     = ps,
    stop("Dataset not recognized")
  )
  
  DATA <- switch(
    dataset,
    "ps_rel"              = ps_rel,
    "psrel_CDHFD"         = prune_samples(sample_data(fullset)$Diet_short  %in% c("CDHFD"),fullset),
    "psrel_ND"            = prune_samples(sample_data(fullset)$Diet_short  %in% c("ND")  ,fullset),
    "psrel_male"          = prune_samples(sample_data(fullset)$Sex         %in% c("male"),fullset),
    "psrel_male_CDHFD"    = prune_samples(sample_data(fullset)$Sex         %in% c("male") & 
                                            sample_data(fullset)$Diet_short     %in% c("CDHFD"),fullset),
    "psrel_female"        = prune_samples(sample_data(fullset)$Sex         %in% c("female"),fullset),
    "psrel_female_CDHFD"  = prune_samples(sample_data(fullset)$Sex         %in% c("female") & 
                                            sample_data(fullset)$Diet_short %in% c("CDHFD"),fullset),
    "ps"                  = ps,
    "ps_CDHFD"            = prune_samples(sample_data(fullset)$Diet_short  %in% c("CDHFD"),fullset),
    "ps_ND"               = prune_samples(sample_data(fullset)$Diet_short  %in% c("ND")  ,fullset),
    "ps_male"             = prune_samples(sample_data(fullset)$Sex         %in% c("male"),fullset),
    "ps_male_CDHFD"       = prune_samples(sample_data(fullset)$Sex         %in% c("male") & 
                                            sample_data(fullset)$Diet_short  %in% c("CDHFD"),fullset),
    "ps_female"           = prune_samples(sample_data(fullset)$Sex         %in% c("female"),fullset),
    "ps_female_CDHFD"     = prune_samples(sample_data(fullset)$Sex         %in% c("female") & 
                                            sample_data(fullset)$Diet_short  %in% c("CDHFD"), fullset),
    stop("Dataset not recognized")
  )
  
  if(!is.null(subsetTAX)){
    DATA <- switch (subsetTAX, 
                    "Firmicutes" =  subset_taxa(DATA, Phylum == "Firmicutes"),
                    "Lactobacillales" =  subset_taxa(DATA, Order == "Lactobacillales"),
                    "Lactobacillus" =  subset_taxa(DATA, Genus  == "Lactobacillus"),
                    "Clostridiales" =  subset_taxa(DATA, Order  == "Clostridiales")
                    
    )
  }
  
  mphyseq <- psmelt(DATA)
  mphyseq <- subset(mphyseq, Abundance > 0)
  
  p<-ggplot(data = mphyseq, 
         mapping = aes_string(x = Category,   y = "Abundance",      color = Color, fill = Color)) +
    geom_violin(fill = NA) +
    geom_point(size = 1, alpha = 0.3, 
               position = position_jitter(width = 0.3)) +
    facet_wrap(facets = Facet, ncol = 3) + 
    scale_y_log10() +
    labs(title = title) +
    theme(legend.position = legend)
  
  ggsave(plot=p,
         filename= paste0("02_GeneratedData/FK49_Microbiome/10Abundance_",nameadd,".png"),
         width = WIDTH, height = HEIGHT, dpi = 300,limitsize= FALSE)
  return(p)
}


plot_abundance(dataset="ps_rel",subsetTAX= NULL, "Microbial Abundance per Phylum All Samples",Category = "Feces",nameadd= "Phylum_full",Facet = "Phylum",WIDTH = 9, HEIGHT= 9)
plot_abundance(dataset="ps_rel",subsetTAX= "Firmicutes", "Microbial Abundance per Phylum All Samples",Category = "T_T",nameadd= "Order_firmi",Facet = "Order",WIDTH = 9, HEIGHT= 9)

plot_abundance(dataset="psrel_CDHFD",subsetTAX= NULL , "Microbial Abundance per family CDHFD Samples",Category = "Treatment",Facet="Genus",nameadd= "CDHFD",
               WIDTH=30,HEIGHT=60)


unique_phylums <- unique(tax_table(ps_rel)[, "Phylum"])
unique_phylums <- unique_phylums[!is.na(unique_phylums)]  # remove NAs
for (phyla in unique_phylums){
  ps_firm <-  subset_taxa(ps_rel, Phylum == phyla)
  p1<-plot_abundance(ps_firm,    
                     title = paste0("Microbial Abundance on ",phyla),
                      Category = "Feces",
                      Facet = "Order",
                      Color = "Order")
  ggsave(plot=p1,
         filename= paste0("02_GeneratedData/FK49_Microbiome/11Abundance_in_Phylum_",phyla,"all.png"),
         width = 9, height = 6, dpi = 300)
 }
rm(unique_phylums)

unique_phylums <- unique(tax_table(only_CDHFD)[, "Phylum"])
unique_phylums <- unique_phylums[!is.na(unique_phylums)]  # remove NAs
for (phyla in unique_phylums){
  ps_firm <-  subset_taxa(only_CDHFD, Phylum == phyla)
  p1<-plot_abundance(ps_firm,    
                     title = paste0("Microbial Abundance on ",phyla),
                     Category = "Sex",
                     Facet = "Order",
                     Color = "Order")
  ggsave(plot=p1,
         filename= paste0("02_GeneratedData/FK49_Microbiome/21Abundance_in_P_",phyla,"CDHFD.png"),
         width = 9, height = 6, dpi = 300)
}
rm(unique_phylums)



unique_phylums <- unique(tax_table(only_male)[, "Phylum"])
unique_phylums <- unique_phylums[!is.na(unique_phylums)]  # remove NAs
for (phyla in unique_phylums){
  ps_firm <-  subset_taxa(only_male, Phylum == phyla)
  p1<-plot_abundance(ps_firm,    
                     title = paste0("Microbial Abundance on ",phyla),
                     Category = "Treatment",
                     Facet = "Order",
                     Color = "Order")
  ggsave(plot=p1,
         filename= paste0("02_GeneratedData/FK49_Microbiome/31Abundance_in_P_",phyla,"male.png"),
         width = 9, height = 6, dpi = 300)
}
rm(unique_phylums)
unique_phylums <- unique(tax_table(only_female)[, "Phylum"])
unique_phylums <- unique_phylums[!is.na(unique_phylums)]  # remove NAs
for (phyla in unique_phylums){
  ps_firm <-  subset_taxa(only_female, Phylum == phyla)
  p1<-plot_abundance(ps_firm,    
                     title = paste0("Microbial Abundance on ",phyla),
                     Category = "Treatment",
                     Facet = "Order",
                     Color = "Order")
  ggsave(plot=p1,
         filename= paste0("02_GeneratedData/FK49_Microbiome/41Abundance_in_P_",phyla,"female.png"),
         width = 9, height = 6, dpi = 300)
}

unique_phylums <- unique(tax_table(only_male)[, "Order"])
unique_phylums <- unique_phylums[!is.na(unique_phylums)]  # remove NAs
for (phyla in unique_phylums){
  ps_firm <-  subset_taxa(only_male, Order == phyla)
  p1<-plot_abundance(ps_firm,    
                     title = paste0("Microbial Abundance on ",phyla),
                     Category = "Treatment",
                     Facet = "Family",
                     Color = "Order")
  ggsave(plot=p1,
         filename= paste0("02_GeneratedData/FK49_Microbiome/31Abundance_in_O_",phyla,"male.png"),
         width = 9, height = 6, dpi = 300)
}
rm(unique_phylums)
unique_phylums <- unique(tax_table(only_male)[, "Family"])
unique_phylums <- unique_phylums[!is.na(unique_phylums)]  # remove NAs

  ps_firm <-  subset_taxa(only_male, Family == "Lactobacillaceae")
  p1<-plot_abundance(ps_firm,    
                     title = paste0("Microbial Abundance on Lactobacillales"),
                     Category = "Treatment",
                     Facet = "Genus",
                     Color = "Family")
  ggsave(plot=p1,
         filename= paste0("02_GeneratedData/FK49_Microbiome/31Abundance_in_O_Lacto_male.png"),
         width = 9, height = 6, dpi = 300)

rm(unique_phylums)


#Generate data.frame with OTUs and metadata
ps_wilcox <- data.frame(t(data.frame(phyloseq::otu_table(ps_rel))))
ps_wilcox$Treatment <- phyloseq::sample_data(ps_rel)$Treatment
#Define functions to pass to map
wilcox_model <- function(df){
  wilcox.test(abund ~ Treatment, data = df)
}
wilcox_pval <- function(df){
  wilcox.test(abund ~ Treatment, data = df)$p.value
}
#Create nested data frames by OTU and loop over each using map 
wilcox_results <- ps_wilcox %>%
  gather(key = OTU, value = abund, -Treatment) %>%
  group_by(OTU) %>%
  nest() %>%
  mutate(wilcox_test = map(data, wilcox_model),
         p_value = map(data, wilcox_pval))                       
#Show results
head(wilcox_results)
wilcox_results <- wilcox_results %>%
  dplyr::select(OTU, p_value) %>%
  unnest(cols = c(p_value))

taxa_info <- data.frame(tax_table(ps_rel))
taxa_info <- taxa_info %>% rownames_to_column(var = "OTU")
#Computing FDR corrected p-values
wilcox_results <- wilcox_results %>%
  full_join(taxa_info) %>%
  arrange(p_value) %>%
  mutate(BH_FDR = p.adjust(p_value, "BH")) %>%
  filter(BH_FDR < 0.05) %>%
  dplyr::select(OTU, p_value, BH_FDR, everything())
print.data.frame(wilcox_results)  

sig_otus <- wilcox_results$OTU[wilcox_results$BH_FDR < 0.05]


## 1. OTU abundance matrix
## ---------------------------

otu_mat <- as(otu_table(ps_rel), "matrix")

if (!taxa_are_rows(ps_rel)) {
  otu_mat <- t(otu_mat)
}

otu_sig <- otu_mat[sig_otus, , drop = FALSE]


## ---------------------------
## 2. Column annotation
## ---------------------------

col_annot <- metadata %>%
  dplyr::select(Treatment, Sex, Diet)

col_annot$Treatment <- factor(
  col_annot$Treatment,
  levels = c("EtOH", "TAM")
)

col_annot$Diet <- factor(
  col_annot$Diet,
  levels = c("ND", "ND_TAM", "CDHFD_3", "CDHFD_7", "CDHFD_11")
)

## Ensure identical sample order

col_annot <- col_annot[colnames(otu_sig), , drop = FALSE]


## ---------------------------
## 3. Remove ND samples
## ---------------------------

keep_samples <- rownames(col_annot)[
  !col_annot$Diet %in% c("ND", "ND_TAM")
]

keep_samples <- intersect(
  keep_samples,
  colnames(otu_sig)
)

otu_sig <- otu_sig[, keep_samples, drop = FALSE]

col_annot <- col_annot[
  keep_samples,
  ,
  drop = FALSE
]


## ---------------------------
## 4. Cluster columns within Treatment
## ---------------------------

ordered_cols <- c()

for(trt in levels(col_annot$Treatment)) {
  
  idx <- which(col_annot$Treatment == trt)
  
  if(length(idx) == 0) next
  
  if(length(idx) == 1) {
    
    ordered_cols <- c(
      ordered_cols,
      rownames(col_annot)[idx]
    )
    
  } else {
    
    hc <- hclust(
      dist(
        t(
          otu_sig[, idx, drop = FALSE]
        )
      ),
      method = "complete"
    )
    
    ordered_cols <- c(
      ordered_cols,
      rownames(col_annot)[idx][hc$order]
    )
  }
}

otu_sig <- otu_sig[, ordered_cols, drop = FALSE]

## ---------------------------
## 5. Taxonomy annotation
## ---------------------------

taxa_info <- as.data.frame(tax_table(ps_rel)) %>%
  tibble::rownames_to_column("OTU") %>%
  dplyr::select(OTU, Genus)

taxa_info$Genus[
  is.na(taxa_info$Genus) |
    taxa_info$Genus == ""
] <- "Unclassified"

row_annot <- taxa_info %>%
  dplyr::filter(OTU %in% rownames(otu_sig)) %>%
  tibble::column_to_rownames("OTU")

row_annot <- row_annot[
  rownames(otu_sig),
  ,  
  drop = FALSE
]


## ---------------------------
## 6. Cluster OTUs within Genus
## ---------------------------

ordered_rows <- c()

for(g in unique(row_annot$Genus)) {
  
  idx <- which(row_annot$Genus == g)
  
  if(length(idx) == 1) {
    
    ordered_rows <- c(
      ordered_rows,
      rownames(otu_sig)[idx]
    )
    
  } else {
    
    hc <- hclust(
      dist(
        otu_sig[idx, , drop = FALSE]
      ),
      method = "complete"
    )
    
    ordered_rows <- c(
      ordered_rows,
      rownames(otu_sig)[idx][hc$order]
    )
  }
}

## Apply row order to BOTH objects

otu_sig <- otu_sig[
  ordered_rows,
  ,
  drop = FALSE
]

otu_sig <- otu_sig[
  ordered_rows,
  ,
  drop = FALSE
]
otu_sig <- otu_sig[
  complete.cases(otu_sig),
  ,
  drop = FALSE
]

otu_sig <- otu_sig[
  apply(otu_sig, 1, sd, na.rm = TRUE) > 0,
  ,
  drop = FALSE
]
col_annot <- col_annot[
  ordered_cols,
  ,
  drop = FALSE
]


row_annot <- row_annot[
  ordered_rows,
  ,
  drop = FALSE
]


## ---------------------------
## 7. Row gaps (between genera)
## ---------------------------

genus_lengths <- rle(
  as.character(row_annot$Genus)
)$lengths

gaps_row <- cumsum(genus_lengths)

gaps_row <- gaps_row[-length(gaps_row)]
## ---------------------------
## 8. Column gaps (between treatments)
## ---------------------------

gaps_col <- c(
  sum(col_annot$Treatment == "EtOH")
)


## ---------------------------
## 9. Colors
## ---------------------------

ann_colors <- list(
  
  Treatment = c(
    EtOH = "#4D4D4DBF",
    TAM  = "#8B0000BF"
  ),
  
  Sex = c(
    male   = "darkblue",
    female = "pink"
  ),
  
  Diet = c(
    ND       = "wheat3",
    ND_TAM   = "gold",
    CDHFD_3  = "purple",
    CDHFD_7  = "purple3",
    CDHFD_11 = "purple4"
  )
)


## ---------------------------
## 10. Heatmap
## ---------------------------

heat_sig<-pheatmap(
  otu_sig,
  
  annotation_col = col_annot,
  annotation_row = row_annot,
  annotation_colors = ann_colors,
  
  scale = "row",
  
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  
  #gaps_row = gaps_row,
  gaps_col = gaps_col,
  
  show_colnames = FALSE,
  
  border_color = "black",
  
  color = colorRampPalette(
    c("navy", "white", "firebrick3")
  )(50),
  
  main = "Significant OTUs grouped by Genus and clustered within Treatment"
)
ggsave(plot=heat_sig,
       filename= paste0("02_GeneratedData/FK49_Microbiome/HeatMAp_OTUs_sig_scaled.png"),
       width = 6, height = 9, dpi = 300,bg = "white")



library(phyloseq)
library(dplyr)
library(tibble)
library(pheatmap)

## ---------------------------
## 1. Bile-acid relevant genera
## ---------------------------

bile_acid_genera <- c(
  "Lactobacillus",
  "Lachnoclostridium",
  "Ruminiclostridium",
  "Ruminiclostridium_5",
  "Ruminiclostridium_6",
  "Ruminiclostridium_9",
  "Oscillibacter",
  "Flavonifractor",
  "Intestinimonas",
  "Roseburia",
  "Butyricicoccus",
  "Bilophila",
  "Enterorhabdus",
  "Olsenella",
  "Alistipes",
  "Bacteroides"
)


## ---------------------------
## 2. OTU abundance matrix
## ---------------------------

otu_mat <- as(otu_table(ps_rel), "matrix")

if (!taxa_are_rows(ps_rel)) {
  otu_mat <- t(otu_mat)
}

## ---------------------------
## 3. Taxonomy
## ---------------------------

taxa_info <- as.data.frame(tax_table(ps_rel)) %>%
  rownames_to_column("OTU") %>%
  select(OTU, Genus)

taxa_info$Genus[
  is.na(taxa_info$Genus) |
    taxa_info$Genus == ""
] <- "Unclassified"

## ---------------------------
## 4. Aggregate to genus level
## ---------------------------

otu_df <- as.data.frame(otu_mat)
otu_df$OTU <- rownames(otu_df)

otu_genus <- otu_df %>%
  left_join(taxa_info, by = "OTU") %>%
  filter(Genus %in% bile_acid_genera) %>%
  select(-OTU) %>%
  group_by(Genus) %>%
  summarise(across(everything(), sum))

otu_sig <- as.data.frame(otu_genus)

rownames(otu_sig) <- otu_sig$Genus
otu_sig$Genus <- NULL

otu_sig <- as.matrix(otu_sig)

## ---------------------------
## 5. Column annotation
## ---------------------------

col_annot <- metadata %>%
  select(Treatment,Diet, Sex,)

col_annot$Treatment <- factor(
  col_annot$Treatment,
  levels = c("EtOH", "TAM")
)

col_annot$Diet <- factor(
  col_annot$Diet,
  levels = c(
    "ND",
    "ND_TAM",
    "CDHFD_3",
    "CDHFD_7",
    "CDHFD_11"
  )
)

col_annot <- col_annot[
  colnames(otu_sig),
  ,
  drop = FALSE
]

## ---------------------------
## 6. Remove ND groups
## ---------------------------

keep_samples <- rownames(col_annot)[
  !col_annot$Diet %in% c("ND", "ND_TAM")
]

otu_sig <- otu_sig[
  ,
  keep_samples,
  drop = FALSE
]

col_annot <- col_annot[
  keep_samples,
  ,
  drop = FALSE
]
otu_sig <- otu_sig[
  complete.cases(otu_sig),
  ,
  drop = FALSE
]

otu_sig <- otu_sig[
  apply(otu_sig, 1, sd, na.rm = TRUE) > 0,
  ,
  drop = FALSE
]
## ---------------------------
## 7. Cluster columns within treatment
## ---------------------------

# ordered_cols <- c()
# 
# for(trt in levels(col_annot$Treatment)) {
#   
#   idx <- which(col_annot$Treatment == trt)
#   
#   if(length(idx) == 0) next
#   
#   if(length(idx) == 1) {
#     
#     ordered_cols <- c(
#       ordered_cols,
#       rownames(col_annot)[idx]
#     )
#     
#   } else {
#     
#     hc <- hclust(
#       dist(
#         t(
#           otu_sig[, idx, drop = FALSE]
#         )
#       ),
#       method = "complete"
#     )
#     
#     ordered_cols <- c(
#       ordered_cols,
#       rownames(col_annot)[idx][hc$order]
#     )
#   }
# }
# 
# otu_sig <- otu_sig[
#   ,
#   ordered_cols,
#   drop = FALSE
# ]
# 
# col_annot <- col_annot[
#   ordered_cols,
#   ,
#   drop = FALSE
# ]

## ---------------------------
## 8. Gap between EtOH and TAM
## ---------------------------

gaps_col <- c(
  sum(col_annot$Treatment == "EtOH")
)

## ---------------------------
## 9. Colors
## ---------------------------

ann_colors <- list(
  
  Treatment = c(
    EtOH = "#4D4D4DBF",
    TAM = "#8B0000BF"
  ),
  
  Sex = c(
    male = "darkblue",
    female = "pink"
  ),
  
  Diet = c(
    ND = "wheat3",
    ND_TAM = "gold",
    CDHFD_3 = "purple",
    CDHFD_7 = "purple3",
    CDHFD_11 = "purple4"
  )
)

## ---------------------------
## 10. Heatmap
## ---------------------------

heat_BA<-pheatmap(
  otu_sig,
  
  annotation_col = col_annot,
  annotation_colors = ann_colors,
  
  scale = "row",
  
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  
  gaps_col = gaps_col,
  
  show_colnames = FALSE,
  
  border_color = "black",
  
  color = colorRampPalette(
    c("navy", "white", "firebrick3")
  )(25),

  main = "Bile acid metabolism-associated genera"
)
ggsave(plot=heat_BA,
       filename= paste0("02_GeneratedData/FK49_Microbiome/HeatMAp_BA_scaled.png"),
       width = 6, height = 7, dpi = 300,bg = "white")
## ---------------------------
## 1. Bile-acid relevant genera
## ---------------------------

BA_deconjugators <- c(
  "Lactobacillus",
  "Bacteroides",
  "Olsenella",
  "Enterorhabdus"
)

BA_2_producers <- c(
  "Lachnoclostridium",
  "Ruminiclostridium",
  "Ruminiclostridium_5",
  "Ruminiclostridium_6",
  "Ruminiclostridium_9"
)

BA_modifiers <- c(
  "Oscillibacter",
  "Flavonifractor",
  "Intestinimonas",
  "Roseburia",
  "Butyricicoccus"
)

BA_associated <- c(
  "Bilophila",
  "Alistipes"
)

all_BA_genera <- c(
  BA_deconjugators,
  BA_2_producers,
  BA_modifiers ,
  BA_associated
)

## ---------------------------
## 2. OTU abundance matrix
## ---------------------------

otu_mat <- as(otu_table(ps_rel), "matrix")

if (!taxa_are_rows(ps_rel)) {
  otu_mat <- t(otu_mat)
}

## ---------------------------
## 3. Taxonomy
## ---------------------------

taxa_info <- as.data.frame(tax_table(ps_rel)) %>%
  rownames_to_column("OTU") %>%
  select(OTU, Genus)

taxa_info$Genus[is.na(taxa_info$Genus) |    taxa_info$Genus == ""] <- "Unclassified"

## ---------------------------
## 4. Aggregate to genus level
## ---------------------------

otu_df <- as.data.frame(otu_mat)
otu_df$OTU <- rownames(otu_df)

otu_genus <- otu_df %>%
  left_join(taxa_info, by = "OTU") %>%
  filter(Genus %in% all_BA_genera) %>%
  select(-OTU) %>%
  group_by(Genus) %>%
  summarise(across(everything(), sum))

otu_sig <- as.data.frame(otu_genus)

rownames(otu_sig) <- otu_sig$Genus
otu_sig$Genus <- NULL

otu_sig <- as.matrix(otu_sig)

## ---------------------------
## 5. Row annotations
## ---------------------------

row_annot <- data.frame(
  Category = case_when(
    rownames(otu_sig) %in% BA_deconjugators ~ "Deconjugators",
    rownames(otu_sig) %in% BA_2_producers ~ "Secondary BA producers",
    rownames(otu_sig) %in% BA_modifiers ~ "BA modifiers",
    rownames(otu_sig) %in% BA_associated ~ "BA associated"
  )
)

rownames(row_annot) <- rownames(otu_sig)

row_annot$Category <- factor(
  row_annot$Category,
  levels = c(
    "Deconjugators",
    "Secondary BA producers",
    "BA modifiers",
    "BA associated"
  )
)

## ---------------------------
## 6. Sort rows by category
## ---------------------------

row_order <- order(row_annot$Category)

otu_sig <- otu_sig[row_order, ]
row_annot <- row_annot[row_order, , drop = FALSE]

## ---------------------------
## 7. Column annotation
## ---------------------------

col_annot <- metadata %>%
  select(Treatment, Sex, Diet)

common_samples <- intersect(
  rownames(col_annot),
  colnames(otu_sig)
)

col_annot <- col_annot[
  common_samples,
  ,
  drop = FALSE
]

otu_sig <- otu_sig[
  ,
  common_samples,
  drop = FALSE
]
rm(common_samples)
col_annot$Treatment <- factor(
  col_annot$Treatment,
  levels = c("EtOH", "TAM")
)

# Sort by Treatment
col_annot <- col_annot %>%
  arrange(Treatment,Sex)

# Reorder abundance matrix columns accordingly
otu_sig <- otu_sig[
  ,
  rownames(col_annot),
  drop = FALSE
]

col_annot$Diet <- factor(
  col_annot$Diet,
  levels = c(
    "ND",
    "ND_TAM",
    "CDHFD_3",
    "CDHFD_7",
    "CDHFD_11"
  )
)

col_annot <- col_annot[
  colnames(otu_sig),
  ,
  drop = FALSE
]

## ---------------------------
## 8. Remove unwanted diets
## ---------------------------

keep_samples <- rownames(col_annot)[
  !col_annot$Diet %in% c(
    "ND",
    "ND_TAM"
    
  )
]

otu_sig <- otu_sig[
  ,
  keep_samples,
  drop = FALSE
]

col_annot <- col_annot[
  keep_samples,
  ,
  drop = FALSE
]

## ---------------------------
## 9. Remove invariant rows
## ---------------------------

otu_sig <- otu_sig[
  complete.cases(otu_sig),
  ,
  drop = FALSE
]

otu_sig <- otu_sig[
  apply(otu_sig, 1, sd) > 0,
  ,
  drop = FALSE
]

row_annot <- row_annot[
  rownames(otu_sig),
  ,
  drop = FALSE
]

## ---------------------------
## 10. Gap positions
## ---------------------------

gaps_col <- c(
  sum(col_annot$Treatment == "EtOH")
)

gaps_row <- c(
  sum(row_annot$Category == "Deconjugators"),
  sum(row_annot$Category %in% c(
    "Deconjugators",
    "Secondary BA producers"
  )),
  sum(row_annot$Category %in% c(
    "Deconjugators",
    "Secondary BA producers",
    "BA modifiers"
  ))
)

## ---------------------------
## 11. Colors
## ---------------------------

ann_colors <- list(
  
  Treatment = c(
    EtOH = "#4D4D4DBF",
    TAM = "#8B0000BF"
  ),
  
  Sex = c(
    male = "darkblue",
    female = "pink"
  ),
  
  Diet = c(
    ND = "wheat3",
    ND_TAM = "gold",
    CDHFD_3 = "purple",
    CDHFD_7 = "purple3",
    CDHFD_11 = "purple4"
  ),
  
  Category = c(
    "Deconjugators" = "#1b9e77",
    "Secondary BA producers" = "#d95f02",
    "BA modifiers" = "#7570b3",
    "BA associated" = "#e7298a"
  )
)

## ---------------------------
## 12. Heatmap
## ---------------------------

heat_BA <- pheatmap(
  otu_sig,
  
  annotation_col = col_annot,
  annotation_row = row_annot,
  
  annotation_colors = ann_colors,
  
  scale = "row",
  
  cluster_rows = FALSE,
  cluster_cols = TRUE,
  
  gaps_col = gaps_col,
  gaps_row = gaps_row,
  
  show_colnames = FALSE,
  
  border_color = "black",
  
  color = colorRampPalette(
    c("navy", "white", "firebrick3")
  )(25),
  
  main = "Bile acid metabolism-associated genera"
)

## ---------------------------
## 13. Save
## ---------------------------

ggsave(
  filename = "02_GeneratedData/FK49_Microbiome/HeatMAp_BA_CDHFD11_clusteredscaled.png",
  plot = heat_BA$gtable,
  width = 9,
  height = 6,
  dpi = 300,
  bg = "white"
)



## Separate heatmaps
## ===========================

plot_ba_heatmap <- function(mat, title, file_name){
  
  if(nrow(mat) == 0) return(NULL)
  
  p <- pheatmap(
    mat,
    
    annotation_col = col_annot,
    annotation_colors = ann_colors,
    
    scale = "row",
    
    cluster_rows = FALSE,
    cluster_cols = TRUE,
    
    show_colnames = FALSE,
    
    border_color = "black",
    
    color = colorRampPalette(
      c("navy", "white", "firebrick3")
    )(50),
    
    main = title
  )
  
  ggsave(
    filename = file_name,
    plot = p$gtable,
    width = 7,
    height = 4,
    dpi = 300,
    bg = "white"
  )
  
  return(p)
}

## ---------------------------
## Deconjugators
## ---------------------------

heat_deconjugators <- plot_ba_heatmap(
  otu_sig[
    rownames(otu_sig) %in% BA_deconjugators,
    ,
    drop = FALSE
  ],
  title = "Primary bile acid deconjugators",
  file_name = "02_GeneratedData/FK49_Microbiome/Heatmap_BA_Deconjugators.png"
)

## ---------------------------
## Secondary BA producers
## ---------------------------

heat_secondary <- plot_ba_heatmap(
  otu_sig[
    rownames(otu_sig) %in% BA_2_producers,
    ,
    drop = FALSE
  ],
  title = "Secondary bile acid producers",
  file_name = "02_GeneratedData/FK49_Microbiome/Heatmap_BA_SecondaryProducers.png"
)

## ---------------------------
## BA modifiers
## ---------------------------

heat_modifiers <- plot_ba_heatmap(
  otu_sig[
    rownames(otu_sig) %in% BA_modifiers,
    ,
    drop = FALSE
  ],
  title = "Bile acid modifiers",
  file_name = "02_GeneratedData/FK49_Microbiome/Heatmap_BA_Modifiers.png"
)