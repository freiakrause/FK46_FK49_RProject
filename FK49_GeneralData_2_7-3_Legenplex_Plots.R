rm(list=ls())
gc()
library(tidyverse)
library(pheatmap)
library(ggfortify)
library(factoextra)
library(patchwork)
library(ggnewscale)
library(ggrepel)
source("FK49_Definitions.R")

ExpID <- "FK49"

# Read data ----------------------------------------------------------------
if(ExpID == "FK49"){
  output_pwd <- PATHS$legendplex$FK49_output
  d1 <- readRDS(file.path(dirname(dirname(output_pwd)), "01_RawData/FK49_Legendplex_clean.Rds"))
  param_list <- PARAMETERS$Legendplex$cytokine_list
  stats <- read.csv2(file.path(output_pwd, "FK49_Legendplex_Statistics.csv"))
  Treatment_levels <- c("Ctrl","TAM")
  Treatment_colors_plot <- Treatment_colors[c("Ctrl","TAM")]
} else if(ExpID == "FK46") {
  output_pwd <- PATHS$legendplex$FK46_output
  d1 <- readRDS(file.path(dirname(dirname(output_pwd)), "01_RawData/FK46_Legendplex_clean.Rds"))
  param_list <- PARAMETERS$Legendplex$cytokine_list
  stats <- read.csv2(file.path(output_pwd, "FK46_Legendplex_Statistics.csv"))
  Treatment_levels <- c("Ctrl","TAM","Tumor")
  Treatment_colors_plot <- c(Treatment_colors[c("Ctrl","TAM")], Tumor="grey")
} else stop("Give me an existing Experiment ID.")

str(d1)
str(stats)

# Dotplot ------------------------------------------------------------------
do_legendplex <- function(inputdata,value,batch="2",sex="both",y_title,path_images,
                          normal_range=NULL,lowlimit=NULL,hilimit=NULL,p_value_override=NULL){
  d <- inputdata %>% filter(complete.cases(.data[[value]]))
  censored_col <- paste0(value,"_censored")
  direction_col <- paste0(value,"_direction")
  d <- d %>% mutate(
    value_numeric=as.numeric(.data[[value]]),
    censored=as.character(.data[[censored_col]]),
    direction=.data[[direction_col]],
    cens_logical=censored=="TRUE",
    censor_status_combined=case_when(
      censored=="TRUE"&direction=="<"~"Below LOD",
      censored=="TRUE"&direction==">"~"Above ULOQ",
      TRUE~"Detected"
    ))
  
  y_vals <- d$value_numeric
  y_max <- max(y_vals,na.rm=TRUE); y_min <- min(y_vals,na.rm=TRUE)
  y_pos <- if(is.finite(y_max)&y_max>y_min) y_max+0.15*(y_max-y_min) else y_max*1.05
  
  p1 <- ggplot(d,aes(x=Treatment,y=.data[[value]]))+
    stat_summary(fun=mean,geom="bar",aes(fill=Treatment,color=Treatment),alpha=.5,width=.75)+
    scale_color_manual(values=Treatment_colors_plot,limits=Treatment_levels)+
    guides(color=guide_legend(order=1,nrow=2,byrow=TRUE))+
    ggnewscale::new_scale_color()+
    stat_summary(fun.data=mean_sdl,fun.args=list(mult=1),geom="errorbar",width=.2,color="black")+
    scale_fill_manual(name="Treatment",values=Treatment_colors_plot,limits=Treatment_levels)+
    guides(fill=guide_legend(order=1,nrow=2,byrow=TRUE))+
    ggnewscale::new_scale_fill()+
    geom_point(aes(shape=Sex,color=censor_status_combined),fill="lightgrey",alpha=.8,
               position=position_jitter(width=.15,height=0),size=5.3,stroke=1.8)+
    scale_shape_manual(name="Sex",values=Sex_shape)+
    scale_color_manual(name="Censoring",
                       values=c("Below LOD"="navyblue","Above ULOQ"="sienna4","Detected"="black"))+
    scale_y_continuous(name=y_title,expand=expansion(mult=c(.05,.15)))+
    theme_minimal()+
    theme(legend.position="bottom",legend.box="vertical",legend.box.just="left",
          legend.title=element_text(size=11,face="bold"),legend.text=element_text(size=9),
          axis.line=element_line(color="black",linewidth=.5),axis.ticks=element_line(color="black",linewidth=.5),
          axis.title=element_text(size=20,face="bold"),axis.title.x=element_blank(),
          axis.text=element_text(size=19,face="bold"),panel.grid=element_blank())+
    guides(shape=guide_legend(title="Sex",order=2,nrow=1,byrow=TRUE),
           color=guide_legend(title="Censoring",order=4,nrow=1,byrow=TRUE),
           fill=guide_legend(title="Group",order=3,nrow=1,byrow=TRUE))
  
  if(!is.null(normal_range)&&all(is.finite(normal_range)))
    p1 <- p1+geom_hline(yintercept=normal_range,linetype="dotted")
  
  list(plot_raw=p1,y_pos=y_pos)
}

plots <- lapply(param_list,function(p){
  stat_row <- stats %>% filter(parameter==p$value)
  res <- do_legendplex(d1,p$value,batch="ALL",sex="both",y_title=p$y_title,
                       path_images=output_pwd,normal_range=p$normal_range%||%NULL,
                       lowlimit=p$lowlimit%||%NULL,hilimit=p$hilimit%||%NULL)
  p_final <- res$plot_raw+
    annotate("text",x=1.5,y=res$y_pos,
             label=paste0("adj p = ",format.pval(stat_row$p_adj,digits=3)),
             size=5.5,fontface="italic")
  fname_val <- gsub("[^[:alnum:]_]","_",p$value)
  ggsave(paste0(ExpID,"_",fname_val,"_Treatment.png"),p_final,path=output_pwd,
         width=4,height=11,dpi=300)
  p_final
})

# Correlation heatmap: animals ---------------------------------------------
d_mat <- d1 %>% select(all_of(PARAMETERS$Legendplex$cytokines)) %>% as.matrix()
rownames(d_mat) <- d1$Animal
d_mat <- t(d_mat)
ann <- data.frame(Treatment=factor(d1$Treatment,levels=Treatment_levels),Sex=d1$Sex)
rownames(ann) <- d1$Animal

cor_spearman <- cor(d_mat,method="spearman",use="pairwise.complete.obs")
p_cor <- pheatmap(cor_spearman,cluster_rows=TRUE,cluster_cols=TRUE,fontsize=10,
                  fontsize_main=7,fontsize_row=9,fontsize_col=9,fontsize_number=8,
                  annotation_col=ann,color=colorRampPalette(c("blue","white","red"))(100),
                  breaks=seq(-1,1,length.out=101),
                  annotation_colors=list(Treatment=Treatment_colors_plot,Sex=Sex_colors),
                  main=paste0("Correlation HeatMap ",ExpID),border_color="black",
                  cellwidth=10,cellheight=10,angle_col=90,gaps_row=3,
                  show_rownames=FALSE,show_colnames=FALSE,display_numbers=FALSE)

ggsave(paste0(ExpID,"_legendplex_Correlation_Animals.png"),p_cor,path=output_pwd,
       width=8,height=5,dpi=300,bg="white")

# Correlation heatmap: cytokines -------------------------------------------
cor_spearman <- cor(t(d_mat),method="spearman",use="pairwise.complete.obs")
p_cor <- pheatmap(cor_spearman,cluster_rows=TRUE,cluster_cols=TRUE,fontsize=10,
                  fontsize_main=7,fontsize_row=9,fontsize_col=9,fontsize_number=8,
                  color=colorRampPalette(c("blue","white","red"))(100),
                  breaks=seq(-1,1,length.out=101),border_color="black",
                  cellwidth=10,cellheight=10,angle_col=90,gaps_row=3,
                  show_rownames=TRUE,show_colnames=TRUE,display_numbers=FALSE)

ggsave(paste0(ExpID,"_legendplex_Correlation_Cytokines.png"),p_cor,path=output_pwd,
       width=8,height=5,dpi=300,bg="white")

# Scaled heatmap ------------------------------------------------------------
cens_cols <- paste0(PARAMETERS$Legendplex$cytokines,"_censored")
direction_cols <- paste0(PARAMETERS$Legendplex$cytokines,"_direction")
cens_mat <- t(d1 %>% select(all_of(cens_cols)) %>% as.matrix())
direction_mat <- t(d1 %>% select(all_of(direction_cols)) %>% as.matrix())
rownames(cens_mat) <- rownames(direction_mat) <- PARAMETERS$Legendplex$cytokines
d_scaled <- t(scale(t(d_mat)))

stats <- stats %>%
  mutate(effect_size_type=factor(effect_size_type,
                                 levels=c("standardized model effect","GMR","NA")))

p_order <- stats %>% arrange(is.na(effect_size_type),effect_size_type,p_adj) %>% pull(parameter)

ann$Treatment <- factor(ann$Treatment,levels=Treatment_levels)
ann$Sex <- factor(ann$Sex,levels=c("female","male"))
column_order <- order(ann$Treatment,ann$Sex)

row_annot <- stats %>%
  select(parameter,p_adj) %>%
  mutate(adj.p=cut(p_adj,breaks=c(-Inf,.0001,.001,.01,.05,Inf),
                   labels=c("<0.0001","0.0001–0.001","0.001–0.01","0.01–0.05",">0.05"),
                   include.lowest=TRUE)) %>%
  select(parameter,adj.p) %>% column_to_rownames("parameter")

row_annot$adj.p <- as.character(row_annot$adj.p)
row_annot$adj.p[is.na(row_annot$adj.p)] <- "NA"
row_annot$adj.p <- factor(row_annot$adj.p,
                          levels=c("<0.0001","0.0001–0.001","0.001–0.01","0.01–0.05",">0.05","NA"))

d_scaled <- d_scaled[p_order,column_order,drop=FALSE]
ann <- ann[column_order,,drop=FALSE]
row_annot <- row_annot[p_order,,drop=FALSE]
cens_mat <- cens_mat[p_order,column_order,drop=FALSE]
direction_mat <- direction_mat[p_order,column_order,drop=FALSE]

below_detection <- cens_mat==TRUE&direction_mat=="<"
above_detection <- cens_mat==TRUE&direction_mat==">"
d_scaled[below_detection] <- -4
d_scaled[above_detection] <- 4

heatmap_colors <- c("#FFF5CC",colorRampPalette(c("#FFE699","orange","red"))(98),"darkred")
heatmap_breaks <- seq(-4,4,length.out=length(heatmap_colors)+1)

p_color_list <- c("<0.0001"="darkgreen","0.0001–0.001"="forestgreen",
                  "0.001–0.01"="green3","0.01–0.05"="green1",
                  ">0.05"="white","NA"="grey80")

d_s <- pheatmap(d_scaled,cluster_rows=FALSE,cluster_cols=FALSE,fontsize=10,
                fontsize_main=7,fontsize_row=9,fontsize_col=9,fontsize_number=8,
                annotation_col=ann,annotation_row=row_annot,color=heatmap_colors,
                breaks=heatmap_breaks,
                annotation_colors=list(Treatment=Treatment_colors_plot,
                                       Sex=Sex_colors,adj.p=p_color_list),
                main=paste0("Scaled HeatMap ",ExpID),border_color="black",
                cellwidth=10,cellheight=10,angle_col=90,gaps_row=0,
                display_numbers=FALSE,number_format="%.3f",show_colnames=FALSE,
                legend_breaks=c(-4,-2,0,2,4),
                legend_labels=c("Below LOD","-2","0","2","Above ULOQ"))

ggsave(paste0(ExpID,"_Legendplex_scaled.png"),d_s,path=output_pwd,
       width=8,height=5,dpi=300,bg="white")

# Effect sizes --------------------------------------------------------------
stats_eff_size <- stats %>% filter(!is.na(effect_size)) %>%
  mutate(parameter=factor(parameter,levels=rev(p_order)))
stats_lm <- stats_eff_size %>% filter(effect_size_type=="standardized model effect")
stats_gmr <- stats_eff_size %>% filter(effect_size_type=="GMR")

p_lm <- ggplot(stats_lm,aes(x=effect_size,y=parameter))+
  geom_vline(xintercept=0,linetype="dashed")+
  geom_errorbar(aes(xmin=effect_CI_low,xmax=effect_CI_high),height=0)+
  geom_point(size=2)+scale_y_discrete(limits=rev(p_order))+
  labs(x="Standardized model-based effect TAM - Ctrl",y=NULL)+theme_classic()

p_gmr <- ggplot(stats_gmr,aes(x=effect_size,y=parameter))+
  geom_vline(xintercept=1,linetype="dashed")+
  geom_errorbar(aes(xmin=effect_CI_low,xmax=effect_CI_high),height=0)+
  geom_point(size=2)+scale_x_log10()+scale_y_discrete(limits=rev(p_order))+
  labs(x="Geometric mean ratio (TAM / Ctrl)",y=NULL)+theme_classic()

ggsave(p_lm,file=paste0(ExpID,"_Legendplex_Effect_size_LM.png"),dpi=300,width=5,height=3,path=output_pwd)
ggsave(p_gmr,file=paste0(ExpID,"_Legendplex_Effect_size_Cens.png"),dpi=300,width=5,height=3,path=output_pwd)

p_effect <- p_lm+p_gmr+patchwork::plot_layout(widths=c(1,1))
ggsave(p_effect,file=paste0(ExpID,"_Legendplex_Effect_size_LM+Cens.png"),
       dpi=300,width=6,height=3,path=output_pwd)

# PCA -----------------------------------------------------------------------
pca_data <- d1 %>% select(Animal,Sex,Treatment,all_of(PARAMETERS$Legendplex$cytokines)) %>% drop_na()
numerical_d <- pca_data %>% select(all_of(PARAMETERS$Legendplex$cytokines))
data.pca <- prcomp(scale(numerical_d))

summary(data.pca)
data.pca$rotation[,1:2]

f1 <- fviz_eig(data.pca,addlabels=TRUE)
f2 <- fviz_pca_var(data.pca,col.var="black")
f3 <- fviz_cos2(data.pca,choice="var",axes=1:2)
f4 <- fviz_pca_var(data.pca,col.var="cos2",
                   gradient.cols=c("black","orange","green"),repel=TRUE)
f5 <- fviz_contrib(data.pca,choice="var",axes=1,top=15,sort.val="desc")

f6 <- autoplot(data.pca,data=pca_data,x=1,y=2,size=3,fill="Treatment",
               color="Treatment",shape="Sex",frame=TRUE,frame.type="norm",
               frame.level=.95)+
  theme_bw()+theme_classic()+ggtitle("Principal Component Analysis")+
  scale_color_manual(values=Treatment_colors_plot,limits=Treatment_levels)+
  scale_fill_manual(values=Treatment_colors_plot,limits=Treatment_levels)+
  scale_shape_manual(values=Sex_shape)+theme(text=element_text(size=20))

ggsave(f1,path=output_pwd,file=paste0(ExpID,"_PCA_ScreePlot.png"),dpi=300,width=4,height=2.5)
ggsave(f2,path=output_pwd,file=paste0(ExpID,"_PCA_VariablePlot.png"),dpi=300,width=5,height=4)
ggsave(f3,path=output_pwd,file=paste0(ExpID,"_PCA_Cos2.png"),dpi=300,width=5,height=4)
ggsave(f4,path=output_pwd,file=paste0(ExpID,"_PCA_Variable_Cos2.png"),dpi=300,width=5,height=4)
ggsave(f5,path=output_pwd,file=paste0(ExpID,"_PCA_Variable_Contribution_PC1.png"),dpi=300,width=5,height=4)
ggsave(f6,path=output_pwd,file=paste0(ExpID,"_PCA_Scores.png"),dpi=300,width=9,height=7)

# Volcano plot ---------------------------------------------------------------
volcano_data <- stats %>%
  mutate(
    FC=mean_tam/mean_ctrl,
    log2FC=log2(FC),
    negLog10FDR=-log10(p_adj),
    direction=case_when(
      p_adj<.05&log2FC<0~"TAM lower",
      p_adj<.05&log2FC>0~"TAM higher",
      TRUE~"Not significant"),
    significant=!is.na(p_adj)&p_adj<.05,
    Metabolite=parameter) %>%
  filter(is.finite(log2FC),is.finite(negLog10FDR))

p_volcano <- ggplot(volcano_data,aes(x=log2FC,y=negLog10FDR))+
  geom_point(aes(fill=direction),alpha=.5,size=3,stroke=.5,
             position=position_jitter(width=.08),shape=21,color="black")+
  scale_fill_manual(values=c("TAM lower"="blue","Not significant"="grey60","TAM higher"="firebrick"))+
  geom_vline(xintercept=c(-.5,.5),linetype="dashed",color="grey80")+
  geom_hline(yintercept=-log10(.05),linetype="dashed",color="grey80")+
  labs(title="Volcano plot - Treatment at TP11",
       x=expression(paste("log"[2]," FC (TAM / Ctrl)")),
       y=expression(paste("-log"[10],"(adj.p.value)")))+
  theme_classic()+
  theme(panel.grid=element_line(color="grey90",linewidth=.1))+
  ggrepel::geom_text_repel(data=volcano_data %>% filter(significant),
                           aes(label=Metabolite),size=4,max.overlaps=20)+
  coord_cartesian(xlim=c(-4,4),ylim=c(0,5))

print(p_volcano)

ggsave(p_volcano,filename=paste0(ExpID,"_Cytokine_volcano.png"),
       width=6,height=9,dpi=300,path=output_pwd)