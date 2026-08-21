rm(list = ls())
gc()
library(dplyr)
library(tidyr)
library(ggplot2)
library('corrr')
library(ggcorrplot)
library("FactoMineR")
library(factoextra)
library(ggrepel)
source("FK49_Definitions.R")
ExpId = "FK49"
if (ExpId=="FK49") {
  meta_cols <-PARAMETERS$BA$meta_cols
  output_pwd <- PATHS$BA$output
  BA<-readRDS(file= paste0(dirname(dirname(output_pwd)),"/01_RawData/FK49_BA_raw.rds"))
}  else if (ExpId == "FK46"){
  print("You dont have data for this Experiment")
} else{
  print("You dont have data for this Experiment")
}
#Check for Cholemia -----
#- increased bile acids already at TP1
#### Plotting Sum BA whole Dataset -----
Cholemia_Check_Plot<-ggplot(data=BA, aes(x=Timepoint, y = SUM,label = Animal))+
  geom_point()+
  theme_classic()+ 
  labs(title = "Sum of all measured BA per animal")+
  ylab("Sum of Serum Bile Acids [umol/L]")+ # for the y axis label
  geom_text(data = subset(BA, SUM> 500),check_overlap = TRUE,nudge_x = +0.3,nudge_y = 0)

ggsave(plot = Cholemia_Check_Plot, width= 4, height= 6, dpi=300, filename="BG1_CholemiaCheck.png", path =paste0(output_pwd,"/QC/" ))

p_spaghetti <- ggplot(BA, aes(x = Timepoint, y = SUM, group = Animal,color=Treatment,label = Animal)) +
  geom_line(alpha = 0.4, size= 0.1) +
  geom_point(size = 1.5)+
  scale_color_manual(values = Treatment_colors)+
  labs( y = "Sum of serum bile acids [umol/L]", x = "Timepoint") +
  theme_classic() +
  theme(strip.text = element_text(size = 10, hjust = 0.8))+
  stat_summary(aes(group = Treatment), fun = mean, geom = "line", size = 1.2,alpha = 0.5)+
  geom_text(data = subset(BA, SUM> 500),check_overlap = TRUE,nudge_x = +0.3,nudge_y = 0)

ggsave(plot= p_spaghetti, filename= paste0("BG1_CholemiaCheck_s.png"), width= 18, height = 9, dpi = 300,path =paste0(output_pwd,"/QC/" ))

## Check Cholemia again wo outlier 199 ------
#199 is absoulutely removed. 
BA<-BA%>%filter(Animal!="iAL199")
### All BA -----
After199_Cholemia_Check_Plot<-ggplot(data=BA, aes(x=Timepoint, y = SUM,label = Animal))+
  geom_point()+
  theme_classic()+ 
  labs(title = "Sum of all measured BA per animal")+
  ylab("Sum of Serum Bile Acids [µM]")+ # for the y axis label
  geom_text(data = subset(BA, SUM> 500),check_overlap = TRUE,nudge_x = +0.3,nudge_y = 0)


ggsave(plot = After199_Cholemia_Check_Plot, width= 4, height= 6, dpi=300, filename="BG2_CholemiaCheck.png", path =paste0(output_pwd,"/QC/" ))

p_spaghetti_after199 <- ggplot(BA, aes(x = Timepoint, y = SUM, group = Animal, color=Treatment,label = Animal)) +
  geom_line(alpha = 0.4, size= 0.1) +
  geom_point(size = 1.5, position = position_jitterdodge(jitter.width = 0.00, dodge.width = 0.03))+
  scale_color_manual(values = Treatment_colors)+
  labs( y = "Sum of serum bile acids [umol/L]", x = "Timepoint") +
  theme_classic() +
  theme(strip.text = element_text(size = 10, hjust = 0.8))+
  stat_summary(aes(group = Treatment), fun = mean, geom = "line", size = 1.2,alpha = 0.5)+
  geom_text(data = subset(BA, SUM> 500),check_overlap = TRUE,nudge_x = +0.3,nudge_y = 0)

ggsave(plot= p_spaghetti_after199,   filename= paste0("BG2_CholemiaCheck_s.png"), width= 18, height = 9, dpi = 300, path =paste0(output_pwd,"/QC/" ))

### 1°BA ----
After199_Cholemia_Check_Plotprim<-ggplot(data=BA, aes(x=Timepoint, y = SUMprim,label = Animal))+
  geom_point()+
  theme_classic()+ 
  labs(title = "Sum of all measured BA per animal")+
  ylab("Sum of 1°Serum Bile Acids [µM]")+ # for the y axis label
  geom_text(data = subset(BA, SUM> 500),check_overlap = TRUE,nudge_x = +0.3,nudge_y = 0)

ggsave(plot = After199_Cholemia_Check_Plotprim, width= 4, height= 6, dpi=300, filename="BG2_CholemiaCheck1.png", path =paste0(output_pwd,"/QC/" ))

p_spaghetti_after199prim <- ggplot(BA, aes(x = Timepoint, y = SUMprim, group = Animal, color=Treatment,label = Animal)) +
  geom_line(alpha = 0.4, size= 0.1) +
  geom_point(size = 1.5, position = position_jitterdodge(jitter.width = 0.00, dodge.width = 0.03))+
  scale_color_manual(values = Treatment_colors)+
  labs( y = "Sum of 1°serum bile acids [umol/L]", x = "Timepoint") +
  theme_classic() +
  theme(strip.text = element_text(size = 10, hjust = 0.8))+
  stat_summary(aes(group = Treatment), fun = mean, geom = "line", size = 1.2,alpha = 0.5)+
  geom_text(data = subset(BA, SUM> 500),check_overlap = TRUE,nudge_x = +0.3,nudge_y = 0)

ggsave(plot= p_spaghetti_after199prim, filename= paste0("BG2_CholemiaCheck1_s.png"), width= 18, height = 9, dpi = 300, path =paste0(output_pwd,"/QC/" ))

### 2°BA ----
After199_Cholemia_Check_Plotsec<-ggplot(data=BA, aes(x=Timepoint, y = SUMsec,label = Animal))+
  geom_point()+
  theme_classic()+ 
  labs(title = "Sum of all measured BA per animal")+
  ylab("Sum of 2°Serum Bile Acids [µM]")+ # for the y axis label
  geom_text(data = subset(BA, SUM> 500),check_overlap = TRUE,nudge_x = +0.3,nudge_y = 0)

ggsave(plot = After199_Cholemia_Check_Plotsec, width= 4, height= 6, dpi=300, filename="BG2_CholemiaCheck2.png", path =paste0(output_pwd,"/QC/" ))

p_spaghetti_after199sec <- ggplot(BA, aes(x = Timepoint, y = SUMsec, group = Animal, color=Treatment,label = Animal)) +
  geom_line(alpha = 0.4, size= 0.1) +
  geom_point(size = 1.5, position = position_jitterdodge(jitter.width = 0.00, dodge.width = 0.03))+
  scale_color_manual(values = Treatment_colors)+
  labs( y = "Sum of 2°serum bile acids [umol/L]", x = "Timepoint") +
  theme_classic() +
  theme(strip.text = element_text(size = 10, hjust = 0.8))+
  stat_summary(aes(group = Treatment), fun = mean, geom = "line", size = 1.2,alpha = 0.5)+
  geom_text(data = subset(BA, SUM> 500),check_overlap = TRUE,nudge_x = +0.3,nudge_y = 0)

ggsave(plot= p_spaghetti_after199sec, filename= paste0("BG2_CholemiaCheck2_s.png"), width= 18, height = 9, dpi = 300, path =paste0(output_pwd,"/QC/" ))

## Check again after unsure animal 186 and 187 are filtered out ------
bad_animals <- BA$Animal[BA$SUM > 500 & BA$Timepoint == "-1"]
BA_filtered <- BA[!BA$Animal %in% bad_animals, ] # showed increased BA in early and late timepoint indicating potential spontaneous cholemia
#Cell Mol Gastroenterol Hepatol. 2021 Dec 6;13(3):875–878. doi: 10.1016/j.jcmgh.2021.11.012

### All BA ----
After_Cholemia_Check_Plot<-ggplot(data=BA_filtered, aes(x=Timepoint, y = SUM,label = Animal))+
  geom_point()+
  theme_classic()+ 
  labs(title = "Sum of all measured BA per animal")+
  ylab("Sum of Serum Bile Acids [µM]") # for the y axis label


ggsave(plot = After_Cholemia_Check_Plot, width= 4, height= 6, dpi=300, filename="BG3_CholemiaCheck.png",path =paste0(output_pwd,"/QC/" ))

p_spaghetti_after <- ggplot(BA_filtered, aes(x = Timepoint, y = SUM, group = Animal, color=Treatment)) +
  geom_line(alpha = 0.4, size= 0.1) +
  geom_point(size = 1.5, position = position_jitterdodge(jitter.width = 0.00, dodge.width = 0.03))+
  scale_color_manual(values = Treatment_colors)+
  labs( y = "Sum of serum bile acids [umol/L]", x = "Timepoint") +
  theme_classic() +
  theme(strip.text = element_text(size = 10, hjust = 0.8))+
  stat_summary(aes(group = Treatment), fun = mean, geom = "line", size = 1.2,alpha = 0.5)

ggsave(plot= p_spaghetti_after, filename= paste0("BG3_CholemiaCheck_s.png"), width= 18, height = 9, dpi = 300, path =paste0(output_pwd,"/QC/" ))

### 1°BA ----
After_Cholemia_Check_Plot1<-ggplot(data=BA_filtered, aes(x=Timepoint, y = SUMprim,label = Animal))+
  geom_point()+
  theme_classic()+ 
  labs(title = "Sum of 1° BA per animal")+
  ylab("Sum of 1°Serum Bile Acids [µM]") # for the y axis label

ggsave(plot = After_Cholemia_Check_Plot1, width= 4, height= 6, dpi=300, filename="BG3_CholemiaCheck1.png", path =paste0(output_pwd,"/QC/" ))

p_spaghetti_after1 <- ggplot(BA_filtered, aes(x = Timepoint, y = SUMprim, group = Animal, color=Treatment)) +
  geom_line(alpha = 0.4, size= 0.1) +
  geom_point(size = 1.5, position = position_jitterdodge(jitter.width = 0.00, dodge.width = 0.03))+
  scale_color_manual(values = Treatment_colors)+
  labs( y = "Sum of 1°serum bile acids [umol/L]", x = "Timepoint") +
  theme_classic() +
  theme(strip.text = element_text(size = 10, hjust = 0.8))+
  stat_summary(aes(group = Treatment), fun = mean, geom = "line", size = 1.2,alpha = 0.5)

ggsave(plot= p_spaghetti_after1, filename= paste0("BG3_CholemiaCheck1_s.png"), width= 18, height = 9, dpi = 300, path =paste0(output_pwd,"/QC/" ))
### 2°BA ----

After_Cholemia_Check_Plot1<-ggplot(data=BA_filtered, aes(x=Timepoint, y = SUMsec,label = Animal))+
  geom_point()+
  theme_classic()+ 
  labs(title = "Sum of 2° BA per animal")+
  ylab("Sum of 2°Serum Bile Acids [µM]") # for the y axis label


ggsave(plot = After_Cholemia_Check_Plot1, width= 4, height= 6, dpi=300, filename="BG3_CholemiaCheck2.png", path =paste0(output_pwd,"/QC/" ))

p_spaghetti_after2 <- ggplot(BA_filtered, aes(x = Timepoint, y = SUMsec, group = Animal, color=Treatment)) +
  geom_line(alpha = 0.4, size= 0.1) +
  geom_point(size = 1.5, position = position_jitterdodge(jitter.width = 0.00, dodge.width = 0.03))+
  scale_color_manual(values = Treatment_colors)+
  labs( y = "Sum of 2°serum bile acids [umol/L]", x = "Timepoint") +
  theme_classic() +
  theme(strip.text = element_text(size = 10, hjust = 0.8))+
  stat_summary(aes(group = Treatment), fun = mean, geom = "line", size = 1.2,alpha = 0.5)

ggsave(plot= p_spaghetti_after2, filename= paste0("BG3_CholemiaCheck2_s.png"), width= 18, height = 9, dpi = 300, path =paste0(output_pwd,"/QC/" ))
## Clean Up after Check and filtering
rm(list= c(ls(pattern = "^p_spaghetti|After")),Cholemia_Check_Plot)
gc()

BA_filtered<-BA_filtered%>%select(-SUM,-SUMprim,-SUMsec)
BA<-BA%>%select(-SUM)
BA<-BA%>%select(-SUMprim,-SUMsec)
saveRDS(BA, file= paste0(dirname(dirname(output_pwd)),"/01_RawData/FK49_BA_preprocessed.rds"))
saveRDS(BA_filtered, file=paste0(dirname(dirname(output_pwd)),"/01_RawData/FK49_BAfiltered_preprocessed.rds"))
rm(list = ls())
gc()