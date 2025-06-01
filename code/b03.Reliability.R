#Script that tests the retest reliability of the neural efficiency
#########################################################
### (A) Installing and loading required packages
#########################################################
if (!require("dplyr")) {
  install.packages("dplyr", dependencies = TRUE)
  library(dplyr)
}
if (!require("tidyverse")) {
  install.packages("tidyverse", dependencies = TRUE)
  library(tidyverse)
}
if (!require("lme4")) {
  install.packages("lme4", dependencies = TRUE)
  library(lme4)
}
if (!require("lmerTest")) {
  install.packages("lmerTest", dependencies = TRUE)
  library(lmerTest)
}
if (!require("ggplot2")) {
  install.packages("ggplot2", dependencies = TRUE)
  library(ggplot2)
}
if (!require("writexl")) {
  install.packages("writexl", dependencies = TRUE)
  library(writexl)
}
#########################################################
### (B) Set paths
#########################################################
if(file.exists("/MyWorkingDirectory/derivatives")){
  datadir1 <-"MyWorkingDirectory/derivatives/stats"
  datadir2 <-"/Volumes/NIHDATA/TAU/Data/derivatives/Cohort1/Efficiency"
  datadir3 <-"/Volumes/NIHDATA/TAU/Data/derivatives/Cohort2/Efficiency"
  figuredir <-"MyWorkingDirectory/derivatives/stats/Figures/Figure_2"
}
#########################################################
### (C) Load the respective files
#########################################################
setwd(datadir1)
TAU<-read.csv("Phenotype.csv")

setwd(datadir2)
Eff11f<-read.csv("Efficiency_TAU1_Ses1_out.norm.nonparametric_full.csv")
Eff11f$Results<-as.numeric(Eff11f$Results)
Eff11p<-read.csv("Efficiency_TAU1_Ses1_out.norm.nonparametric_partial.csv")
Eff11p$Results<-as.numeric(Eff11p$Results)
Eff12f<-read.csv("Efficiency_TAU1_Ses2_out.norm.nonparametric_full.csv")
Eff12f$Results<-as.numeric(Eff12f$Results)
Eff12p<-read.csv("Efficiency_TAU1_Ses2_out.norm.nonparametric_partial.csv")
Eff12p$Results<-as.numeric(Eff12p$Results)

setwd(datadir3)
Eff21f<-read.csv("Efficiency_TAU2_Ses1_out.norm.nonparametric_full.csv")
Eff21f$Results<-as.numeric(Eff21f$Results)
Eff21p<-read.csv("Efficiency_TAU2_Ses1_out.norm.nonparametric_partial.csv")
Eff21p$Results<-as.numeric(Eff21p$Results)
Eff22f<-read.csv("Efficiency_TAU2_Ses2_out.norm.nonparametric_full.csv")
Eff22f$Results<-as.numeric(Eff22f$Results)
Eff22p<-read.csv("Efficiency_TAU2_Ses2_out.norm.nonparametric_partial.csv")
Eff22p$Results<-as.numeric(Eff22p$Results)

data11<-full_join(Eff11f, Eff11p, by = "ID")
colnames(data11) <- c("ID", "S1_Eff_full_glob", "S1_Eff_part_glob")
data12<-full_join(Eff12f, Eff12p, by = "ID")
colnames(data12) <- c("ID", "S2_Eff_full_glob", "S2_Eff_part_glob")
data21<-full_join(Eff21f, Eff21p, by = "ID")
colnames(data21) <- c("ID", "S1_Eff_full_glob", "S1_Eff_part_glob")
data22<-full_join(Eff22f, Eff22p, by = "ID")
colnames(data22) <- c("ID", "S2_Eff_full_glob", "S2_Eff_part_glob")

data1<-left_join(data11,data12, by = "ID")
data2<-left_join(data21,data22, by = "ID")

newdata<-rbind(data1,data2)
newdata$ID <- gsub("sub-s", "", newdata$ID)

TAU$ID <- as.character(TAU$ID)
TAU<-left_join(newdata,TAU,by="ID")
rm(newdata, data1, data11, data12, data2, data22, data21, Eff11f, Eff11p, Eff12f, Eff12p, Eff21f, Eff21p, Eff22f, Eff22p)

TAU1 <- TAU %>% filter(Cohort == 1)
TAU2 <- TAU %>% filter(Cohort == -1)

Results <- data.frame(
  Model = character(),  # Empty character column for model names
  ICC = numeric(),      # Empty numeric column for ICC values
  CI_low = numeric(),   # Empty numeric column for lower confidence interval
  CI_high = numeric()   # Empty numeric column for upper confidence interval
)

#####################################################################################
# Calculate Reliability of Neural Efficiency
#####################################################################################
# (1) Total Sample
#####################################################################################
# (1.1) For Efficiency derived from the full correlations
#create dataframe that has the correct format and contains only the relevant variables
TAUshortfull<-TAU[,c("ID","S1_Eff_full_glob","S2_Eff_full_glob","Days.Btw.Scans","DX","Cohort","Scanner")]
TAUshortfull <- TAUshortfull[complete.cases(TAUshortfull[,3]),]
TAUshortfull <- TAUshortfull[(TAUshortfull[,5]) == "HV",]
cor.test(TAUshortfull$S1_Eff_full_glob,TAUshortfull$S2_Eff_full_glob,use = "complete.obs",method=c("spearman"))
TAUshortfulllong<-gather(TAUshortfull,Session,Efficiency,S1_Eff_full_glob,S2_Eff_full_glob)
TAUshortfulllong<-TAUshortfulllong[order(TAUshortfulllong$ID),]
TAUshortfullmodel <- lmer(Efficiency ~ Session + Days.Btw.Scans + (1 | ID) , data = TAUshortfulllong)
icc_model<-performance::icc(TAUshortfullmodel, ci=TRUE, ci_level = 0.95)
Results[2,1]<-"Total Sample Full Correlation, with GR"
Results[2,2]<-round(icc_model[1,1],3)
Results[2,3]<-round(icc_model[2,1],3)
Results[2,4]<-round(icc_model[3,1],3)

#####################################################################################
# (1.2) For Efficiency derived from the partial correlations
#create dataframe that has the correct format and contains only the relevant variables
TAUshortpart<-TAU[,c("ID","S1_Eff_part_glob","S2_Eff_part_glob","Days.Btw.Scans","DX","Cohort","Scanner")]
TAUshortpart <- TAUshortpart[complete.cases(TAUshortpart[,3]),]
TAUshortpart <- TAUshortpart[(TAUshortpart[,5]) == "HV",]
cor.test(TAUshortpart$S1_Eff_part_glob,TAUshortpart$S2_Eff_part_glob,use = "complete.obs",method=c("spearman"))
TAUshortpartlong<-gather(TAUshortpart,Session,Efficiency,S1_Eff_part_glob,S2_Eff_part_glob)
TAUshortpartlong<-TAUshortpartlong[order(TAUshortpartlong$ID),]
TAUshortpartmodel <- lmer(Efficiency ~ Session + Days.Btw.Scans + (1 | ID) , data = TAUshortpartlong)
icc_model <- performance::icc(TAUshortpartmodel, ci=TRUE, ci_level = 0.95)
Results[6,1]<-"Total Sample Partial Correlation, with GR"
Results[6,2]<-round(icc_model[1,1],3)
Results[6,3]<-round(icc_model[2,1],3)
Results[6,4]<-round(icc_model[3,1],3)

#Make Figure 2A
figure2A <- ggplot(data=TAUshortpart,aes(x=S1_Eff_part_glob,y=S2_Eff_part_glob,color=as.character(Cohort)))+
  geom_point(size  = 6)+
  geom_smooth(method = lm, se = F, col = "black",size = 1,alpha = .8)+
  scale_color_manual(name = "Cohort",values = c("-1" = "snow3", "1" = "black"))+
  theme_classic()+
  theme(aspect.ratio=1)+
  ylim(.1,.25)+
  xlim(.1,.25)+
  theme(plot.title = element_text(size=24,face="bold"),
        axis.title.x = element_text(face="bold", size=20),
        axis.title.y = element_text(face="bold", size=20),
        axis.text.x = element_text(size=16), axis.text.y = element_text(size=16))+
  labs(x="Neural Efficiency | Timepoint 1", y="Neural Efficiency | Timepoint 2", title = "Retest Reliability, global signal regression")

################################################################################
#(1.2.1) for the first cohort only
TAU1shortpart<-TAU1[,c("ID","S1_Eff_part_glob","S2_Eff_part_glob","Days.Btw.Scans","DX","Cohort","Scanner")]
TAU1shortpart <- TAU1shortpart[complete.cases(TAU1shortpart[,3]),]
TAU1shortpart <- TAU1shortpart[(TAU1shortpart[,5]) == "HV",]
cor.test(TAU1shortpart$S1_Eff_part_glob,TAU1shortpart$S2_Eff_part_glob,use = "complete.obs",method=c("spearman"))
TAU1shortpartlong<-gather(TAU1shortpart,Session,Efficiency,S1_Eff_part_glob,S2_Eff_part_glob)
TAU1shortpartlong<-TAU1shortpartlong[order(TAU1shortpartlong$ID),]
TAU1shortpartmodel <- lmer(Efficiency ~ Session + Days.Btw.Scans + (1 | ID) , data = TAU1shortpartlong)
icc_model <-performance::icc(TAU1shortpartmodel, ci=TRUE, ci_level = 0.95)
Results[7,1]<-"Cohort 1 Partial Correlation, with GR"
Results[7,2]<-round(icc_model[1,1],3)
Results[7,3]<-round(icc_model[2,1],3)
Results[7,4]<-round(icc_model[3,1],3)

################################################################################
#(1.2.2) for the second cohort only
TAU2shortpart<-TAU2[,c("ID","S1_Eff_part_glob","S2_Eff_part_glob","Days.Btw.Scans","DX","Cohort","Scanner")]
TAU2shortpart <- TAU2shortpart[complete.cases(TAU2shortpart[,3]),]
TAU2shortpart <- TAU2shortpart[(TAU2shortpart[,5]) == "HV",]
cor.test(TAU2shortpart$S1_Eff_part_glob,TAU2shortpart$S2_Eff_part_glob,use = "complete.obs",method=c("spearman"))
TAU2shortpartlong<-gather(TAU2shortpart,Session,Efficiency,S1_Eff_part_glob,S2_Eff_part_glob)
TAU2shortpartlong<-TAU2shortpartlong[order(TAU2shortpartlong$ID),]
TAU2shortpartmodel <- lmer(Efficiency ~ Session + Days.Btw.Scans + (1 | ID) , data = TAU2shortpartlong)
icc_model<-performance::icc(TAU2shortpartmodel, ci=TRUE, ci_level = 0.95)
Results[8,1]<-"Cohort 2 Partial Correlation, with GR"
Results[8,2]<-round(icc_model[1,1],3)
Results[8,3]<-round(icc_model[2,1],3)
Results[8,4]<-round(icc_model[3,1],3)

######################################
# Make a plot
MyPlot <- Results[c(2, 6, 7, 8), ]
MyPlot[1,1]<-"Full"
MyPlot[2,1]<-"Partial"
MyPlot[2,3]<-0.424
MyPlot[3,1]<-"Partial Cohort 1"
MyPlot[3,3]<-0.256
MyPlot[4,1]<-"Partial Cohort 2"
MyPlot[4,3]<-0.359
figure2B<-ggplot(MyPlot, aes(x = Model, y = ICC)) +
  geom_point(size = 4, color = "limegreen") +  # Plot mean points
  geom_errorbar(aes(ymin = CI_low, ymax = CI_high), width = 0.2) +  # Error bars
  theme_classic() +
  labs(x = "Type of Correlation + Cohort", y = "ICC") +
  theme(axis.title = element_text(size = 22, face = "bold"),
        text = element_text(size = 18))+
  theme(axis.text.x = element_text(angle = 35, hjust = 1))+ 
  coord_fixed(ratio = 3) 

################################################################################
# Write Results
setwd(figuredir)
ggsave("Figure2B.svg", plot = figure2B, device = "svg")
ggsave("Figure2B.png", plot = figure2B, device = "png")
ggsave("Figure2A.svg", plot = figure2A, device = "svg")
ggsave("Figure2A.png", plot = figure2A, device = "png")
write_xlsx(Results, "ICC_Overview.xlsx")

setwd(datadir)
write.csv(TAU, "CombinedData.csv", row.names = FALSE, quote = FALSE)