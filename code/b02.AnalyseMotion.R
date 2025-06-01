### This script investigates associations between anxiety and motion
#########################################################

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
if (!require("car")) {
  install.packages("car", dependencies = TRUE)
  library(car)
}
if (!require("ppcor")) {
  install.packages("ppcor", dependencies = TRUE)
  library(ppcor)
}
#########################################################
### (B) Set paths
#########################################################
if(file.exists("/MyWorkingDirectory/derivatives/")){
  datadir <- "/MyWorkingDirectory/derivatives/stats"
  figuredir1<-"/MyWorkingDirectory/derivatives/stats/Figures/Figure_S1"
  figuredir2<-"/MyWorkingDirectory/derivatives/stats/Figures/Figure_S2"
  }

################################################################################
### (C) Load the respective files
################################################################################
setwd(datadir)
TAU<-read.csv("Phenotype.csv")
TAU1<-read.csv("Motion_TAU1_Ses1.csv")
TAU2<-read.csv("Motion_TAU2_Ses1.csv")

Motion<-rbind(TAU1,TAU2)
Motion<-left_join(TAU,Motion,by="ID")
Motion$task.framewise_displacement <- rowMeans(Motion[, c("ses1.framewise_displacement", "ses2.framewise_displacement")], na.rm = TRUE)
Motion$task.dvars <- rowMeans(Motion[, c("ses1.dvars", "ses2.dvars")], na.rm = TRUE)
Motion$task.trans_X <- rowMeans(Motion[, c("ses1.trans_X", "ses2.trans_X")], na.rm = TRUE)
Motion$task.trans_Y <- rowMeans(Motion[, c("ses1.trans_Y", "ses2.trans_Y")], na.rm = TRUE)
Motion$task.trans_Z <- rowMeans(Motion[, c("ses1.trans_Z", "ses2.trans_Z")], na.rm = TRUE)
Motion$task.rot_X <- rowMeans(Motion[, c("ses1.rot_X", "ses2.rot_X")], na.rm = TRUE)
Motion$task.rot_Y <- rowMeans(Motion[, c("ses1.rot_Y", "ses2.rot_Y")], na.rm = TRUE)
Motion$task.rot_Z <- rowMeans(Motion[, c("ses1.rot_Z", "ses2.rot_Z")], na.rm = TRUE)

################################################################################
### (D) Relationship with anxiety
################################################################################
rest<-Motion[, c("SCARED.BASELINE.PARENT-RATING","rest.framewise_displacement", "rest.dvars", "rest.trans_X", "rest.trans_Y", "rest.trans_Z", "rest.rot_X", "rest.rot_Y", "rest.rot_Z")]
rest<-na.omit(rest)
cor.test(rest$SCARED.BASELINE.PARENT-RATING,rest$rest.rot_X,method="spearman")
cor.test(rest$SCARED.BASELINE.PARENT-RATING,rest$rest.rot_Y,method="spearman")
cor.test(rest$SCARED.BASELINE.PARENT-RATING,rest$rest.rot_Z,method="spearman")
cor.test(rest$SCARED.BASELINE.PARENT-RATING,rest$rest.trans_X,method="spearman")
cor.test(rest$SCARED.BASELINE.PARENT-RATING,rest$rest.trans_Y,method="spearman")
cor.test(rest$SCARED.BASELINE.PARENT-RATING,rest$rest.trans_Z,method="spearman")

task<-Motion[, c("SCARED.BASELINE.PARENT-RATING","task.framewise_displacement", "task.dvars", "task.trans_X", "task.trans_Y", "task.trans_Z", "task.rot_X", "task.rot_Y", "task.rot_Z")]
task<-na.omit(task)
cor.test(task$SCARED.BASELINE.PARENT-RATING,task$task.rot_X,method="spearman")
cor.test(task$SCARED.BASELINE.PARENT-RATING,task$task.rot_Y,method="spearman")
cor.test(task$SCARED.BASELINE.PARENT-RATING,task$task.rot_Z,method="spearman")
cor.test(task$SCARED.BASELINE.PARENT-RATING,task$task.trans_X,method="spearman")
cor.test(task$SCARED.BASELINE.PARENT-RATING,task$task.trans_Y,method="spearman")
cor.test(task$SCARED.BASELINE.PARENT-RATING,task$task.trans_Z,method="spearman")

#rest rot, rest trans, task rot, task trans
pvalues<-c(0.09711,0.9689,0.7321,0.9302,0.04638,0.9929,0.3717,0.6647,0.265,0.1684,0.1403,0.08291)
fdr_p <- p.adjust(pvalues, method = "fdr")
print(fdr_p)

#############################################################################
#Making Figure for Supplement
#############################################################################
#Resting-state
transx<-Motion[,c("SCARED.BASELINE.PARENT-RATING","rest.trans_X","DX")]
transx$type <- "trans_x"
colnames(transx)<-c("SCARED","Value","DX","Type")
transy<-Motion[,c("SCARED.BASELINE.PARENT-RATING","rest.trans_Y","DX")]
transy$type <- "trans_y"
colnames(transy)<-c("SCARED","Value","DX","Type")
transz<-Motion[,c("SCARED.BASELINE.PARENT-RATING","rest.trans_Z","DX")]
transz$type <- "trans_z"
colnames(transz)<-c("SCARED","Value","DX","Type")
rotx<-Motion[,c("SCARED.BASELINE.PARENT-RATING","rest.rot_X","DX")]
rotx$type <- "rot_x"
colnames(rotx)<-c("SCARED","Value","DX","Type")
roty<-Motion[,c("SCARED.BASELINE.PARENT-RATING","rest.rot_Y","DX")]
roty$type <- "rot_y"
colnames(roty)<-c("SCARED","Value","DX","Type")
rotz<-Motion[,c("SCARED.BASELINE.PARENT-RATING","rest.rot_Z","DX")]
rotz$type <- "rot_z"
colnames(rotz)<-c("SCARED","Value","DX","Type")

plotrest<-rbind(transx,transy,transz,rotx,roty, rotz)
plotrest<-na.omit(plotrest)
FigureS1<-ggplot(plotrest, aes(x = SCARED, y = Value,color=as.character(DX))) +
  geom_point(alpha = 0.6,size  = 2) +
  scale_color_manual(name = "DX",values = c("HV" = "#00bfc4ff", "ANX" = "#f8766dff"))+
  geom_smooth(method = "lm", color = "black", se = FALSE) + # Add regression line
  facet_wrap(~Type, scales = "free") + 
  labs(title = "Associations between motion during rest and anxiety",
       x = "Parent-rated anxiety (SCARED)", y = "Motion") +
  theme_classic()+
  ylim(-1, 1)+
  theme(plot.title = element_text(size=18,face="bold"),
        axis.title.x = element_text(face="bold", size=16),
        axis.title.y = element_text(face="bold", size=16),
        axis.text.x = element_text(size=14), axis.text.y = element_text(size=14))
print(FigureS1)
setwd(figuredir1)
ggsave("FigureS1.svg", plot = FigureS1, device = "svg")
ggsave("FigureS1.png", plot = FigureS1,device = "png")

#############################################################################
#Task
transx<-Motion[,c("SCARED.BASELINE.PARENT-RATING","task.trans_X","DX")]
transx$type <- "trans_x"
colnames(transx)<-c("SCARED","Value","DX","Type")
transy<-Motion[,c("SCARED.BASELINE.PARENT-RATING","task.trans_Y","DX")]
transy$type <- "trans_y"
colnames(transy)<-c("SCARED","Value","DX","Type")
transz<-Motion[,c("SCARED.BASELINE.PARENT-RATING","task.trans_Z","DX")]
transz$type <- "trans_z"
colnames(transz)<-c("SCARED","Value","DX","Type")
rotx<-Motion[,c("SCARED.BASELINE.PARENT-RATING","task.rot_X","DX")]
rotx$type <- "rot_x"
colnames(rotx)<-c("SCARED","Value","DX","Type")
roty<-Motion[,c("SCARED.BASELINE.PARENT-RATING","task.rot_Y","DX")]
roty$type <- "rot_y"
colnames(roty)<-c("SCARED","Value","DX","Type")
rotz<-Motion[,c("SCARED.BASELINE.PARENT-RATING","task.rot_Z","DX")]
rotz$type <- "rot_z"
colnames(rotz)<-c("SCARED","Value","DX","Type")

plottask<-rbind(transx,transy,transz,rotx,roty, rotz)
plottask<-na.omit(plottask)
FigureS2<-ggplot(plottask, aes(x = SCARED, y = Value,color=as.character(DX))) +
  geom_point(alpha = 0.6,size  = 2) +
  scale_color_manual(name = "DX",values = c("HV" = "#00bfc4ff", "ANX" = "#f8766dff"))+
  geom_smooth(method = "lm", color = "black", se = FALSE) + # Add regression line
  facet_wrap(~Type, scales = "free") + 
  labs(title = "Associations between motion during task and anxiety",
       x = "Parent-rated anxiety (SCARED)", y = "Motion") +
  theme_classic()+
  ylim(-1, 1)+
  theme(plot.title = element_text(size=18,face="bold"),
        axis.title.x = element_text(face="bold", size=16),
        axis.title.y = element_text(face="bold", size=16),
        axis.text.x = element_text(size=14), axis.text.y = element_text(size=14))
print(FigureS2)
setwd(figuredir2)
ggsave("FigureS2.svg", plot = FigureS2, device = "svg")
ggsave("FigureS2.png", plot = FigureS2,device = "png")

#############################################################################
### (E) Effect of diagnostic status
#############################################################################
DX.rest.frame<-aov (rest.framewise_displacement ~ DX, data=Motion)
Anova(DX.rest.frame, type="III")
DX.rest.dvars<-aov (rest.dvars ~ DX, data=Motion)
Anova(DX.rest.dvars, type="III")
DX.rest.tx<-aov (rest.trans_X ~ DX, data=Motion)
Anova(DX.rest.tx, type="III")
DX.rest.ty<-aov (rest.trans_Y ~ DX, data=Motion)
Anova(DX.rest.ty, type="III")
DX.rest.tz<-aov (rest.trans_Z ~ DX, data=Motion)
Anova(DX.rest.tz, type="III")
DX.rest.rx<-aov (rest.rot_X ~ DX, data=Motion)
Anova(DX.rest.rx, type="III")
DX.rest.ry<-aov (rest.rot_Y ~ DX, data=Motion)
Anova(DX.rest.ry, type="III")
DX.rest.rz<-aov (rest.rot_Z ~ DX, data=Motion)
Anova(DX.rest.rz, type="III")

DX.task.frame<-aov (task.framewise_displacement ~ DX, data=Motion)
Anova(DX.task.frame, type="III")
DX.task.dvars<-aov (task.dvars ~ DX, data=Motion)
Anova(DX.task.dvars, type="III")
DX.task.tx<-aov (task.trans_X ~ DX, data=Motion)
Anova(DX.task.tx, type="III")
DX.task.ty<-aov (task.trans_Y ~ DX, data=Motion)
Anova(DX.task.ty, type="III")
DX.task.tz<-aov (task.trans_Z ~ DX, data=Motion)
Anova(DX.task.tz, type="III")
DX.task.rx<-aov (task.rot_X ~ DX, data=Motion)
Anova(DX.task.rx, type="III")
DX.task.ry<-aov (task.rot_Y ~ DX, data=Motion)
Anova(DX.task.ry, type="III")
DX.task.rz<-aov (task.rot_Z ~ DX, data=Motion)
Anova(DX.task.rz, type="III")

pvaluesdx<-c(0.7923,0.0444935,0.9022,0.40156,0.85913,0.4977,0.9444,0.3284,0.4566,0.14356,0.2072,0.54053)
fdr_p_dx <- p.adjust(pvaluesdx, method = "fdr")
print(fdr_p_dx)
