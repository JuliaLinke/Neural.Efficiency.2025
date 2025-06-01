#How does neural efficiency relate to EX-CBT?
#########################################################
### (A) Installing and loading required packAGEs
#########################################################
if (!require("dplyr")) install.packages("dplyr", dependencies = TRUE)
if (!require("lme4")) install.packages("lme4", dependencies = TRUE)
if (!require("lmerTest")) install.packages("lmerTest", dependencies = TRUE)
if (!require("ggplot2")) install.packages("ggplot2", dependencies = TRUE)
if (!require("sjPlot")) install.packages("sjPlot", dependencies = TRUE)

library(dplyr)
library(lme4)
library(lmerTest)
library(ggplot2)
library(sjPlot)

#########################################################
### (B) Set paths
#########################################################
if(file.exists("/MyWorkingDirectory/stats")){
  datadir <-"/MyWorkingDirectory/stats"
  figuredir <-"/MyWorkingDirectory/stats/Figures/Figure_4"
}

#########################################################
### (C) Load the respective files
#########################################################
setwd(datadir)
TAU<-read.csv("CombinedData.csv")

#########################################################
### (D) Determine whether neural efficiency changes during EX-CBT
#########################################################

##############################
###Reorganizing the data 
##############################
ChangeTX<-TAU[,c("PARTICIPANT.ID","NEURAL.EFFICIENCY_BASELINE_PARTIAL.CORRELATIONS","NEURAL.EFFICIENCY_T2_PARTIAL.CORRELATIONS","ABMT","COHORT","SCANNER","Days.Btw.Scans","AGE","SCARED.BASELINE.PARENT.RATING","PARS.BASELINE.CLINICIAN.RATING","WASI_FULL_2_IQ","SEX","KSADS.MAIN.DIAGNOSIS")]
T2<-ChangeTX[,c("PARTICIPANT.ID","AGE","NEURAL.EFFICIENCY_T2_PARTIAL.CORRELATIONS","COHORT","SCARED.BASELINE.PARENT.RATING","SCANNER","WASI_FULL_2_IQ","SEX","KSADS.MAIN.DIAGNOSIS","ABMT")]
colnames(T2)[3] <- "NeuralEff"
T2$Timepoint <- '2'
T2<-T2[!is.na(T2$NeuralEff),]

T1<-ChangeTX[,c("PARTICIPANT.ID","AGE","NEURAL.EFFICIENCY_BASELINE_PARTIAL.CORRELATIONS","NEURAL.EFFICIENCY_T2_PARTIAL.CORRELATIONS","COHORT","SCARED.BASELINE.PARENT.RATING","SCANNER","WASI_FULL_2_IQ","SEX","KSADS.MAIN.DIAGNOSIS","ABMT")]
colnames(T1)[3] <- "NeuralEff"
T1$Timepoint <- '1'
T1<-T1[!is.na(T1$NEURAL.EFFICIENCY_T2_PARTIAL.CORRELATIONS),]
T1<-T1[,c("PARTICIPANT.ID","AGE","NeuralEff","COHORT","SCARED.BASELINE.PARENT.RATING","SCANNER","WASI_FULL_2_IQ","SEX","KSADS.MAIN.DIAGNOSIS","ABMT","Timepoint")]
table(T1$KSADS.MAIN.DIAGNOSIS)

Change <- rbind(T1,T2)
Change <- Change[order(Change$PARTICIPANT.ID),] 
Change$Timepoint <- as.numeric(Change$Timepoint)
Change$ABMT <- ifelse(Change$KSADS.MAIN.DIAGNOSIS == "HV", -1, Change$ABMT)
Change<-Change[complete.cases(Change[1:10]),]

##############################
###Linear mixed effects models
##############################
change1 <- lmer(NeuralEff ~ COHORT + SCANNER + ABMT + AGE + WASI_FULL_2_IQ + (1 | PARTICIPANT.ID), data=Change, REML=F)
summary(change1)
rand(change1)

change2 <- lmer(NeuralEff ~ KSADS.MAIN.DIAGNOSIS + COHORT + SCANNER + ABMT + AGE + WASI_FULL_2_IQ + (1 | PARTICIPANT.ID), data=Change, REML=F)
summary(change2)
rand(change2)

change3 <- lmer(NeuralEff ~ Timepoint*KSADS.MAIN.DIAGNOSIS + COHORT + SCANNER + ABMT + AGE + WASI_FULL_2_IQ + (1 | PARTICIPANT.ID), data=Change, REML=F)
summary(change3)
plot_model(change3, type = "int")
rand(change3)
anova(change1,change2,change3)
tab_model(change1,change2,change3)

#Exploratory
#change4 <- lmer(NeuralEff ~ Timepoint*WASI_FULL_2_IQ + KSADS.MAIN.DIAGNOSIS + COHORT + SCANNER + ABMT + AGE + (1 | PARTICIPANT.ID), data=Change, REML=F)
#summary(change4)

#########################################################
### (E) Determine predictors of TX response
#########################################################
PredTX<-subset(TAU, KSADS.MAIN.DIAGNOSIS=="ANX")
PredTX <- PredTX[!is.na(PredTX$Perc.Improved), ]
PredTX<-PredTX[,c("PARTICIPANT.ID","AGE","NEURAL.EFFICIENCY_BASELINE_PARTIAL.CORRELATIONS","NEURAL.EFFICIENCY_T2_PARTIAL.CORRELATIONS","ABMT","COHORT","SCANNER","Days.Btw.Scans","SCARED.BASELINE.PARENT.RATING","PARS.BASELINE.CLINICIAN.RATING","PARS.MID-TREATMENT.CLINICIAN.RATING","PARS.POST-TREATMENT.CLINICIAN.RATING","WASI_FULL_2_IQ","SEX")]

##############################
###Reorganizing the data for a in depth analysis of TX response
##############################
PARS0<-PredTX[,c("PARTICIPANT.ID","AGE","NEURAL.EFFICIENCY_BASELINE_PARTIAL.CORRELATIONS","ABMT","COHORT","SCARED.BASELINE.PARENT.RATING","PARS.BASELINE.CLINICIAN.RATING","PARS.BASELINE.CLINICIAN.RATING","SCANNER","WASI_FULL_2_IQ","SEX")]
colnames(PARS0)[8] <- "PARS"
PARS0$Timepoint <- '1'
PARS0<-PARS0[!is.na(PARS0$PARS),]

PARS4<-PredTX[,c("PARTICIPANT.ID","AGE","NEURAL.EFFICIENCY_BASELINE_PARTIAL.CORRELATIONS","ABMT","COHORT","SCARED.BASELINE.PARENT.RATING","PARS.BASELINE.CLINICIAN.RATING","PARS.MID-TREATMENT.CLINICIAN.RATING","SCANNER","WASI_FULL_2_IQ","SEX")]
colnames(PARS4)[8] <- "PARS"
PARS4$Timepoint <- '2'
PARS4<-PARS4[!is.na(PARS4$PARS),]

PARS8<-PredTX[,c("PARTICIPANT.ID","AGE","NEURAL.EFFICIENCY_BASELINE_PARTIAL.CORRELATIONS","ABMT","COHORT","SCARED.BASELINE.PARENT.RATING","PARS.BASELINE.CLINICIAN.RATING","PARS.POST-TREATMENT.CLINICIAN.RATING","SCANNER","WASI_FULL_2_IQ","SEX")]
colnames(PARS8)[8] <- "PARS"
PARS8$Timepoint <- '3'
PARS8<-PARS8[!is.na(PARS8$PARS),]

PARSlong <- rbind(PARS0,PARS4,PARS8)
PARSlong <- PARSlong[order(PARSlong$PARTICIPANT.ID),] 
PARSlong$Timepoint <- as.numeric(PARSlong$Timepoint)
PARSlong<-PARSlong[complete.cases(PARSlong[1:11]),]
PARSlong$PARS <- as.numeric(PARSlong$PARS)

##############################
###Run the linear mixed effect models
##############################
m0 <- lmer(PARS ~ Timepoint + (1 | PARTICIPANT.ID), data=PARSlong, REML=F)
summary(m0)
rand(m0)

m1 <- lmer(PARS ~ Timepoint + SCANNER + COHORT + ABMT + PARS.BASELINE.CLINICIAN.RATING   + AGE + WASI_FULL_2_IQ + SEX + (1 | PARTICIPANT.ID), data=PARSlong, REML=F)
summary(m1)
rand(m1)

m2 <- lmer(PARS ~ Timepoint + NEURAL.EFFICIENCY_BASELINE_PARTIAL.CORRELATIONS + SCANNER + COHORT + ABMT + PARS.BASELINE.CLINICIAN.RATING   + AGE + WASI_FULL_2_IQ + SEX + (1 | PARTICIPANT.ID), data=PARSlong, REML=F)
summary(m2)
rand(m2)

m3 <- lmer(PARS ~ Timepoint*NEURAL.EFFICIENCY_BASELINE_PARTIAL.CORRELATIONS + SCANNER + COHORT + ABMT + PARS.BASELINE.CLINICIAN.RATING   + AGE + WASI_FULL_2_IQ + SEX + (1 | PARTICIPANT.ID), data=PARSlong, REML=F)
summary(m3)
rand(m3)

figure4B<-plot_model(m2, type = "int", terms = c("Timepoint", "PARS")) +
  theme_classic() +  # Apply a clean theme
  labs(x = "Timepoints", 
       y = "Clinician-rated anxiety",
       color = "PARTICIPANT.ID",
       title = NULL) +  # Remove title
  theme(
    text = element_text(size = 22),   # Adjust text size
    axis.title = element_text(face = "bold"),  # Make axis labels bold
    plot.title = element_blank()  # Remove any automatically added title
  )+
 # scale_y_continuous(breaks = c(5, 10, 15))+
  scale_x_continuous(breaks = c(1, 2, 3))+  # Set x-axis breaks to 1, 2, 3
  coord_fixed(ratio = 0.15)
figure4B

setwd(figuredir)
tab_model(m0,m1,m2,m3,file="TableS5.doc")
table<-data.frame(anova(m0,m1,m2,m3))
write_csv(table, "TableS6.csv")

##############################
###Additionally: Calculate TX response as symptom reduction (>20% symptom reduction = responders, Liebert et al., 2015, J Child Adolesc Psychopharm)
##############################
PredTX$PARS.MID-TREATMENT.CLINICIAN.RATING <- as.numeric(as.character(PredTX$PARS.MID-TREATMENT.CLINICIAN.RATING))
PredTX$PARS.POST-TREATMENT.CLINICIAN.RATING <- as.numeric(as.character(PredTX$PARS.POST-TREATMENT.CLINICIAN.RATING))
PredTX$SymptReduct <- ((PredTX$PARS.BASELINE.CLINICIAN.RATING - PredTX$PARS.POST-TREATMENT.CLINICIAN.RATING)/PredTX$PARS.BASELINE.CLINICIAN.RATING)*-100
PredTX$SymptReduct.Bin[PredTX$SymptReduct >= 20] <- 1
PredTX$SymptReduct.Bin[PredTX$SymptReduct < 20] <- -1

##############################
###Calculate reliable change index (RCI) Jacobs & Truax, 1991, J Consult Clin Psychol
###PARS retest-reliability=.55; SD at baseline = 3.21; Group RUoPPAS, 2002, JAACAP
##############################
se <- 3.21*sqrt(1-.55)
sdiff <- sqrt(2*se^2)
PredTX$RCI <- as.numeric((PredTX$PARS.POST-TREATMENT.CLINICIAN.RATING - PredTX$PARS.BASELINE.CLINICIAN.RATING)/sdiff)
PredTX$RCIBin[PredTX$RCI < -1.96] <- 1
PredTX$RCIBin[PredTX$RCI >= -1.96] <- -1

##############################
###Correlation between neural efficiency and two measures of TX response
##############################
cor.test(PredTX$NEURAL.EFFICIENCY_BASELINE_PARTIAL.CORRELATIONS,PredTX$SymptReduct,use = "complete.obs",method=c("spearman"))
cor.test(PredTX$NEURAL.EFFICIENCY_BASELINE_PARTIAL.CORRELATIONS,PredTX$RCI,use = "complete.obs",method=c("spearman"))

##############################
###Make figures for the paper
##############################
figure4A<-ggplot(data = PARSlong, aes(x = as.factor(Timepoint), y = PARS, group = NEURAL.EFFICIENCY_BASELINE_PARTIAL.CORRELATIONS, color = NEURAL.EFFICIENCY_BASELINE_PARTIAL.CORRELATIONS))+ 
  geom_line()+ 
  #facet_grPARTICIPANT.ID(. ~ COHORT)+ 
  #stat_summary(aes(group = 1), geom = "point", fun.y = mean,shape = 17, size = 3)+
  theme_classic()+
  theme(aspect.ratio=1)+
  #scale_colour_gradientn(colours=rainbow(2))
  scale_colour_gradientn(colours = c("red2", "snow1","deepskyblue4"), limits = c(0.10, 0.25), breaks = seq(0.10,0.25, by=0.05))+
  labs(x = "Timing During CBT", y = "Clinician-Rated Anxiety")+
  theme(panel.background = element_rect(fill = "white"))+
  theme(panel.grPARTICIPANT.ID.major=element_line(colour="white"),panel.grPARTICIPANT.ID.minor=element_line(colour="white"))+
  theme(panel.grPARTICIPANT.ID.major=element_line(colour="white"),panel.grPARTICIPANT.ID.minor=element_line(colour="white"))+
  labs(color = "Neural Efficiency")+
  scale_x_discrete(expand = c(0.05, 0)) +
  scale_y_continuous(limits = c(0,22.5), breaks = seq(0,22.5,5))+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(axis.title = element_text(face = "bold", size = 22))+
  theme(axis.text.x = element_text(size = 18))+ 
  theme(axis.text.y = element_text(size = 18))+ 
  theme(title = element_text(size = 22))+
  theme(legend.title = element_text(size = 18))+
  theme(legend.text = element_text(size = 18))+
  theme(legend.key.size = unit(0.5, "cm"))
figure4A

figure4C<-ggplot(data=PredTX,aes(x=NEURAL.EFFICIENCY_BASELINE_PARTIAL.CORRELATIONS,y=as.numeric(SymptReduct),color=factor(COHORT)))+
  geom_point(size  = 4)+
  geom_smooth(method = lm, se = F, col = "black",size = 1,alpha = .8)+
  theme_classic()+
  theme(aspect.ratio=1)+
  ylim(-100,50)+
  xlim(.1,.25)+
  theme(plot.title = element_text(size=18,face="bold"),
        axis.title.x = element_text(face="bold", size=22),
        axis.title.y = element_text(face="bold", size=22),
        axis.text.x = element_text(size=18, angle = 20), axis.text.y = element_text(size=18))+
  labs(x="Neural efficiency", y="% Symptom Reduction")+
  scale_color_manual(labels = c("COHORT2", "COHORT1"), values = c("snow4", "black"), aesthetics = c("colour","fill"))
figure4C

setwd(figuredir)
ggsave("figure_4A.svg", figure4A, wPARTICIPANT.IDth = 8, height = 8,dpi=300, device = "svg")
ggsave("figure_4A.png", figure4A, wPARTICIPANT.IDth = 8, height = 8,dpi=300, device = "png")
ggsave("figure_4B.svg", figure4B, wPARTICIPANT.IDth = 8, height = 8,dpi=300, device = "svg")
ggsave("figure_4B.png", figure4B, wPARTICIPANT.IDth = 8, height = 8,dpi=300, device = "png")
ggsave("figure_4C.svg", figure4C, wPARTICIPANT.IDth = 8, height = 8,dpi=300, device = "svg")
ggsave("figure_4C.png", figure4C, wPARTICIPANT.IDth = 8, height = 8,dpi=300, device = "png")

################################################################################
### Extra Analyses exploring ABMT effect on neural efficiency
################################################################################
ABMT<- Change[Change$KSADS.MAIN.DIAGNOSIS == 'ANX',]
ABMT<- ABMT[!is.na(ABMT$ABMT),]

NE.ABMT.base <- lmer(NeuralEff ~ COHORT + SCANNER + AGE  + WASI_FULL_2_IQ + (1 | PARTICIPANT.ID), data=ABMT, REML=F)
summary(NE.ABMT.base)
rand(NE.ABMT.base)

NE.ABMT.effect1 <- lmer(NeuralEff ~ Timepoint*ABMT + COHORT + SCANNER  + AGE + WASI_FULL_2_IQ + (1 | PARTICIPANT.ID), data=ABMT, REML=F)
summary(NE.ABMT.effect1)
rand(NE.ABMT.effect1)

#recode ABMT into 4 conditions (although that is essentially confounded with the COHORT variable)
ABMT <- ABMT %>%
  mutate(ABMT4 = case_when(
    COHORT == 1  & ABMT == 1  ~ "COHORT1.Active",
    COHORT == 1 & ABMT == -1 ~ "COHORT1.Sham",
    COHORT == -1  & ABMT == 1 ~ "COHORT2.Active",
    COHORT == -1 & ABMT == -1  ~ "COHORT2.Sham",
    TRUE              ~ NA_character_
  ))
NE.ABMT.effect2 <- lmer(NeuralEff ~ Timepoint*ABMT4 + SCANNER + AGE + WASI_FULL_2_IQ + (1 | PARTICIPANT.ID), data=ABMT, REML=F)
summary(NE.ABMT.effect2)
rand(NE.ABMT.effect2)


################################################################################
### Extra Analyses exploring ABMT effect on PARS
################################################################################
PARS.ABMT.base <- lmer(PARS ~ COHORT + SCANNER + PARS.BASELINE.CLINICIAN.RATING + AGE  + WASI_FULL_2_IQ + SEX + (1 | PARTICIPANT.ID), data=PARSlong, REML=F)
summary(PARS.ABMT.base)
rand(PARS.ABMT.base)

NE.ABMT.effect1 <- lmer(PARS ~ Timepoint*ABMT + COHORT + SCANNER  + AGE + WASI_FULL_2_IQ + PARS.BASELINE.CLINICIAN.RATING + SEX + (1 | PARTICIPANT.ID), data=PARSlong, REML=F)
summary(NE.ABMT.effect1)
rand(NE.ABMT.effect1)

NE.ABMT.Int.effect1 <- lmer(PARS ~ Timepoint*ABMT*NEURAL.EFFICIENCY_BASELINE_PARTIAL.CORRELATIONS + COHORT + SCANNER  + AGE + WASI_FULL_2_IQ + PARS.BASELINE.CLINICIAN.RATING + SEX + (1 | PARTICIPANT.ID), data=PARSlong, REML=F)
summary(NE.ABMT.Int.effect1)
rand(NE.ABMT.Int.effect1)

setwd(figuredir)
tab_model(m0,m1,m3,NE.ABMT.Int.effect1,file="TableS5.doc")
table2<-data.frame(anova(m0,m1,m3,NE.ABMT.Int.effect1))
write_csv(table2, "TableS6.csv")