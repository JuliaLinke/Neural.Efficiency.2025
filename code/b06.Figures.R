#Script that mostly creates figures and does some analyses requested by the reviewers
#########################################################
### (A) Installing and loading required packages
#########################################################
if (!require("dplyr")) install.packages("dplyr", dependencies = TRUE)
if (!require("ggplot2")) install.packages("ggplot2", dependencies = TRUE)
if (!require("car")) install.packages("car", dependencies = TRUE)
if (!require("effsize")) install.packages("effsize", dependencies = TRUE)
if (!require("ppcor")) install.packages("ppcor", dependencies = TRUE)
if (!require("cocor")) install.packages("cocor", dependencies = TRUE)

library(dplyr)
library(ggplot2)
library(car)
library(effsize)
library(ppcor)
library(cocor)

################################################################################
### (B) Set paths
################################################################################
if(file.exists("/MyWorkingDirectory/derivatives")){
  datadir <-"/MyWorkingDirectory/derivatives/stats"
  figuredir1<-"/MyWorkingDirectory/derivatives/stats/Figures/Figure_3"
  figuredir2<-"/MyWorkingDirectory/derivatives/stats/Figures/Figure_S4"
  figuredir3<-"/MyWorkingDirectory/derivatives/stats/Figures/Figure_S5"
}

if (!dir.exists(figuredir1)) {
  dir.create(figuredir1, recursive = TRUE)
}
if (!dir.exists(figuredir2)) {
  dir.create(figuredir2, recursive = TRUE)
}
if (!dir.exists(figuredir3)) {
  dir.create(figuredir3, recursive = TRUE)
}

################################################################################
### (C) Load the respective files
################################################################################
setwd(datadir)
Mydata<-read.csv("CombinedData.csv")

################################################################################
### (D) Figure drift rate & neural efficiency 
### Only in Supplement
################################################################################
### (D.1) Drift Rate | no GR
figureS5A<-ggplot(data=Mydata,aes(x=NEURAL.EFFICIENCY_BASELINE_PARTIAL.CORRELATIONS,y=as.numeric(DOT.PROBE.BASELINE.DDM.DRIFT.RATE),color=as.character(COHORT)))+
  geom_point(size  = 4)+
  geom_smooth(method = lm, se = F, col = "black",size = 1,alpha = .8)+
  scale_color_manual(name = "Cohort",values = c("-1" = "snow3", "1" = "black"))+
  theme_classic()+
  theme(aspect.ratio=1)+
  #ylim(.,.3)+
  #xlim(.05,.3)+
  theme(axis.title.x = element_text(face="bold", size=18),
        axis.title.y = element_text(face="bold", size=18),
        axis.text.x = element_text(size=16), axis.text.y = element_text(size=16))+
  labs(x="Neural Efficiency | Timepoint 1", y="DDM | Drifte Rate")

setwd(figuredir3)
ggsave("FigureS5A.svg", plot = figureS5A, device = "svg")
ggsave("FigureS5A.png", plot = figureS5A,device = "png")


################################################################################
# (E) Diagnostic Group
################################################################################
DX.GR<-aov (NEURAL.EFFICIENCY_BASELINE_PARTIAL.CORRELATIONS ~ KSADS.MAIN.DIAGNOSIS + COHORT, data=Mydata)
Anova(DX.GR, type="III")
DX.GR.Race<-aov (NEURAL.EFFICIENCY_BASELINE_PARTIAL.CORRELATIONS ~ KSADS.MAIN.DIAGNOSIS + COHORT + RACE.WHITE + RACE.BLACK + ETHNICITY, data=Mydata)
Anova(DX.GR.Race, type="III")
DX.GR.All<-aov (NEURAL.EFFICIENCY_BASELINE_PARTIAL.CORRELATIONS ~ KSADS.MAIN.DIAGNOSIS + COHORT + RACE.WHITE + RACE.BLACK + ETHNICITY + AGE, data=Mydata)
Anova(DX.GR.All, type="III")

figure3A<-ggplot(Mydata,aes(x=KSADS.MAIN.DIAGNOSIS,y=NEURAL.EFFICIENCY_BASELINE_PARTIAL.CORRELATIONS, fill=KSADS.MAIN.DIAGNOSIS)) + 
  geom_violin(width=1)+
  geom_boxplot(width=0.1, color="grey",alpha=0.2)+
  theme_classic()+
  theme(aspect.ratio=1)+
  theme(axis.title.x = element_text(face="bold", size=0),
        axis.title.y = element_text(face="bold", size=18),
        axis.text.x = element_text(size=16), axis.text.y = element_text(size=16))+
  labs(y="Neural Efficiency | Timepoint 1")

setwd(figuredir1)
ggsave("Figure3A.svg", plot = figure3A, width = 8, height = 8, device = "svg")
ggsave("Figure3A.png", plot = figure3A, width = 8, height = 8, device = "png")

HV<- Mydata %>% filter(KSADS.MAIN.DIAGNOSIS == "HV")
HV<-HV$NEURAL.EFFICIENCY_BASELINE_PARTIAL.CORRELATIONS
ANX<- Mydata %>% filter(KSADS.MAIN.DIAGNOSIS == "ANX")
ANX<-ANX$NEURAL.EFFICIENCY_BASELINE_PARTIAL.CORRELATIONS
cohen.d(HV, ANX, paired = FALSE)

################################################################################
# (F) Anxiety as Dimension
################################################################################
# (F.1) Parent-rated anxiety
figure3B<-ggplot(data=Mydata,aes(x=NEURAL.EFFICIENCY_BASELINE_PARTIAL.CORRELATIONS,y=as.numeric(SCARED.BASELINE.PARENT.RATING),color=as.character(COHORT)))+
  geom_point(size  = 4)+
  geom_smooth(method = lm, se = F, col = "black",size = 1,alpha = .8)+
  scale_color_manual(name = "Cohort",values = c("-1" = "snow3", "1" = "black"))+
  theme_classic()+
  theme(aspect.ratio=1)+
  #ylim(-.06,.06)+
  #xlim(.1,.3)+
  theme(axis.title.x = element_text(face="bold", size=18),
        axis.title.y = element_text(face="bold", size=18),
        axis.text.x = element_text(size=16), axis.text.y = element_text(size=16))+
  labs(x="Neural efficiency | Timepoint 1", y="parent-rated anxiety (SCARED)")
figure3B

ggsave("Figure3B.svg", plot = figure3A, width = 8, height = 8, device = "svg")
ggsave("Figure3B.png", plot = figure3A, width = 8, height = 8, device = "png")

###############################################################################
# (G) Impairment 
ANX <- Mydata[(Mydata[,"KSADS.MAIN.DIAGNOSIS"]) == "ANX",]
CGAS.ANX.GR<-ANX[, c("NEURAL.EFFICIENCY_BASELINE_PARTIAL.CORRELATIONS","CGAS_BASELINE_SCORE","COHORT","SCANNER","AGE")]
CGAS.ANX.GR<-na.omit(CGAS.ANX.GR)
CGAS.ANX.GR$CGAS <- (100-CGAS.ANX.GR$CGAS_BASELINE_SCORE)
pcor(CGAS.ANX.GR)

figure3C<-ggplot(data=CGAS.ANX.GR,aes(x=NEURAL.EFFICIENCY_BASELINE_PARTIAL.CORRELATIONS,y=as.numeric(CGAS_BASELINE_SCORE),color=as.character(COHORT)))+
  geom_point(size  = 4)+
  geom_smooth(method = lm, se = F, col = "black",size = 1,alpha = .8)+
  scale_color_manual(name = "Cohort",values = c("-1" = "snow3", "1" = "black"))+
  theme_classic()+
  theme(aspect.ratio=1)+
  theme(axis.title.x = element_text(face="bold", size=18),
        axis.title.y = element_text(face="bold", size=18),
        axis.text.x = element_text(size=16), axis.text.y = element_text(size=16))+
  labs(x="Neural efficiency | Timepoint 1", y="Clinican-rated severity \n(recoded CGAS, higher values = more severe)")

ggsave("Figure3C.svg", plot = figure3C, width = 8, height = 8, device = "svg")
ggsave("Figure3C.png", plot = figure3C, width = 8, height = 8, device = "png")

###############################################################################
# (H) Race
t.test(NEURAL.EFFICIENCY_BASELINE_PARTIAL.CORRELATIONS ~ RACE.WHITE, data = Mydata)
t.test(NEURAL.EFFICIENCY_BASELINE_PARTIAL.CORRELATIONS ~ RACE.BLACK, data = Mydata)
t.test(NEURAL.EFFICIENCY_BASELINE_PARTIAL.CORRELATIONS ~ RACE.ASIAN, data = Mydata)
t.test(NEURAL.EFFICIENCY_BASELINE_PARTIAL.CORRELATIONS ~ RACE.MULTIPLE, data = Mydata)
t.test(NEURAL.EFFICIENCY_BASELINE_PARTIAL.CORRELATIONS ~ ETHNICITY, data = Mydata)
t.test(SCARED.BASELINE.PARENT.RATING ~ RACE.WHITE, data = Mydata) #white participants more anxious
t.test(SCARED.BASELINE.PARENT.RATING ~ RACE.BLACK, data = Mydata) #black participants less anxious
t.test(SCARED.BASELINE.PARENT.RATING ~ RACE.ASIAN, data = Mydata)
t.test(SCARED.BASELINE.PARENT.RATING ~ RACE.MULTIPLE, data = Mydata)
t.test(SCARED.BASELINE.PARENT.RATING ~ ETHNICITY, data = Mydata)

Mydata %>% group_by(RACE.WHITE) %>% summarise(correlation = cor(NEURAL.EFFICIENCY_BASELINE_PARTIAL.CORRELATIONS, as.numeric(SCARED.BASELINE.PARENT.RATING), use = "complete.obs"))
table(Mydata$RACE.WHITE)
cocor.indep.groups(-0.216,-0.207,78,128)
figureS4A <- ggplot(data=Mydata, aes(x=NEURAL.EFFICIENCY_BASELINE_PARTIAL.CORRELATIONS, y=as.numeric(SCARED.BASELINE.PARENT.RATING), color=as.character(RACE.WHITE))) +
  geom_point(size = 4) +
  geom_smooth(method = lm, se = FALSE, aes(color = as.character(RACE.WHITE)), size = 1, alpha = 0.8) + # Different regression lines
  scale_color_manual(name = "Cohort",values = c("-1" = "#1f77b4", "1" = "#ff7f0e"))+
  theme_classic() +
  theme(aspect.ratio = 1) +
  theme(axis.title.x = element_text(face="bold", size=18),
        axis.title.y = element_text(face="bold", size=18),
        axis.text.x = element_text(size=16), axis.text.y = element_text(size=16)) +
  labs(x = "Neural efficiency | Timepoint 1", y = "Parent-rated anxiety (SCARED)")
figureS4A

Mydata %>% group_by(RACE.BLACK) %>% summarise(correlation = cor(NEURAL.EFFICIENCY_BASELINE_PARTIAL.CORRELATIONS, as.numeric(SCARED.BASELINE.PARENT.RATING), use = "complete.obs"))
table(Mydata$RACE.BLACK)
cocor.indep.groups(-0.200,-0.538,173,33)
figureS4B <- ggplot(data=Mydata, aes(x=NEURAL.EFFICIENCY_BASELINE_PARTIAL.CORRELATIONS, y=as.numeric(SCARED.BASELINE.PARENT.RATING), color=as.character(RACE.BLACK))) +
  geom_point(size = 4) +
  geom_smooth(method = lm, se = FALSE, aes(color = as.character(RACE.BLACK)), size = 1, alpha = 0.8) + # Different regression lines
  scale_color_manual(name = "Cohort",values = c("-1" = "#17becf", "1" = "#e377c2"))+
  theme_classic() +
  theme(aspect.ratio = 1) +
  theme(axis.title.x = element_text(face="bold", size=18),
        axis.title.y = element_text(face="bold", size=18),
        axis.text.x = element_text(size=16), axis.text.y = element_text(size=16)) +
  labs(x = "Neural efficiency | Timepoint 1", y = "Parent-rated anxiety (SCARED)")
figureS4B

Mydata %>% group_by(ETHNICTY) %>% summarise(correlation = cor(NEURAL.EFFICIENCY_BASELINE_PARTIAL.CORRELATIONS, as.numeric(SCARED.BASELINE.PARENT.RATING), use = "complete.obs"))
table(Mydata$ETHNICTY)
cocor.indep.groups(-0.187,-0.176,173,33)
figureS4C <- ggplot(data=Mydata, aes(x=NEURAL.EFFICIENCY_BASELINE_PARTIAL.CORRELATIONS, y=as.numeric(SCARED.BASELINE.PARENT.RATING), color=as.character(ETHNICTY))) +
  geom_point(size = 4) +
  geom_smooth(method = lm, se = FALSE, aes(color = as.character(ETHNICTY)), size = 1, alpha = 0.8) + # Different regression lines
  scale_color_manual(name = "Cohort",values = c("-1" = "#2ca02c", "1" = "#ffbf00"))+
  theme_classic() +
  theme(aspect.ratio = 1) +
  theme(axis.title.x = element_text(face="bold", size=18),
        axis.title.y = element_text(face="bold", size=18),
        axis.text.x = element_text(size=16), axis.text.y = element_text(size=16)) +
  labs(x = "Neural efficiency | Timepoint 1", y = "Parent-rated anxiety (SCARED)")
figureS4C

pvalues.race<-c(0.9486,0.0441,0.9542)
fdr_p_race <- p.adjust(pvalues.race, method = "fdr")
print(fdr_p_race)

setwd(figuredir2)
ggsave("FigureS4A.svg", plot = figureS4A, device = "svg")
ggsave("FigureS4A.png", plot = figureS4A,device = "png")
ggsave("FigureS4B.svg", plot = figureS4B, device = "svg")
ggsave("FigureS4B.png", plot = figureS4B,device = "png")
ggsave("FigureS4C.svg", plot = figureS4C, device = "svg")
ggsave("FigureS4C.png", plot = figureS4C,device = "png")