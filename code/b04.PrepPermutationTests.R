#Script that prepares the input for the permutation tests in Palm
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
if (!require("writexl")) {
  install.packages("writexl", dependencies = TRUE)
  library(writexl)
}

#########################################################
### (B) Set paths
#########################################################
if(file.exists("/MyWorkingDirectory/derivatives/Stats")){
  datadir1 <-"/MyWorkingDirectory/derivatives/Stats"
  datadir2 <-"/MyWorkingDirectory/derivatives/Cohort1/Efficiency"
  datadir3 <-"/MyWorkingDirectory/derivatives/Cohort2/Efficiency"
  outdir <-"/MyWorkingDirectory/derivatives/Stats/Palm"
  listdir <-"/MyWorkingDirectory/lists"
  covbatdir <-"/MyWorkingDirectory/derivatives/CovBat"
}
#########################################################
### (C) Load the respective files
#########################################################
setwd(datadir1)
Pheno <- read.csv('CombinedData.csv', header = TRUE)

#########################################################
#Create the subset of participants with available DDM parameters
DDM <- Pheno %>% filter(!is.na(DOT.PROBE.BASELINE.DDM.DRIFT.RATE)) %>% select(ID = PARTICIPANT.ID)
DDM$ID <- as.numeric(gsub("sub-s0", "", DDM$ID))

setwd(listdir)
write.table(DDM$ID, file = "T1.DDMsample.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)

#########################################################
#Compile contribution to neural efficiency
setwd(datadir2)
Contrib11p<-read.csv("Contrib_Efficiency_TAU1_Ses1_out.norm.nonparametric_partial.csv", header = FALSE)
setwd(datadir3)
Contrib21p<-read.csv("Contrib_Efficiency_TAU2_Ses1_out.norm.nonparametric_partial.csv", header = FALSE)

Contrib<-rbind(Contrib11p,Contrib21p)
ID<-as.data.frame(Pheno$ID)
Contrib<-cbind(ID,Contrib)
colnames(Contrib)[1] <- "ID"

#########################################################
#Make Input for all participants
setwd(outdir)

AllEff <- Pheno$NEURAL.EFFICIENCY_BASELINE_PARTIAL.CORRELATIONS
write.table(AllEff, "Input_All_Eff.csv", sep = ",", row.names = FALSE, col.names = FALSE, quote = FALSE)

AllContrib <- Contrib[, -1]
AllContrib <- data.frame(lapply(AllContrib, as.numeric))
setwd(outdir)
write.table(AllContrib, "Input_All_Contrib.csv", sep = ",", row.names = FALSE, col.names = FALSE, quote = FALSE)

#########################################################
#Make Input for DDM analysis in the supplementary materials

DDMEff <- left_join(DDM,Pheno, by="ID")
DDMEff <- DDMEff$NEURAL.EFFICIENCY_BASELINE_PARTIAL.CORRELATIONS
write.table(DDMEff, "Input_DDM_Eff.csv", sep = ",", row.names = FALSE, col.names = FALSE, quote = FALSE)

DDMContrib <- left_join(DDM,Contrib, by="ID")
DDMContrib <- DDMContrib[, -1]
DDMContrib <- data.frame(lapply(DDMContrib, as.numeric))
write.table(DDMContrib, "Input_DDM_Contrib.csv", sep = ",", row.names = FALSE, col.names = FALSE, quote = FALSE)

#########################################################
#Make Input for comparison between rest and dot probe in the supplement
setwd(covbatdir)
Cohort1Rest<-read.csv("TAU_11_rest_out.norm.nonparametric.partial.csv", header = FALSE)
Cohort1Rest<-as.data.frame(t(Cohort1Rest))
Cohort1Task<-read.csv("TAU_11_task_out.norm.nonparametric.partial.csv", header = FALSE)
Cohort1Task<-as.data.frame(t(Cohort1Task))
Cohort2Rest<-read.csv("TAU_21_rest_out.norm.nonparametric.partial.csv", header = FALSE)
Cohort2Rest<-as.data.frame(t(Cohort2Rest))
Cohort2Task<-read.csv("TAU_21_task_out.norm.nonparametric.partial.csv", header = FALSE)
Cohort2Task<-as.data.frame(t(Cohort2Task))

Rest<-rbind(Cohort1Rest,Cohort2Rest)
Task<-rbind(Cohort1Task,Cohort2Task)
Comp<-rbind(Rest,Task)
setwd(outdir)
write.table(Comp, "Input_Comp_Rest_TAU.csv", sep = ",", row.names = FALSE, col.names = FALSE, quote = FALSE)
