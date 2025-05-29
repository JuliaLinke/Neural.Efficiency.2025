### Harmonizing connectivity matrices for efficiency analysis

#########################################################
### (A) Installing and loading required packages
#########################################################
if (!require("dplyr")) {
  install.packages("dplyr", dependencies = TRUE)
  library(dplyr)
}

if (!require("devtools")) {
  install.packages("devtools", dependencies = TRUE)
  library(devtools)
}

if (!require("CovBat")) {
  devtools::install_github("andy1764/CovBat_Harmonization/R")
  library(CovBat)
}

#########################################################
### (B) Set paths
#########################################################
if(file.exists("/MyWorkingDirectory/derivatives/")){
datadir <-"/MyWorkingDirectory/derivatives"
datadir1 <-"/MyWorkingDirectory/derivatives/Cohort1/Efficiency"
datadir2 <-"/MyWorkingDirectory/derivatives/Cohort2/Efficiency"
}

#########################################################
# (C) Load the respective files
#########################################################
setwd(datadir)
data<-read.csv("Phenotype.csv")

### Read in connectivity based on full correlations
setwd(datadir1)
data11_full<-read.csv("Input_Covbat_TAU1_Ses1_full.csv",header = FALSE)
data12_full<-read.csv("Input_Covbat_TAU1_Ses2_full.csv",header = FALSE)

setwd(datadir2)
data21_full<-read.csv("Input_Covbat_TAU2_Ses1_full.csv",header = FALSE)
data22_full<-read.csv("Input_Covbat_TAU2_Ses2_full.csv",header = FALSE)

data_full<-cbind(data11_full,data12_full,data21_full,data22_full)
rm(data11_full,data12_full,data21_full,data22_full)

### Read in connectivity based on partial correlations
setwd(datadir1)
data11_partial<-read.csv("Input_Combat_TAU1_Ses1_partial.csv",header = FALSE)
data12_partial<-read.csv("Input_Combat_TAU1_Ses2_partial.csv",header = FALSE)

data21_partial<-read.csv("Input_Combat_TAU2_Ses1_partial.csv",header = FALSE)
data22_partial<-read.csv("Input_Combat_TAU2_Ses2_partial.csv",header = FALSE)

data_partial<-cbind(data11_partial,data12_partial,data21_partial,data22_partial)
rm(data11_partial,data12_partial,data21_partial,data22_partial)

#########################################################
# (D) Create additional vectors
#########################################################
data <- data %>%
  mutate(CovBatVar = case_when(
    Cohort == 1 & Scanner == 1  ~ 1,
    Cohort == 1 & Scanner == -1 ~ 2,
    Cohort == -1 & Scanner == 1 ~ 3,
    Cohort == -1 & Scanner == -1 ~ 4
  ))

data_11 <- data %>% filter(Cohort == 1)
data_11 <- data_11 %>% mutate(DX = recode(DX, "HV" = -1, "ANX" = 1))
data_21 <- data %>% filter(Cohort == -1)
data_21 <- data_21 %>% mutate(DX = recode(DX, "HV" = -1, "ANX" = 1))
data_2 <- data %>% filter(!is.na(Ses2.Sim.Part))
data_12 <- data_2 %>% filter(Cohort == 1)
data_12 <- data_12 %>% mutate(DX = recode(DX, "HV" = -1, "ANX" = 1))
data_22 <- data_2 %>% filter(Cohort == -1)
data_22 <- data_22 %>% mutate(DX = recode(DX, "HV" = -1, "ANX" = 1))

Covb1<-t(data_11$CovBatVar)
Covb2<-t(data_12$CovBatVar) 
Covb3<-t(data_21$CovBatVar)
Covb4<-t(data_22$CovBatVar)
Covbatvector <- c(Covb1, Covb2, Covb3, Covb4)
rm(Covb1, Covb2, Covb3, Covb4, data_11, data_12, data_2, data_21, data_22)

##############################################################
#Now do the actual thing
##############################################################
norm.nonparametric.partial   <- CovBat::covbat(data_partial, Covbatvector, parametric=FALSE)
norm.nonparametric.full      <- CovBat::covbat(data_full, Covbatvector, parametric=FALSE)

#Break them apart again

out.norm.nonparametric.partial <- as.data.frame(norm.nonparametric.partial$dat.covbat, row.names = NULL)
data_11_out.norm.nonparametric.partial <- out.norm.nonparametric.partial %>% select(1:86)
data_11_rest_out.norm.nonparametric.partial <- data_11_out.norm.nonparametric.partial %>% slice(1:(n()/2))       
data_11_task_out.norm.nonparametric.partial <- data_11_out.norm.nonparametric.partial %>% slice(((n()/2) + 1):n()) 
data_12_out.norm.nonparametric.partial <- out.norm.nonparametric.partial %>% select(87:(86+34))
data_12_rest_out.norm.nonparametric.partial <- data_12_out.norm.nonparametric.partial %>% slice(1:(n()/2))       
data_12_task_out.norm.nonparametric.partial <- data_12_out.norm.nonparametric.partial %>% slice(((n()/2) + 1):n()) 
data_21_out.norm.nonparametric.partial <- out.norm.nonparametric.partial %>% select((86+34+1):(86+34+120))
data_21_rest_out.norm.nonparametric.partial <- data_21_out.norm.nonparametric.partial %>% slice(1:(n()/2))       
data_21_task_out.norm.nonparametric.partial <- data_21_out.norm.nonparametric.partial %>% slice(((n()/2) + 1):n())
data_22_out.norm.nonparametric.partial <- out.norm.nonparametric.partial %>% select((86+34+120+1):(86+34+120+56))
data_22_rest_out.norm.nonparametric.partial <- data_22_out.norm.nonparametric.partial %>% slice(1:(n()/2))       
data_22_task_out.norm.nonparametric.partial <- data_22_out.norm.nonparametric.partial %>% slice(((n()/2) + 1):n()) 


out.norm.nonparametric.full <- as.data.frame(norm.nonparametric.full$dat.covbat, row.names = NULL)
data_11_out.norm.nonparametric.full <- out.norm.nonparametric.full %>% select(1:86)
data_11_rest_out.norm.nonparametric.full <- data_11_out.norm.nonparametric.full %>% slice(1:(n()/2))       
data_11_task_out.norm.nonparametric.full <- data_11_out.norm.nonparametric.full %>% slice(((n()/2) + 1):n())  
data_12_out.norm.nonparametric.full <- out.norm.nonparametric.full %>% select(87:(86+34))
data_12_rest_out.norm.nonparametric.full <- data_12_out.norm.nonparametric.full %>% slice(1:(n()/2))       
data_12_task_out.norm.nonparametric.full <- data_12_out.norm.nonparametric.full %>% slice(((n()/2) + 1):n()) 
data_21_out.norm.nonparametric.full <- out.norm.nonparametric.full %>% select((86+34+1):(86+34+120))
data_21_rest_out.norm.nonparametric.full <- data_21_out.norm.nonparametric.full %>% slice(1:(n()/2))       
data_21_task_out.norm.nonparametric.full <- data_21_out.norm.nonparametric.full %>% slice(((n()/2) + 1):n()) 
data_22_out.norm.nonparametric.full <- out.norm.nonparametric.full %>% select((86+34+120+1):(86+34+120+56))
data_22_rest_out.norm.nonparametric.full <- data_22_out.norm.nonparametric.full %>% slice(1:(n()/2))       
data_22_task_out.norm.nonparametric.full <- data_22_out.norm.nonparametric.full %>% slice(((n()/2) + 1):n()) 

#Write output

setwd(datadir1)
write.table(data_11_rest_out.norm.nonparametric.partial, quote=FALSE, sep = ",", file="data_11_rest_out.norm.nonparametric.partial.csv", row.names=FALSE, col.names = FALSE)
write.table(data_12_rest_out.norm.nonparametric.partial, quote=FALSE, sep = ",", file="data_12_rest_out.norm.nonparametric.partial.csv", row.names=FALSE, col.names = FALSE)
write.table(data_11_rest_out.norm.nonparametric.full, quote=FALSE, sep = ",", file="data_11_rest_out.norm.nonparametric.full.csv", row.names=FALSE, col.names = FALSE)
write.table(data_12_rest_out.norm.nonparametric.full, quote=FALSE, sep = ",", file="data_12_rest_out.norm.nonparametric.full.csv", row.names=FALSE, col.names = FALSE)
write.table(data_11_task_out.norm.nonparametric.partial, quote=FALSE, sep = ",", file="data_11_task_out.norm.nonparametric.partial.csv", row.names=FALSE, col.names = FALSE)
write.table(data_12_task_out.norm.nonparametric.partial, quote=FALSE, sep = ",", file="data_12_task_out.norm.nonparametric.partial.csv", row.names=FALSE, col.names = FALSE)
write.table(data_11_task_out.norm.nonparametric.full, quote=FALSE, sep = ",", file="data_11_task_out.norm.nonparametric.full.csv", row.names=FALSE, col.names = FALSE)
write.table(data_12_task_out.norm.nonparametric.full, quote=FALSE, sep = ",", file="data_12_task_out.norm.nonparametric.full.csv", row.names=FALSE, col.names = FALSE)

setwd(datadir2)
write.table(data_21_rest_out.norm.nonparametric.partial, quote=FALSE, sep = ",", file="data_21_rest_out.norm.nonparametric.partial.csv", row.names=FALSE, col.names = FALSE)
write.table(data_22_rest_out.norm.nonparametric.partial, quote=FALSE, sep = ",", file="data_22_rest_out.norm.nonparametric.partial.csv", row.names=FALSE, col.names = FALSE)
write.table(data_21_rest_out.norm.nonparametric.full, quote=FALSE, sep = ",", file="data_21_rest_out.norm.nonparametric.full.csv", row.names=FALSE, col.names = FALSE)
write.table(data_22_rest_out.norm.nonparametric.full, quote=FALSE, sep = ",", file="data_22_rest_out.norm.nonparametric.full.csv", row.names=FALSE, col.names = FALSE)
write.table(data_21_task_out.norm.nonparametric.partial, quote=FALSE, sep = ",", file="data_21_task_out.norm.nonparametric.partial.csv", row.names=FALSE, col.names = FALSE)
write.table(data_22_task_out.norm.nonparametric.partial, quote=FALSE, sep = ",", file="data_22_task_out.norm.nonparametric.partial.csv", row.names=FALSE, col.names = FALSE)
write.table(data_21_task_out.norm.nonparametric.full, quote=FALSE, sep = ",", file="data_21_task_out.norm.nonparametric.full.csv", row.names=FALSE, col.names = FALSE)
write.table(data_22_task_out.norm.nonparametric.full, quote=FALSE, sep = ",", file="data_22_task_out.norm.nonparametric.full.csv", row.names=FALSE, col.names = FALSE)

