### Prepare analysis of motion parameters
#########################################################
Cohort <- as.character(2)    # which cohort? (1 or 2)
Timepoint <- 1               # which timepoint? (1 or 2)

#########################################################
### (A) Required libraries
#########################################################
if (!require("dplyr")) install.packages("dplyr", dependencies = TRUE)
if (!require("tidyverse")) install.packages("tidyverse", dependencies = TRUE)
library(dplyr)
library(tidyverse)

#########################################################
### (B) Set paths
#########################################################
if (file.exists("/MyWorkingDirectory/derivatives/")) {
  datadir <- paste0("/MyWorkingDirectory/derivatives/TAU", Cohort, "/RegressNuissance")
  listdir <- "/MyWorkingDirectory/lists"
  outdir <- "/MyWorkingDirectory/stats"
  figuredir <- "/MyWorkingDirectory/stats/Figures"
} else {
  stop("Directory does not exist.")
}

#########################################################
### (C) Read data
#########################################################
setwd(listdir)
sublist <- readLines(paste0("TAU", Cohort, "_T", Timepoint, ".txt"))
sublist <- as.list(sublist)

tasklist <- list("rest", paste0("TAU", Cohort, "_run-1"), paste0("TAU", Cohort, "_run-2"))

Results <- data.frame(matrix(ncol = 1 + 8 * 3, nrow = length(sublist)))
colnames(Results) <- c("ID",
                       paste0(rep(c("rest", "ses1", "ses2"), each = 8), 
                              c(".framewise_displacement", ".dvars", ".trans_X", ".trans_Y", ".trans_Z", ".rot_X", ".rot_Y", ".rot_Z")))

counter <- 0

for (subj in sublist) {
  setwd(datadir)
  counter <- counter + 1
  Results[counter, 1] <- subj
  
  for (task in tasklist) {
    file_path <- paste0(subj, "_ses-", Timepoint, "_task-", task, "_desc-confounds_timeseries.tsv")
    mydata <- read_tsv(file_path)
    mydata <- data.frame(lapply(mydata, as.numeric))
    
    offset <- switch(task,
                     "rest" = 2,
                     paste0("TAU", Cohort, "_run-1") = 10,
                     paste0("TAU", Cohort, "_run-2") = 18)
    
    Results[counter, offset:(offset + 7)] <- sapply(
      c("framewise_displacement", "dvars", "trans_x", "trans_y", "trans_z", "rot_x", "rot_y", "rot_z"),
      function(col) mean(mydata[[col]], na.rm = TRUE)
    )
  }
}

Results$ID <- as.numeric(gsub("sub-s", "", Results$ID))

setwd(outdir)
write.table(Results, paste0("Motion_TAU", Cohort, "_Ses", Timepoint, ".csv"),
            sep = ",", row.names = FALSE, col.names = TRUE, quote = FALSE)