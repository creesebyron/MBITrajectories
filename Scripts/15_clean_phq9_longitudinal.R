library(tidyverse)

# Load data

mh <- read.csv("./Data/PROTECT Source Data/PID PROT_UK_DA_091 PROTECT UK Mental Health FINAL v1.0.csv")


mh <- mh %>%
  mutate(DateCompleted = as.Date(DateCompleted, "%d/%m/%Y")) %>%
  arrange(ParticipantId , DateCompleted)

# Keep only relevant variables 

colnames(mh)

phq <- mh %>% select(ParticipantId, DateCompleted, Timepoint, MH6:MH14)

colnames(phq)

# Convert all PHQ9 to numeric (this will mean that NULL/DA are coded as missing), change 99 to missing

phq <- phq %>%
  mutate(across(MH6:MH14, as.numeric)) %>%
  mutate(across(MH6:MH14, ~ if_else(.x == 99, NA_real_, .x)))

lapply(phq, class)

# get missingness

lapply(phq[, 4:12], function(x) table(x, useNA = "ifany"))


# Make PHQ total score.  PHQ is on scale 1-4, change so it is on the original scale of 0-3: subtract 9 from total score

phq <- phq %>%
  mutate(PHQ9Tot = rowSums(across(MH6:MH14), na.rm = FALSE), 
         PHQCat = case_when(
           PHQ9Tot < 14 ~ 0,
           PHQ9Tot >= 14 & PHQ9Tot <= 18 ~ 1,
           PHQ9Tot >= 19 & PHQ9Tot <= 23 ~ 2,
           PHQ9Tot >= 24 ~ 3
         )) %>%
  mutate(PHQ9Tot = PHQ9Tot-9)

# Remove rows where PHQ9Tot is missing

length(unique(phq$ParticipantId)) # N=29586 

phq_nomiss <- phq %>% filter(!is.na(PHQ9Tot))

length(unique(phq_nomiss$ParticipantId)) # N=29362 

length(unique(phq_nomiss$ParticipantId)) - length(unique(phq$ParticipantId)) 

# 224 dropped due to missing

head(phq_nomiss)

# Save 

write_rds(phq_nomiss, "./Data/Interim/phq_long_clean.rds")




