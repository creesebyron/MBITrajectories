library(tidyverse)

# Load data

iqs <- read.csv("./Data/PROTECT Source Data/PID TEC PROT_UK_DA_090 IQCode January 2024 v1.0 QR.csv")

# Longitudinal IQCODE data

iqs <- iqs %>%
  mutate(DateCompleted = as.Date(DateCompleted, "%d/%m/%Y")) %>%  # ensure date column is in Date format
  arrange(ParticipantId, DateCompleted)

# Convert all to numeric (this will mean that DA are coded as missing)

colnames(iqs)

iqs <- iqs %>%
  mutate(across(c(IQCODES_1:IQCODES_16), ~ as.numeric(as.character(.))))

# Create IQCODE average score, then make cognitive impairment variable 

iqs <- iqs %>% mutate(IQCODEMeanS = rowMeans(across(IQCODES_1:IQCODES_16), na.rm = FALSE),
                      SCI3.3 = case_when(
                        IQCODEMeanS  <3.3 ~ 0,
                        IQCODEMeanS >3.3 ~ 1,
                        TRUE ~ NA_real_))

iqs <- iqs %>% mutate(IQCODEMeanS = rowMeans(across(IQCODES_1:IQCODES_16), na.rm = FALSE),
                      SCI3.6 = case_when(
                        IQCODEMeanS  <3.6 ~ 0,
                        IQCODEMeanS >3.6 ~ 1,
                        TRUE ~ NA_real_))

# Rename DateCompleted to avoid confusion with other source data files

iqs <- iqs %>%
  rename(DateCompletedIQS = DateCompleted)

# Create a total visits variable

iqs <- iqs %>%
  group_by(ParticipantId) %>%
  mutate(TotalVisitsIQS = n()) %>%
  ungroup()

summary(iqs$TotalVisitsIQS)
table(iqs$TotalVisitsIQS)

# Create a visit number variable by assigning sequential visit numbers for each participant

iqs <- iqs %>%
  group_by(ParticipantId) %>%
  mutate(VisitNumIQS = row_number()) %>%
  ungroup()

# create Index Date (earliest of DateCompleted)

iqs <- iqs %>%
  group_by(ParticipantId) %>%
  mutate(IndexDateIQS = DateCompletedIQS[which.min(VisitNumIQS)]) %>%
  ungroup()

# Save (only cols needed for analysis, reordered to logical order)

colnames(iqs)

write_rds(iqs[, c(1,4,23:27)], "./Data/Interim/iqs_clean.rds")





