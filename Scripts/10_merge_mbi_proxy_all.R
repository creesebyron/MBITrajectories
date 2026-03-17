library(tidyverse)
library(data.table)

# Read in data for MBI, IQCODE, Demographics, Medical History, Mental Health, Cognition and APOE

mbi <- readRDS("./Data/Interim/mbi_proxy_clean.rds")
iqs <- readRDS("./Data/Interim/iqs_clean.rds")
dem <- readRDS("./Data/Interim/demogs_clean.rds")
mnh <- readRDS("./Data/Interim/mentalhealth_clean.rds")
mdh <- readRDS("./Data/Interim/medhist_clean.rds")
ape <- readRDS("./Data/Interim/apoe_clean.rds")

# MBI is the index file, so first MBI assessment is the index date for participants

# IQCODE observation closest to IndexDate in MBI data needs to be the index IQCODE date.  This is necessary as MBI was introduced after IQCODE so it is not always the case that the first MBI = the first IQCODE

# Subset mbi to one index date per person (e.g., earliest)

length(unique(mbi$ParticipantId))

mbi_index <- mbi %>%
  group_by(ParticipantId) %>%
  slice_min(VisitNumProxy , with_ties = FALSE) %>%
  ungroup()

# Join mbi index date to all iqs rows, compute date difference
iqs_with_index <- iqs %>%
  left_join(mbi_index %>% select(ParticipantId, IndexDateMBIProxy), by = "ParticipantId") %>%
  mutate(date_diff = as.numeric(DateCompletedIQS - IndexDateMBIProxy),
         abs_diff = abs(date_diff))

# Keep only iqs rows within ±6 months (365.25 days)
nrow(iqs_with_index)
length(unique(iqs_with_index$ParticipantId))

iqs_within_6m <- iqs_with_index %>%
  filter(abs_diff <= 182.63)

nrow(iqs_within_6m)
length(unique(iqs_within_6m$ParticipantId))

# Find the closest row within that window
closest_iqs <- iqs_within_6m %>%
  group_by(ParticipantId) %>%
  slice_min(order_by = abs_diff, n = 1, with_ties = FALSE) %>%
  select(ParticipantId, closest_iqs_date = DateCompletedIQS)

# Keep all subsequent iqs rows from the closest date onward
iqs_subsequent <- iqs %>%
  left_join(closest_iqs, by = "ParticipantId") %>%
  filter(!is.na(closest_iqs_date), DateCompletedIQS >= closest_iqs_date)

# Fuzzy join

# Create interval for mbis (±6 months around DateCompletedMBIS: two new variables one 6 months earlier; one 6 months later)

mbi <- as.data.table(mbi)

mbi[, start := DateCompletedMBII %m-% months(6)]  # 6 months prior
mbi[, end := DateCompletedMBII %m+% months(6)]    # 6 months after

# Create point interval for iqs_subsequent (start = end = DateCompletedIQS)

iqs_subsequent <- as.data.table(iqs_subsequent)

iqs_subsequent[, start := DateCompletedIQS]
iqs_subsequent[, end := DateCompletedIQS]

# Set keys

setkey(mbi, ParticipantId, start, end)
setkey(iqs_subsequent, ParticipantId, start, end)

# Perform fuzzy join using foverlaps (only keep matches within window and same ID)

nrow_beforeiqs <- length(unique(mbi$ParticipantId))

iqsmbi <- foverlaps(iqs_subsequent, mbi,
                    by.x = c("ParticipantId", "start", "end"),
                    by.y = c("ParticipantId", "start", "end"),
                    type = "within",
                    nomatch = 0L)

# The fuzzy join can introduce duplicates, remove them, keeping the earliest observation

iqsmbi_unique <- iqsmbi %>%
  group_by(ParticipantId, DateCompletedMBII) %>%
  slice_min(order_by = DateCompletedIQS, n = 1, with_ties = FALSE) %>%
  ungroup()

length(unique(iqsmbi_unique$ParticipantId))

# 17833 index IQCODE not within 6m of MBI-C, or no MBI-C at all


# Compute absolute time difference (in days)

iqsmbi_unique <- iqsmbi_unique %>% 
  mutate(date_diff = as.numeric(DateCompletedIQS - DateCompletedMBII),
         abs_diff = abs(date_diff))

# Confirm that all IQ and MBI observations are within 6 months

summary(iqsmbi_unique$abs_diff)

# Merge the IQCODE and MBI data with the baseline demographic, medical, mental health

clinicalmerge <- iqsmbi_unique %>%
  left_join(dem, by = "ParticipantId") %>%           
  left_join(mdh, by = "ParticipantId") %>%
  left_join(mnh, by = "ParticipantId")


# Drop any missing data across MBI, IQCODE, PHQ9, BMI, EDU, SEX

clinicalmerge_nomiss <- clinicalmerge %>%
  filter(if_all(c(MBIITotal , SCI3.3, PHQCat, BMIcat, EDU, SEX), ~ !is.na(.)))

n_clinicalmerge <- length(unique(clinicalmerge$ParticipantId))
n_clinicalmerge_nomiss <- length(unique(clinicalmerge_nomiss$ParticipantId)) # check how many participants in full merged after filtering out missingness on MBI and IQCODE

n_clinicalmerge - n_clinicalmerge_nomiss # 1216 cases dropped 

# check there are no missing data across key covariates

misscheck <- clinicalmerge_nomiss %>%
  filter(if_all(c(MBIITotal , SCI3.3, PHQCat, BMIcat, EDU, SEX), ~ is.na(.)))

nrow(misscheck) # no missing


# Determine time between each visit

# Compute time differences between visits for each participant
clinicalmerge_nomiss <- clinicalmerge_nomiss %>%
  arrange(ParticipantId, DateCompletedMBII) %>%
  group_by(ParticipantId) %>%
  mutate(DaysSinceLast = as.numeric(DateCompletedMBII - lag(DateCompletedMBII))) %>%
  ungroup()

avg_gap_overall <- clinicalmerge_nomiss %>%
  summarise(mean_gap_days = mean(DaysSinceLast, na.rm = TRUE),
            sd_gap_days = sd(DaysSinceLast, na.rm = TRUE))

ggplot(clinicalmerge_nomiss, aes(x = DaysSinceLast)) +
  geom_histogram(binwidth = 30, fill = "skyblue", color = "white") +
  labs(title = "Distribution of Time Between Visits",
       x = "Days Between Visits",
       y = "Count") +
  theme_minimal()



# Create an age at visit variable

clinicalmerge_nomiss <- as.data.frame(clinicalmerge_nomiss %>% 
                                        mutate(
                                          AgeVisit = interval(DOBapprox, DateCompletedMBII) %/% years(1)
                                        ))

# Create Age at baseline variable

clinicalmerge_nomiss <- as.data.frame(clinicalmerge_nomiss %>% group_by(ParticipantId) %>%
                                        mutate(AgeBL = AgeVisit[which.min(DateCompletedMBII)]))

summary(clinicalmerge_nomiss$AgeBL)


# Merge in APOE data

fullmerge <- clinicalmerge_nomiss %>%
  inner_join(ape, join_by("TECID" == "TEC_ID"))

# drop the  missing genotypes

fullmerge <- fullmerge %>%
  filter(!is.na(APOEbin))

n_clinicalmerge_nomiss <- length(unique(clinicalmerge_nomiss$ParticipantId)) # check how many participants in full merged 
n_fullmerge <- length(unique(fullmerge$ParticipantId))

n_clinicalmerge_nomiss - n_fullmerge # 8171 with no genetic data

# keep only observations with an Index date in the 12 months to 31 Jan 2017 - wave 1

length(unique(fullmerge$ParticipantId))

fullmerge_wave1 <- fullmerge %>%
  filter(between(IndexDateMBIProxy, as.Date("2016-02-01"), as.Date("2017-01-31")))

length(unique(fullmerge_wave1$ParticipantId))

summary(fullmerge_wave1$IndexDateMBIProxy)

# Drop out cases with only 1 visit completed

# Create a total visits variable

fullmerge_wave1 <- fullmerge_wave1 %>%
  group_by(ParticipantId) %>%
  mutate(TotalVisits = n()) %>%
  ungroup()

summary(fullmerge_wave1$TotalVisits)

fullmerge_wave1_2visit <- fullmerge_wave1 %>%
  filter(TotalVisits>1 )

n_fullmerge_wave1_2visit <- length(unique(fullmerge_wave1_2visit$ParticipantId))
n_fullmerge_wave1 <- length(unique(fullmerge_wave1$ParticipantId))

n_fullmerge_wave1_2visit - n_fullmerge_wave1 # 545 excluded

# Drop out <50 at baseline

fullmerge_wave1_2visit50 <- as.data.frame(fullmerge_wave1_2visit %>% 
                                            filter(AgeBL >= 50))

summary(fullmerge_wave1_2visit50$AgeBL)

n_fullmerge_wave1_2visit50 <- length(unique(fullmerge_wave1_2visit50$ParticipantId))

n_fullmerge_wave1_2visit - n_fullmerge_wave1_2visit50 # 1 excluded



# Create average duration of follow up variable 

fullmerge_wave1_2visit50$FollowUpTime <- abs(time_length(difftime(fullmerge_wave1_2visit50$DateCompletedMBII, fullmerge_wave1_2visit50$IndexDateMBIProxy), "months"))

fullmerge_wave1_2visit50<- as.data.frame(fullmerge_wave1_2visit50 %>%
                                           group_by(ParticipantId) %>%
                                           mutate(FollowUpMax = FollowUpTime[which.max(FollowUpTime)]))

summary(fullmerge_wave1_2visit50$FollowUpMax)

# Create a binary indicator variable for dropping out prior to completing the per protocol 7 years from Jan 2017 (allow a six month window, so 78 months)

fullmerge_wave1_2visit50 <- fullmerge_wave1_2visit50 %>% mutate(
  DropOut = case_when(
    FollowUpMax  >= 78 ~ 0,
    FollowUpMax <78 ~ 1
  ))


# Create a grand mean centred Age Visit

grand_mean_age <- mean(fullmerge_wave1_2visit50$AgeVisit)

fullmerge_wave1_2visit50$AgeC <- fullmerge_wave1_2visit50$AgeVisit - grand_mean_age

hist(fullmerge_wave1_2visit50$AgeC)

length(unique(fullmerge_wave1_2visit50$ParticipantId))

# Save 

write_rds(fullmerge_wave1_2visit50, "./Data/Analytical Cohort/AnalyticalData_MBIProxy.rds")



































