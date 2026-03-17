library(tidyverse)
library(data.table)


# Read in data for MBI, IQCODE, Demographics, Medical History, Mental Health, Cognition and APOE

mbi <- readRDS("./Data/Analytical Cohort/AnalyticalData_MBISelf.rds") # read use the cleaned and merged MBI data as a start
phq9 <- readRDS("./Data/Interim/phq_long_clean.rds")


# MBI is the index file, so first MBI assessment is the index date for participants

# IQCODE observation closest to IndexDate in MBI data needs to be the index IQCODE date.  This is necessary as MBI was introduced after IQCODE so it is not always the case that the first MBI = the first IQCODE

# Subset mbi to one index date per person (e.g., earliest)

index <- mbi %>%
  group_by(ParticipantId) %>%
  slice_min(VisitNumSelf , with_ties = FALSE) %>%
  ungroup()

# Join mbi index date to all phq rows, compute date difference with IQS
phq9_with_index <- phq9 %>%
  left_join(index %>% select(ParticipantId, IndexDateMBISelf), by = "ParticipantId") %>%
  mutate(date_diff = as.numeric(DateCompleted - IndexDateMBISelf),
         abs_diff = abs(date_diff))

# Keep only phq rows within ±6 months of the index date (365.25 days)
phq9_within_6m <- phq9_with_index %>%
  filter(abs_diff <= 182.63)

# Find the closest row within that window
closest_phq <- phq9_within_6m %>%
  group_by(ParticipantId) %>%
  slice_min(order_by = abs_diff, n = 1, with_ties = FALSE) %>%
  select(ParticipantId, closest_phq_date = DateCompleted)

# Keep all subsequent iqs rows from the closest date onward
phq_subsequent <- phq9 %>%
  left_join(closest_phq, by = "ParticipantId") %>%
  filter(!is.na(closest_phq_date), DateCompleted >= closest_phq_date)

# Fuzzy join

# Create interval for iqs (±6 months around DateCompletedIQS: two new variables one 6 months earlier; one 6 months later)

# drop the previous start and end variables (created in 08_merge_mbi_all) and all those other variables that depend on time (these will be created later)

mbi <- mbi %>%
  select(-all_of(c("start","end", "i.start", "i.end", "AgeVisit", "AgeC", "TotalVisits", "FollowUpTime", "FollowUpMax", "DropOut")))

mbi <- as.data.table(mbi)

mbi[, start := DateCompletedIQS  %m-% months(6)]  # 6 months prior
mbi[, end := DateCompletedIQS  %m+% months(6)]    # 6 months after

# Create point interval for cog_subsequaent (start = end = DateCompleted)

phq_subsequent <- as.data.table(phq_subsequent)

phq_subsequent[, start := DateCompleted]
phq_subsequent[, end := DateCompleted]

# Set keys for the cognitive data and IQS date

setkey(mbi, ParticipantId, start, end)
setkey(phq_subsequent, ParticipantId, start, end)

# Perform fuzzy join using foverlaps (only keep matches within window and same ID)

nrow_beforephq <- length(unique(mbi$ParticipantId))

iqsphq <- foverlaps(phq_subsequent, mbi,
                    by.x = c("ParticipantId", "start", "end"),
                    by.y = c("ParticipantId", "start", "end"),
                    type = "within",
                    nomatch = 0L)

# The fuzzy join can introduce duplicates, remove them, keeping the earliest observation

iqsphq_unique <- iqsphq %>%
  group_by(ParticipantId, DateCompletedIQS) %>%
  slice_min(order_by = DateCompleted, n = 1, with_ties = FALSE) %>%
  ungroup()

nrow_afterphq <- length(unique(iqsphq_unique$ParticipantId))

nrow_beforephq - nrow_afterphq # 116 cases dropped


# Compute absolute time difference (in days)

iqsphq_unique <- iqsphq_unique %>% 
  mutate(date_diff = as.numeric(DateCompletedIQS - DateCompleted),
         abs_diff = abs(date_diff))

# Confirm that all IQ and PHQ observations are within 6 months

summary(iqsphq_unique$abs_diff)

# Determine time between each visit

# Compute time differences between visits for each participant (generates NA for all baselines as there is no lag there)
forplot <- iqsphq_unique %>%
  arrange(ParticipantId, DateCompleted) %>%
  group_by(ParticipantId) %>%
  mutate(DaysSinceLast = as.numeric(DateCompleted - lag(DateCompleted))) %>%
  ungroup()

avg_gap_overall <- iqsphq_unique %>%
  summarise(mean_gap_days = mean(DaysSinceLast, na.rm = TRUE),
            sd_gap_days = sd(DaysSinceLast, na.rm = TRUE))

ggplot(forplot, aes(x = DaysSinceLast)) +
  geom_histogram(binwidth = 30, fill = "skyblue", color = "white") +
  labs(title = "Distribution of Time Between Visits",
       x = "Days Between Visits",
       y = "Count") +
  theme_minimal()


# Create an age at visit variable

iqsphq_unique <- as.data.frame(iqsphq_unique %>% 
                                        mutate(
                                          AgeVisit = interval(DOBapprox, DateCompleted) %/% years(1)
                                        ))

# Create Age at baseline variable

iqsphq_unique <- as.data.frame(iqsphq_unique %>% group_by(ParticipantId) %>%
                                        mutate(AgeBL = AgeVisit[which.min(DateCompleted)]))


summary(iqsphq_unique)

# Check all in wave 1

summary(iqsphq_unique$IndexDateMBISelf)

# Check total Visits

iqsphq_unique <- iqsphq_unique %>%
  group_by(ParticipantId) %>%
  mutate(TotalVisits = n()) %>%
  ungroup()

summary(iqsphq_unique$TotalVisits)

# Drop out total visits <=1

iqsphq_unique_2visit <- iqsphq_unique %>%
  filter(TotalVisits>1 )

summary(iqsphq_unique_2visit$TotalVisits)

n_iqsphq_unique_2visit <- length(unique(iqsphq_unique_2visit$ParticipantId))
n_iqsphq_unique <- length(unique(iqsphq_unique$ParticipantId))

n_iqsphq_unique_2visit - n_iqsphq_unique # 85 excluded

# Check age at BL >50

summary(iqsphq_unique_2visit$AgeBL) # min is >50 so OK


# Create average duration of follow up variable 

iqsphq_unique_2visit$FollowUpTime <- abs(time_length(difftime(iqsphq_unique_2visit$DateCompleted, iqsphq_unique_2visit$IndexDateMBISelf), "months"))

iqsphq_unique_2visit<- as.data.frame(iqsphq_unique_2visit %>%
                                       group_by(ParticipantId) %>%
                                       mutate(FollowUpMax = FollowUpTime[which.max(FollowUpTime)]))

summary(iqsphq_unique_2visit$FollowUpMax)

# Drop case where follow up is < 6 months

niqsphq_unique_2visit <- length(unique(iqsphq_unique_2visit$ParticipantId))

iqsphq_unique_2visit <- iqsphq_unique_2visit %>%
  filter(FollowUpMax>6)

niqsphq_unique_2visit_6mfilter <- length(unique(iqsphq_unique_2visit$ParticipantId))

niqsphq_unique_2visit_6mfilter - niqsphq_unique_2visit# 0 excluded


summary(iqsphq_unique_2visit$FollowUpMax)

# Create a binary indicator variable for dropping out prior to completing the per protocol 7 years from Jan 2017 (allow a six month window, so 78 months)

iqsphq_unique_2visit <- iqsphq_unique_2visit %>% mutate(
  DropOut = case_when(
    FollowUpMax  >= 81 ~ 0,
    FollowUpMax <81 ~ 1
  ))


# Create a grand mean centred Age Visit

grand_mean_age <- mean(iqsphq_unique_2visit$AgeVisit)

iqsphq_unique_2visit$AgeC <- iqsphq_unique_2visit$AgeVisit - grand_mean_age

hist(iqsphq_unique_2visit$AgeC)

# Save 

write_rds(iqsphq_unique_2visit, "./Data/Analytical Cohort/AnalyticalData_PHQ.rds")



