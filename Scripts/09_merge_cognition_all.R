library(tidyverse)
library(data.table)


# Read in data for MBI, IQCODE, Demographics, Medical History, Mental Health, Cognition and APOE

mbi <- readRDS("./Data/Analytical Cohort/AnalyticalData_MBISelf.rds") # read use the cleaned and merged MBI data as a start
cog <- readRDS("./Data/Interim/cognition_clean.rds")


# MBI is the index file, so first MBI assessment is the index date for participants

# IQCODE observation closest to IndexDate in MBI data needs to be the index IQCODE date.  This is necessary as MBI was introduced after IQCODE so it is not always the case that the first MBI = the first IQCODE

# Subset mbi to one index date per person (e.g., earliest)

index <- mbi %>%
  group_by(ParticipantId) %>%
  slice_min(VisitNumSelf , with_ties = FALSE) %>%
  ungroup()

# Join mbi index date to all cog rows, compute date difference with IQS
cog_with_index <- cog %>%
  left_join(index %>% select(ParticipantId, IndexDateMBISelf), by = "ParticipantId") %>%
  mutate(date_diff = as.numeric(DateCompleted - IndexDateMBISelf),
         abs_diff = abs(date_diff))

# Keep only cog rows within ±6 months of the index date (365.25 days)
cog_within_6m <- cog_with_index %>%
  filter(abs_diff <= 182.63)

# Find the closest row within that window
closest_cog <- cog_within_6m %>%
  group_by(ParticipantId) %>%
  slice_min(order_by = abs_diff, n = 1, with_ties = FALSE) %>%
  select(ParticipantId, closest_cog_date = DateCompleted)

# Keep all subsequent iqs rows from the closest date onward
cog_subsequent <- cog %>%
  left_join(closest_cog, by = "ParticipantId") %>%
  filter(!is.na(closest_cog_date), DateCompleted >= closest_cog_date)

# Fuzzy join

# Create interval for iqs (±6 months around DateCompletedIQS: two new variables one 6 months earlier; one 6 months later)

# drop the previous start and end variables (created in 08_merge_mbi_all) and all those other variables that depend on time (these will be created later)

mbi <- mbi %>%
  select(-all_of(c("start","end", "i.start", "i.end", "AgeVisit", "AgeC", "TotalVisits", "FollowUpTime", "FollowUpMax", "DropOut")))

mbi <- as.data.table(mbi)

mbi[, start := DateCompletedIQS  %m-% months(6)]  # 6 months prior
mbi[, end := DateCompletedIQS  %m+% months(6)]    # 6 months after

# Create point interval for cog_subsequaent (start = end = DateCompleted)

cog_subsequent <- as.data.table(cog_subsequent)

cog_subsequent[, start := DateCompleted]
cog_subsequent[, end := DateCompleted]

# Set keys for the cognitive data and IQS date

setkey(mbi, ParticipantId, start, end)
setkey(cog_subsequent, ParticipantId, start, end)

# Perform fuzzy join using foverlaps (only keep matches within window and same ID)

nrow_beforecog <- length(unique(mbi$ParticipantId))

iqscog <- foverlaps(cog_subsequent, mbi,
                    by.x = c("ParticipantId", "start", "end"),
                    by.y = c("ParticipantId", "start", "end"),
                    type = "within",
                    nomatch = 0L)

# The fuzzy join can introduce duplicates, remove them, keeping the earliest observation

iqscog_unique <- iqscog %>%
  group_by(ParticipantId, DateCompletedIQS) %>%
  slice_min(order_by = DateCompleted, n = 1, with_ties = FALSE) %>%
  ungroup()

nrow_aftercog <- length(unique(iqscog_unique$ParticipantId))

nrow_beforecog - nrow_aftercog # 125 cases dropped


# Compute absolute time difference (in days)

iqscog_unique <- iqscog_unique %>% 
  mutate(date_diff = as.numeric(DateCompletedIQS - DateCompleted),
         abs_diff = abs(date_diff))

# Confirm that all IQ and Cog observations are within 6 months

summary(iqscog_unique$abs_diff)

# Determine time between each visit

# Compute time differences between visits for each participant (generates NA for all baselines as there is no lag there)
forplot <- iqscog_unique %>%
  arrange(ParticipantId, DateCompleted) %>%
  group_by(ParticipantId) %>%
  mutate(DaysSinceLast = as.numeric(DateCompleted - lag(DateCompleted))) %>%
  ungroup()

avg_gap_overall <- iqscog_unique %>%
  summarise(mean_gap_days = mean(DaysSinceLast, na.rm = TRUE),
            sd_gap_days = sd(DaysSinceLast, na.rm = TRUE))

ggplot(forplot, aes(x = DaysSinceLast)) +
  geom_histogram(binwidth = 30, fill = "skyblue", color = "white") +
  labs(title = "Distribution of Time Between Visits",
       x = "Days Between Visits",
       y = "Count") +
  theme_minimal()


# Create an age at visit variable

iqscog_unique <- as.data.frame(iqscog_unique %>% 
                                    mutate(
                                      AgeVisit = interval(DOBapprox, DateCompleted) %/% years(1)
                                    ))

# Create Age at baseline variable

iqscog_unique <- as.data.frame(iqscog_unique %>% group_by(ParticipantId) %>%
                                        mutate(AgeBL = AgeVisit[which.min(DateCompleted)]))

summary(iqscog_unique$AgeBL)

# Check all in wave 1

summary(iqscog_unique$IndexDateMBISelf)

# Check total Visits

iqscog_unique <- iqscog_unique %>%
  group_by(ParticipantId) %>%
  mutate(TotalVisits = n()) %>%
  ungroup()

summary(iqscog_unique$TotalVisits)

# Drop out total visits <=1

iqscog_unique_2visit <- iqscog_unique %>%
  filter(TotalVisits>1 )

summary(iqscog_unique_2visit$TotalVisits)

n_iqscog_unique_2visit <- length(unique(iqscog_unique_2visit$ParticipantId))
n_iqscog_unique <- length(unique(iqscog_unique$ParticipantId))

n_iqscog_unique_2visit - n_iqscog_unique # 680 excluded

# Check age at BL >50

summary(iqscog_unique_2visit$AgeBL) # min is >50 so OK


# Create average duration of follow up variable 

iqscog_unique_2visit$FollowUpTime <- abs(time_length(difftime(iqscog_unique_2visit$DateCompleted, iqscog_unique_2visit$IndexDateMBISelf), "months"))

iqscog_unique_2visit<- as.data.frame(iqscog_unique_2visit %>%
                                           group_by(ParticipantId) %>%
                                           mutate(FollowUpMax = FollowUpTime[which.max(FollowUpTime)]))

summary(iqscog_unique_2visit$FollowUpMax)

# Drop case where follow up is < 6 months

niqscog_unique_2visit <- length(unique(iqscog_unique_2visit$ParticipantId))

iqscog_unique_2visit <- iqscog_unique_2visit %>%
                            filter(FollowUpMax>6)

niqscog_unique_2visit_6mfilter <- length(unique(iqscog_unique_2visit$ParticipantId))

niqscog_unique_2visit_6mfilter- niqscog_unique_2visit# 6 excluded


summary(iqscog_unique_2visit$FollowUpMax)

# Create a binary indicator variable for dropping out prior to completing the per protocol follow up

iqscog_unique_2visit <- iqscog_unique_2visit %>% mutate(
  DropOut = case_when(
    FollowUpMax  >= 60 ~ 0,
    FollowUpMax <60 ~ 1
  ))


# Create a grand mean centred Age Visit

grand_mean_age <- mean(iqscog_unique_2visit$AgeVisit)

iqscog_unique_2visit$AgeC <- iqscog_unique_2visit$AgeVisit - grand_mean_age

hist(iqscog_unique_2visit$AgeC)

# Save 

write_rds(iqscog_unique_2visit, "./Data/Analytical Cohort/AnalyticalData_Cognitive.rds")








