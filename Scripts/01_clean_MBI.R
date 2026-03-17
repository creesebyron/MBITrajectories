library(tidyverse)
library(ggplot2)

# Do MBI self and MBI proxy separately to maximise sample size (fewer people nominated a study partner)

### MBI SELF ####

# Load data

mbi_self <- read.csv("./Data/PROTECT Source Data/PID TEC PROT_UK_DA_090 MBI Self January 2024 v1.0 QR.csv")

# Convert all to numeric (this will mean that DA are coded as missing)

colnames(mbi_self)

mbi_self <- mbi_self %>%
  mutate(across(c(MBIS_1:MBIS_38), ~ as.numeric(as.character(.))))

# Recode 99 to missing

mbi_self <- mbi_self %>%
  mutate(across(MBIS_1:MBIS_38, ~ if_else(.x == 99, NA_real_, .x)))

# Convert DateCompleted to Date format and sort by ParticipantId and DateCompleted

mbi_self <- mbi_self %>%
  mutate(DateCompleted = as.Date(DateCompleted, "%d/%m/%Y")) %>%
  arrange(ParticipantId, DateCompleted)

# Create a total visits variable

mbi_self <- mbi_self %>%
  group_by(ParticipantId) %>%
  mutate(TotalVisitsMBIS = n()) %>%
  ungroup()

summary(mbi_self$TotalVisitsMBIS)
table(mbi_self$TotalVisitsMBIS)

# Create a visit number variable by assigning sequential visit numbers for each participant

mbi_self <- mbi_self %>%
  group_by(ParticipantId) %>%
  mutate(VisitNumSelf = row_number()) %>%
  ungroup()

# create Index Date (earliest of DateCompleted)

mbi_self <- mbi_self %>%
  group_by(ParticipantId) %>%
  mutate(IndexDateMBISelf = DateCompleted[which.min(VisitNumSelf)]) %>%
  ungroup()

head(mbi_self)

# Visualise waves of enrollment dates

# Draw histogram with monthly breaks
ggplot(mbi_self, aes(x = IndexDateMBISelf)) +
  geom_histogram(binwidth = 10, color = "black", fill = "skyblue") +
  labs(
    title = "Histogram of Registration Dates",
    x = "Date Completed",
    y = "Frequency"
  ) +
  theme_minimal()


# Calculate MBI total score and domain total scores

mbi_self <- mbi_self %>%
  mutate(mbiSTotal = rowSums(across(MBIS_1:MBIS_38), na.rm = FALSE),   # Total mbi score
         mbiSTotD1 = rowSums(across(MBIS_1:MBIS_6), na.rm = FALSE),    # Domain 1
         mbiSTotD2 = rowSums(across(MBIS_8:MBIS_13), na.rm = FALSE),   # Domain 2
         mbiSTotD3 = rowSums(across(MBIS_15:MBIS_26), na.rm = FALSE),   # Domain 3
         mbiSTotD4 = rowSums(across(MBIS_28:MBIS_32), na.rm = FALSE),   # Domain 4
         mbiSTotD5 = rowSums(across(MBIS_34:MBIS_38), na.rm = FALSE))   # Domain 5

# Rename DateCompleted to avoid confusion with other source data files

mbi_self <- mbi_self %>%
  rename(DateCompletedMBIS = DateCompleted)

# Save (only cols needed for analysis, reordered to logical order)

colnames(mbi_self)

write_rds(mbi_self[, c(1,2,4,40:48)], "./Data/Interim/mbi_clean.rds")

### MBI PROXY ###


# Load data

mbi_proxy <- read.csv("./Data/PROTECT Source Data/PID TEC PROT_UK_DA_090 MBI Informant January 2024 v1.0 QR.csv")

# Convert all to numeric (this will mean that DA are coded as missing)

colnames(mbi_proxy)

mbi_proxy <- mbi_proxy %>%
  mutate(across(c(MBII_1:MBII_38), ~ as.numeric(as.character(.))))

# Recode 99 to missing

mbi_proxy <- mbi_proxy %>%
  mutate(across(MBII_1:MBII_38, ~ if_else(.x == 99, NA_real_, .x)))

# Convert DateCompleted to Date format and sort by ParticipantId and DateCompleted

mbi_proxy <- mbi_proxy %>%
  mutate(DateCompleted = as.Date(DateCompleted, "%d/%m/%Y")) %>%
  arrange(ParticipantId, DateCompleted)

# Create a total visits variable

mbi_proxy <- mbi_proxy %>%
  group_by(ParticipantId) %>%
  mutate(TotalVisitsMBII = n()) %>%
  ungroup()

summary(mbi_proxy$TotalVisitsMBII)
table(mbi_proxy$TotalVisitsMBII)

# Create a visit number variable by assigning sequential visit numbers for each participant

mbi_proxy <- mbi_proxy %>%
  group_by(ParticipantId) %>%
  mutate(VisitNumProxy = row_number()) %>%
  ungroup()

# create Index Date (earliest of DateCompleted)

mbi_proxy <- mbi_proxy %>%
  group_by(ParticipantId) %>%
  mutate(IndexDateMBIProxy = DateCompleted[which.min(VisitNumProxy)]) %>%
  ungroup()

head(mbi_proxy)

# Visualise waves of enrollment dates

# Draw histogram with monthly breaks
ggplot(mbi_proxy, aes(x = IndexDateMBIProxy)) +
  geom_histogram(binwidth = 10, color = "black", fill = "skyblue") +
  labs(
    title = "Histogram of Registration Dates",
    x = "Date Completed",
    y = "Frequency"
  ) +
  theme_minimal()


# Calculate MBI total score and domain total scores

mbi_proxy <- mbi_proxy %>%
  mutate(MBIITotal = rowSums(across(MBII_1:MBII_38), na.rm = FALSE),   # Total mbi score
         MBIITotD1 = rowSums(across(MBII_1:MBII_6), na.rm = FALSE),    # Domain 1
         MBIITotD2 = rowSums(across(MBII_8:MBII_13), na.rm = FALSE),   # Domain 2
         MBIITotD3 = rowSums(across(MBII_15:MBII_26), na.rm = FALSE),   # Domain 3
         MBIITotD4 = rowSums(across(MBII_28:MBII_32), na.rm = FALSE),   # Domain 4
         MBIITotD5 = rowSums(across(MBII_34:MBII_38), na.rm = FALSE))   # Domain 5

# Rename DateCompleted to avoid confusion with other source data files

mbi_proxy <- mbi_proxy %>%
  rename(DateCompletedMBII = DateCompleted)

# Save (only cols needed for analysis, reordered to logical order)

colnames(mbi_proxy)

write_rds(mbi_proxy[, c(1,2,5,41:49)], "./Data/Interim/mbi_proxy_clean.rds")

