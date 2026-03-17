library(tidyverse)

# Load data

mental <- read.csv("./Data/PROTECT Source Data/PID PROT_UK_DA_091 PROTECT UK Mental Health FINAL v1.0.csv")

# Mental health baseline data only needed (these will be merged in with Index MBI data from wave 1, non-matches will be dropped as they will be from later waves)

mental_bl <- mental %>%
  mutate(DateCompleted = as.Date(DateCompleted, "%d/%m/%Y")) %>%
  arrange(ParticipantId , DateCompleted) %>%
  group_by(ParticipantId ) %>%
  slice(1) %>%
  ungroup()

# Convert all PHQ9 to numeric (this will mean that NULL/DA are coded as missing), change 99 to missing

mental_bl <- mental_bl %>%
  mutate(across(MH6:MH14, as.numeric)) %>%
  mutate(across(MH6:MH14, ~ if_else(.x == 99, NA_real_, .x)))

# Recode and convert to factor.  Treat NULL as No, due to optional selection of questions (i.e., those with no symptoms more likely to have left both blank)

# Vector of target variable names

vars_to_recode <- c("MH3_DEP","MH3_MANIA","MH3_ANXIETY","MH3_SOCIAL","MH3_AGRO","MH3_PANIC","MH3_OCD","MH4_ANOREXIA","MH4_BULIMIA","MH4_OVEREATING","MH4_SCHIZL"     
  ,"MH4_PSYCHOSIS","MH4_PERSONALITY","MH4_AUTISM","MH4_ADHD")

mental_bl <- mental_bl %>%
  mutate(across(all_of(vars_to_recode), ~ factor(
    case_when(
      .x == "0" ~ 0,
      .x == "1" ~ 1,
      .x == "NULL" ~ 0,
      TRUE ~ NA_real_
    ),
    levels = c(0, 1),
    labels = c("No", "Yes")
  )))

# Create PHQ9 Total score, then create depression variable (in PROTECT this is coded 1-4 instead of 0-3 so do normal PHQ cut offs +9)

mental_bl <- mental_bl %>%
  mutate(PHQ9Tot = rowSums(across(MH6:MH14), na.rm = FALSE), 
    PHQCat = case_when(
    PHQ9Tot < 14 ~ 0,
    PHQ9Tot >= 14 & PHQ9Tot <= 18 ~ 1,
    PHQ9Tot >= 19 & PHQ9Tot <= 23 ~ 2,
    PHQ9Tot >= 24 ~ 3
  ))



# Test relationship between PHQ9 and mental health var.  If highly correlated, just include PHQ9 to avoid multicolinearity

# List of mental health and medical variables
mh_vars <- c("MH3_DEP","MH3_MANIA","MH3_ANXIETY","MH3_SOCIAL","MH3_AGRO","MH3_PANIC","MH3_OCD","MH4_ANOREXIA","MH4_BULIMIA","MH4_OVEREATING","MH4_SCHIZL"     
             ,"MH4_PSYCHOSIS","MH4_PERSONALITY","MH4_AUTISM","MH4_ADHD")


# Pivot to long format
long_datmh <- mental_bl %>%
  select(PHQCat, all_of(mh_vars)) %>%
  pivot_longer(cols = all_of(mh_vars), names_to = "Condition", values_to = "Level") %>%
  filter(!is.na(PHQCat), !is.na(Level))

# Summarise counts and percentages
summary_table <- long_datmh %>%
  group_by(Condition, PHQCat, Level) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(Condition, Level) %>%
  mutate(pct = round(100 * n / sum(n), 1),
         cell = paste0(n, " (", pct, "%)")) %>%
  ungroup() %>%
  select(Condition, Level, PHQCat, cell) %>%
  pivot_wider(names_from = PHQCat, values_from = cell)

# View result
print(summary_table, n = Inf)

# Save (only cols needed for analysis, reordered to logical order)


colnames(mental_bl) 

write_rds(mental_bl[,c(1,7:27,118,119)], "./Data/Interim/mentalhealth_clean.rds")



