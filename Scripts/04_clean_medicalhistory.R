library(tidyverse)

# Load data

medhist <- read.csv("./Data/PROTECT Source Data/PID TEC PROT_UK_DA_090 Medical History January 2024 v1.0 QR.csv")

#  Medical history baseline data only needed (these will be merged in with Index MBI data from wave 1, non-matches will be dropped as they will be from later waves)

medhist_bl <- medhist %>%
  mutate(DateCompleted = as.Date(DateCompleted, "%d/%m/%Y")) %>%
  arrange(ParticipantId , DateCompleted) %>%
  group_by(ParticipantId ) %>%
  slice(1) %>%
  ungroup()

# Convert all diagnoses to factor, recode NULL to missing

# Vector of target variable names

vars_to_recode <- c(
  "DIOG_HBP", "DIOG_STROKE", "DIOG_HEART", "DIOG_DIAB", "DIOG_MCI", "DIOG_PD", 
  "DIOG_HC", "DIOG_HP", "DIOG_HYP", "DIOG_ARTH", "DIOG_HD", "DIOG_CC", "DIOG_CFR", 
  "DIOG_OSTEO", "DIOG_ASTHMA", "DIOG_EPILEPSY", "DIOG_MND", "DIOG_MS", "DIOG_PAGD", 
  "DIOG_DVT", "DIOG_HIV", "DIOG_AIDS", "DIOG_HEPC"
)

# Recode and convert to factor

medhist_bl <- medhist_bl %>%
  mutate(across(all_of(vars_to_recode), ~ factor(
    case_when(
      .x == "0" ~ 0,
      .x == "1" ~ 1,
      .x == "NULL" ~ 2,
      TRUE ~ NA_real_
    ),
    levels = c(0, 1, 2),
    labels = c("No", "Yes", "Missing")
  )))

# Convert height and weight to numeric

medhist_bl <- medhist_bl %>%
  mutate(across(c(WEIGHT:HEIGHT), ~ as.numeric(as.character(.))))

# Calculate BMI, then convert to levels 

medhist_bl <- medhist_bl %>%
  mutate(
    BMI = WEIGHT / ((HEIGHT / 100) ^ 2)
  )

medhist_bl <- medhist_bl %>%
  mutate(BMIcat = case_when(
    BMI < 18.5 ~ "underweight",
    BMI >= 18.5 & BMI <= 24.9 ~ "healthy range",
    BMI >= 25 & BMI <= 29.9 ~ "overweight",
    BMI >= 30  ~ "obese/severe obese",
    TRUE ~ NA_character_  # for missing or invalid BMI
  ),
  BMIcat = factor(BMIcat, 
                  levels = c("underweight", "healthy range", "overweight", "obese/severe obese")),
  BMIcat = relevel(BMIcat, ref = "healthy range"))


# Drop any cases with missing BMICat

nrow(medhist_bl)

medhist_bl_nomiss <- medhist_bl %>% filter(!is.na(BMIcat))

nrow(medhist_bl_nomiss)

table(medhist_bl_nomiss$BMIcat, useNA = "ifany")

# Save (only chronic disease cols needed for analysis, reordered to logical order)

colnames(medhist_bl) 

write_rds(medhist_bl[,c(1,8:31,49,50)], "./Data/Interim/medhist_clean.rds")


