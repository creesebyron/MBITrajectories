library(glmmTMB)
library(tidyverse)
library(ggplot2)
library(purrr)
library(broom.mixed)
library(DHARMa)
library(dplyr)
library(tidyr)
library(stringr)

# Load data

dat <- readRDS("./Data/Analytical Cohort/AnalyticalData_MBISelf.rds")

# Cross tabs of PHQ cat score against all mental health measures to determine covariates (leave out mental health conditions if PHQ9 is strongly associated).  Use baseline data

# Use baseline data

datbl <- dat %>% filter(VisitNumSelf==1)

# List of mental health and medical variables
mh_vars <- c(
  "MH3_DEP", "MH3_MANIA", "MH3_ANXIETY", "MH3_SOCIAL", "MH3_AGRO",
  "MH3_PANIC", "MH3_OCD", "MH4_SCHIZL", "MH4_PSYCHOSIS", "MH4_AUTISM", "MH4_ADHD"
)

med_vars <- c(
  "DIOG_HBP", "DIOG_STROKE", "DIOG_HEART", "DIOG_DIAB", "DIOG_MCI", "DIOG_PD",
  "DIOG_HC", "DIOG_HP", "DIOG_HYP", "DIOG_ARTH", "DIOG_HD", "DIOG_CC",
  "DIOG_CFR", "DIOG_OSTEO", "DIOG_ASTHMA", "DIOG_EPILEPSY", "DIOG_MND", "DIOG_MS",
  "DIOG_PAGD", "DIOG_DVT", "DIOG_HIV", "DIOG_AIDS", "DIOG_HEPC", "BMIcat"
)

# Pivot to long format
long_dat <- datbl %>%
  select(PHQCat, all_of(mh_vars)) %>%
  pivot_longer(cols = all_of(mh_vars), names_to = "Condition", values_to = "Level") %>%
  filter(!is.na(PHQCat), !is.na(Level))

# Summarise counts and percentages
summary_table <- long_dat %>%
  group_by(Condition, PHQCat, Level) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(Condition, PHQCat) %>%
  mutate(pct = round(100 * n / sum(n), 1),
         cell = paste0(n, " (", pct, "%)")) %>%
  ungroup() %>%
  select(Condition, Level, PHQCat, cell) %>%
  pivot_wider(names_from = PHQCat, values_from = cell)

# View result
print(summary_table, n = Inf)


# Pivot to long format and summarise mean/SD of mbiSTotal
summary_tablemh <- datbl %>%
  select(mbiSTotal, all_of(mh_vars)) %>%
  pivot_longer(cols = all_of(mh_vars), names_to = "Condition", values_to = "Level") %>%
  filter(!is.na(Level), !is.na(mbiSTotal)) %>%
  group_by(Condition, Level) %>%
  summarise(
    Mean_MBIStotal = round(mean(mbiSTotal, na.rm = TRUE), 2),
    SD_MBIStotal   = round(sd(mbiSTotal, na.rm = TRUE), 2),
    .groups = "drop"
  )

# View the result
print(summary_tablemh, n = Inf)

# Pivot to long format and summarise mean/SD of mbiSTotal
summary_tablemed <- datbl %>%
  select(mbiSTotal, all_of(med_vars)) %>%
  pivot_longer(cols = all_of(med_vars), names_to = "Condition", values_to = "Level") %>%
  filter(!is.na(Level), !is.na(mbiSTotal)) %>%
  group_by(Condition, Level) %>%
  summarise(
    Mean_MBIStotal = round(mean(mbiSTotal, na.rm = TRUE), 2),
    SD_MBIStotal   = round(sd(mbiSTotal, na.rm = TRUE), 2),
    .groups = "drop"
  )

# View the result
print(summary_tablemed, n = Inf)




# Specify formulae

formula_negbin_min <- as.formula(paste0(outcome_var, " ~ poly(AgeC, 2)*SCI3.3 + SEX + EDU + DropOut + PHQCat + APOEbin + (1 | ParticipantId)"))

formula_negbin_all <- as.formula(paste0("mbiSTotal", " ~ poly(AgeC, 2)*SCI3.3 + SEX + EDU + EMPLOY + ETH + DropOut + APOEbin +
                                        DIOG_HBP+DIOG_STROKE + DIOG_HEART +DIOG_DIAB+ DIOG_MCI+ DIOG_PD +DIOG_HC +DIOG_HP+DIOG_HYP+ DIOG_ARTH+ DIOG_HD+ DIOG_CC +
                                        DIOG_CFR+ DIOG_OSTEO+ DIOG_ASTHMA+ DIOG_EPILEPSY+ DIOG_MND +DIOG_MS+ DIOG_PAGD +DIOG_DVT +DIOG_HIV +DIOG_AIDS+ DIOG_HEPC +BMIcat + 
                                        MH3_DEP +MH3_MANIA+ MH3_ANXIETY +MH3_SOCIAL+ MH3_AGRO +MH3_PANIC+ MH3_OCD +
                                        MH4_SCHIZL + MH4_PSYCHOSIS+ MH4_AUTISM+ MH4_ADHD + (1 | ParticipantId)"))

formula_zi <- as.formula(paste0(outcome_var, " ~ poly(AgeC, 2)*SCI3.3 + SEX + EDU + DropOut + APOEbin +  (1 | ParticipantId)"))


## SELF RATED MBI-C ##

# Overall then subset by APOE status

negbin_MBISTot <- glmmTMB(mbiSTotal  ~ poly(AgeC, 2)*SCI3.3 + SEX + EDU + DropOut + APOEbin + (1|ParticipantId), data = dat, family=nbinom2 )
summary(negbin_MBISTot )

negbin1_MBISTot <- glmmTMB(formula_negbin_all, data = dat, family=nbinom1 )
summary(negbin1_MBISTot )

zi_MBISTot <- glmmTMB(mbiSTotal  ~ poly(AgeC, 2)*SCI3.3 + SEX + EDU + DropOut + APOEbin + (1|ParticipantId), ziformula = ~SCI3.3+APOEbin,data = dat, family=nbinom2 )
summary(zi_MBISTot )


# Subset APOE carriers and non-carriers

dat_apoe0 <- dat %>% filter(APOEbin==0)
dat_apoe1 <- dat %>% filter(APOEbin==1)

negbin_MBISTot <- glmmTMB(mbiSTotal  ~ poly(AgeC, 2)*SCI3.3 + SEX + EDU + DropOut + (1|ParticipantId), data = dat_apoe0, family=nbinom2 )
summary(negbin_MBISTot )

negbin_MBISTot <- glmmTMB(mbiSTotal  ~ poly(AgeC, 2)*SCI3.3 + SEX + EDU + DropOut + (1|ParticipantId), data = dat_apoe1, family=nbinom2 )
summary(negbin_MBISTot )

# List of outcome variables

outcomes <- c("mbiSTotal", "mbiSTotD1", "mbiSTotD2", "mbiSTotD3", "mbiSTotD4", "mbiSTotD5")

# Function to fit model and extract tidy summary

fit_and_summarise <- function(outcome_var) {
  formula <- as.formula(paste0(outcome_var, " ~ poly(AgeC, 2)*SCI3.3 + SEX + EDU + DropOut + APOEbin + (1 | ParticipantId)"))
  
  model <- glmmTMB(formula, data = dat, family = nbinom2)
  
  tidy(model, effects = "fixed", conf.int = TRUE, exponentiate = TRUE) %>%
    mutate(outcome = outcome_var)
}

# Apply function to each outcome and bind into a single dataframe
results_negbin <- map_dfr(outcomes, fit_and_summarise)

# Save results

write.csv(results_negbin, "./Outputs/MBISTotal_WholeCohort.csv")




