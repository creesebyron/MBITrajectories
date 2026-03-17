library(tidyverse)
library(psych)

# Read in data

mbi <- read_rds("./Data/Analytical Cohort/AnalyticalData_MBISelf.rds")
mbii <- read_rds("./Data/Analytical Cohort/AnalyticalData_MBIProxy.rds")

# Cross tabs of PHQ cat score against all mental health measures to determine covariates (leave out mental health conditions if PHQ9 is strongly associated).  Use baseline data

# Use baseline data

datbl <- mbi %>% filter(VisitNumSelf==1)
datblproxy <- mbii %>% filter(VisitNumProxy==1)

# Table of key demographics 

dflist <- c("datbl", "datblproxy")

demogstab <- list(lapply(dflist, function(x) describeBy(get(x)$AgeBL)),
             lapply(dflist, function(x) describeBy(get(x)$FollowUpMax)),
             lapply(dflist, function(x) table(get(x)$PHQCat)),
             lapply(dflist, function(x) table(get(x)$SEX)),
             lapply(dflist, function(x) table(get(x)$EDU)),
             lapply(dflist, function(x) table(get(x)$DropOut)),
             lapply(dflist, function(x) table(get(x)$BMIcat)),
             lapply(dflist, function(x) table(get(x)$EMPLOY)),
             lapply(dflist, function(x) table(get(x)$ETH)),
             lapply(dflist, function(x) table(get(x)$APOEbin))
             )

names(demogstab) <- c("Age", "Follow up duration", "PHQ9", "Sex", "Education level", "Drop Out", "BMI", "Employment status", "Ethnicity", "APOE carrier")

write_rds(demogstab, "./Outputs/Demographics_Clinical_self_proxy.rds")


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

long_datproxy <- datblproxy %>%
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

summary_tableproxy <- long_datproxy %>%
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
print(summary_tableproxy, n = Inf)

write_rds(summary_table, "./Outputs/SelfTableMentalHealthByPHQ.rds")
write_rds(summary_tableproxy, "./Outputs/ProxyTableMentalHealthByPHQ.rds")

# Start from long_dat
plot_data <- long_dat %>%
  group_by(Condition, PHQCat, Level) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(Condition, PHQCat) %>%
  mutate(pct = 100 * n / sum(n)) %>%
  ungroup()

phqplot <- ggplot(plot_data, aes(x = PHQCat, y = pct, fill = Level)) +
  geom_col(position = "fill") +
  facet_wrap(~ Condition, scales = "free_y") +
  labs(
    x = "PHQ-9 Category",
    y = "Proportion",
    fill = "Condition Present?",
    title = "Distribution of Mental Health Conditions Across PHQ-9 Categories"
  ) +
  theme_minimal() +
  scale_y_continuous(labels = scales::percent_format())

ggsave("./Plots/phqplot_a5_highres.tiff", plot = phqplot, width = 8.27, height = 5.83, dpi = 600, units = "in", compression = "lzw", bg = "white")

# Start from long_dat
plot_dataproxy <- long_datproxy %>%
  group_by(Condition, PHQCat, Level) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(Condition, PHQCat) %>%
  mutate(pct = 100 * n / sum(n)) %>%
  ungroup()

phqplotproxy <- ggplot(plot_dataproxy, aes(x = PHQCat, y = pct, fill = Level)) +
  geom_col(position = "fill") +
  facet_wrap(~ Condition, scales = "free_y") +
  labs(
    x = "PHQ-9 Category",
    y = "Proportion",
    fill = "Condition Present?",
    title = "Distribution of Mental Health Conditions Across PHQ-9 Categories"
  ) +
  theme_minimal() +
  scale_y_continuous(labels = scales::percent_format())

ggsave("./Plots/phqplotproxy_a5_highres.tiff", plot = phqplotproxy, width = 8.27, height = 5.83, dpi = 600, units = "in", compression = "lzw", bg = "white")


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

write_rds(summary_tablemh, "./Outputs/SelfTableMentalHealthMBIS.rds")


mhplot <- ggplot(summary_tablemh, aes(x = Level, y = Mean_MBIStotal, fill = Level)) +
  geom_col() +
  facet_wrap(~Condition) +
  theme_minimal()

ggsave("./Plots/mhplot_a5_highres.tiff", plot = mhplot, width = 8.27, height = 5.83, dpi = 600, units = "in", compression = "lzw", bg = "white")

# Pivot to long format and summarise mean/SD of MBIITotal
summary_tablemhproxy <- datblproxy %>%
  select(MBIITotal, all_of(mh_vars)) %>%
  pivot_longer(cols = all_of(mh_vars), names_to = "Condition", values_to = "Level") %>%
  filter(!is.na(Level), !is.na(MBIITotal)) %>%
  group_by(Condition, Level) %>%
  summarise(
    Mean_MBIITotal = round(mean(MBIITotal, na.rm = TRUE), 2),
    SD_MBIITotal   = round(sd(MBIITotal, na.rm = TRUE), 2),
    .groups = "drop"
  )

# View the result
print(summary_tablemhproxy, n = Inf)


write_rds(summary_tablemhproxy, "./Outputs/ProxyTableMentalHealthMBII.rds")

mhplotproxy <- ggplot(summary_tablemhproxy, aes(x = Level, y = Mean_MBIITotal, fill = Level)) +
  geom_col() +
  facet_wrap(~Condition) +
  theme_minimal()

ggsave("./Plots/mhplotproxy_a5_highres.tiff", plot = mhplotproxy, width = 8.27, height = 5.83, dpi = 600, units = "in", compression = "lzw", bg = "white")


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

summary_tablemedproxy <- datblproxy %>%
  select(MBIITotal, all_of(med_vars)) %>%
  pivot_longer(cols = all_of(med_vars), names_to = "Condition", values_to = "Level") %>%
  filter(!is.na(Level), !is.na(MBIITotal)) %>%
  group_by(Condition, Level) %>%
  summarise(
    Mean_MBIITotal = round(mean(MBIITotal, na.rm = TRUE), 2),
    SD_MBIITotal   = round(sd(MBIITotal, na.rm = TRUE), 2),
    .groups = "drop"
  )


# View the result
print(summary_tablemed, n = Inf)
print(summary_tablemedproxy, n = Inf)

write_rds(summary_tablemed, "./Outputs/SelfTableMedicalMBIS.rds")
write_rds(summary_tablemedproxy, "./Outputs/ProxyTableMedicalMBII.rds")


