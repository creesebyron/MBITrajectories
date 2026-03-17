library(tidyverse)
library(DescTools)

# Load data

cog <- read.csv("./Data/PROTECT Source Data/PID TEC PROT_UK_DA_090 CTP1 FINAL V1.0.csv")

# Keep only relevant variables 

colnames(cog)

cog <- cog %>%
   select(1,2,4,6,7,8,11,14,17)     

head(cog)

# Check class

lapply(cog, class)
lapply(cog, summary)

# Convert columns 6 to 9 to numeric (this will change "NULL" to NA)

cog <- cog %>%
  mutate(across(6:9, as.numeric))

lapply(cog, class)
lapply(cog, summary)

# Convert 'DateCompleted' to Date format 
cog <- cog %>%
  mutate(DateCompleted = as.Date(DateCompleted, "%d/%m/%Y"))  

lapply(cog, class)

# New variable which is the count how many repeats were done at each time point

cog <- cog %>%
            group_by(ParticipantId, Timepoint) %>%
            mutate(NumRepeats = max(RepeatNumber, na.rm = TRUE)) %>%
            ungroup()

# take the mean of the three repeats, ignoring any missing values

cog <- cog %>% arrange(ParticipantId, Timepoint, RepeatNumber) %>%
            group_by(ParticipantId, Timepoint) %>%
              mutate(VERBmean = mean(VERBTOT, na.rm = TRUE),
                     SPANmean = mean(SPANTOT, na.rm = TRUE),
                     PALmean = mean(PALTOT, na.rm = TRUE),
                     SEARCHmean = mean(SEARCHTOT, na.rm = TRUE)) %>%
                     ungroup()

head(cog)

# keep only repeat number 1 (now that the mean of all three has been created, only one observation per time point per person should be present)

cog <- as.data.frame(cog %>%
                       filter(RepeatNumber==1))

head(cog)

# covert to z scores

cog$SPANTOT_z  <- scale(cog$SPANmean , center=TRUE, scale = TRUE)
cog$PALTOT_z   <- scale(cog$PALmean , center=TRUE, scale = TRUE)
cog$VERBTOT_z  <- scale(cog$VERBmean , center=TRUE, scale = TRUE)
cog$SEARCHTOT_z <- scale(cog$SEARCHmean, center=TRUE, scale = TRUE)

# Evaluate outliers at each time point

ggplot(cog, aes(x = VERBmean)) +
  geom_histogram(bins = 30, fill = "steelblue", color = "white") +
  facet_wrap(~Timepoint, scales = "free") +
  theme_minimal() +
  labs(title = "Histograms by and Timepoint",
       x = "Value", y = "Count")

ggplot(cog, aes(x = PALmean)) +
  geom_histogram(bins = 30, fill = "steelblue", color = "white") +
  facet_wrap(~Timepoint, scales = "free") +
  theme_minimal() +
  labs(title = "Histograms by and Timepoint",
       x = "Value", y = "Count")

ggplot(cog, aes(x = SPANmean)) +
  geom_histogram(bins = 30, fill = "steelblue", color = "white") +
  facet_wrap(~Timepoint, scales = "free") +
  theme_minimal() +
  labs(title = "Histograms by and Timepoint",
       x = "Value", y = "Count")

ggplot(cog, aes(x = SEARCHmean)) +
  geom_histogram(bins = 30, fill = "steelblue", color = "white") +
  facet_wrap(~Timepoint, scales = "free") +
  theme_minimal() +
  labs(title = "Histograms by and Timepoint",
       x = "Value", y = "Count")

# Store row count before filtering

nrow_before <- nrow(cog)

# Filter to keep only rows with no missing values in selected columns
cog <- cog %>%
        filter(if_all(VERBmean:SEARCHmean, ~ !is.na(.x)))

# Store row count after filtering
nrow_after <- nrow(cog)

# Calculate how many were dropped
rows_dropped <- nrow_before - nrow_after # removes ~5% of observations

# G factor (1st unrotated PC of 4 cognitive tests)

# subset to leave only cognitive scores at each time point

colnames(cog)

gfactor <- cog %>%
              select(15:18)

g1pca <- princomp(gfactor, cor = TRUE)

summary(g1pca)
loadings(g1pca)

cog <- cbind(cog, g1pca$scores[,1])

#rename pca columns

colnames(cog)[colnames(cog)=="g1pca$scores[, 1]"] <- "G"

head(cog)
nrow(cog)

# Check distribution

qqnorm(cog$G)
qqline(cog$G, col = "red", lwd = 2)

# Save (keeping only relevant variables)

write_rds(cog[, c(1,2,3,10,19)], "./Data/Interim/cognition_clean.rds")






