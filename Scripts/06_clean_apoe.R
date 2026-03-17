library(tidyverse)

# Load data

apoe <- read.csv("./Data/PROTECT Source Data/PROTECT_APOE4.csv")

table(apoe$Genotype)

#recode apoes 

apoe <- apoe %>%
  mutate(
    APOE = case_when(
      Genotype == "CC" ~ 2,         # Recode 'CC' to 2
      Genotype %in% c("CT", "TC") ~ 1,  # Recode 'CT' and 'TC' to 1
      Genotype == "TT" ~ 0,         # Recode 'TT' to 0
      TRUE ~ NA_real_                # Handle any unexpected values, set as NA
    )
  )

table(apoe$APOE)

#convert APOE to factor

apoe$APOE <- as.factor(apoe$APOE)

#create APOEbin variable (carrier of at least one e4 allele)

apoe$APOEbin <- ifelse(apoe$APOE==0,0,1)

table(apoe$APOEbin,apoe$APOE)

# Save (only cols needed for analysis)

colnames(apoe) 

write_rds(apoe[,c(1,8:9)], "./Data/Interim/apoe_clean.rds")




