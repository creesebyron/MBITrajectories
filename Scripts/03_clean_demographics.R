library(tidyverse)


# Load data

demogs <- read.csv("./Data/PROTECT Source Data/PID TEC PROT_UKDA_090 Demographics January 2024 v2.0 QR.csv")

# Demographics baseline data only needed (these will be merged in with Index MBI data from wave 1, non-matches will be dropped as they will be from later waves)

demog_bl <- demogs %>%
  mutate(DateCompleted = as.Date(DateCompleted, "%d/%m/%Y")) %>%
  arrange(ParticipantId , DateCompleted) %>%
  group_by(ParticipantId ) %>%
  slice(1) %>%  # Keep the earliest baseline observation
  ungroup()

# recode sex, ethnicity, employment status and education variables to factors.  Collapse EDU into three level factor

# Sex

demog_bl <- demog_bl %>%
  mutate(SEX = factor(SEX, levels = c(1, 2), labels = c("Male", "Female")),
         SEX = relevel(SEX, ref = "Female"))

# Education

demog_bl <- demog_bl %>%
  mutate(
    EDU = factor(EDU, 
                 levels = 1:6,
                 labels = c("Secondary", "College", "Vocational", "Degree", "Masters", "PhD")),
    EDU = fct_collapse(EDU,
                       "Post-secondary" = c("Secondary", "College", "Vocational"),
                       "PostGraduate" = c("Masters", "PhD")),
    EDU = fct_relevel(EDU, 
                      "PostGraduate", "Degree", "Post-secondary" )
  )

# Ethnicity

demog_bl <- demog_bl %>%
  mutate(
    ETH = factor(ETH, 
                 levels = 1:19,
                 labels = c("White: English / Welsh / Scottish / Northern Irish / British",
                            "White: Irish",
                            "White: Gypsy or Irish Traveller",
                            "White: European",
                            "White: Non-European",
                            "Mixed: White and Black Caribbean",
                            "Mixed: White and Black African",
                            "Mixed: White and Asian",
                            "Mixed: Any other Mixed / Multiple ethnic background",
                            "Asian / Asian British: Indian",
                            "Asian / Asian British: Pakistani",
                            "Asian / Asian British: Bangladeshi",
                            "Asian / Asian British: Chinese",
                            "Asian / Asian British: Any other Asian background",
                            "Black / African / Caribbean / Black British: African",
                            "Black / African / Caribbean / Black British: Caribbean",
                            "Any other Black / African / Caribbean background",
                            "Other ethnic group: Arab",
                            "Any other ethnic group")),
    ETH = fct_collapse(ETH,
                       "White" = c("White: English / Welsh / Scottish / Northern Irish / British", 
                                   "White: Irish", 
                                   "White: Gypsy or Irish Traveller",
                                   "White: European",
                                   "White: Non-European"),
                       "Mixed" = c("Mixed: White and Black Caribbean", 
                                   "Mixed: White and Black African",
                                   "Mixed: White and Asian",
                                   "Mixed: Any other Mixed / Multiple ethnic background"),
                       "Asian / Asian British" = c("Asian / Asian British: Indian",
                                                   "Asian / Asian British: Pakistani",
                                                   "Asian / Asian British: Bangladeshi",
                                                   "Asian / Asian British: Chinese",
                                                   "Asian / Asian British: Any other Asian background"),
                       "Black / African / Caribbean / Black British" = c("Black / African / Caribbean / Black British: African",
                                                                         "Black / African / Caribbean / Black British: Caribbean",
                                                                         "Any other Black / African / Caribbean background"),
                       "Other" = c("Other ethnic group: Arab", "Any other ethnic group")
  ))


# create an approx DOB variable

demog_bl <- as.data.frame(demog_bl %>% group_by(ParticipantId) %>%
                                        mutate(
                                          # Compute the VisitYear from the earliest visit date
                                          VisitYear = year(DateCompleted),
                                          DOBapprox = as.Date(paste0(VisitYear - AGE , "-06-01"))))  # Approximate birthdate as June 1st of the year
                                          

# Save (only cols needed for analysis, reordered to logical order)

colnames(demog_bl) 

write_rds(demog_bl[,c(1,15,6,7,9,11)], "./Data/Interim/demogs_clean.rds")

