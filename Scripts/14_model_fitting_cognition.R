library(tidyverse)
library(glmmTMB)
library(DHARMa)

# Load data

cogn <- readRDS("./Data/Analytical Cohort/AnalyticalData_Cognitive.rds")


# Define formula parts 

covars <- "SEX + EDU +  PHQCat + BMIcat + DropOut + NumRepeats"
random_effect <- "(1 | ParticipantId)"
interaction_base <- "poly(AgeC, %s)*SCI3.3 + poly(AgeC, %s)*APOEbin"

#  Create formulas 
make_formula <- function(order = 2, outcome = "G") {
  interaction_term <- sprintf(interaction_base, order, order)
  rhs <- paste(interaction_term, covars, random_effect, sep = " + ")
  as.formula(paste(outcome, "~", rhs))
}

# Fit and evaluate models
fit_model <- function(formula, data) {
  model <- glmmTMB(formula, data, family = gaussian)
  sim <- simulateResiduals(model)
  list(
    model = model,
    summary = summary(model),
    residuals = sim
  )
}

# Fit  model, quadratic
results_g_quad <- fit_model(make_formula(order = 2), cogn)
# Fit  model, cubic
results_g_cub <- fit_model(make_formula(order = 3), cogn)


# Compare models 
AIC(results_g_cub$model, results_g_quad$model)
BIC(results_g_cub$model, results_g_quad$model)
anova(results_g_cub$model, results_g_quad$model)

# Interaction term

form_interaction_quad <- as.formula(paste(
  "G ~ poly(AgeC, 2)*SCI3.3*APOEbin +", covars, "+", random_effect
))

results_int_cog <- fit_model(form_interaction_quad, cogn)

anova(results_g_quad$model, results_int_cog$model)

# Three way interaction term does not improve model fit, stick with the quadratic term model

# Explore removing outliers +/- 3SD

describeBy(cogn$G)

threshold <- 3*sd(cogn$G)

cogn3sd <- cogn %>% filter(abs(G) < threshold)
# number dropped:
length(unique(cogn3sd$ParticipantId)) - length(unique(cogn$ParticipantId))

qqnorm(cogn3sd$G)
qqline(cogn3sd$G)

# Fit NB model, quadratic
results_g_quad3sd <- fit_model(make_formula(order = 2), cogn3sd)

results_g_quad3sd$summary
results_g_quad$summary

plot(results_g_quad3sd$residuals)

# No major difference after removing outliers, but does improve the distribution

# write out results

write_list <- list(
  results_g_cub = results_g_cub,
  results_g_quad = results_g_quad
)

lapply(names(write_list), function(name) {
  saveRDS(write_list[[name]], file = paste0("./Outputs/Cognitive/", name, ".rds"))
})

