library(tidyverse)
library(glmmTMB)
library(DHARMa)

# Load data

phq <- read_rds("./Data/Analytical Cohort/AnalyticalData_PHQ.rds")


#### PHQ9 #####
##################

# Define common formula parts 
covars <- "SEX + EDU + BMIcat + DropOut"
random_effect <- "(1 | ParticipantId)"
interaction_base <- function(order, sci_var = "SCI3.3") {
  sprintf("poly(AgeC, %s)*%s", order, sci_var)
}


#  Create formulas 
make_formula <- function(order = 2, sci_var = "SCI3.3", outcome = "PHQ9Tot", apoe = TRUE) {
  interaction_term <- interaction_base(order, sci_var)
  rhs <- paste(interaction_term, covars, if (apoe) "APOEbin", random_effect, sep = " + ")
  as.formula(paste(outcome, "~", rhs))
}

make_ziformula <- function(order = 2, sci_var = "SCI3.3", apoe = TRUE) {
  zi_terms <- c(interaction_base(order, sci_var), covars)
  if (apoe) zi_terms <- c(zi_terms, "APOEbin")
  as.formula(paste("~", paste(zi_terms, collapse = " + ")))
}

# Fit and evaluate models
fit_model <- function(formula, phq, ziformula = ~1, control= glmmTMBControl()) {
  model <- glmmTMB(formula, data = phq, family = nbinom2, ziformula = ziformula)
  sim <- simulateResiduals(model)
  list(
    model = model,
    summary = summary(model),
    residuals = sim,
    outliers = testOutliers(sim, type = "bootstrap"),
    zi_test = testZeroInflation(sim),
    dispersion = testDispersion(sim)
  )
}

# Fit models 
form_linear <- make_formula(order = 1, sci_var = "SCI3.3")
form_quad <- make_formula(order = 2, sci_var = "SCI3.3")
form_cub <- make_formula(order = 3, sci_var = "SCI3.3")
ziform_linear <- make_ziformula(order = 1, sci_var = "SCI3.3")
ziform_quad <- make_ziformula(order = 2, sci_var = "SCI3.3")
ziform_cub <- make_ziformula(order = 3, sci_var = "SCI3.3")


# Fit reduced model with linear term 
results_zi_lin_min <- fit_model(form_linear, phq, ziformula = ~1) # convergence failure 

# Fit full model with linear term 
results_zi_lin_full <- fit_model(form_linear, phq, ziformula = ziform_linear)

# Fit reduced quadratic ZI model 
results_zi_quad_min <- fit_model(form_quad, phq, ziformula = ~1)

# Fit full quadratic ZI model 
results_zi_quad_full <- fit_model(form_quad, phq, ziformula = ziform_quad) #  convergence failure 

# Fit reduced ZI model
results_zi_cub_min <- fit_model(form_cub, phq, ziformula = ~1)

# Fit cubic ZI model
results_zi_cub_full <- fit_model(form_cub, phq, ziformula = ziform_cub) # convergence failure 


# Compare models 
AIC(results_zi_lin_min$model, results_zi_lin_full$model, results_zi_quad_min$model, results_zi_quad_full$model, results_zi_cub_min$model)
BIC(results_zi_lin_min$model, results_zi_lin_full$model, results_zi_quad_min$model, results_zi_quad_full$model, results_zi_cub_min$model)
anova(results_zi_lin_min$model, results_zi_lin_full$model, results_zi_quad_min$model, results_zi_quad_full$model, results_zi_cub_min$model)

# Use quadratic model without the ZI

# Interaction model 
form_interaction <- as.formula(paste(
  "mbiSTotal ~ poly(AgeC, 2)*SCI3.3*APOEbin +", covars, "+", random_effect
))

results_int <- fit_model(form_interaction, phq)

# Compare the interaction model with the no interaction model 

AIC(results_int$model, results_zi_quad_min$model)
BIC(results_int$model, results_zi_quad_min$model)
anova(results_int$model, results_zi_quad_min$model)

# Use the interaction term model for plotting 

write_list <- list(
  results_zi_lin_min = results_zi_lin_min,
  results_zi_lin_full = results_zi_lin_full,
  results_zi_quad_min = results_zi_quad_min,
  results_zi_cub_full = results_zi_cub_full,
  results_zi_quad_full = results_zi_quad_full,
  results_zi_cub_min = results_zi_cub_min,
  results_int = results_int
)

lapply(names(write_list), function(name) {
  saveRDS(write_list[[name]], file = paste0("./Outputs/Self/Main analysis/", name, "_phq.rds"))
})



