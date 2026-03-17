library(tidyverse)
library(glmmTMB)
library(DHARMa)

# Load data

mbi <- read_rds("./Data/Analytical Cohort/AnalyticalData_MBISelf.rds")


#### MBI SELF #####
##################

# Define common formula parts 
covars <- "SEX + EDU + PHQCat + BMIcat + DropOut"
covarsnoPHQ <- "SEX + EDU + BMIcat + DropOut"
random_effect <- "(1 | ParticipantId)"
interaction_base <- function(order, sci_var = "SCI3.3") {
  sprintf("poly(AgeC, %s)*%s", order, sci_var)
}


#  Create formulas 
make_formula <- function(order = 2, sci_var = "SCI3.3", outcome = "mbiSTotal", apoe = TRUE) {
  interaction_term <- interaction_base(order, sci_var)
  rhs <- paste(interaction_term, covarsnoPHQ, if (apoe) "APOEbin", random_effect, sep = " + ")
  as.formula(paste(outcome, "~", rhs))
}

make_ziformula <- function(order = 2, sci_var = "SCI3.3", apoe = TRUE) {
  zi_terms <- c(interaction_base(order, sci_var), covars)
  if (apoe) zi_terms <- c(zi_terms, "APOEbin")
  as.formula(paste("~", paste(zi_terms, collapse = " + ")))
}

# Fit and evaluate models
fit_model <- function(formula, mbi, ziformula = ~1, control= glmmTMBControl()) {
  model <- glmmTMB(formula, data = mbi, family = nbinom2, ziformula = ziformula)
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
form_quad <- make_formula(order = 2, sci_var = "SCI3.3")
form_cub <- make_formula(order = 3, sci_var = "SCI3.3")
ziform_cub <- make_ziformula(order = 3, sci_var = "SCI3.3")
ziform_quad <- make_ziformula(order = 2, sci_var = "SCI3.3")

# Fit NB model
results_nb_quad <- fit_model(form_quad, mbi)
# Fit ZI model
results_nb_cub <- fit_model(form_cub, mbi)
# Fit quad ZI model
results_zi_quad_full <- fit_model(form_quad, mbi, ziformula = ziform_quad)
# Fit reduced ZI model
results_zi_cub_min <- fit_model(form_cub, mbi, ziformula = ~1)
# Fit cubic ZI model
results_zi_cub_full <- fit_model(form_cub, mbi, ziformula = ziform_cub)


# Compare models 
AIC(results_nb_quad$model, results_nb_cub$model, results_zi_cub_min$model, results_zi_cub_full$model, results_zi_quad_full$model)
BIC(results_nb_quad$model, results_nb_cub$model, results_zi_cub_min$model, results_zi_cub_full$model, results_zi_quad_full$model)
anova(results_nb_quad$model, results_nb_cub$model, results_zi_cub_min$model, results_zi_cub_full$model, results_zi_quad_full$model)

# Use cubic term model

# Interaction model 
form_interaction <- as.formula(paste(
  "mbiSTotal ~ poly(AgeC, 3)*SCI3.3*APOEbin +", covars, "+", random_effect
))
ziform_interaction <- as.formula(paste(
  "~ poly(AgeC, 3)*SCI3.3*APOEbin +", covars
))

results_int <- fit_model(form_interaction, mbi, ziformula = ziform_interaction)


# Compare the interaction model with the no interaction model 

AIC(results_int$model, results_zi_cub_full$model)

# Use the interaction term model for plotting 

write_list <- list(
  results_nb_quad = results_nb_quad,
  results_nb_cub = results_nb_cub,
  results_zi_cub_min = results_zi_cub_min,
  results_zi_cub_full = results_zi_cub_full,
  results_zi_quad_full = results_zi_quad_full,
  results_int = results_int
)

lapply(names(write_list), function(name) {
  saveRDS(write_list[[name]], file = paste0("./Outputs/Self/Main analysis/", name, "_self.rds"))
})


##### EXPLORATORY ANALYSIS #####
###############################

# Use the final model (i.e., zi with cubic term) from above 

# AGE >=65

mbi65 <- mbi %>% filter(AgeBL >=65)

# Re-centre Age

grand_mean_age <- mean(mbi65$AgeVisit)

mbi65$AgeC <- mbi65$AgeVisit - grand_mean_age

hist(mbi65$AgeC)

# Fit model

results_zi_cub_full_65 <- fit_model(make_formula(order = 3, apoe = TRUE), mbi65, ziformula =  make_ziformula(order = 3, apoe = TRUE))
results_zi_quad_full_65 <- fit_model(make_formula(order = 2, apoe = TRUE), mbi65, ziformula =  make_ziformula(order = 2, apoe = TRUE))
results_zi_lin_full_65 <- fit_model(make_formula(order = 1, apoe = TRUE), mbi65, ziformula =  make_ziformula(order = 1, apoe = TRUE))

# Compare models. Linear did not converge. 
AIC(results_zi_cub_full_65$model, results_zi_quad_full_65$model, results_zi_lin_full_65$model)
BIC(results_zi_cub_full_65$model, results_zi_quad_full_65$model, results_zi_lin_full_65$model)
anova(results_zi_cub_full_65$model, results_zi_quad_full_65$model, results_zi_lin_full_65$model)

#  There is no benefit to the cubic or quadratic models so use the linear

# Stratified by APOE 
mbi65_apoe_list <- list(
  apoe0 = mbi65 %>% filter(APOEbin == 0),
  apoe1 = mbi65 %>% filter(APOEbin == 1)
)

strat_models_quad_full_65 <- lapply(mbi65_apoe_list, function(df) {
  fit_model(make_formula(order = 1, apoe = FALSE), df, ziformula = make_ziformula(order = 1, apoe = FALSE))
})


# Stratification did not converge.  Use the full model with no interaction

# Write out results

write_list <- list(
  results_zi_lin_full_65 = results_zi_lin_full_65
)

lapply(names(write_list), function(name) {
  saveRDS(write_list[[name]], file = paste0("./Outputs/Self/Age 65/", name, "_self.rds"))
})

# SCI 3.6

# Fit models 

results_zi_cub_full_sci3.6 <- fit_model(make_formula(order = 3,sci_var = "SCI3.6", apoe = TRUE), mbi, ziformula =  make_ziformula(order = 3,sci_var = "SCI3.6", apoe = TRUE))
results_zi_quad_full_sci3.6 <- fit_model(make_formula(order = 2,sci_var = "SCI3.6", apoe = TRUE), mbi, ziformula =  make_ziformula(order = 2,sci_var = "SCI3.6", apoe = TRUE))
results_zi_lin_full_sci3.6 <- fit_model(make_formula(order = 1,sci_var = "SCI3.6", apoe = TRUE), mbi, ziformula =  make_ziformula(order = 1,sci_var = "SCI3.6", apoe = TRUE))

# Compare models. 
AIC(results_zi_cub_full_sci3.6$model, results_zi_quad_full_sci3.6$model,results_zi_lin_full_sci3.6$model)
BIC(results_zi_cub_full_sci3.6$model, results_zi_quad_full_sci3.6$model,results_zi_lin_full_sci3.6$model)
anova(results_zi_cub_full_sci3.6$model, results_zi_quad_full_sci3.6$model,results_zi_lin_full_sci3.6$model)

# Use cubic term model

# Stratified by APOE 
mbi_apoe_list <- list(
  apoe0 = mbi %>% filter(APOEbin == 0),
  apoe1 = mbi %>% filter(APOEbin == 1)
)

form_strat_sci3.6 <- make_formula(order = 3, sci_var = "SCI3.6", apoe = FALSE)
ziform_strat <- make_ziformula(order = 3, apoe = FALSE)

strat_models_sci3.6 <- lapply(mbi_apoe_list, function(df) {
  fit_model(form_strat_sci3.6, df, ziformula = make_ziformula(order = 3,sci_var = "SCI3.6"))
})

# Interaction model 
form_interaction_sci3.6 <- as.formula(paste(
  "mbiSTotal ~ poly(AgeC, 3)*SCI3.6*APOEbin +", covars, "+", random_effect
))

results_int_sci3.6 <- fit_model(form_interaction_sci3.6, mbi, ziformula = make_ziformula(order = 3,sci_var = "SCI3.6", apoe = TRUE))

# Interaction model fails to converge. Use the full model with no interaction term 


# Write out results. Use interaction term model 

write_list <- list(
  results_zi_cub_full_sci3.6 = results_zi_cub_full_sci3.6
)

lapply(names(write_list), function(name) {
  saveRDS(write_list[[name]], file = paste0("./Outputs/Self/SCI 3.6/", name, "_self.rds"))
})


# MBI DOMAINS

domain_list <- paste0("mbiSTotD", 1:5)


results_domains_self_full <- lapply(domain_list, function(outcome_var) {
  # Build formulas for each outcome
  form_interaction <- as.formula(paste(
    outcome_var, "~ poly(AgeC, 3)*SCI3.3*APOEbin +", covars, "+", random_effect
  ))
  ziform_interaction <- as.formula(paste(
    "~ poly(AgeC, 3)*SCI3.3*APOEbin +", covars
  ))
  # Fit model and return results
  fit_model(formula = form_interaction, mbi = mbi, ziformula = ziform_interaction)
})    


# models did not converge, so simplify by setting ZI term to constant

results_domains_self_min <- lapply(domain_list, function(outcome_var) {
  # Build formulas for each outcome
  form_interaction <- as.formula(paste(
    outcome_var, "~ poly(AgeC, 3)*SCI3.3*APOEbin +", covars, "+", random_effect
  ))
  # Fit model and return results
  fit_model(formula = form_interaction, mbi = mbi, ziformula = ~1)
})   

# Name the list by outcome
names(results_domains_self_min) <- domain_list
names(results_domains_self_full) <- domain_list


# Write out results. 

write_list <- list(
  results_domains_self_min = results_domains_self_min,
  results_domains_self_full = results_domains_self_full
  
)

lapply(names(write_list), function(name) {
  saveRDS(write_list[[name]], file = paste0("./Outputs/Self/Domains/", name, "_self.rds"))
})





