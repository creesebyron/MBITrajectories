library(tidyverse)
library(glmmTMB)
library(DHARMa)

# Load data

mbi <- read_rds("./Data/Analytical Cohort/AnalyticalData_MBIProxy.rds")


#### MBI PROXY #####
##################

# Define common formula parts 
covars <- "SEX + EDU + PHQCat + BMIcat + DropOut"
random_effect <- "(1 | ParticipantId)"
interaction_base <- function(order, sci_var = "SCI3.3") {
  sprintf("poly(AgeC, %s)*%s", order, sci_var)
}

#  Create formulas 
make_formula <- function(order = 2, outcome = "MBIITotal", sci_var = "SCI3.3" , apoe = TRUE) {
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
fit_model <- function(formula, mbi, ziformula = ~1, control = glmmTMBControl()) {
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

# Fit models, start with cubic zi model as per self analysis

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

# Quadratic model ZI model is the best fit for the data

# Interaction model 
form_interaction <- as.formula(paste(
  "MBIITotal ~ poly(AgeC, 2)*SCI3.3*APOEbin +", covars, "+", random_effect
))
ziform_interaction <- as.formula(paste(
  "~ poly(AgeC, 2)*SCI3.3*APOEbin +", covars
))

results_int <- fit_model(form_interaction, mbi, ziformula = ziform_interaction)

# Compare the interaction model with the no interaction model 

AIC(results_int$model, results_zi_quad_full$model)

# There is no statistical interaction for proxy so just present the full model, fit indices also suggest this one is the better fit

# Run without three way interaction and with two, two-way interactions (SCI and APOE)

form_interaction_2 <- as.formula(paste(
  "MBIITotal ~ poly(AgeC, 2)*SCI3.3 + poly(AgeC, 2)*APOEbin +", covars, "+", random_effect))
ziform_interaction_2 <- as.formula(paste(
  "~ poly(AgeC, 2)*SCI3.3 + poly(AgeC, 2)*APOEbin +", covars
))

results_zi_quad_full_2twoway <- fit_model(form_interaction_2, mbi, ziformula = ziform_interaction_2)


# Compare the model with the SCI interaction only with the two-way interaction with SCI and APOE

AIC(results_zi_quad_full$model, results_int2$model)
BIC(results_zi_quad_full$model, results_int2$model)
anova(results_zi_quad_full$model, results_int2$model)

# There is no additional benefit to the model of adding APOE*Time interaction

write_list <- list(
  results_nb_quad = results_nb_quad,
  results_nb_cub = results_nb_cub,
  results_zi_cub_min = results_zi_cub_min,
  results_zi_cub_full = results_zi_cub_full,
  results_zi_quad_full = results_zi_quad_full,
  results_zi_quad_full_2twoway = results_zi_quad_full_2twoway
)

lapply(names(write_list), function(name) {
  saveRDS(write_list[[name]], file = paste0("./Outputs/Proxy/Main Analysis/", name, "_proxy.rds"))
})

saveRDS()

##### EXPLORATORY ANALYSIS #####
###############################

# Use the final model 

# AGE >=65

mbi65 <- mbi %>% filter(AgeBL >=65)

# Re-centre Age

grand_mean_age <- mean(mbi65$AgeVisit)

mbi65$AgeC <- mbi65$AgeVisit - grand_mean_age

hist(mbi65$AgeC)


results_zi_cub_full_65 <- fit_model(make_formula(order = 3, apoe = TRUE), mbi65, ziformula =  make_ziformula(order = 3, apoe = TRUE))
results_zi_quad_full_65 <- fit_model(make_formula(order = 2, apoe = TRUE), mbi65, ziformula =  make_ziformula(order = 2, apoe = TRUE))
results_zi_lin_full_65 <- fit_model(make_formula(order = 1, apoe = TRUE), mbi65, ziformula =  make_ziformula(order = 1, apoe = TRUE))

# Compare models 
AIC(results_zi_cub_full_65$model, results_zi_quad_full_65$model,results_zi_lin_full_65$model)
BIC(results_zi_cub_full_65$model, results_zi_quad_full_65$model,results_zi_lin_full_65$model)
anova(results_zi_cub_full_65$model, results_zi_quad_full_65$model,results_zi_lin_full_65$model)

# Linear model is the best fit

# Interaction term 

form_interaction_mbi65 <- as.formula(paste(
  "MBIITotal  ~ poly(AgeC, 1)*SCI3.6*APOEbin +", covars, "+", random_effect
))


results_int_mbi65 <- fit_model(form_interaction_mbi65, mbi65, ziformula = make_ziformula(order = 1,sci_var = "SCI3.6", apoe = TRUE))

# Compare models 
AIC(results_int_mbi65$model,  results_zi_lin_full_65$model)
BIC(results_int_mbi65$model,  results_zi_lin_full_65$model)
anova(results_int_mbi65$model,results_zi_lin_full_65$model)

# Interaction does not improve fit.  Use full model linear term no interaction 

# Write out results

write_list <- list(
  results_zi_lin_full_65 = results_zi_lin_full_65,
  results_zi_quad_full_65 = results_zi_quad_full_65,
  results_zi_cub_full_65 = results_zi_cub_full_65,
  results_int_mbi65 = results_int_mbi65
  )

lapply(names(write_list), function(name) {
  saveRDS(write_list[[name]], file = paste0("./Outputs/Proxy/Age 65/", name, "_proxy.rds"))
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

# Use quadratic term model

# Interaction model 
form_interaction_sci3.6 <- as.formula(paste(
  "MBIITotal  ~ poly(AgeC, 2)*SCI3.6*APOEbin +", covars, "+", random_effect
))


results_int_sci3.6 <- fit_model(form_interaction_sci3.6, mbi, ziformula = make_ziformula(order = 2,sci_var = "SCI3.6", apoe = TRUE))

# Compare models. 
AIC(results_int_sci3.6$model, results_zi_quad_full_sci3.6$model)
BIC(results_int_sci3.6$model, results_zi_quad_full_sci3.6$model)
anova(results_int_sci3.6$model, results_zi_quad_full_sci3.6$model)

# There is no benefit to the interaction term.  Use the full quadratic model no interaction term


# Write out results. 

write_list <- list(
  results_zi_cub_full_sci3.6 = results_zi_cub_full_sci3.6,
  results_zi_quad_full_sci3.6 = results_zi_quad_full_sci3.6,
  results_zi_lin_full_sci3.6 = results_zi_lin_full_sci3.6,
  results_int_sci3.6 = results_int_sci3.6
)

lapply(names(write_list), function(name) {
  saveRDS(write_list[[name]], file = paste0("./Outputs/Proxy/SCI 3.6/", name, "_proxy.rds"))
})


# MBI DOMAINS, set zi formula to ~1 so that models run for all domains

# MBI DOMAINS

domain_list <- paste0("MBIITotD", 1:5)

covars_domains <- "APOEbin + SEX + EDU + PHQCat + BMIcat + DropOut"

# Mirror the main analysis - i.e., just model main effect of APOE, no trajectory 

results_domains_proxy_full <- lapply(domain_list, function(outcome_var) {
  # Build formulas for each outcome
  form_interaction <- as.formula(paste(
    outcome_var, "~ poly(AgeC, 2)*SCI3.3 +", covars_domains, "+", random_effect
  ))
  ziform_interaction <- as.formula(paste(
    "~ poly(AgeC, 2)*SCI3.3 +", covars_domains
  ))
  # Fit model and return results
  fit_model(formula = form_interaction, mbi = mbi, ziformula = ziform_interaction)
})    


# models did not converge, so simplify by setting ZI term to constant

results_domains_proxy_min <- lapply(domain_list, function(outcome_var) {
  # Build formulas for each outcome
  form_interaction <- as.formula(paste(
    outcome_var, "~ poly(AgeC, 2)*SCI3.3 +", covars_domains, "+", random_effect
  ))
  # Fit model and return results
  fit_model(formula = form_interaction, mbi = mbi, ziformula = ~1)
})   

# Name the list by outcome
names(results_domains_proxy_min) <- domain_list
names(results_domains_proxy_full) <- domain_list


# Write out results. Use interaction term model 

write_list <- list(
  results_domains_proxy_min = results_domains_proxy_min,
  results_domains_proxy_full = results_domains_proxy_full
  
)

lapply(names(write_list), function(name) {
  saveRDS(write_list[[name]], file = paste0("./Outputs/Proxy/Domains/", name, "_proxy.rds"))
})


