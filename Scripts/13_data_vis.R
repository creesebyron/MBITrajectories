library(ggplot2)
library(tidyverse)

# Load model output and data

mod_results_self <- readRDS("./Outputs/Self/Main analysis/results_int_self.rds")
mod_results_proxy <- readRDS("./Outputs/Proxy/Main Analysis/results_zi_quad_full_proxy.rds")
mod_results_cog <- readRDS("./Outputs/Cognitive/results_g_quad.rds")
mod_results_domains_self <- readRDS("./Outputs/Self/Domains/results_domains_self_min_self.rds")
mod_results_domains_proxy <- readRDS("./Outputs/Proxy/Domains/results_domains_proxy_min_proxy.rds")
mod_results_65_self <- readRDS("./Outputs/Self/Age 65/results_zi_lin_full_65_self.rds")
mod_results_65_proxy <-readRDS("./Outputs/Proxy/Age 65/results_zi_lin_full_65_proxy.rds")
mod_results_sci3.6_self <- readRDS("./Outputs/Self/SCI 3.6/results_zi_cub_full_sci3.6_self.rds")
mod_results_sci3.6_proxy <-readRDS("./Outputs/Proxy/SCI 3.6/results_zi_quad_full_sci3.6_proxy.rds")
mod_results_phq <- readRDS("./Outputs/Self/Main analysis/results_int_phq.rds")
mbi <- readRDS("./Data/Analytical Cohort/AnalyticalData_MBISelf.rds")
mbii <- readRDS("./Data/Analytical Cohort/AnalyticalData_MBIProxy.rds")
cog <- readRDS("./Data/Analytical Cohort/AnalyticalData_Cognitive.rds")
phq <- readRDS("./Data/Analytical Cohort/AnalyticalData_PHQ.rds")


#### GET DESCRIPTIVES ####
##########################


##### SELF #####
################

# Visualise the data

# Get modal levels
get_mode <- function(x) {
  names(which.max(table(x)))
}

summary_stats_mbi <- mbi %>%
  summarize(
    SEX    = get_mode(SEX),
    EDU    = get_mode(EDU),
    BMIcat = get_mode(BMIcat),
    PHQCat = mean(PHQCat),
    DropOut = mean(DropOut, na.rm = TRUE)
  )



newdata <- expand.grid(
  AgeC = seq(-15, 15, length.out = 100),
  SCI3.3 = c(0, 1),
  APOEbin = c(0, 1),
  SEX = summary_stats_mbi$SEX,
  EDU = summary_stats_mbi$EDU,
  PHQCat = summary_stats_mbi$PHQCat,
  BMIcat = summary_stats_mbi$BMIcat,
  DropOut = summary_stats_mbi$DropOut,
  ParticipantId = NA
)

# Add raw age

newdata$Age <- 65 + newdata$AgeC 


# Plot negative binomial

pred_out <- predict(mod_results_self$model,
                    newdata = newdata,
                    type = "response",
                    se.fit = TRUE,
                    allow.new.levels = TRUE)

# Add confidence intervals

newdata$predicted <- pred_out$fit
newdata$se <- pred_out$se.fit
newdata$LCI <- newdata$predicted - 1.96 * newdata$se
newdata$UCI <- newdata$predicted + 1.96 * newdata$se

newdata <- newdata %>%
  mutate(SCI3.3 = factor(SCI3.3,
                          levels = 0:1,
                          labels = c("No SCI", "SCI")))

mbis_plot <- ggplot(newdata, aes(x = Age, y = predicted,
                                 color = factor(APOEbin))) +
  geom_line(linewidth = 1) +
  geom_ribbon(aes(ymin = LCI, ymax = UCI, fill = factor(APOEbin), group = interaction(APOEbin, SCI3.3)),
              alpha = 0.2, color = NA) +
  labs(x = "Age (Years)",
       y = "Predicted (Expected) MBI-C Total Score (Self Rated)",
       color = "APOE Status",
       fill = "APOE Status",
       linetype = "SCI Group") +
  theme_minimal() +
  scale_color_manual(
    values = c("gray30", "steelblue"),
    labels = c("APOE ε4 non-carrier", "APOE ε4 carrier")
  ) +
  scale_fill_manual(
    values = c("gray30", "steelblue"),
    labels = c("APOE ε4 non-carrier", "APOE ε4 carrier")
  ) +
  facet_grid(. ~ SCI3.3    ) +  # facets side by side (1 row)
  theme(legend.position = "right")

# Plot zero-inflated probabilities

newdata_zi <- expand.grid(
  AgeC = seq(-15, 15, length.out = 100),
  SCI3.3 = c(0, 1),
  APOEbin = c(0, 1),
  SEX = summary_stats_mbi$SEX,
  EDU = summary_stats_mbi$EDU,
  PHQCat = summary_stats_mbi$PHQCat,
  BMIcat = summary_stats_mbi$BMIcat,
  DropOut = summary_stats_mbi$DropOut,
  ParticipantId = NA
)

# Add raw age

newdata_zi$Age <- 65 + newdata_zi$AgeC 

pred_out_zi <- predict(mod_results_self$model,
                    newdata = newdata_zi,
                    type = "zprob",
                    se.fit = TRUE,
                    allow.new.levels = TRUE)

# Add confidence intervals

newdata_zi$zprob <- pred_out_zi$fit
newdata_zi$se <- pred_out_zi$se.fit
newdata_zi$LCI <- plogis(qlogis(newdata_zi$zprob)) - 1.96 * newdata_zi$se
newdata_zi$UCI <- plogis(qlogis(newdata_zi$zprob)) + 1.96 * newdata_zi$se

# make SCI factor

newdata_zi <- newdata_zi %>%
  mutate(SCI3.3 = factor(SCI3.3,
                         levels = 0:1,
                         labels = c("No SCI", "SCI")))

mbis_plot_zi <- ggplot(newdata_zi, aes(x = Age, y = zprob,
                                 color = factor(APOEbin))) +
  geom_line(linewidth = 1) +
  geom_ribbon(aes(ymin = LCI, ymax = UCI, fill = factor(APOEbin), group = interaction(APOEbin, SCI3.3)),
              alpha = 0.2, color = NA) +
  labs(x = "Age (Years)",
       y = "Probability of Structural Zero, MBI-C (Self Rated)",
       color = "APOE Status",
       fill = "APOE Status",
       linetype = "SCI Group") +
  theme_minimal() +
  scale_color_manual(
    values = c("gray30", "steelblue"),
    labels = c("APOE ε4 non-carrier", "APOE ε4 carrier")
  ) +
  scale_fill_manual(
    values = c("gray30", "steelblue"),
    labels = c("APOE ε4 non-carrier", "APOE ε4 carrier")
  ) +
  facet_grid(. ~ SCI3.3) +  # facets side by side (1 row)
  scale_y_continuous(limits = c(-.05, 1)) + # change limits to 0, 1 for probability
  theme(legend.position = "right")



# Save plots

ggsave("./Plots/mbis_plot_a5_highres.tiff", plot = mbis_plot, width = 8.27, height = 5.83, dpi = 600, units = "in", compression = "lzw", bg = "white")
ggsave("./Plots/mbis_plot_zi_a5_highres.tiff", plot = mbis_plot_zi, width = 8.27, height = 5.83, dpi = 600, units = "in", compression = "lzw", bg = "white")

newdata <- expand.grid(
  AgeC = seq(-15, 15, length.out = 100),
  SCI3.3 = c(0, 1),
  APOEbin = c(0, 1),
  SEX = summary_stats_mbi$SEX,
  EDU = summary_stats_mbi$EDU,
  PHQCat = summary_stats_mbi$PHQCat,
  BMIcat = summary_stats_mbi$BMIcat,
  DropOut = summary_stats_mbi$DropOut,
  ParticipantId = NA
)

# Plot domain trajectories


pred_outD1S <- predict(mod_results_domains_self[[1]]$model,
                    newdata = newdata,
                    type = "response",
                    se.fit = TRUE,
                    allow.new.levels = TRUE)

# Add confidence intervals

newdataD1S <- newdata

newdataD1S$predicted <- pred_outD1S$fit
newdataD1S$se <- pred_outD1S$se.fit
newdataD1S$LCI <- newdataD1S$predicted - 1.96 * newdataD1S$se
newdataD1S$UCI <- newdataD1S$predicted + 1.96 * newdataD1S$se

pred_outD2S <- predict(mod_results_domains_self[[2]]$model,
                       newdata = newdata,
                       type = "response",
                       se.fit = TRUE,
                       allow.new.levels = TRUE)

# Add confidence intervals

newdataD2S <- newdata

newdataD2S$predicted <- pred_outD2S$fit
newdataD2S$se <- pred_outD2S$se.fit
newdataD2S$LCI <- newdataD2S$predicted - 1.96 * newdataD2S$se
newdataD2S$UCI <- newdataD2S$predicted + 1.96 * newdataD2S$se

pred_outD3S <- predict(mod_results_domains_self[[3]]$model,
                       newdata = newdata,
                       type = "response",
                       se.fit = TRUE,
                       allow.new.levels = TRUE)

# Add confidence intervals

newdataD3S <- newdata

newdataD3S$predicted <- pred_outD3S$fit
newdataD3S$se <- pred_outD3S$se.fit
newdataD3S$LCI <- newdataD3S$predicted - 1.96 * newdataD3S$se
newdataD3S$UCI <- newdataD3S$predicted + 1.96 * newdataD3S$se

# Domains 4 and 5 did not converge so do not plot

# Bind all domain dfs

newdataDs <- bind_rows(newdataD1S, newdataD2S, newdataD3S)

# Assign Domain names

newdataDs$Domain <- c(rep(1, 400), rep(2, 400), rep(3, 400))
newdataDs <- newdataDs %>%
  mutate(Domain = factor(Domain,
                         levels = 1:3,
                         labels = c("Apathy", "Affect", "Impulse")))
  
newdataDs <- newdataDs %>%
  mutate(SCI3.3 = factor(SCI3.3,
                         levels = 0:1,
                         labels = c("No SCI", "SCI")))

# convert to raw age

newdataDs$Age <- 65 + newdataDs$AgeC 

# Plot with facets
mbisDs_plot <- ggplot(newdataDs, aes(
  x = Age, y = predicted,
  color = Domain,
  linetype = factor(APOEbin)
)) +
  geom_line(linewidth = 1) +
  geom_ribbon(
    aes(ymin = LCI, ymax = UCI,
        fill = Domain,
        group = interaction(Domain, APOEbin)),
    alpha = 0.2, color = NA
  ) +
  labs(
    x = "Age (Years)",
    y = "Predicted (Expected) MBI-C Total Score (Self Rated)",
    color = "Domain",
    fill = "Domain",
    linetype = "APOE Status"
  ) +
  theme_minimal() +
  scale_color_manual(values = c(
    "Apathy"    = "gray30",
    "Affect"    = "steelblue",
    "Impulse"   = "darkorange2"
  )) +
  scale_fill_manual(values = c(
    "Apathy"    = "gray30",
    "Affect"    = "steelblue",
    "Impulse"   = "darkorange2"
  )) +
  scale_linetype_manual(values = c(
    "0" = "solid", "1"= "dotted"),
    labels = c(
    "APOE ε4 non-carrier", "APOE ε4 carrier")) +
  scale_y_continuous(limits = c(0, 3)) +
  facet_grid(. ~ SCI3.3    ) +  # facets side by side (1 row)
  theme(legend.position = "right")


# Save plot

ggsave("./Plots/mbis_plot_domains_a5_highres.tiff", plot = mbisDs_plot, width = 8.27, height = 5.83, dpi = 600, units = "in", compression = "lzw", bg = "white")



##### PROXY #####
#################


# Visualise the data


# Visualise the data

summary_stats_mbii <- mbii %>%
  summarize(
    SEX    = get_mode(SEX),
    EDU    = get_mode(EDU),
    APOEbin   = get_mode(APOEbin),
    BMIcat = get_mode(BMIcat),
    PHQCat = mean(PHQCat, na.rm = TRUE),
    DropOut = mean(DropOut, na.rm = TRUE)
  )


newdata_proxy <- expand.grid(
  AgeC = seq(-15, 15, length.out = 100),
  SCI3.3 = c(0, 1),
    APOEbin = c(0, 1),
  SEX = summary_stats_mbii$SEX,
  EDU = summary_stats_mbii$EDU,
  PHQCat = summary_stats_mbii$PHQCat,
  BMIcat = summary_stats_mbii$BMIcat,
  DropOut = summary_stats_mbii$DropOut,
  ParticipantId = NA
)

# Add in raw age

newdata_proxy$Age <- 65 + newdata_proxy$AgeC

# plot negative binomial

pred_out <- predict(mod_results_proxy$model,
                    newdata = newdata_proxy,
                    type = "response",
                    se.fit = TRUE,
                    allow.new.levels = TRUE)

# Add confidence intervals

newdata_proxy$predicted <- pred_out$fit
newdata_proxy$se <- pred_out$se.fit
newdata_proxy$LCI <- newdata_proxy$predicted - 1.96 * newdata_proxy$se
newdata_proxy$UCI <- newdata_proxy$predicted + 1.96 * newdata_proxy$se

newdata_proxy <- newdata_proxy %>%
  mutate(SCI3.3 = factor(SCI3.3,
                          levels = 0:1,
                          labels = c("No SCI", "SCI")))

mbii_plot <- ggplot(newdata_proxy, aes(x = Age, y = predicted,
                                 color = factor(APOEbin))) +
  geom_line(linewidth = 1) +
  geom_ribbon(aes(ymin = LCI, ymax = UCI, fill = factor(APOEbin), group = interaction(APOEbin, SCI3.3)),
              alpha = 0.2, color = NA) +
  labs(x = "Age (Years)",
       y = "Predicted (Expected) MBI-C Total Score (Proxy Rated)",
       color = "APOE Status",
       fill = "APOE Status",
       linetype = "SCI Group") +
  theme_minimal() +
  scale_color_manual(
    values = c("gray30", "steelblue"),
    labels = c("APOE ε4 non-carrier", "APOE ε4 carrier")
  ) +
  scale_fill_manual(
    values = c("gray30", "steelblue"),
    labels = c("APOE ε4 non-carrier", "APOE ε4 carrier")
  ) +
  facet_grid(. ~ SCI3.3    ) +  # facets side by side (1 row)
  theme(legend.position = "right")

# Plot zero-inflated probabilities

newdata_zi <- expand.grid(
  AgeC = seq(-15, 15, length.out = 100),
  SCI3.3 = c(0,1),
  APOEbin = c(0, 1),
  SEX = summary_stats_mbii$SEX,
  EDU = summary_stats_mbii$EDU,
  PHQCat = summary_stats_mbii$PHQCat,
  BMIcat = summary_stats_mbii$BMIcat,
  DropOut = summary_stats_mbii$DropOut,
  ParticipantId = NA
)

# Add raw age

newdata_zi$Age <- 65 + newdata_zi$AgeC 

pred_out_zi <- predict(mod_results_proxy$model,
                       newdata = newdata_zi,
                       type = "zprob",
                       se.fit = TRUE,
                       allow.new.levels = TRUE)

# Add confidence intervals

newdata_zi$zprob <- pred_out_zi$fit
newdata_zi$se <- pred_out_zi$se.fit
newdata_zi$LCI <- plogis(qlogis(newdata_zi$zprob)) - 1.96 * newdata_zi$se
newdata_zi$UCI <- plogis(qlogis(newdata_zi$zprob)) + 1.96 * newdata_zi$se

# make SCI factor

newdata_zi <- newdata_zi %>%
  mutate(SCI3.3 = factor(SCI3.3,
                         levels = 0:1,
                         labels = c("No SCI", "SCI")))

mbis_plot_zi <- ggplot(newdata_zi, aes(x = Age, y = zprob,
                                       color = factor(APOEbin))) +
  geom_line(linewidth = 1) +
  geom_ribbon(aes(ymin = LCI, ymax = UCI, fill = factor(APOEbin), group = interaction(APOEbin, SCI3.3)),
              alpha = 0.2, color = NA) +
  labs(x = "Age (Years)",
       y = "Probability of Structural Zero, MBI-C (Proxy Rated)",
       color = "APOE Status",
       fill = "APOE Status") +
  theme_minimal() +
  scale_color_manual(
    values = c("gray30", "steelblue"),
    labels = c("APOE ε4 non-carrier", "APOE ε4 carrier")
  ) +
  scale_fill_manual(
    values = c("gray30", "steelblue"),
    labels = c("APOE ε4 non-carrier", "APOE ε4 carrier")
  ) +
  facet_grid(. ~ SCI3.3)+
  scale_y_continuous(limits = c(-.05, 1)) + # change limits to 0, 1 for probability
  theme(legend.position = "right")


# Save plots

ggsave("./Plots/mbii_plot_a5_highres.tiff", plot = mbii_plot, width = 8.27, height = 5.83, dpi = 600, units = "in", compression = "lzw", bg = "white")
ggsave("./Plots/mbii_plot_zi_a5_highres.tiff", plot = mbis_plot_zi, width = 8.27, height = 5.83, dpi = 600, units = "in", compression = "lzw", bg = "white")


# Plot domain trajectories

newdata_proxy <- expand.grid(
  AgeC = seq(-15, 15, length.out = 100),
  SCI3.3 = c(0, 1),
  APOEbin = c(0, 1),
  SEX = summary_stats_mbii$SEX,
  EDU = summary_stats_mbii$EDU,
  PHQCat = summary_stats_mbii$PHQCat,
  BMIcat = summary_stats_mbii$BMIcat,
  DropOut = summary_stats_mbii$DropOut,
  ParticipantId = NA
)

pred_outD1I <- predict(mod_results_domains_proxy[[1]]$model,
                       newdata = newdata_proxy,
                       type = "response",
                       se.fit = TRUE,
                       allow.new.levels = TRUE)

# Add confidence intervals

newdataD1I <- newdata_proxy

newdataD1I$predicted <- pred_outD1I$fit
newdataD1I$se <- pred_outD1I$se.fit
newdataD1I$LCI <- newdataD1I$predicted - 1.96 * newdataD1I$se
newdataD1I$UCI <- newdataD1I$predicted + 1.96 * newdataD1I$se

pred_outD2I <- predict(mod_results_domains_proxy[[2]]$model,
                       newdata = newdata_proxy,
                       type = "response",
                       se.fit = TRUE,
                       allow.new.levels = TRUE)

# Add confidence intervals

newdataD2I <- newdata_proxy

newdataD2I$predicted <- pred_outD2I$fit
newdataD2I$se <- pred_outD2I$se.fit
newdataD2I$LCI <- newdataD2I$predicted - 1.96 * newdataD2I$se
newdataD2I$UCI <- newdataD2I$predicted + 1.96 * newdataD2I$se

pred_outD3I <- predict(mod_results_domains_proxy[[3]]$model,
                       newdata = newdata_proxy,
                       type = "response",
                       se.fit = TRUE,
                       allow.new.levels = TRUE)

# Add confidence intervals

newdataD3I <- newdata_proxy

newdataD3I$predicted <- pred_outD3I$fit
newdataD3I$se <- pred_outD3I$se.fit
newdataD3I$LCI <- newdataD3I$predicted - 1.96 * newdataD3I$se
newdataD3I$UCI <- newdataD3I$predicted + 1.96 * newdataD3I$se


pred_outD4I <- predict(mod_results_domains_proxy[[4]]$model,
                       newdata = newdata_proxy,
                       type = "response",
                       se.fit = TRUE,
                       allow.new.levels = TRUE)

# Add confidence intervals

newdataD4I <- newdata_proxy

newdataD4I$predicted <- pred_outD4I$fit
newdataD4I$se <- pred_outD4I$se.fit
newdataD4I$LCI <- newdataD4I$predicted - 1.96 * newdataD4I$se
newdataD4I$UCI <- newdataD4I$predicted + 1.96 * newdataD4I$se


pred_outD5I <- predict(mod_results_domains_proxy[[5]]$model,
                       newdata = newdata_proxy,
                       type = "response",
                       se.fit = TRUE,
                       allow.new.levels = TRUE)

# Add confidence intervals

newdataD5I <- newdata_proxy

newdataD5I$predicted <- pred_outD5I$fit
newdataD5I$se <- pred_outD5I$se.fit
newdataD5I$LCI <- newdataD5I$predicted - 1.96 * newdataD5I$se
newdataD5I$UCI <- newdataD5I$predicted + 1.96 * newdataD5I$se

# Bind all domain dfs

newdataDI <- bind_rows(newdataD1I, newdataD2I, newdataD3I, newdataD4I, newdataD5I)

# Assign Domain names

newdataDI$Domain <- c(rep(1, 400), rep(2, 400), rep(3, 400), rep(4, 400),  rep(5, 400))
newdataDI <- newdataDI %>%
  mutate(Domain = factor(Domain,
                         levels = 1:5,
                         labels = c("Apathy", "Affect", "Impulse", "Social", "Psychosis")))

# Make SCI a factor

newdataDI <- newdataDI %>%
            mutate(SCI3.3 = factor(SCI3.3,
                                   levels = 0:1,
                                   labels = c("No SCI", "SCI")))

# add in raw age

newdataDI$Age <- 65 + newdataDI$AgeC 

# Plot

mbiDI_plot <- ggplot(newdataDI, aes(
  x = Age, y = predicted,
  color = Domain,
  linetype = factor(APOEbin)
)) +
  geom_line(linewidth = 1) +
  geom_ribbon(
    aes(ymin = LCI, ymax = UCI,
        fill = Domain,
        group = interaction(Domain, factor(APOEbin))),
    alpha = 0.2, color = NA
  ) +
  labs(
    x = "Age (Years)",
    y = "Predicted (Expected) MBI-C Total Score (Proxy Rated)",
    color = "Domain",
    fill = "Domain",
    linetype = "APOE Status"
  ) +
  theme_minimal() +
  scale_color_manual(values = c(
    "Apathy"    = "gray30",
    "Affect"    = "steelblue",
    "Impulse"   = "darkorange2",
    "Social"    = "seagreen3",
    "Psychosis" = "mediumpurple3"
  )) +
  scale_fill_manual(values = c(
    "Apathy"    = "gray30",
    "Affect"    = "steelblue",
    "Impulse"   = "darkorange2",
    "Social"    = "seagreen3",
    "Psychosis" = "mediumpurple3"
  )) +
  scale_linetype_manual(values = c(
    "0" = "solid", "1"= "dotted"),
    labels = c(
      "APOE ε4 non-carrier", "APOE ε4 carrier")) +
  scale_y_continuous(limits = c(0, 3)) +
  facet_grid(. ~SCI3.3)+
theme(legend.position = "right")



# Save plot

ggsave("./Plots/mbii_plot_domains_a5_highres.tiff", plot = mbiDI_plot, width = 8.27, height = 5.83, dpi = 600, units = "in", compression = "lzw", bg = "white")




##### COGNITION #####
#####################


# Visualise the data

# Get modal levels

summary_stats_cog <- cog %>%
  summarize(
    SEX    = get_mode(SEX),
    EDU    = get_mode(EDU),
    BMIcat = get_mode(BMIcat),
    PHQCat = mean(PHQCat, na.rm = TRUE),
    DropOut = mean(DropOut, na.rm = TRUE),
    NumRepeats = mean(NumRepeats, na.rm = TRUE)
  )



newdata_cog <- expand.grid(
  AgeC = seq(-15, 15, length.out = 100),
  SCI3.3 = c(0, 1),
  APOEbin = c(0, 1),
  SEX = summary_stats_cog$SEX,
  EDU = summary_stats_cog$EDU,
  PHQCat = summary_stats_cog$PHQCat,
  BMIcat = summary_stats_cog$BMIcat,
  DropOut = summary_stats_cog$DropOut,
  NumRepeats = summary_stats_cog$NumRepeats,
  ParticipantId = NA
)

# Add in raw age

newdata_cog$Age <- 65 + newdata_cog$AgeC

pred_out_cog <- predict(mod_results_cog$model,
                    newdata = newdata_cog,
                    type = "response",
                    se.fit = TRUE,
                    allow.new.levels = TRUE)

# Add confidence intervals

newdata_cog$predicted <- pred_out_cog$fit
newdata_cog$se <- pred_out_cog$se.fit
newdata_cog$LCI <- newdata_cog$predicted - 1.96 * newdata_cog$se
newdata_cog$UCI <- newdata_cog$predicted + 1.96 * newdata_cog$se


cog_plot <- ggplot(newdata_cog, aes(x = Age, y = predicted,
                                 color = factor(SCI3.3),
                                 linetype = factor(APOEbin))) +
  geom_line(linewidth = 1) +
  geom_ribbon(aes(ymin = LCI, ymax = UCI, fill = factor(SCI3.3), group = interaction(SCI3.3, APOEbin)),
              alpha = 0.2, color = NA) +
  labs(x = "Age (Years)",
       y = "Predicted (Expected) Cognitive Composite Score",
       color = "SCI Group",
       fill = "SCI Group",
       linetype = "APOE ε4 carrier") +
  theme_minimal() +
  scale_color_manual(
    values = c("gray30", "steelblue"),
    labels = c("No SCI", "SCI")
  ) +
  scale_fill_manual(
    values = c("gray30", "steelblue"),
    labels = c("No SCI", "SCI")
  ) +
  theme(legend.position = "right")

# Save plot

ggsave("./Plots/cog_plot_a5_highres.tiff", plot = cog_plot, width = 8.27, height = 5.83, dpi = 600, units = "in", compression = "lzw", bg = "white")

# Dual plot of cognition and MBI

# Plot lines by SCI 

newdata_dualmbi <- expand.grid(
  AgeC = seq(-15, 15, length.out = 100),
  SCI3.3 = 0,
  APOEbin = 0,
  SEX = summary_stats_mbi$SEX,
  EDU = summary_stats_mbi$EDU,
  PHQCat = summary_stats_mbi$PHQCat,
  BMIcat = summary_stats_mbi$BMIcat,
  DropOut = summary_stats_mbi$DropOut,
  ParticipantId = NA
)

# Add in raw age

newdata_dualmbi$Age <- 65 + newdata_dualmbi$AgeC

pred_out_dualmbi <- predict(mod_results_self$model,
                    newdata = newdata_dualmbi,
                    type = "response",
                    se.fit = TRUE,
                    allow.new.levels = TRUE)

# Add confidence intervals

newdata_dualmbi$predicted <- pred_out_dualmbi$fit
newdata_dualmbi$se <- pred_out_dualmbi$se.fit
newdata_dualmbi$LCI <- newdata_dualmbi$predicted - 1.96 * newdata_dualmbi$se
newdata_dualmbi$UCI <- newdata_dualmbi$predicted + 1.96 * newdata_dualmbi$se


newdata_dualcog <- expand.grid(
  AgeC = seq(-15, 15, length.out = 100),
  SCI3.3 = 0,
  APOEbin = 0,
  SEX = summary_stats_cog$SEX,
  EDU = summary_stats_cog$EDU,
  PHQCat = summary_stats_cog$PHQCat,
  BMIcat = summary_stats_cog$BMIcat,
  DropOut = summary_stats_cog$DropOut,
  NumRepeats = summary_stats_cog$NumRepeats,
  ParticipantId = NA
)
     
# Add in raw age

newdata_dualcog$Age <- 65 + newdata_dualcog$AgeC

pred_out_dualcog <- predict(mod_results_cog$model,
                            newdata = newdata_dualcog,
                            type = "response",
                            se.fit = TRUE,
                            allow.new.levels = TRUE)

# Add confidence intervals

newdata_dualcog$predicted <- pred_out_dualcog$fit
newdata_dualcog$se <- pred_out_dualcog$se.fit
newdata_dualcog$LCI <- newdata_dualcog$predicted - 1.96 * newdata_dualcog$se
newdata_dualcog$UCI <- newdata_dualcog$predicted + 1.96 * newdata_dualcog$se

# Make sure SCI is a factor
newdata_dualmbi$SCI3.3 <- factor(newdata_dualmbi$SCI3.3)
newdata_dualcog$SCI3.3 <- factor(newdata_dualcog$SCI3.3)

# Get the range of predicted values
range_mbi <- range(newdata_dualmbi$predicted, na.rm = TRUE)
range_cog <- range(newdata_dualcog$predicted, na.rm = TRUE)

# Define transformation functions
scale_factor <- diff(range_mbi) / diff(range_cog)
offset <- range_mbi[1] - range_cog[1] * scale_factor

rescale_cog <- function(x) x * scale_factor + offset
inv_rescale_cog <- function(x) (x - offset) / scale_factor

# Create merged dataset

# MBI-C data: no transformation needed
mbi_df <- newdata_dualmbi %>%
  transmute(
    AgeC,
    Age,
    SCI3.3,
    APOEbin,
    pred = predicted,
    lci = LCI,
    uci = UCI,
    plot_type = "MBI-C"
  )

# Cognition data: rescale predicted values
cog_df <- newdata_dualcog %>%
  transmute(
    AgeC,
    Age,
    SCI3.3,
    APOEbin,
    pred = rescale_cog(predicted),
    lci = rescale_cog(LCI),
    uci = rescale_cog(UCI),
    plot_type = "Cognition"
  )

# Merge 
plot_df <- bind_rows(mbi_df, cog_df)

plot_df <- plot_df %>%
  mutate(SCI3.3 = factor(SCI3.3,
                          levels = 0:1,
                          labels = c("No SCI", "SCI")))

plot_df <- plot_df %>%
  mutate(APOEbin = factor(APOEbin,
                         levels = 0:1,
                         labels = c("APOE ε4 non-carrier", "APOE ε4 carrier")))

# Plot
dualplot <- ggplot(plot_df, aes(x = Age, y = pred, group = interaction(plot_type))) +
  geom_ribbon(aes(ymin = lci, ymax = uci), alpha = 0.2) +
   geom_line(aes(color = plot_type), size = 1) +
   scale_y_continuous(
    name = "Predicted (Expected) MBI-C Self Total",
    sec.axis = sec_axis(~inv_rescale_cog(.), name = "Predicted (Expected) Cognitive Composite")
  ) +
     scale_color_manual(
    values = c("MBI-C" = "grey30", "Cognition" = "steelblue"),
    labels = c("Cognitive Composite", "MBI-C"),
    name = "Outcome"
  ) +
    scale_fill_manual(
    values = c("MBI-C" = "grey30", "Cognition" = "steelblue"),
    labels = c("Cognitive Composite","MBI-C" ),
    name = "Outcome"
  ) +
    theme_minimal(base_size = 14) +
  labs(x = "Age (Years)") +
  theme(
    axis.title.y.left = element_text(color = "grey30"),
    axis.title.y.right = element_text(color = "steelblue"),
    legend.position = "top"
  )+
  theme(legend.position = "right")



# Save plot

ggsave("./Plots/dual_plot_a5_highres.tiff", plot = dualplot, width = 8.27, height = 5.83, dpi = 600, units = "in", compression = "lzw", bg = "white")



##### PHQ9 #####
################

# Visualise the data

# Get modal levels
get_mode <- function(x) {
  names(which.max(table(x)))
}

summary_stats_phq <- phq %>%
  summarize(
    SEX    = get_mode(SEX),
    EDU    = get_mode(EDU),
    BMIcat = get_mode(BMIcat),
    DropOut = mean(DropOut, na.rm = TRUE)
  )



newdata_phq <- expand.grid(
  AgeC = seq(-15, 15, length.out = 100),
  SCI3.3 = c(0, 1),
  APOEbin = c(0, 1),
  SEX = summary_stats_phq$SEX,
  EDU = summary_stats_phq$EDU,
  BMIcat = summary_stats_phq$BMIcat,
  DropOut = summary_stats_phq$DropOut,
  ParticipantId = NA
)

# Add in raw age

newdata_phq$Age <- 65+ newdata_phq$AgeC

# Plot negative binomial

pred_out <- predict(mod_results_phq$model,
                    newdata = newdata_phq,
                    type = "response",
                    se.fit = TRUE,
                    allow.new.levels = TRUE)

# Add confidence intervals

newdata_phq$predicted <- pred_out$fit
newdata_phq$se <- pred_out$se.fit
newdata_phq$LCI <- newdata_phq$predicted - 1.96 * newdata_phq$se
newdata_phq$UCI <- newdata_phq$predicted + 1.96 * newdata_phq$se

# make sci factor

newdata_phq <- newdata_phq %>% 
                mutate(SCI3.3 =  factor(SCI3.3, 
                                        levels = 0:1,
                                        labels = c("No SCI", "SCI")))

phq_plot <- ggplot(newdata_phq, aes(x = Age, y = predicted,
                                 color = factor(APOEbin))) +
  geom_line(linewidth = 1) +
  geom_ribbon(aes(ymin = LCI, ymax = UCI, fill = factor(APOEbin), group = interaction(APOEbin, SCI3.3)),
              alpha = 0.2, color = NA) +
  labs(x = "Age (Years)",
       y = "Predicted (Expected) PHQ-9 Total Score",
       color = "APOE Status",
       fill = "APOE Status",
       linetype = "APOE Status") +
  theme_minimal() +
  scale_color_manual(
    values = c("gray30", "steelblue"),
    labels = c("APOE ε4 non-carrier", "APOE ε4 carrier")
  ) +
  scale_fill_manual(
    values = c("gray30", "steelblue"),
    labels = c("APOE ε4 non-carrier", "APOE ε4 carrier")
  ) +
  theme(legend.position = "right")+
  facet_grid(. ~ SCI3.3)


# Save plots

ggsave("./Plots/phq_plot_a5_highres.tiff", plot = phq_plot, width = 8.27, height = 5.83, dpi = 600, units = "in", compression = "lzw", bg = "white")

