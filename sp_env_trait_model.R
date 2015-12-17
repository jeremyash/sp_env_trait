# Code to develop a hierarchical model (linear mixed model) following the framework developed by Pollock, L. J., Morris, W. K., & Vesk, P. A. (2012). The role of functional traits in species distributions revealed through a hierarchical model. Ecography, 35(8), 716–725. http://doi.org/10.1111/j.1600-0587.2011.07085.x.  Code below is modeled after what Pollock et al. present in their appendix, with attempts at replicating the figures presented in their paper.



#----------------
# load libraries
#----------------

library(plyr)
library(dplyr)
library(ggplot2)
library(tidyr)
library(stringr)
library(readr)
library(reshape2)
library(lme4)
library(broom)
library(gridExtra)
library(GGally)
#------------------------------------------------------------------------------


#----------------
# functions
#----------------

fixef_models_dat <- function(LMM) {
# function to extract fixed effects in tidy format for model LMM
  lmer_effects <- tidy(LMM, effects = "fixed") %>%
    mutate(overlap_zero = ifelse(xor(estimate - 1.96*std.error > 0, estimate + 1.96*std.error < 0), "NO", "YES"))
}



ranef_models_dat <- function(LMM) {
# function that exttracts random effects for each species from the model (LMM), calculates the standard normal quantile, standard error, confidence interval and whehter the confidence interval overlaps zero

  ranef_dat <- lme4::ranef(LMM, condVar=TRUE)
  colnames(ranef_dat[[1]])[1] <- "intercept"

  sp_names <- rownames(ranef_dat[[1]])
  var_names <- colnames(ranef_dat[[1]])

  ranef_cond_mode <- ranef_dat[[1]] %>%
   mutate(sp = sp_names) %>%
   gather(ranef_term, ranef_term_val, -sp, convert=TRUE)

  ranef_cond_mode_snq <-ranef_cond_mode %>%
   group_by(ranef_term) %>%
   arrange(ranef_term_val) %>%
   mutate(ranef_term_snq = qnorm(ppoints(ranef_term_val))) %>%
   ungroup() 

  #extract variances and calculate se 
  ranef_dat_var <- attr(ranef_dat[[1]], which="postVar")

  dimnames(ranef_dat_var) <- list(var_names, var_names, sp_names)


  se_df <- adply(ranef_dat_var, c(1,2)) %>%
    filter(X1 == X2) %>%
    mutate(ranef_term = var_names) %>%
    dplyr::select(-X1, -X2) %>%
    gather(sp, ranef_term_var, -ranef_term, convert=TRUE) %>%
    mutate(ranef_term_se = sqrt(ranef_term_var)) %>%
    dplyr::select(-ranef_term_var)

    
  ranef_all_dat <- left_join(ranef_cond_mode_snq, se_df, by = c('ranef_term', 'sp')) %>%
    mutate(ranef_term_loci = ranef_term_val-1.96*ranef_term_se, ranef_term_upci = ranef_term_val+1.96*ranef_term_se) %>%
    mutate(overlap_zero = ifelse(xor(ranef_term_loci > 0, ranef_term_upci < 0), "NO", "YES"))
}





ranef_dotplot <- function(DAT) {
  # function that plots the random effects for each species, colored yellow when the confidence interval does not overlap zero.  The plots are ordered by standard normal quantiles given the number of species

  plot_name <- unique(DAT$ranef_term)
  
  ggplot(DAT) + 
    geom_errorbar(aes(x = ranef_term_snq, ymin=ranef_term_val-1.96*ranef_term_se, ymax=ranef_term_val+1.96*ranef_term_se), color="gray60") +
    geom_point(aes(x = ranef_term_snq, y = ranef_term_val, color = overlap_zero)) +
    scale_color_manual(values = c("darkgoldenrod3", "steelblue4")) +
    labs(y = plot_name, x = "Standard Normal Quantiles") +
    geom_hline(yintercept = 0, size = 0.4, linetype = "dashed") +
    scale_y_continuous(limits = c(round(model_3_limits$min_val), round(model_3_limits$max_val))) +
    coord_flip() 
  }


#------------------------------------------------------------------------------

#----------------
# data
#----------------

dat_raw <- read.csv("../data/sp_envi_traits.csv", stringsAsFactors=FALSE) # data with change in species abundance at each site as rows (delta_quads), environmental variables and traits in columns
#------------------------------------------------------------------------------


#----------------
# model fitting
#----------------

env_trait_model <- lmer(delta_quads ~ bio4 + bio17 + bio9 + Sand_prop + bio15 + Ca_ppm + precip_spr + #environment variables
  ldmc:bio4 + ldmc:bio17 + ldmc:bio9 + ldmc:Sand_prop + ldmc:bio15 + ldmc:Ca_ppm + ldmc:precip_spr + #correlation of responses to traits 
  sla:bio4 + sla:bio17 + sla:bio9 + sla:Sand_prop + sla:bio15 + sla:Ca_ppm + sla:precip_spr +
  l_thick:bio4 + l_thick:bio17 + l_thick:bio9 + l_thick:Sand_prop + l_thick:bio15 + l_thick:Ca_ppm + l_thick:precip_spr +
  veg_height_log:bio4 + veg_height_log:bio17 + veg_height_log:bio9 + veg_height_log:Sand_prop + veg_height_log:bio15 + veg_height_log:Ca_ppm + veg_height_log:precip_spr + 
  (1 +  bio4 + bio17 + bio9 + Sand_prop + bio15 + Ca_ppm + precip_spr | sp), #response to enviro variables allowed to vary by taxa
  data=dat, control=lmerControl(optCtrl=list(maxfun=1e5))))

#------------------------------------------------------------------------------

#----------------
# plots
#----------------
env_trait_model_ranef <- ranef_models_dat(env_trait_model)

#find limits for plotting
env_trait_model_limits <- env_trait_model_ranef %>%
 summarise(min_val = min(ranef_term_loci), max_val = max(ranef_term_upci))

env_trait_model_ranef_plots <- dlply(env_trait_model_ranef, .(ranef_term), function(x) ranef_dotplot(x))


env_trait_model_ranef_plots[[1]]
# ...model_3_ranef_plots[[n]]

#------------------------------------------------------------------------------


#----------------
# fixed effects
#----------------

model_3_fixef <- fixef_models_dat(model_3)

#------------------------------------------------------------------------------

# working on generalizing remaining code to generate other figures from Pollock et al.
