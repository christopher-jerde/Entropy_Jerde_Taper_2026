#Empirical resampling of cutthroat trout data
# CL Jerde
# Modified 12 Nov 2026


#Clear the environment
rm(list=ls())

#libraries
library(tidyverse)
library(here)
library(janitor)
library(car)
library(glmmTMB)
library(lterdatasampler)
library(MuMIn)
library(rstatix)

#####################################
data_raw <- and_vertebrates |> clean_names()
data_clean <- data_raw |>
  filter(species == "Cutthroat trout") |>
  select(section, weight_g) |>
  drop_na() |>
  transform(section = as.factor(section))

alpha  <- 0.05
n_reps <- 1000
n_tot_values <- c(10, 26, 50, 76, 100, 250, 500, 750, 1000, 2500, 5000, 7500, 10000)

# Master list to hold one data frame per n_tot
all_results <- vector("list", length(n_tot_values))
names(all_results) <- paste0("n_tot_", n_tot_values)

for (j in seq_along(n_tot_values)) {
  
  n_tot <- n_tot_values[j]
  message("Running n_tot = ", n_tot)
  
  results_list <- vector("list", n_reps)
  
  for (i in seq_len(n_reps)) {
    
    rep <- data_clean |> group_by(section) |>
      slice_sample(n = n_tot / 2) |> ungroup(section)
    
    norm_out <- rep |> group_by(section) |> shapiro_test(weight_g)
    cc_norm  <- norm_out |> filter(section == "CC") |> mutate(result = ifelse(p < alpha, "reject", "fail")) |> pull(result)
    og_norm  <- norm_out |> filter(section == "OG") |> mutate(result = ifelse(p < alpha, "reject", "fail")) |> pull(result)
    
    var_test <- leveneTest(weight_g ~ section, data = rep)
    var_out  <- ifelse(var_test[["Pr(>F)"]][1] < alpha, "reject", "fail")
    
    students     <- t.test(weight_g ~ section, data = rep, equal.var = TRUE)
    students_out <- ifelse(students$p.value < alpha, "reject", "fail")
    
    welchs     <- t.test(weight_g ~ section, data = rep, equal.var = FALSE)
    welchs_out <- ifelse(welchs$p.value < alpha, "reject", "fail")
    
    model.1.norm  <- glmmTMB(weight_g ~ 1, data = rep)
    model.1.gamma <- glmmTMB(weight_g ~ 1, family = Gamma(link = "log"), data = rep)
    model.2.norm  <- glmmTMB(weight_g ~ section, data = rep)
    model.2.gamma <- glmmTMB(weight_g ~ section, family = Gamma(link = "log"), data = rep)
    model.3.norm  <- glmmTMB(weight_g ~ 1, dispformula = ~ section, data = rep)
    model.3.gamma <- glmmTMB(weight_g ~ 1, family = Gamma(link = "log"), dispformula = ~ section, data = rep)
    model.4.norm  <- glmmTMB(weight_g ~ section, dispformula = ~ section, data = rep)
    model.4.gamma <- glmmTMB(weight_g ~ section, family = Gamma(link = "log"), dispformula = ~ section, data = rep)
    
    ms_out <- as.data.frame(model.sel(
      model.1.norm, model.1.gamma,
      model.2.norm, model.2.gamma,
      model.3.norm, model.3.gamma,
      model.4.norm, model.4.gamma,
      rank = "BIC"
    ))
    
    best   <- ms_out |> filter(delta == 0)
    m_name <- rownames(best)[1]  # take top model only, guards against ties
    
    temp <- c(cc_norm, og_norm, var_out, students_out, welchs_out, m_name)
    results_list[[i]] <- temp
  }
  
  # Build data frame for this n_tot and tag each row with n_tot and replicate number
  results_df <- as.data.frame(do.call(rbind, results_list))
  colnames(results_df) <- c("cc_norm", "og_norm", "var_out", "students_out", "welchs_out", "best_model")
  results_df$n_tot     <- n_tot
  results_df$replicate <- seq_len(n_reps)
  
  # Store in master list
  all_results[[j]] <- results_df
}

# Combine everything into one tidy data frame
final_df <- bind_rows(all_results)


