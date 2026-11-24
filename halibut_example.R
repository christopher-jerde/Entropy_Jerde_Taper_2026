# 2 July 2025 (updated 23 November 2026)
# Models  used in evidence paper (Halibut example)
# CL Jerde
# #########################################################

#clear the environment
rm(list = ls()) 


library(here) #localized file paths
library(tidyverse) #compendium of useful packages
library(janitor) #cleans data
library(car) # needed for using the Levene's test
library(glmmTMB) #allows for linear and non-linear mixed effects models
library(MuMIn) # AIC and BIC table maker

#####################################
#import raw data and format variable names
data_raw<-read.csv(here("data", "trawl_data.csv")) |> clean_names()

#clean data
data_clean<-data_raw |>
  select(trawl_led, plasma_cortisol) |>
  drop_na() |> 
  transform(trawl_led = as.factor(trawl_led))

#summary statistics
summary_table<-data_clean |>
  group_by(trawl_led) |> 
  summarize(obs=n(),mean=mean(plasma_cortisol),median=median(plasma_cortisol),variance=var(plasma_cortisol))

#simple plot to visually inspect data
pc_box <- ggplot(data_clean, aes(x = trawl_led, y = plasma_cortisol)) +
  geom_boxplot(fill = "lightgray", color = "black", outlier.shape = NA) +
  geom_jitter(color = "black", shape = 16, width = 0.15, size = 2, alpha = 0.7) +
  labs(x = "Led (off=0, on=1)", y = "Plasma cortisol (ng/mL)") +
  theme_bw()


#####################################
#Classical approaches
#####################################
# Student's t test
t.test(plasma_cortisol ~ trawl_led, data=data_clean, var.equal=TRUE)

#Normality test
data_clean |>
  group_by(trawl_led) |>
  summarise(shapiro_p = shapiro.test(plasma_cortisol)$p.value) # normality testing

leveneTest(plasma_cortisol ~ trawl_led, data = data_clean) 

#Welch's T-test
t.test(plasma_cortisol ~ trawl_led, data=data_clean, var.equal=FALSE)


#####################################
#multi-model approach
#consider two error distributions; normal and gamma 
#####################################

#same mean, same variance or dispersion
#####################################

#normal
model.1.norm <- glmmTMB( plasma_cortisol ~ 1, data = data_clean)
summary(model.1.norm)

#gamma
model.1.gamma <- glmmTMB( plasma_cortisol ~ 1,family = Gamma(link = "log"), data = data_clean)
summary(model.1.gamma)

#different mean, same variance or dispersion
#normal
model.2.norm <- glmmTMB( plasma_cortisol ~ trawl_led, data = data_clean)
summary(model.2.norm)

#gamma
model.2.gamma <- glmmTMB( plasma_cortisol ~ trawl_led,family = Gamma(link = "log"), data = data_clean)
summary(model.2.gamma)

#same mean, different variance or dispersion
#####################################

#normal
model.3.norm <- glmmTMB(
  plasma_cortisol ~ 1,
  dispformula = ~ trawl_led,       # model heteroscedasticity
  data = data_clean
)
summary(model.3.norm)

#gamma
model.3.gamma <- glmmTMB(
  plasma_cortisol ~ 1,
  family = Gamma(link = "log"),
  dispformula = ~ trawl_led,       # model heteroscedasticity
  data = data_clean
)
summary(model.3.gamma)


#different mean, different variance or dispersion
#####################################
#normal
model.4.norm <- glmmTMB(
  plasma_cortisol ~ trawl_led,
  dispformula = ~ trawl_led,       # model heteroscedasticity
  data = data_clean
)
summary(model.4.norm)

#gamma
model.4.gamma <- glmmTMB(
  plasma_cortisol ~ trawl_led,
  family = Gamma(link = "log"),
  dispformula = ~ trawl_led,       # model heteroscedasticity
  data = data_clean
)
summary(model.4.gamma)


#look at the ICs and deltaICs


ms_out_bic_halibut<-model.sel(model.1.norm, model.1.gamma, model.2.norm, model.2.gamma,
                                model.3.norm, model.3.gamma, model.4.norm,model.4.gamma,
                                rank = "BIC") 


summary(model.3.norm)

