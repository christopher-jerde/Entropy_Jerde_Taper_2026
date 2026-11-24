# 3 July 2025 (updated 23 November 2026)
# Models used in evidence paper (Cutthroat example)
# CL Jerde
# #########################################################

#clear the environment
rm(list = ls()) 


library(here) #localized file paths
library(tidyverse) #compendium of useful packages
library(janitor) #cleans data
library(car) # needed for using the Levene's test
library(glmmTMB) #allows for linear and non-linear mixed effects models
library(lterdatasampler) #data and idea used from https://rossyndicate.github.io/Intro_R_Data_Analysis/
library(MuMIn) # AIC and BIC table maker

#####################################
#get data from lterdatasampler data and format variable names
data_raw<-and_vertebrates |> clean_names() 

#clean data
data_clean<-data_raw |>
  filter(species=="Cutthroat trout") |>
  select(section, weight_g) |>
  drop_na() |> 
  transform(section = as.factor(section))

#summary statistics
summary_table<-data_clean |>
  group_by(section) |> 
  summarize(obs=n(),mean=mean(weight_g),median=median(weight_g),variance=var(weight_g))

#simple plot to visually inspect data
fish_box <- ggplot(data_clean, aes(x = section, y = weight_g)) +
  geom_boxplot(fill = "lightgray", color = "black", outlier.shape = NA) +
  geom_jitter(color = "gray", shape = 1, width = 0.25, size = 2, alpha = 0.2) +
  labs(x = "Forest type", y = "Fish weight (g)") +
  theme_bw()


#####################################
#Classical approaches
#####################################

#Student's t-test
t.test(weight_g ~ section, data=data_clean, var.equal=TRUE)

trout_QQ<- ggplot(data_clean, aes(sample = weight_g)) +
  stat_qq() +
  stat_qq_line() +
  facet_wrap(~ section) +
  theme_minimal() +
  labs(title = "QQ‑Plot of Trout Weights by Forest Section",
       x = "Theoretical Quantiles",
       y = "Fish weight (g)")

#test for equal variance
leveneTest(weight_g ~ section, data = data_clean) 

# Welch's T-test (for completeness)
t.test(weight_g ~ section, data=data_clean, var.equal=FALSE)


#####################################
#multi-model approach
#consider two error distributions; normal and gamma 
#####################################

#same mean, same variance
#####################################

#normal
model.1.norm <- glmmTMB(weight_g ~ 1, data = data_clean)
summary(model.1.norm)

#gamma
model.1.gamma <- glmmTMB( weight_g ~ 1,family = Gamma(link = "log"), data = data_clean)
summary(model.1.gamma)


#different mean, same variance
#normal
model.2.norm <- glmmTMB( weight_g ~ section, data = data_clean)
summary(model.2.norm)

#gamma
model.2.gamma <- glmmTMB( weight_g ~ section,family = Gamma(link = "log"), data = data_clean)
summary(model.2.gamma)

#same mean, different variance or dispersion parameter
#####################################

#normal
model.3.norm <- glmmTMB(
  weight_g ~ 1,
  dispformula = ~ section,       # model heteroscedasticity
  data = data_clean
)
summary(model.3.norm)

#gamma
model.3.gamma <- glmmTMB(
  weight_g ~ 1,
  family = Gamma(link = "log"),
  dispformula = ~ section,       # model heteroscedasticity
  data = data_clean
)
summary(model.3.gamma)


#different mean, different variance or dispersion parameter
#####################################
#normal
model.4.norm <- glmmTMB(
  weight_g ~ section,
  dispformula = ~ section,       # model heteroscedasticity
  data = data_clean
)
summary(model.4.norm)

#gamma
model.4.gamma <- glmmTMB(
  weight_g ~ section,
  family = Gamma(link = "log"),
  dispformula = ~ section,       # model heteroscedasticity
  data = data_clean)
summary(model.4.gamma)



#look at the ICs and deltaICs

ms_out_bic_cutthroat<-model.sel(model.1.norm, model.1.gamma, model.2.norm, model.2.gamma,
                                model.3.norm, model.3.gamma, model.4.norm,model.4.gamma,
          rank = "BIC") 

#best model as selected by delta_BIC is gamma distributed error distributions with the different mean but same dispersion 
summary(model.2.gamma)

