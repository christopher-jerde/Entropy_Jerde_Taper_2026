#simulation in glmmTMB for cutthroat data set to get bootstrapped evidence profile
#30 Sept 2025 Modified 5 Jan 2026 
#CL Jerde
########################

#clear the environment
rm(list = ls()) #best practice when working with many scripts


library(here) #localized file paths
library(tidyverse) #compendium of useful packages
library(janitor) #cleans data
library(glmmTMB) #allows for linear and non-linear mixed effects models
library(lterdatasampler) #data and idea used from https://rossyndicate.github.io/Intro_R_Data_Analysis/
library(MuMIn) # AIC and BIC table maker

#cutthroat data
#get data from lterdatasampler data and format variable names
data_raw_cutthroat<-and_vertebrates |> clean_names() 

#clean data
data_clean_cutthroat<-data_raw_cutthroat |>
  filter(species=="Cutthroat trout") |>
  select(section, weight_g) |>
  drop_na() |> 
  transform(section = as.factor(section))

#single bootstrap replicate function
single_boot_BIC<-function(data) {
  rep<- data |> group_by(section) |>
    sample_frac(replace=TRUE)

  
  model.2.norm <- glmmTMB( weight_g ~ section, data = rep)
  model.2.gamma <- glmmTMB( weight_g ~ section,family = Gamma(link = "log"), data = rep)
  model.4.norm <- glmmTMB(weight_g ~ section,dispformula = ~ section,data = rep)
  model.4.gamma <- glmmTMB(weight_g ~ section,family = Gamma(link = "log"),dispformula = ~ section,data = rep)

  ms_out<-as.data.frame(model.sel(model.2.norm, model.2.gamma,model.4.norm,model.4.gamma,
                                  rank = "BIC"))
  best<- ms_out |> filter(delta==0)
  m_name<-rownames(best)
  beta<-if(m_name %in% c("model.1.norm","model.1.gamma","model.3.norm","model.3.gamma")){0}
  else{
    if(m_name %in% c("model.2.norm", "model.4.norm")){as.numeric(fixef(get(m_name))$cond[2])}
    else{exp(as.numeric(fixef(get(m_name))$cond[1]) + as.numeric(fixef(get(m_name))$cond[2])) - exp(as.numeric(fixef(get(m_name))$cond[1]))
    }}
  
  return(list(best_model = m_name, beta_1=beta))
}


n_boot=1000
results <- data.frame(model = character(), est = numeric())

# For BIC
set.seed(1345) # for pub data.
for(i in 1:n_boot) {
  output<-single_boot_BIC(data_clean_cutthroat)
  results<-rbind(results, data.frame(model=output$best_model,est=output$beta_1))
}

results # provides the best selected model and the estimate of beta_1





# organize and save results
write.csv(results,here("data","BIC_profile_with_est_1000_cutthroat_MOD_2026.csv"))































