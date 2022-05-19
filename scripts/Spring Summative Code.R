# Loaded packages in 
library(tidyverse)
library(readxl)
library(GGally)
library(emmeans)
library(rstatix)
library(kableExtra)
library(performance)
library(patchwork)
library(here)


# Importing the sheets via a read_excel path
f0lifespan <- (read_excel(path = "Data/elegans.xlsx", sheet = "lifespan_f0", na = "NA"))
f1lifespan <- (read_excel(path = "Data/elegans.xlsx", sheet = "f1_lifespan", na = "NA"))
f0reproduction <- (read_excel(path = "Data/elegans.xlsx", sheet = "reproduction_f0", na = "NA"))
f1reproduction <- (read_excel(path = "Data/elegans.xlsx", sheet = "reproduction_f1", na = "NA"))

# Lowercased everything so the data fits the same format
f0lifespan$rnai <- tolower(f0lifespan$rnai)
f0lifespan$treatment <- tolower(f0lifespan$treatment)
f1lifespan$parental_rnai <- tolower(f1lifespan$parental_rnai)
f1lifespan$parental_treatment <- tolower(f1lifespan$parental_treatment)
f1lifespan$treatment <- tolower(f1lifespan$treatment) 

# Calculating the longevity for f0 
f0lifespan$longevity <- as.numeric (difftime
                                     (f0lifespan$death_date, f0lifespan$set_up_date,
                                       units="days"))
# Calculating the longevity for f1 
f1lifespan$longevity <- as.numeric (difftime
                                     (f1lifespan$death_date, f1lifespan$set_up_date,
                                       units="days"))

# Using GGally to look at whole plots of the datasets 
GGally::ggpairs(f0lifespan)
GGally::ggpairs(f1lifespan)
GGally::ggpairs(f0reproduction)
GGally::ggpairs(f1reproduction)

# Created a summary for rnai and offspring for f0
 f0reproduction_summary <- f0reproduction %>%
  group_by(rnai) %>%
  summarise(mean = mean(offspring),
            sd=sd(offspring))
 
# Created a table for mean and sd of f0 offspring based on rnai treatment
f0reproduction_summary %>%
  kbl(caption="The mean and sd offspring from f0 when treated with ev or raga rnai") %>% 
  kable_styling(bootstrap_options = "striped", full_width = T, position = "left")

#### Testing and working with different models 
        # MODEL 1 
# Gene knockdown AND light/dark treatment 
# F0-lifespan based on THEIR rnai gene and light/dark treatment 
     #  F0 - longevity, rnai, treatment 

# Visualising the data with a box plot 
q2 <- ggplot(f0lifespan, aes(x=rnai, y=longevity, fill=treatment)) +
  geom_boxplot() +
  labs('title' = 'Treatment and genes on longevity',
       y = 'Longevity',
       x = 'RNAi treatment', 
       fill = "Treatment") +
  scale_fill_manual(values = alpha(c("#a84632", "#8ba832"),.3)) 
  
# Creating a linear model 
f0lifespanls1 <- lm(longevity ~ rnai + treatment + rnai + rnai:treatment, data = f0lifespan)

# Assumption checking 
# Doing a performance check to look for normality in the model
performance::check_model(f0lifespanls1)
performance::check_model(f0lifespanls1, check="homogeneity") 
performance::check_model(f0lifespanls1, check=c("normality","qq"))
performance::check_model(f0lifespanls1, check="outliers")
# Doesn't look normal 

# Data Transformations - using BoxCox
#MASS::boxcox(f0lifespanls1)
# Shows parameter = (0-0.5) 

# Using sqrt to transform data 
f0lifespanls1 <- lm(sqrt(longevity) ~ rnai + treatment + rnai + rnai:treatment, data = f0lifespan)
# Using performance check to check for normality 
performance::check_model(f0lifespanls1, check=c("homogeneity", "qq"))

# Using log to transform data 
f0lifespanls1 <- lm(log(longevity) ~ rnai + treatment + rnai + rnai:treatment, data = f0lifespan)
# Using performance check to check for normality 
performance::check_model(f0lifespanls1, check=c("homogeneity", "qq"))

# Found that - Using residual deviance and degrees of freedom 619/394 = 1.6

# using glm to transform data 
f0lifespanls1 <- glm(formula = longevity ~ rnai + treatment + rnai:treatment,
                     family = quasipoisson(), data = f0lifespan)
# Using performance check to check for normality 
performance::check_model(f0lifespanls1, check=c("homogeneity", "qq"))

# Using summary and broom tidy to give a summary of f0lifespanls1 results 
summary(f0lifespanls1)
broom::tidy(f0lifespanls1) 

# model with exp values 
broom::tidy(f0lifespanls1, 
            exponentiate=T, 
            conf.int=T)

# Using drop1 function to see if rnai:treatment interaction term should be kept  
drop1(f0lifespanls1, test = "F")
# Do not keep interaction effect, there is no statistical significance 

# Model without the interaction effect
f0lifespanls1 <- glm(formula = longevity ~ rnai + treatment,
                     family = quasipoisson(), data = f0lifespan)
# Using performance check to check for normality 
performance::check_model(f0lifespanls1, check=c("homogeneity", "qq"))

# Using emmeans for means value for RNAi and treatment
means <- emmeans::emmeans(f0lifespanls1, specs = ~rnai)
means <- emmeans::emmeans(f0lifespanls1, specs = ~treatment)

# Calculating actual values
exp(1.67)
exp(1.81)
exp(1.74)
exp(2.20)
exp(2.32)
exp(2.15)
exp(2.26)
exp(2.27)
exp(2.38)
exp(-24.85)
exp(2.67) 
exp (2.78)
# 14- 16 days 
exp(0.05)
exp(0.19)
# 1 - 1.21 days 
exp(-1.13)
exp (-0.97)
# 0.32 - 0.38 days 


# Using summary and broom::tidy for analysis of new model 
summary(f0lifespanls1)
broom::tidy(f0lifespanls1)

# Making a table of f0lifespan for the write-up based on their rnai gene and light/dark treatment 
f0lifespanls1table <-  f0lifespanls1 %>% broom::tidy(conf.int = T) %>% 
  select(-`std.error`) %>% 
  mutate_if(is.numeric, round, 2) %>% 
  kbl(col.names = c("Predictors",
                    "Estimates",
                    "Z-value",
                    "P",
                    "Lower 95% CI",
                    "Upper 95% CI"),
      caption = "FO Lifepan - Dark/Light Conditions and RNAi gene treatment", 
      booktabs = TRUE) %>% 
  kable_styling(full_width = FALSE, font_size=16, latex_options = c("striped", "hold_position"))



f0lifespanls1table




# MODEL 2
#  F0 longevity based on whether they were in light/dark
 # F0 - longevity, treatment 

# Visualising the data 
q1 <- ggplot(f0lifespan, aes(x=treatment, y=longevity, fill=treatment))+ 
geom_boxplot() + 
labs('title' = ' Treatment on longevity',
       y = 'Longevity',
       x = 'Treatment', 
       fill = "Treatment") +
  scale_fill_manual(values = alpha(c("#a84632", "#8ba832"),.3)) 

# patchwork to create joint plots of the two
q1 + q2

# Saving the plot via ggsave 
ggsave("figures/treatmentslongevity.png", plot= q1+q2 , dpi=900, width = 7, height = 5)

# Creating linear model 
f0lifespanls2 <- lm(longevity ~ treatment, data = f0lifespan)

# performance function - checked for normality
performance::check_model(f0lifespanls2)
performance::check_model(f0lifespanls2, check=c("normality","qq")) 
# Doesn't look great

# Did a transformation test
#MASS::boxcox(f0lifespanls2)
# 0-05

# Using sqrt to transform data 
f0lifespanls2 <- lm(sqrt(longevity) ~  treatment, data = f0lifespan)
# performance function - checked for normality
performance::check_model(f0lifespanls2, check=c("homogeneity", "qq"))
# Looks OK 

# using log to transform data 
f0lifespanls2 <- lm(log(longevity) ~  treatment, data = f0lifespan)
# performance function - checked for normality
performance::check_model(f0lifespanls2, check=c("homogeneity", "qq"))
# sqrt looks better 

# using glm with quasi poisson to transform data 
f0lifespanls2 <- glm(formula = longevity ~ treatment,
                     family = quasipoisson(), data = f0lifespan)
# performance function - checked for normality
performance::check_model(f0lifespanls2, check=c("homogeneity", "qq"))
# Looks the best 

# Using broom::tidy and summary function 
broom::tidy(f0lifespanls2)
summary(f0lifespanls2)

# created a table using KableExtra for the final model choice 
f0lifespanls2table <- 
  f0lifespanls2 %>% broom::tidy(conf.int = T) %>% 
  select(-`std.error`) %>% 
  mutate_if(is.numeric, round, 2) %>% 
  kbl(col.names = c("Predictors",
                    "Estimates",
                    "Z-value",
                    "P",
                    "Lower 95% CI",
                    "Upper 95% CI"),
      caption = "FO Lifepan - Dark/Light conditions", 
      booktabs = TRUE) %>% 
  kable_styling(full_width = FALSE, font_size=16, latex_options = c("striped", "hold_position"))

f0lifespanls2table


# Looked at emmeans data for amount of offspring f0 generation have vs rnai treatment
f0offspringmeans <- emmeans::emmeans(f0lifespanls2, specs = ~treatment)
f0offspringmeans
# Visualised the amount of offspring FO have based on rnai treatment using emmeans data
f0offspringmeans %>% 
  as_tibble() %>% 
  ggplot(aes(x=treatment, 
             y=emmean))+
  geom_pointrange(aes(
    ymin=lower.CL, 
    ymax=upper.CL))


# MODEL 3

# Visualising the data 
f0reproductionls2plot <-  ggplot(f0reproduction, aes(x=rnai, y=offspring, fill=treatment))+
  geom_boxplot()+
  scale_fill_manual(values = alpha(c("#a84632", "#8ba832"),.3))+
  labs( x = "RNAi Treatment",
        y = "Offspring",
        fill ="Treatment")

# Creating a linear model 
f0reproductionls1 <- lm(offspring ~ rnai + treatment + rnai:treatment, data = f0reproduction)

# Using performance check to check for normality 
performance::check_model(f0reproductionls1)
performance::check_model(f0reproductionls1, check=c("normality","qq"))
# Doesn't look good 

# Transforming data using BoxCox
#MASS::boxcox(f0reproductionls1)
# 0-0.5

# Using sqrt to transform data 
f0reproductionls1 <- lm(sqrt(offspring) ~ rnai + treatment + rnai:treatment, data = f0reproduction)
# Using performance check to check for normality 
performance::check_model(f0reproductionls1, check=c("homogeneity", "qq"))

# Using log to transform data 
f0reproductionls1 <- lm(log(offspring) ~ rnai + treatment +  rnai:treatment, data = f0reproduction)
# Using performance check to check for normality 
performance::check_model(f0reproductionls1, check=c("homogeneity", "qq"))

#  glm test 
f0reproductionls1 <- glm(formula = offspring ~ rnai + treatment 
                      + rnai:treatment,
                      family = quasipoisson(), data = f0reproduction)
# Using performance check to check for normality 
performance::check_model(f0reproductionls1, check=c("homogeneity", "qq"))

#  sqrt was best to transform data 

# Using drop1 to see if interaction term should be kept
drop1(f0reproductionls1, test = "F")
# No significant difference 
# Can drop interaction effect from model 

# New linear model without interaction effect  
f0reproductionls2 <- lm(sqrt(offspring) ~ rnai + treatment, data = f0reproduction)
performance::check_model(f0reproductionls2, check=c("homogeneity", "qq"))

# Using broom::tidy and summary to create a summary
broom::tidy(f0reproductionls2)
summary(f0reproductionls2)


# Calculating mean values of treatment and rnai using emmeans function
f0offspringmeans2 <- emmeans::emmeans(f0reproductionls2, specs = ~treatment)
f0offspringmeans2
f0offspringmeans21 <- emmeans::emmeans(f0reproductionls2, specs = ~rnai)
f0offspringmeans21



# Creating a table for the write-up
f0reproductionls2table <- 
  f0reproductionls2 %>% broom::tidy(conf.int = T) %>% 
  select(-`std.error`) %>% 
  mutate_if(is.numeric, round, 2) %>% 
  kbl(col.names = c("Predictors",
                    "Estimates",
                    "Z-value",
                    "P",
                    "Lower 95% CI",
                    "Upper 95% CI"),
      caption = "F0 Reproduction - Dark/Light conditions and RNAi gene treatment", 
      booktabs = TRUE) %>% 
  kable_styling(full_width = FALSE, font_size=16, latex_options = c("striped", "hold_position"))


f0reproductionls2table




# MODEL 4 
#  F1-generation longevity based on their parent's rnai
  # f1 - longevity, parent's rnai

# Visualise the data 
f1lifespanls1plot <- ggplot(f1lifespan, aes(x=parental_rnai, y=longevity, fill= parental_rnai)) + 
  geom_violin()+
   labs( x = "Parental RNAi Treatment", 
         y = "Longevity",
         fill = "Treatment")+
  scale_fill_manual(values = alpha(c("#7732a8", "#32a854"),.3)) 





# creating a linear model
f1lifespanls1 <- lm(longevity ~ parental_rnai, data = f1lifespan )

# Checking the data 
performance::check_model(f1lifespanls1)
# Doesn't look good 

# Using performance check to check for normality 
performance::check_model(f1lifespanls1, check = c("normality", "qq"))


# BoxCox to transform data 
#MASS::boxcox(f1lifespanls1)

# Using sqrt to transform data 
f1lifespanls1 <- lm(sqrt(longevity) ~ parental_rnai, data = f1lifespan)
# Using performance check to check for normality 
performance::check_model(f1lifespanls1, check=c("homogeneity", "qq")) 

# Using log to transform data 
f1lifespanls1 <- lm(log(longevity) ~ parental_rnai, data = f1lifespan)
# Using performance check to check for normality 
performance::check_model(f1lifespanls1, check=c("homogeneity", "qq")) 


# poisson glm
f1lifespanls1<- glm(formula = longevity ~ parental_rnai,
                      family = quasipoisson(), data = f1lifespan)
# Using performance check to check for normality 
performance::check_model(f1lifespanls1, check=c("homogeneity", "qq")) 
# Use Poisson quasi glm 

# Summary of the models
summary(f1lifespanls1)
broom::tidy(f1lifespanls1)



# Doing an anova 
anova(f1lifespanls1)

# looking for emmeans 
f1lifespan1 <- emmeans::emmeans(f1lifespanls1, specs = ~parental_rnai)
f1lifespan1

# making a table of f1lifespan from parental rnai 
f1lifespanls1table <- 
  f1lifespanls1 %>% broom::tidy(conf.int = T) %>% 
  select(-`std.error`) %>% 
  mutate_if(is.numeric, round, 2) %>% 
  kbl(col.names = c("Predictors",
                    "Estimates",
                    "Z-value",
                    "P",
                    "Lower 95% CI",
                    "Upper 95% CI"),
      caption = "How parental RNAi treatment influences longevity", 
      booktabs = TRUE) %>% 
  kable_styling(full_width = FALSE, font_size=16, latex_options = c("striped", "hold_position"))
# ISSUE WITH TABLE 

f1lifespanls1table




# MODEL 5
# F1 - generation longevity based on their parent's rnai and parent's treatment
# F1 - longevity, parent's rnai, parent's treatment 

# Visualising the data 
f1lifespanls3plot <- ggplot(data=f1lifespan, aes(x = tolower(parental_treatment), y = longevity)) +
  geom_boxplot(aes(fill = parental_rnai),
               alpha = 0.2, 
               width = 0.5, 
               outlier.shape=NA)+
  geom_jitter(aes(colour = longevity),
              width=0.2)+
  labs('title' = 'Effect of parental treatment on offspring longevity',
       x = 'Parental Treatment',
       y = 'Longevity of offspring', 
       fill = "Parental RNAi")



# Creating a model 
f1lifespanls2 <-  lm(longevity ~ parental_rnai + parental_treatment + parental_rnai:parental_treatment, data = f1lifespan)
summary(f1lifespanls2)
# Assumption checking 
performance::check_model(f1lifespanls2)
performance::check_model(f1lifespanls2, check=c("normality", "qq"))
# Doesn't look normal 

# Transformation with BoxCox 
# run this, pick a transformation and retest the model fit
#MASS::boxcox(f1lifespanls2)

# transforming data with sqrt 
f1lifespanls2 <- lm(sqrt(longevity) ~ parental_rnai + parental_treatment + parental_treatment:parental_rnai, data = f1lifespan)
performance::check_model(f1lifespanls2, check=c("normality", "qq"))

f1lifespanls2 <- lm(log(longevity) ~ parental_rnai + parental_treatment + parental_treatment:parental_rnai, data = f1lifespan)
performance::check_model(f1lifespanls2, check=c("normality", "qq"))
# Better - the best

# Doing a glm test 
f1lifespanls2<- glm(longevity ~ parental_rnai + parental_treatment 
                     + parental_treatment:parental_rnai, 
              family = gaussian(link = "identity"),
              data = f1lifespan)
# Assumption checking 
performance::check_model(f1lifespanls2, check=c("normality", "qq"))



summary(f1lifespanls2)
# glm exactly the same as lm
# Still doesn't look okay 


#  keep interaction term? 
drop1(f1lifespanls2, test = "F")
#  Don't keep interaction term, no significance 


# New model - use log 
f1lifespanls3 <- lm(log(longevity) ~ parental_rnai + parental_treatment, data = f1lifespan)
# Assumption checking                      
performance::check_model(f1lifespanls3, check=c("homogeneity", "qq"))


# Using broom::tidy as a summary function 
broom::tidy(f1lifespanls3)
summary(f1lifespanls3)

# Making a table for f1 lifespan based on parental rnai and parental treatment 
f1lifespanls3table <- 
  f1lifespanls3 %>% broom::tidy(conf.int = T) %>% 
  select(-`std.error`) %>% 
  mutate_if(is.numeric, round, 2) %>% 
  kbl(col.names = c("Predictors",
                    "Estimates",
                    "Z-value",
                    "P",
                    "Lower 95% CI",
                    "Upper 95% CI"),
      caption = "Longevity of F1 - based on parental RNAi and parental treatment", 
      booktabs = TRUE) %>% 
kable_styling(full_width = FALSE, font_size=12, latex_options = c("striped", "hold_position"))

f1lifespanls3table



# MODEL 6 
# F1 offsprings vs their parent's rnai treatment 
#  F1 offsprings, parent's rnai 

# Visualising the data 
f1reproductionls1plot <- ggplot(f1reproduction, aes(x=parental_rnai, y=offsprings, fill=parental_rnai))+
  geom_boxplot()+
  scale_fill_manual(values = alpha(c("#a84632", "#8ba832"),.3))+
  labs( x = "Parental (FO)  RNAi Treatment",
        y = " F1 Offspring",
        fill ="Parental (FO)  RNAi Treatment",
        title = "F1 offspring based off their parents RNAi treatment")




  
  
# Creating the model 
f1reproductionls1 <- lm(offsprings ~ parental_rnai, data = f1reproduction )

# Assumption Checking 
performance::check_model(f1reproductionls1)
performance::check_model(f1reproductionls1, check=c("homogeneity", "qq"))
# Homogenity not okay 

# Transformation with BoxCox 
# run this, pick a transformation and retest the model fit
#MASS::boxcox(f1reproductionls1)

# use glm 
f1reproductionls1 <- glm(offsprings ~ parental_rnai, 
              family = gaussian(link = "identity"),
              data = f1reproduction)
# Assumption Checking 
performance::check_model(f1reproductionls1, check=c("homogeneity", "qq"))

# looking for emmeans 
f1reproduction1 <- emmeans::emmeans(f1reproductionls1, specs = ~parental_rnai)
f1reproduction1

# Making a summary
broom::tidy(f1reproductionls1)
summary(f1reproductionls1)

# Making a table of f1 reproduction based on parental rnai 
f1reproductionls1table <- 
  f1reproductionls1 %>% broom::tidy(conf.int = T) %>% 
  select(-`std.error`) %>% 
  mutate_if(is.numeric, round, 2) %>% 
  kbl(col.names = c("Predictors",
                    "Estimates",
                    "Z-value",
                    "P",
                    "Lower 95% CI",
                    "Upper 95% CI"),
      caption = "How parental RNAi influences reproduction", 
      booktabs = TRUE) %>% 
  kable_styling(full_width = FALSE, font_size=12, latex_options = c("striped", "hold_position"))



f1reproductionls1table



# MODEL 7
# f1 longevity based on their treatment and their parents treatment 
 # f1 - longevity, parental treatment and their own treatment, with interaction effect  

# Visualise the data 
f1lifespanls5plot <- ggplot(f1lifespan, aes(x=parental_treatment, y=longevity, fill=treatment))+
  geom_boxplot()+
  scale_fill_manual(values = alpha(c("#a84632", "#8ba832"),.3))+
  labs( x = "Parental Treatment",
        y = "Longevity",
        fill ="Own Treatment")

# created a linear model 
f1lifespanls4 <- lm(longevity ~ parental_treatment + treatment + parental_treatment:treatment, data = f1lifespan)
# Assumption Checking 
performance::check_model(f1lifespanls4, check=c("homogeneity", "qq"))

# BoxCox to transform data 
#MASS::boxcox(f1lifespanls4)


# Using sqrt 
f1lifespanls4<- lm(sqrt(longevity) ~ treatment + parental_treatment 
                     + parental_treatment:treatment, data = f1lifespan)
# Assumption Checking 
performance::check_model(f1lifespanls4, check=c("homogeneity", "qq"))

# Using log 
f1lifespanls4 <- lm(log(longevity) ~ treatment + parental_treatment 
                     + parental_treatment:treatment, data = f1lifespan)
# Using performance check to check for normality 
performance::check_model(f1lifespanls4, check=c("homogeneity", "qq"))

# Doing a glm 
f1lifespanls4 <- glm(longevity ~ treatment + parental_treatment 
                      + parental_treatment:parental_rnai, 
                         family = poisson(link = "identity"),
                         data = f1lifespan)
# Using performance check to check for normality 
performance::check_model(f1lifespanls4, check=c("homogeneity", "qq"))

# exponentiate coefficients 
exp(coef(f1lifespanls4)[1]) 
exp(coef(f1lifespanls4)[2])  
exp(coef(f1lifespanls4)[3]) 
exp(coef(f1lifespanls4)[4]) 

# Tidying models (removing log transformation)
broom::tidy(f1lifespanls4, 
            exponentiate=T, 
            conf.int=T)

# fixed mean variance model - chisquared 
drop1(f1lifespanls4, test = "Chisq")

# quasipoisson - overdistribution 
f1lifespanls4 <- glm(longevity ~ treatment + parental_treatment + parental_treatment:treatment,
                      data=f1lifespan, family=quasipoisson(link="log"))
# Using performance check to check for normality 
performance::check_model(f1lifespanls4, check=c("homogeneity", "qq"))

# Summary of model 
summary(f1lifespanls4)
broom::tidy(f1lifespanls4)

# Using drop1 - keep interaction effect?
drop1(f1lifespanls4, 
      test = "F")
# no longer significant 
 

# new model without interaction effect 
f1lifespanls5 <- glm(longevity ~ treatment + parental_treatment,
                      data=f1lifespan, family=quasipoisson(link="log"))
# Using performance check to check for normality 
performance::check_model(f1lifespanls5, check=c("homogeneity", "qq"))

# fixed mean variance model - chisquared
drop1(f1lifespanls5, test = "Chisq")

 
# A summary of the model 
broom::tidy(f1lifespanls5)
summary(f1lifespanls5)

# Looked at emmeans data for amount of offspring f0 generation have vs rnai treatment
f1offspringmeans5 <- emmeans::emmeans(f1lifespanls5, specs = ~treatment)
f1offspringmeans5

# Making a table of f1 lifespan based on their own treatment and parental treatment
f1lifespanls5table <-  
  f1lifespanls5 %>% broom::tidy(conf.int = T) %>% 
  select(-`std.error`) %>% 
  mutate_if(is.numeric, round, 2) %>% 
  kbl(col.names = c("Predictors",
                    "Estimates",
                    "Z-value",
                    "P",
                    "Lower 95% CI",
                    "Upper 95% CI"),
      caption = "Own treatment and parental treatment influence on longevity", 
      booktabs = TRUE) %>% 
  kable_styling(full_width = FALSE, font_size=12, latex_options = c("striped", "hold_position"))


f1lifespanls5table



# Testing to see if replicate/plate values are the same throughout 

## Test for f1 reproduction based on treatment and replicate
f1reproduction %>% group_by(treatment, replicate) %>% summarise(n = n()) %>% view()
# visualising the data 
f1reproduction %>% 
  ggplot(aes(x = as.factor(replicate), y = offsprings, fill = treatment))+geom_boxplot()
# making a linear model 
f1reproductionls5 <- lm(offsprings ~ treatment + as.factor(replicate)
                      + treatment:as.factor(replicate), data = f1reproduction)
# summarising the data 
summary(f1reproductionls5)
broom::tidy(f1reproductionls5)
# using drop1 to check for interaction effect significance 
drop1(f1reproductionls5, test = "F")
# linear model performance check for normality 
performance::check_model(f1reproductionls5, 
                         check = c("homogeneity",
                                   "qq"))

## Test for f0 reproduction based on treatment and replicate
f0reproduction %>% group_by(treatment, replicate) %>% summarise(n = n()) %>% view()
# visualising the data
f0reproduction %>% 
ggplot(aes(x = as.factor(replicate), y = offspring, fill = treatment))+geom_boxplot()
# making a linear model 
f0reproductionls4 <- lm(offspring ~ treatment + as.factor(replicate)
                      + treatment:as.factor(replicate), data = f0reproduction)
# checking for significance of the interaction effect 
drop1(f0reproductionls4, test = "F")
# dropped interaction effect model
f0reproductionls4 <- lm(offspring ~ treatment + as.factor(replicate), data = f0reproduction)
# linear model performance check for normality 
performance::check_model(f0reproductionls4)
# summarising the data 
summary(f0reproductionls4)
broom::tidy(f0reproductionls4)

## Test for f0 longevity based on treatment and replicate

f0lifespan %>% group_by(treatment, plate) %>% summarise(n = n()) %>% view()
# visualising the data
f0lifespan %>% 
  ggplot(aes(x = as.factor(plate), y = longevity, fill = treatment))+geom_boxplot()
# making a linear model 
f0lifespanls4 <- lm(longevity ~ treatment + as.factor(plate)
                      + treatment:as.factor(plate), data = f0lifespan)
# summarising the data 
summary(f0lifespanls4)
broom::tidy(f0lifespanls4)
# checking for significance of the interaction effect 
drop1(f0lifespanls4, test = "F")
# linear model performance check for normality 
performance::check_model(f0lifespanls4, 
                         check = c("homogeneity",
                                   "qq"))



## Test on f1 longevity based on treatment and replicate
f1lifespan %>% group_by(treatment, plate) %>% summarise(n = n()) %>% view()
# visualising the data
f1lifespan %>% 
  ggplot(aes(x = as.factor(plate), y = longevity, fill = treatment))+geom_boxplot()
# making a linear model 
f1longevityls4 <- lm(longevity ~ treatment + as.factor(plate)
                       + treatment:as.factor(plate), data = f1lifespan)
# Using summary function to look at data  
summary(f1longevityls4)
broom::tidy(f1longevityls4)
# testing for significance of the interaction effect 
drop1(f1longevityls4, test = "F")
# using performance to check for normality 
performance::check_model(f1longevityls4, 
                         check = c("homogeneity",
                                   "qq"))
## Test for f1 longevity based on parental treatment and replicate
f1lifespan %>% group_by(parental_treatment, plate) %>% summarise(n = n()) %>% view()
# visualising the data 
f1lifespan %>% 
  ggplot(aes(x = as.factor(plate), y = longevity, fill = parental_treatment))+geom_boxplot()
# making a linear model
f1longevityls5 <- lm(longevity ~ parental_treatment + as.factor(plate)
                     + parental_treatment:as.factor(plate), data = f1lifespan)
# summarising the data 
summary(f1longevityls5)
broom::tidy(f1longevityls5)
drop1(f1longevityls5, test = "F")
performance::check_model(f1longevityls5, 
                         check = c("homogeneity",
                                   "qq"))


