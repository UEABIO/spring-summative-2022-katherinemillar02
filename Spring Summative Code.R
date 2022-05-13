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
                                                    
# Looking for a mean and SD for f0 lifespan depending on rnai treatment
f0lifespan_summary <- f0lifespan %>%
  group_by(rnai) %>%
  summarise(mean = mean(longevity),
            sd=sd(longevity))

# Making a table of f0 lifespan (mean and SD) depending on rnai treatment 
f0lifespan_summary %>%
  kbl(caption="The mean and sd offspring from f0 when treated with ev or raga rnai") %>% 
  kable_styling(bootstrap_options = "striped", full_width = T, position = "left")


# Created a summary for rnai treatment and offspring for f0
 f0reproduction_summary <- f0reproduction %>%
  group_by(rnai) %>%
  summarise(mean = mean(offspring),
            sd=sd(offspring))

# Looking at the mean offspring (and SD value) F0-generation have depending on rnai treatment
f0reproduction_summary %>%
  kbl(caption="The mean and sd offspring from f0 when treated with ev or raga rnai") %>% 
  kable_styling(bootstrap_options = "striped", full_width = T, position = "left")


### Testing and working with different models 

# MODEL 1 
# Gene knockdown AND light/dark treatment 
# F0-lifespan based on THEIR rnai gene and light/dark treatment 
     #  F0 - longevity, rnai, treatment 

# Visualising the data with a box plot 
ggplot(f0lifespan, aes(x=rnai, y=longevity, fill=treatment))+
  geom_boxplot()+
  labs('title' = ' F0 - Effect of treatment and genes on longevity',
       y = 'Longevity',
       x = 'RNAi treatment') 
  

# Creating a linear model 
f0lifespanls1 <- lm(longevity ~ rnai + treatment + rnai + rnai:treatment, data = f0lifespan)

# Using summary and broom tidy to give a summary of f0lifespanls1 results 
summary(f0lifespanls1)
broom::tidy(f0lifespanls1) 

# Using drop1 function to see if rnai:treatment interaction term should be kept  
drop1(f0lifespanls1, test = "F")
# Keep interaction term - there is a significant difference 
# Use this as the final model 

# Assumption checking 
# Doing a performance check to look for normality in the model
performance::check_model(f0lifespanls1)

performance::check_model(f0lifespanls1, check="homogeneity") 
# Doesn't look normal 

performance::check_model(f0lifespanls1, check=c("normality","qq"))
# Doesn't look normal 

performance::check_model(f0lifespanls1, check="outliers")
# 

# Coefficients 
coef(f0lifespanls1)

# Data Transformations - 
MASS::boxcox(f0lifespanls1)

# 

# Making a table of f0lifespan for the write-up based on their rnai gene and light/dark treatment 
f0lifespanls1table <- 
  f0lifespanls1 %>% broom::tidy(conf.int = T) %>% 
  select(-`std.error`) %>% 
  mutate_if(is.numeric, round, 2) %>% 
  kbl(col.names = c("Predictors",
                    "Estimates",
                    "Z-value",
                    "P",
                    "Lower 95% CI",
                    "Upper 95% CI"),
      caption = "Model 2", 
      booktabs = TRUE) %>% 
  kable_styling(full_width = FALSE, font_size=16)

f0lifespanls1table


# MODEL 2
# Gene knockdown and treatment (with no interaction effect)
# Confidence Intervals - for paired T-test for longevity with rnai and light/dark treatment 
f0lifespanls2 <- lm(longevity ~ rnai + factor(treatment), data = f0lifespan) %>% 
  broom::tidy(., conf.int=T) %>% 
  slice(1:2)
# MAYBE DO NOT KEEP THIS AS INTERACTION EFFECT SEEMS SIGNIFICANT 
# MAYBE DO NOT KEEP THIS AS YOU HAVE ALREADY INCLUDED ONE WITH BOTH SO MAY NOT NEED 



# MODEL 3 
#  F0 longevity based on whether they were in light/dark
 # F0 - longevity, treatment 

# Visualising the data 
ggplot(f0lifespan, aes(x=treatment, y=longevity, fill=treatment))+ 
  geom_boxplot()

# Creating linear model 
f0lifespanls3 <- lm(longevity ~ treatment, data = f0lifespan)

# Using broom::tidy and summary function 
broom::tidy(f0lifespanls3)
summary(f0lifespanls3)

# performance function - checked for normality
performance::check_model(f0lifespanls3)

performance::check_model(f0lifespanls3, check=c("normality","qq"))
# Doesn't look great

performance::check_model(f0lifespanls3, check="homogenity")
#  ? 


# Did a transformation test
MASS::boxcox(f0lifespanls3)

# Doing an anova test with original data 
anova(f0lifespanls3)

 lm(longevity ~ treatment, data = f0lifespan) %>%
  broom::tidy(., conf.int=T) %>% 
  slice(1:2)

f0lifespanls3table <- 
  f0lifespanls3 %>% broom::tidy(conf.int = T) %>% 
  select(-`std.error`) %>% 
  mutate_if(is.numeric, round, 2) %>% 
  kbl(col.names = c("Predictors",
                    "Estimates",
                    "Z-value",
                    "P",
                    "Lower 95% CI",
                    "Upper 95% CI"),
      caption = "Model 1", 
      booktabs = TRUE) %>% 
  kable_styling(full_width = FALSE, font_size=16)

f0lifespanls3table



# Looked at emmeans data for amount of offspring f0 generation have vs rnai treatment
f0offspringmeans <- emmeans::emmeans(f0lifespanls3, specs = ~treatment)
f0offspringmeans
# Visualised the amount of offspring FO have based on rnai treatment using emmeans data
f0offspringmeans %>% 
  as_tibble() %>% 
  ggplot(aes(x=treatment, 
             y=emmean))+
  geom_pointrange(aes(
    ymin=lower.CL, 
    ymax=upper.CL))




# MODEL 4 
# F0-generation offspring against rnai gene and treatment 
# F0 - offspring, rnai, treatment 

# Visualising the data 
ggplot(f0reproduction, aes(x=rnai, y=offspring, fill=treatment))+
  geom_boxplot()
 
# Creating a linear model 
f0reproductionls2 <- lm(offspring ~ rnai + treatment + rnai:treatment, data = f0reproduction)

# Using drop1 to see if interaction term should be kept
drop1(f0reproductionls2, test = "F")
#  


# Testing with confidence intervals 
lm(offspring ~ rnai + treatment + rnai:treatment, data = f0reproduction)%>%
  broom::tidy(., conf.int=T) %>% 
  slice(1:2)

# Using broom::tidy and summary to create a summary
broom::tidy(f0reproductionls2)
summary(f0reproductionls2)

# Using performance check to check for normality 
performance::check_model(f0reproductionls2)

performance::check_model(f0reproductionls2, check=c("normality","qq"))
# No 


# Transforming data 
MASS::boxcox(f0reproductionls2)



# Creating a table for the write-up
f0reproductionls1table <- 
  f0reproductionls1 %>% broom::tidy(conf.int = T) %>% 
  select(-`std.error`) %>% 
  mutate_if(is.numeric, round, 2) %>% 
  kbl(col.names = c("Predictors",
                    "Estimates",
                    "Z-value",
                    "P",
                    "Lower 95% CI",
                    "Upper 95% CI"),
      caption = "Model 3", 
      booktabs = TRUE) %>% 
  kable_styling(full_width = FALSE, font_size=16)

f0reproductionls1table




# MODEL 5 
#  F1-generation longevity based on their parent's rnai
  # f1 - longevity, parent's rnai

# Visualise the data 
ggplot(f1lifespan, aes(x=parental_rnai, y=longevity, fill=longevity))+
  geom_boxplot()

# creating a linear model
f1longevityls1 <- lm(longevity ~ parental_rnai, data = f1lifespan )

lm(longevity ~ parental_rnai, data = f1lifespan ) %>%
  broom::tidy(., conf.int=T) %>% 
  slice(1:2)

# Checking the data 
performance::check_model(f1longevityls1)
# Doesn't look good 

performance::check_model(f1longevityls1, check = "homogenity")



# Doing an anova 
anova(f1longevityls1)



# MODEL 6 
# F1 - generation longevity based on their parent's rnai and parent's treatment
# F1 - longevity, parent's rnai, parent's treatment 

# Visualising the data 
ggplot(data=f1lifespan, aes(x = tolower(parental_treatment), y = longevity)) +
  geom_boxplot(aes(fill = parental_rnai),
               alpha = 0.2, 
               width = 0.5, 
               outlier.shape=NA)+
  geom_jitter(aes(colour = longevity),
              width=0.2)+
  theme(legend.position = "none")+
  labs('title' = 'Effect of parental treatment on offspring longevity',
       x = 'Parental Treatment',
       y = 'Longevity of offspring')

# Creating a model 
f1longevityls2 <-  lm(longevity ~ parental_rnai + parental_treatment + parental_rnai:parental_treatment, data = f1lifespan)

lm(longevity ~ parental_rnai + parental_treatment + parental_rnai:parental_treatment, data = f1lifespan ) %>%
  broom::tidy(., conf.int=T) %>% 
  slice(1:2)

#  keep interaction term? 
drop1(f1longevityls2, test = "F")
#  Don't keep interaction term, no significance 

# New model
f1longevityls6 <-  lm(longevity ~ parental_rnai + parental_treatment, data = f1lifespan)

lm(longevity ~ parental_rnai + parental_treatment, data = f1lifespan ) %>%
  broom::tidy(., conf.int=T) %>% 
  slice(1:2)

# Using broom::tidy as a summary function 
broom::tidy(f1longevityls6)
summary(f1longevityls6)

# Assumption checking 
performance::check_model(f1longevityls6)
# Looks normal at first glance 


f1longevityls6table <- 
  f1longevityls6 %>% broom::tidy(conf.int = T) %>% 
  select(-`std.error`) %>% 
  mutate_if(is.numeric, round, 2) %>% 
  kbl(col.names = c("Predictors",
                    "Estimates",
                    "Z-value",
                    "P",
                    "Lower 95% CI",
                    "Upper 95% CI"),
      caption = "Model 4", 
      booktabs = TRUE) %>% 
  kable_styling(full_width = FALSE, font_size=16)

f1longevityls6table





# MODEL 7 
# F1 offsprings vs their parent's rnai treatment 
#  F1 offsprings, parent's rnai 

# Visualising the data 
 ggplot(f1reproduction, aes(x=parental_rnai, y=offsprings, fill=parental_rnai))+
  geom_boxplot()+
  labs('title' = "F1 offspring and their parent's rnai"
       "subtitle" = "the amount of offspring F1 generation have and what rnai gene their parents had",
       y = 'Offspring',
       x = 'RNAi treatment')

# Creating the model 
f1reproductionls4 <- lm(offsprings ~ parental_rnai, data = f1reproduction )

lm(offsprings ~ parental_rnai, data = f1reproduction )%>%
  broom::tidy(., conf.int=T) %>% 
  slice(1:2)

# Assumption Checking 
performance::check_model(f1reproductionls4)

# Homogenity not okay 

# Transformation with BoxCox 
# run this, pick a transformation and retest the model fit
MASS::boxcox(f1reproductionls4)

f1reproductionls4table <- 
  f1reproductionls4 %>% broom::tidy(conf.int = T) %>% 
  select(-`std.error`) %>% 
  mutate_if(is.numeric, round, 2) %>% 
  kbl(col.names = c("Predictors",
                    "Estimates",
                    "Z-value",
                    "P",
                    "Lower 95% CI",
                    "Upper 95% CI"),
      caption = "Model 4", 
      booktabs = TRUE) %>% 
  kable_styling(full_width = FALSE, font_size=16)

f1reproductionls4table



# MODEL 8 
# f1 longevity based on their treatment and their parents treatment 
 # f1 - longevity, parental treatment and their own treatment, with interaction effect  

# Visualise the data 
ggplot(f1lifespan, aes(x=parental_treatment, y=longevity, fill=treatment))+
  geom_boxplot()

f1longevityls3 <- lm(longevity ~ parental_treatment + treatment + parental_treatment:treatment, data = f1lifespan)


# Summary of model 
broom::tidy(f1longevityls3)
summary(f1longevityls3)

# Using drop1 - keep interaction effect?
drop1(f1longevityls3, test = "F")
#

# Assumption checking 
performance::check_model(f1longevityls3, check="homogeneity")
performance::check_model(f1longevityls3, check="normality")


# Box cox - Transformation 











# Testing to see if replicate/plate values are the same throughout 

## Test for f1 reproduction based on treatment and replicate
f1reproduction %>% group_by(treatment, replicate) %>% summarise(n = n()) %>% view()

f1reproduction %>% 
  ggplot(aes(x = as.factor(replicate), y = offsprings, fill = treatment))+geom_boxplot()

f1reproductionls5 <- lm(offsprings ~ treatment + as.factor(replicate)
                      + treatment:as.factor(replicate), data = f1reproduction)

summary(f1reproductionls5)
broom::tidy(f1reproductionls5)
drop1(f1reproductionls5, test = "F")
performance::check_model(f1reproductionls5, 
                         check = c("homogeneity",
                                   "qq"))

## Test for f0 reproduction based on treatment and replicate
f0reproduction %>% group_by(treatment, replicate) %>% summarise(n = n()) %>% view()
f0reproduction %>% 
  
ggplot(aes(x = as.factor(replicate), y = offspring, fill = treatment))+geom_boxplot()

f0reproductionls3 <- lm(offspring ~ treatment + as.factor(replicate)
                      + treatment:as.factor(replicate), data = f0reproduction)

f0reproductionls4 <- lm(offspring ~ treatment + as.factor(replicate), data = f0reproduction)

drop1(f0reproductionls4, test = "F")

performance::check_model(f0reproductionls4)

f0reproductionls4table <- 
  f0reproductionls4 %>% broom::tidy(conf.int = T) %>% 
  select(-`std.error`) %>% 
  mutate_if(is.numeric, round, 2) %>% 
  kbl(col.names = c("Predictors",
                    "Estimates",
                    "Z-value",
                    "P",
                    "Lower 95% CI",
                    "Upper 95% CI"),
      caption = "Model 4", 
      booktabs = TRUE) %>% 
  kable_styling(full_width = FALSE, font_size=16)

f0reproductionls4table

summary(f0reproductionls3)
broom::tidy(f0reproductionls3)


drop1(f0reproductionls3, test = "F")

performance::check_model(f0reproductionls3, 
                         check = c("homogeneity",
                                   "qq"))

## Test for f0 longevity based on treatment and replicate

f0lifespan %>% group_by(treatment, plate) %>% summarise(n = n()) %>% view()

f0lifespan %>% 
  ggplot(aes(x = as.factor(plate), y = longevity, fill = treatment))+geom_boxplot()

f0lifespanls4 <- lm(longevity ~ treatment + as.factor(plate)
                      + treatment:as.factor(plate), data = f0lifespan)
summary(f0lifespanls4)

broom::tidy(f0lifespanls4)

drop1(f0lifespanls4, test = "F")

performance::check_model(f0lifespanls4, 
                         check = c("homogeneity",
                                   "qq"))

## Test on f1 longevity based on treatment and replicate
f1lifespan %>% group_by(treatment, plate) %>% summarise(n = n()) %>% view()
f1lifespan %>% 
  ggplot(aes(x = as.factor(plate), y = longevity, fill = treatment))+geom_boxplot()

f1longevityls4 <- lm(longevity ~ treatment + as.factor(plate)
                       + treatment:as.factor(plate), data = f1lifespan)
summary(f1longevityls4)
broom::tidy(f1longevityls4)
drop1(f1longevityls4, test = "F")
performance::check_model(f1longevityls4, 
                         check = c("homogeneity",
                                   "qq"))


## Test for f1 longevity based on parental treatment and replicate
f1lifespan %>% group_by(parental_treatment, plate) %>% summarise(n = n()) %>% view()
f1lifespan %>% 
  ggplot(aes(x = as.factor(plate), y = longevity, fill = parental_treatment))+geom_boxplot()


f1longevityls5 <- lm(longevity ~ parental_treatment + as.factor(plate)
                     + parental_treatment:as.factor(plate), data = f1lifespan)

summary(f1longevityls5)
broom::tidy(f1longevityls5)
drop1(f1longevityls5, test = "F")
performance::check_model(f1longevityls5, 
                         check = c("homogeneity",
                                   "qq"))


