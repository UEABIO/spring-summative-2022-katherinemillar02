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


# Importing the sheets
f0lifespan <- (read_excel(path = "Data/elegans.xlsx", sheet = "lifespan_f0", na = "NA"))
f1lifespan <- (read_excel(path = "Data/elegans.xlsx", sheet = "f1_lifespan", na = "NA"))
f0reproduction <- (read_excel(path = "Data/elegans.xlsx", sheet = "reproduction_f0", na = "NA"))
f1reproduction <- (read_excel(path = "Data/elegans.xlsx", sheet = "reproduction_f1", na = "NA"))





# lowercased everything 
f0lifespan$rnai <- tolower(f0lifespan$rnai)
f0lifespan$treatment <- tolower(f0lifespan$treatment)
f1lifespan$parental_rnai <- tolower(f1lifespan$parental_rnai)
f1lifespan$parental_treatment <- tolower(f1lifespan$parental_treatment)
f1lifespan$treatment <- tolower(f1lifespan$treatment) 


# Calculating the longevity for f0 and f1 
f0lifespan$longevity <- as.numeric (difftime
                                     (f0lifespan$death_date, f0lifespan$set_up_date,
                                       units="days"))

f1lifespan$longevity <- as.numeric (difftime
                                     (f1lifespan$death_date, f1lifespan$set_up_date,
                                       units="days"))
                                                    

#### LOOKING AT LIFESPAN OF FO 
# Looking for a mean and SD for f0 lifespan depending on rnai treatment
f0lifespan_summary <- f0lifespan %>%
  group_by(rnai) %>%
  summarise(mean = mean(longevity),
            sd=sd(longevity))



# Making a table of f0 lifespan (mean and SD) depending on rnai treatment 
f0lifespan_summary %>%
  kbl(caption="The mean and sd offspring from f0 when treated with ev or raga rnai") %>% 
  kable_styling(bootstrap_options = "striped", full_width = T, position = "left")


## FIGURES/ VISUALING DATA for F0 LIFESPAN
# Created a Box plot showing the effect of treatment and genes on longevity
f0lplot <- ggplot(f0lifespan, aes(x=rnai, y=longevity, fill=treatment))+
  geom_boxplot()+
  labs('title' = ' F0 - Effect of treatment and genes on longevity',
       y = 'Longevity',
       x = 'RNAi treatment') 

f0lplot



#### REPRODUCTION OF F0

# Created a summary for rnai treatment and offspring for f0
 f0reproduction_summary <- f0reproduction %>%
  group_by(rnai) %>%
  summarise(mean = mean(offspring),
            sd=sd(offspring))

# Looking at the mean offspring (and SD value) f0 gen have depending on rnai treatment
f0reproduction_summary %>%
  kbl(caption="The mean and sd offspring from f0 when treated with ev or raga rnai") %>% 
  kable_styling(bootstrap_options = "striped", full_width = T, position = "left")



#### LIFESPAN OF F1

# FIGURES/ VISUALISING DATA FOR F1 LIFESPAN
# Created a boxplot comparing longevity of offspring and parental treatment 

ggplot(data=f1lifespan, aes(x = tolower(parental_treatment), y = longevity)) +
  geom_boxplot(aes(fill = longevity),
               alpha = 0.2, 
               width = 0.5, 
               outlier.shape=NA)+
  geom_jitter(aes(colour = longevity),
              width=0.2)+
  theme(legend.position = "none")+
  labs('title' = 'Effect of parental treatment on
       offspring longevity',
       x = 'Parental Treatment',
       y = 'Longevity')

f1lplot <- ggplot(f1lifespan, aes(x=parental_rnai, y=longevity, fill=treatment))+ 
  geom_boxplot()+
  labs('title' = ' F1 - Effect of treatment and genes on longevity',
       y = 'Longevity',
       x = 'RNAi treatment') 


#### REPRODUCTION OF F1

# FIGURES/ VISUALISING DATA FOR F1 REPRODUCTION 
# Created a plot looking at amount of OFFSPRING the f1 generation can have vs what treatments they had 
ggplot(f1reproduction, aes(x=parental_rnai, y=offsprings, fill=treatment))+
  geom_boxplot()+
  labs('title' = 'Effect of treatment and genes on amount of offspring for f1',
       y = 'Offspring',
       x = 'RNAi treatment') 

##### COMBINING
# Working with patchwork - looking at longevity of f0 and f1
f0lplot + f1lplot



# summary of models have worked with 

#  model 1 - f0 lifespan based on THEIR rnai gene and treatment 
     #  f0 - longevity, rnai, treatment 
f0ls1 <- lm(longevity ~ rnai + treatment + rnai + rnai:treatment, data = f0lifespan)
broom::tidy(f0ls1)
summary(f0ls1)
f0tidymodel <- broom::tidy(f0ls1) 
f0tidymodel



library(kableExtra)

model2table <- 
  f0ls1 %>% broom::tidy(conf.int = T) %>% 
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
  
model2table



# CI for paired T test for longevity with rnai and light/dark treatment 
lm(longevity ~ rnai + factor(treatment), data = f0lifespan) %>% 
  broom::tidy(., conf.int=T) %>% 
  slice(1:2)
f0tidymodel[[2,2]] / f0tidymodel[[2,3]] # this is the same value as rnai raga with ev statistic

performance::check_model(f0ls1)
performance::check_model(f0ls1, check="homogeneity")

#  keep interaction term? 
drop1(f0ls1, test = "F")

 # model 2 - f0 longevity based on whether they were in light/dark
             # f0 - longevity, treatment 
f0longevityandtreatment <- lm(longevity ~ treatment, data = f0lifespan )

f0longevityandtreatment <- lm(longevity ~ treatment, data = f0lifespan) %>%
  broom::tidy(., conf.int=T) %>% 
  slice(1:2)
anova(f0longevityandtreatment)
f0longevityandtreatment

model1table <- 
  f0longevityandtreatment %>% broom::tidy(conf.int = T) %>% 
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

model1table

# performance function - checked for normality?? 
performance::check_model(f0longevityandtreatment, check=c("normality","qq"))
# Looked at emmeans data for amount of offspring f0 generation have vs rnai treatment
f0means <- emmeans::emmeans(f0longevityandtreatment, specs = ~treatment)
f0means
# Look at emmeans visually 
f0means %>% 
  as_tibble() %>% 
  ggplot(aes(x=treatment, 
             y=emmean))+
  geom_pointrange(aes(
    ymin=lower.CL, 
    ymax=upper.CL))

broom::tidy(f0longevityandtreatment)


performance::check_model(f0longevityandtreatment, check="homogeneity")
performance::check_model(f0longevityandtreatment)

# model 3 - f0 amount of offspring based on rnai gene and treatment 
                      # f0 - offspring, rnai, treatment 
f0rnaiosmodel <- lm(offspring ~ rnai + factor(treatment), data = f0reproduction) 
f0rnaiosmodel2 <- lm(offspring ~ rnai + treatment, data = f0reproduction)

%>%
  broom::tidy(., conf.int=T) %>% 
  slice(1:2)
broom::tidy(f0rnaiosmodel2)
summary(f0rnaiosmodel)

broom::tidy(f0rnaiosmodel) 
performance::check_model(f0rnaiosmodel)

model3table <- 
  f0rnaiosmodel %>% broom::tidy(conf.int = T) %>% 
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

model3table

 # model 4 - f1 longevity based on their parent's rnai and treatment
                  # f1 - longevity, parent's rnai, parent's treatment 
f1longptreatmodel <- lm(longevity ~ parental_rnai + factor(parental_rnai), data = f1lifespan )

lm(longevity ~ parental_rnai + factor(parental_rnai), data = f1lifespan ) %>%
  broom::tidy(., conf.int=T) %>% 
  slice(1:2)


anova(f1longptreatmodel)

f1model <-  lm(longevity ~ parental_rnai + parental_treatment, data = f1lifespan)
summary(f1model)

lm(longevity ~ parental_rnai + parental_treatment, data = f1lifespan ) %>%
  broom::tidy(., conf.int=T) %>% 
  slice(1:2)

model4table <- 
  f1model %>% broom::tidy(conf.int = T) %>% 
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

model4table

#  keep interaction term? 
drop1(f1model, test = "F")
#  Don't keep interaction term? 

 # model 5 - f1 offspring vs their parents rnai and treatment 
            # f1 - offspring, parent's rnai, parent's treatment
f1reproductionls1 <- lm(offsprings ~ parental_rnai + parental_treatment + 
                          parental_rnai:parental_treatment, data = f1reproduction)
f1reproductionls1
performance::check_model(f1reproductionls1)

f1osptreatment <- lm(offsprings ~ parental_rnai + factor(parental_rnai), data = f1reproduction )
anova(f1osptreatment)

lm(offsprings ~ parental_rnai + factor(parental_rnai), data = f1reproduction )%>%
  broom::tidy(., conf.int=T) %>% 
  slice(1:2)




performance::check_model(f1reproductionls1)

#  keep interaction term? 
drop1(f1reproductionls1, test = "F")
# Don't keep interaction term? 

# model 6 - f1 longevity based on their treatment and their parents treatment 
       # f1 - longevity, parental treatment and treatment 
f1treatptreatment <- lm(longevity ~ parental_treatment + factor(treatment), data = f1lifespan)
anova(f1treatptreatment)

broom::tidy(f1treatptreatment)

summary(f1treatptreatment)

performance::check_model(f1treatptreatment, check="homogeneity")
performance::check_model(f1treatptreatment, check="normality")





# testing to see if replicate/plate values are the same throughout 

f1reproduction %>% group_by(treatment, replicate) %>% summarise(n = n()) %>% view()
f1reproduction %>% 
  ggplot(aes(x = as.factor(replicate), y = offsprings, fill = treatment))+geom_boxplot()
f1replicatetest <- lm(offsprings ~ treatment + as.factor(replicate)
                      + treatment:as.factor(replicate), data = f1reproduction)
summary(f1replicatetest)
broom::tidy(f1replicatetest)
drop1(f1replicatetest, test = "F")
performance::check_model(f1replicatetest, 
                         check = c("homogeneity",
                                   "qq"))


f0reproduction %>% group_by(treatment, replicate) %>% summarise(n = n()) %>% view()
f0reproduction %>% 
  ggplot(aes(x = as.factor(replicate), y = offspring, fill = treatment))+geom_boxplot()

f0replicatetest <- lm(offspring ~ treatment + as.factor(replicate)
                      + treatment:as.factor(replicate), data = f0reproduction)
summary(f0replicatetest)
broom::tidy(f0replicatetest)
drop1(f0replicatetest, test = "F")
performance::check_model(f0replicatetest, 
                         check = c("homogeneity",
                                   "qq"))


f0lifespan %>% group_by(treatment, plate) %>% summarise(n = n()) %>% view()
f0lifespan %>% 
  ggplot(aes(x = as.factor(plate), y = longevity, fill = treatment))+geom_boxplot()
f0replicatetest2 <- lm(longevity ~ treatment + as.factor(plate)
                      + treatment:as.factor(plate), data = f0lifespan)
summary(f0replicatetest2)
broom::tidy(f0replicatetest2)
drop1(f0replicatetest2, test = "F")
performance::check_model(f0replicatetest2, 
                         check = c("homogeneity",
                                   "qq"))



f1lifespan %>% group_by(treatment, plate) %>% summarise(n = n()) %>% view()
f1lifespan %>% 
  ggplot(aes(x = as.factor(plate), y = longevity, fill = treatment))+geom_boxplot()

f1replicatetest2 <- lm(longevity ~ treatment + as.factor(plate)
                       + treatment:as.factor(plate), data = f1lifespan)
summary(f1replicatetest2)
broom::tidy(f1replicatetest2)
drop1(f1replicatetest2, test = "F")
performance::check_model(f1replicatetest2, 
                         check = c("homogeneity",
                                   "qq"))



replicate <- read_csv("Data/replicatedata.csv")
