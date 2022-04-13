library(tidyverse)
library(readxl)
library(reshap)
library(ggplot2)
library(dplyr)
library(GGally)
library(emmeans)
library(kableExtra)
library(performance)

# Importing the sheets
f0lifespan <- (read_excel(path = "Data/elegans.xlsx", sheet = "lifespan_f0", na = "NA"))
f1lifespan <- (read_excel(path = "Data/elegans.xlsx", sheet = "f1_lifespan", na = "NA"))
f0reproduction <- (read_excel(path = "Data/elegans.xlsx", sheet = "reproduction_f0", na = "NA"))
f1reproduction <- (read_excel(path = "Data/elegans.xlsx", sheet = "reproduction_f1", na = "NA"))

# Lowercase the lifespan data so it can be more easily grouped 
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
                                                    

#### LIFESPAN OF F0 
# Created a data frame for f0 lifespan
f0l <- na.omit(f0lifespan)
rnai=f0l$rnai
treatment=f0l$treatment
longevity=f0l$longevity
dataf0l=data.frame(rnai, treatment ,  longevity)

# Created a Box plot showing the effect of treatment and genes on longevity
ggplot(dataf0l, aes(x=rnai, y=longevity, fill=treatment))+ 
  geom_boxplot()+
  labs('title' = 'Effect of treatment and genes on longevity',
       y = 'Longevity',
       x = 'RNAi treatment') 

# Looking for a mean and SD for f0 lifespan depending on rnai treatment
f0lifespan_summary1 <- f0l %>%
  group_by(rnai) %>%
  summarise(mean = mean(longevity),
            sd=sd(longevity))

# Making a table of f0 lifespan (mean and SD) depending on rnai treatment 
f0lifespan_summary1 %>%
  kbl(caption="The mean and sd offspring from f0 when treated with ev or raga rnai") %>% 
  kable_styling(bootstrap_options = "striped", full_width = T, position = "left")

# Looking at models for longevity of f0 and what treatment they had
lsmodel4 <- lm(longevity ~ treatment + factor(treatment), data = f0l )
anova(lsmodel4)
lsmodel4





#### REPRODUCTION OF F0

# Data frame for f0 reproduction 
f0r <- na.omit(f0reproduction)
rnaif0r=f0r$rnai
treatmentf0r=f0r$treatment
offspringf0r=f0r$offspring
dataf0r=data.frame(rnaif0r, treatmentf0r ,  offspringf0r)

# Created f0r to emit NA values completley
f0reproduction_summary1 <- f0r %>%
  group_by(rnai) %>%
  summarise(mean = mean(offspring),
            sd=sd(offspring))

# Looking at the mean offspring (and SD value) f0 gen have depending on rnai treatment
f0reproduction_summary1 %>%
  kbl(caption="The mean and sd offspring from f0 when treated with ev or raga rnai") %>% 
  kable_styling(bootstrap_options = "striped", full_width = T, position = "left")

## Experimenting with looking at models for f0 generation  
# Model for offspring that f0 generation have vs their rnai treatment
lm(offspring ~ rnai + factor(rnai), data = f0r) %>%
  broom::tidy(., conf.int=T) %>% 
  slice(1:2)
# Model for offspring that f0 generation have vs their rnai treatmenr
lsmodel2 <- lm(offspring ~ rnai, data = f0r)
lsmodel2
anova(lsmodel2)





#### LIFESPAN OF F1 
# Looking at models for longevity of offspring from f0 vs parental rnai treatment
lm(longevity ~ parental_rnai + factor(parental_rnai), data = f1lifespan ) %>%
  broom::tidy(., conf.int=T) %>% 
  slice(1:2)

# Looking at f1 longevity vs parental rnai treatment 
lsmodel3 <- lm(longevity ~ parental_rnai + factor(parental_rnai), data = f1lifespan )
anova(lsmodel3)

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




#### REPRODUCTION OF F1
# Looking at models for offspring for f1 vs what rnai treatment their parents had 
lsmodel5 <- lm(offsprings ~ parental_rnai + factor(parental_rnai), data = f1reproduction )
anova(lsmodel5)

lm(offsprings ~ parental_rnai + factor(parental_rnai), data = f1reproduction )%>%
broom::tidy(., conf.int=T) %>% 
  slice(1:2)

# Created a plot looking at amount of offspring the f1 generation can have vs what treatments they had 
ggplot(f1reproduction, aes(x=parental_rnai, y=offsprings, fill=treatment))+ 
  geom_boxplot()+
  labs('title' = 'Effect of treatment and genes on amount of offspring for f1',
       y = 'Offspring',
       x = 'RNAi treatment') 

# "full model" of f1 reproduction 
f1reproductionls1 <- lm(offsprings ~ parental_rnai + parental_treatment + 
                          parental_rnai:parental_treatment, data = f1reproduction)
f1reproductionls1

performance::check_model(f1reproductionls1)
