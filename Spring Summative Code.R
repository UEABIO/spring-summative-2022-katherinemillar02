library(tidyverse)
library(readxl)
library(reshap)
library(ggplot2)
library(dplyr)
library(GGally)
library(emmeans)
library(kableExtra)

# Importing the sheets
f0lifespan <- (read_excel(path = "Data/elegans.xlsx", sheet = "lifespan_f0", na = "NA"))
f1lifespan <- (read_excel(path = "Data/elegans.xlsx", sheet = "f1_lifespan", na = "NA"))
f0reproduction <- (read_excel(path = "Data/elegans.xlsx", sheet = "reproduction_f0", na = "NA"))
f1reproduction <- (read_excel(path = "Data/elegans.xlsx", sheet = "reproduction_f1", na = "NA"))

# Lowercase the data so it can be more easily grouped 
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
                                                    
# Created a boxplot comparing longevity of offspring and parental treatment 
ggplot(data=f1lifespan, aes(x = tolower(parental_treatment), y = longevity)) +
  geom_boxplot(aes(fill = longevity),
               alpha = 0.2, 
               width = 0.5, 
               outlier.shape=NA)+
  geom_jitter(aes(colour = longevity),
              width=0.2)+
  theme(legend.position = "none")+
  labs('title' = 'Effect of parental treatment on offspring longevity',
       x = 'Parental Treatment',
       y = 'Longevity')

# Created data frames
f0l <- na.omit(f0lifespan)
rnai=f0l$rnai
treatment=f0l$treatment
longevity=f0l$longevity
data=data.frame(rnai, treatment ,  longevity)

# Created a Box plot showing the effect of treatment and genes on longevity
ggplot(data, aes(x=rnai, y=longevity, fill=treatment))+ 
  geom_boxplot()+
  labs('title' = 'Effect of treatment and genes on longevity',
       y = 'Longevity',
       x = 'RNAi treatment') 

# The parents lifespan
f0lifespan_summary1 <- f0l %>%
  group_by(rnai) %>%
  summarise(mean = mean(longevity),
            sd=sd(longevity))

f0lifespan_summary1 %>%
  kbl(caption="The mean and sd offspring from f0 when treated with ev or raga rnai") %>% 
  kable_styling(bootstrap_options = "striped", full_width = T, position = "left")

# Trying to look at rnai treatment and offspring in a table
#f1reproduction_wide <- f1reproduction %>%
#pivot_wider(names_from = parental_rnai, values_from = offsprings, id_cols = day) %>%
# mutate(difference=ev-raga)

  #f0reproduction_wide <- f0reproduction %>%
# pivot_wider(names_from = rnai, values_from = offspring, id_cols = offspring) %>%
# mutate(difference=raga-ev)

# Data frame for f0 reproduction 
f0r <- na.omit(f0reproduction)
rnaif0r=f0r$rnai
treatmentf0r=f0r$treatment
offspringf0r=f0r$offspring
data=data.frame(rnaif0r, treatmentf0r ,  offspringf0r)

# Trying to look at rnai treatment and offspring in a table

# Created f0r to emit NA values completley
f0reproduction_summary1 <- f0r %>%
  group_by(rnai) %>%
  summarise(mean = mean(offspring),
            sd=sd(offspring))

f0reproduction_summary1 %>%
  kbl(caption="The mean and sd offspring from f0 when treated with ev or raga rnai") %>% 
  kable_styling(bootstrap_options = "striped", full_width = T, position = "left")

# Beginning to look at models 

# Offspring
lsmodel0 <- lm(formula = offspring~1, data = f0r)
lsmodel0
broom::tidy(lsmodel0)

# Looking at statistics for offspring vs rnai treatment
lm(offspring ~ rnai + factor(rnai), data = f0r) %>%
  broom::tidy(., conf.int=T) %>% 
  slice(1:2)

lsmodel2 <- lm(offspring ~ rnai + factor(rnai), data = f0r)
anova(lsmodel2)

# Looking at statistics for longevity of offspring vs parental rnai treatment
lm(longevity ~ parental_rnai + factor(parental_rnai), data = f1lifespan ) %>%
  broom::tidy(., conf.int=T) %>% 
  slice(1:2)

lsmodel3 <- lm(longevity ~ parental_rnai + factor(parental_rnai), data = f1lifespan )
anova(lsmodel3)

str(f0l)
str(f0r)
str(f1reproduction)
str(f1lifespan)

# Looking at models for longevity of f0 and what treatment they had
lsmodel4 <- lm(longevity ~ treatment + factor(treatment), data = f0l )
anova(lsmodel4)
lsmodel4

# Looking at models for offspring vs what rnai treatment their parents had 
lsmodel5 <- lm(offsprings ~ parental_rnai + factor(parental_rnai), data = f1reproduction )
anova(lsmodel5)

lm(offsprings ~ parental_rnai + factor(parental_rnai), data = f1reproduction )%>%
broom::tidy(., conf.int=T) %>% 
  slice(1:2)




