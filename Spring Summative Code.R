library(tidyverse)
library(readxl)
library(reshap)
library(ggplot2)
library(dplyr)

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

# Create a data frame and boxplot comparing longevity with treatment and genes
f0l <- na.omit(f0lifespan) 
rnai=f0l$rnai
treatment=f0l$treatment
longevity=f0l$longevity
data=data.frame(rnai, treatment ,  longevity)
ggplot(data, aes(x=rnai, y=longevity, fill=treatment)) + 
  geom_boxplot()+
  labs('title' = 'Effect of treatment and genes on longevity',
       y = 'Longevity',
       x = 'RNAi treatment')
