library(tidyverse)
library(readxl)

excel_sheets("Data/elegans.xlsx")


readxl::read_excel

celegans <- read.csv("Data/elegans_data.csv")
reproductionf0 <- read.csv("Data/reproduction_f0.csv")
lifespanf0 <- read.csv("Data/lifespan_f0.csv")
# f0 = parental generation 
reproductionf1 <- read.csv("Data/reproduction_f1.csv")
# f1 = offspring 

head(celegans)
head(reproductionf0)
colnames(celegans)
summary(celegans)

celegans %>% 
  duplicated() %>% 
  sum()

celegans %>%
summarise(min=min(death_date, na.rm=TRUE), 
          max=max(death_date, na.rm=TRUE))

celegans%>% 
  distinct(parental_rnai)

celegans%>%
  distinct(worm)

celegans%>%
  is.na()
sum()

celegans %>% 
  ggplot(aes(x=parental_rnai,
             y=death_date))+
  geom_boxplot()+
  coord_flip()



