library(tidyverse)


celegans <- read.csv("Data/elegans_data.csv")

glimpse(celegans)
head(celegans)
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



