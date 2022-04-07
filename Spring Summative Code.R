library(tidyverse)
library(readxl)
library(kableExtra)

#Read Data in
elegansdata <- read_excel("Data/elegans.xlsx")
readxl::read_excel("Data/elegans.xlsx")
excel_sheets("Data/elegans.xlsx")

#The different sheets 
read_excel(path = "Data/elegans.xlsx", sheet = "reproduction_f0")
read_excel(path = "Data/elegans.xlsx", sheet = "reproduction_f1")
read_excel(path = "Data/elegans.xlsx", sheet = "lifespan_f0")
read_excel(path = "Data/elegans.xlsx", sheet = "f1_lifespan")

#Checking the data
summary(read_excel(path = "Data/elegans.xlsx", sheet = "reproduction_f1"))
head(read_excel(path = "Data/elegans.xlsx", sheet = "lifespan_f0"))
head(read_excel(path = "Data/elegans.xlsx", sheet = "f1_lifespan"))
head(read_excel(path = "Data/elegans.xlsx", sheet = "reproduction_f0"))
glimpse(read_excel(path = "Data/elegans.xlsx", sheet = "f1_lifespan"))
summary(read_excel(path = "Data/elegans.xlsx", sheet = "f1_lifespan"))
head(read_excel(path = "Data/elegans.xlsx", sheet = "f1_lifespan"))

# Visualising the data 
(read_excel(path = "Data/elegans.xlsx", sheet = "f1_lifespan")) %>% 
  ggplot(aes(x=treatment,
             y=set_up_date))+
  geom_point()

read_excel(path = "Data/elegans.xlsx", sheet = "reproduction_f0") %>% 
  ggplot(aes(x=treatment,
             y=offspring))+
  geom_col()

read_excel(path = "Data/elegans.xlsx", sheet = "reproduction_f1") %>%
ggplot(aes(x= parental_treatment,
           y=offsprings))+
  geom_col()

# Attempt to look at days living for f1/f0 c.elegans
(read_excel(path = "Data/elegans.xlsx", sheet = "f1_lifespan"))

f1lifespan_wide <- (read_excel(path = "Data/elegans.xlsx", sheet = "f1_lifespan")) %>% 
  pivot_wider(names_from=worm, values_from=worm, id_cols=worm) %>% 
  mutate(difference=death_date-set_up_date)

difftime(set_up_date, death_date, days)
     





