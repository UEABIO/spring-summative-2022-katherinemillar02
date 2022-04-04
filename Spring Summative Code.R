library(tidyverse)
library(readxl)

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


# Visualising the data 
(read_excel(path = "Data/elegans.xlsx", sheet = "f1_lifespan")) %>% 
  ggplot(aes(x=treatment,
             y=plate))+
  geom_point()
     



