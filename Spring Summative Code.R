library(tidyverse)
library(readxl)

elegansdata <- read_excel("Data/elegans.xlsx")
readxl::read_excel("Data/elegans.xlsx")
excel_sheets("Data/elegans.xlsx")
read_excel(path = "Data/elegans.xlsx", sheet = "reproduction_f0")
read_excel(path = "Data/elegans.xlsx", sheet = "reproduction_f1")

head(read_excel(path = "Data/elegans.xlsx", sheet = "reproduction_f1"))
head(read_excel(path = "Data/elegans.xlsx", sheet = "lifespan_f0"))








