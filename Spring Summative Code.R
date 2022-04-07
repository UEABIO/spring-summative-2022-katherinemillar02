library(tidyverse)
library(readxl)


# Importing the sheets
f0lifespan <- (read_excel(path = "Data/elegans.xlsx", sheet = "lifespan_f0", na = "NA"))
f1lifespan <- (read_excel(path = "Data/elegans.xlsx", sheet = "f1_lifespan", na = "NA"))
f0reproduction <- (read_excel(path = "Data/elegans.xlsx", sheet = "reproduction_f0", na = "NA"))
f1reproduction <- (read_excel(path = "Data/elegans.xlsx", sheet = "reproduction_f1", na = "NA"))


# Calculating the longevity for f0/1
f0lifespan$longevity <- as.numeric (difftime
                                     (f0lifespan$death_date,f0lifespan$set_up_date,
                                       units="days"))

f1lifespan$longevity <- as.numeric (difftime
                                     (f1lifespan$death_date,f1lifespan$set_up_date,
                                       units="days"))






                             
                                                                                    


