library(readxl)
library(tidyverse)

d <- read_excel("./Processed Data/WatershedProductivity.xlsx", sheet = "ggplot")

ggplot(d) +
  geom_line(aes(x = Date, y = NDVI, col = Watershed))
