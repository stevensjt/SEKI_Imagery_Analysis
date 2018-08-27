

library(tidyr)
library(dplyr)
library(readxl)
library(ggplot2)

####1. Forestry analyses####
####1a. Data processing####
trees <- 
  read_xlsx("./Raw Data/SugarloafForestryPlotDataCombined.xlsx", 
            sheet = "Treelist")
subplots <-
  read_xlsx("./Raw Data/SugarloafForestryPlotDataCombined.xlsx", 
                     sheet = "SubPlotlist")
regen <- trees[grep("s",trees$DBH),]
colnames(regen)[5] <- "Size"
trees <- trees[-grep("s",trees$DBH),]
trees$DBH <- as.numeric(trees$DBH)
trees$Count[is.na(trees$Count)] <- 1

####1b. Structure calculations####

trees_summ <-
  trees %>%
  group_by(Year, Subplot) %>%
  summarise(n_trees = sum(Count),
            ba = sum(pi*((DBH/2)^2)*0.0001*Count)
            )

#Standardize by ha
#Plots are 0.2 ac (0.0809371 ha) large
trees_summ[,c(3:4)] = trees_summ[,c(3:4)]/0.0809371
#Units now trees/ha and m2/ha

change = subplots[subplots$Year==2017,c(2:11)]
change$dens_change = trees_summ[trees_summ$Year==2017,"n_trees"] -
  trees_summ[trees_summ$Year==1970,"n_trees"]
change$ba_change = trees_summ[trees_summ$Year==2017,"ba"] -
  trees_summ[trees_summ$Year==1970,"ba"]

ggplot(change) +
  geom_bar(aes(x=Burned,y=dens_change),stat="summary",fun.y = "mean")
ggplot(change) +
  geom_bar(aes(x=Burned,y=ba_change),stat="summary",fun.y = "mean")
ggplot(change) +
  geom_bar(aes(x=Last_burned,y=dens_change),stat="summary",fun.y = "mean")
