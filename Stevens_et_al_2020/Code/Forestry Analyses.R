##This is a component of analysis code for Stevens et al. 2020: Forest vegetation change and its impacts on soil water following 47 years of managed wildfire. Ecosystems, In Press.
##Contact: Jens Stevens, stevensjt@gmail.com
##Forestry analyses: This code produces Figures 3 and 4 in the manuscript.

####0. Load libraries ####
library(tidyr) #for gather(); v 0.7.2
library(readxl) #For readxl; v 1.0.0
library(dplyr) #For between, others; v 0.7.4
require(devtools)
install_version("ggplot2", version = "2.2.1", repos = "http://cran.us.r-project.org") #skip updates
library(ggplot2) #for ggplot; v 2.2.1
library(gridExtra) #for gird.arrange(); v 2.3
library(lme4) #for lmer()
library(pbkrtest)
GetME_PVals=function(m){
  require(pbkrtest) #Kenward-Rodger approximation
  # get the KR-approximated degrees of freedom
  df.KR <- get_ddf_Lb(m, fixef(m))
  coefs <- data.frame(coef(summary(m)))
  # get p-values from the t-distribution using the t-values and approximated degrees of freedom
  coefs$p.KR <- round(2 * (1 - pt(abs(coefs$t.value), df.KR)),6); #coefs #Everything is significantly different from UB
  coefs$df = df.KR
  return(coefs)
}

####1. Data processing####
trees <- #Read tree list data
  as.data.frame(read_xlsx("./Data/SugarloafForestryPlotData.xlsx", 
            sheet = "Treelist"))
subplots <- #Read plot-level data
  as.data.frame(read_xlsx("./Data/SugarloafForestryPlotData.xlsx", 
                     sheet = "SubPlotlist"))
trees <- trees[-which(trees$Subplot==79.4),] #likely mis-location of plot
subplots <- subplots[-which(subplots$Subplot==79.4),] #likely mis-location of plot
regen <- #Create a separate dataset for regeneration
  trees[grep("s",trees$DBH),]
colnames(regen)[5] <- "Size"
trees <- #Remove regeneration from the tree list data
  trees[-grep("s",trees$DBH),]
trees$DBH <- as.numeric(trees$DBH)
trees <- #For comparison, both years have trees >3 inches (7.6 cm)
  trees[-which(between(trees$DBH,0.1,7.6)),] #Remove a few 2017 trees below cutoff.
trees$Count[is.na(trees$Count)] <- 1 #2017 trees were all measured so each row gets count=1

####2. Structure calculations, data prep####
#Calculate # trees per plot and ba per plot:
trees_summ <- 
  trees %>%
  group_by(Year, Subplot) %>%
  summarise(n_trees = sum(Count),
            ba = sum(pi*((DBH/2)^2)*0.0001*Count), #in m2/subpot
            ba_ABCO = sum(pi*((DBH[Species=="ABCO"]/2)^2)*0.0001 * 
                            Count[Species=="ABCO"]), #in m2
            ba_ABMA = sum(pi*((DBH[Species=="ABMA"]/2)^2)*0.0001 * 
                            Count[Species=="ABMA"]), #in m2
            ba_PICO = sum(pi*((DBH[Species=="PICO"]/2)^2)*0.0001 * 
                            Count[Species=="PICO"]), #in m2
            ba_PIJE = sum(pi*((DBH[Species=="PIJE"]/2)^2)*0.0001 * 
                            Count[Species=="PIJE"]) #in m2
            
            )
trees_summ <- as.data.frame(trees_summ)

trees_summ_15.2 <- #Trees > 15.2 cm (6 in)
  trees[trees$DBH>15.2,] %>%
  group_by(Year, Subplot) %>%
  summarise(n_trees = sum(Count),
            ba = sum(pi*((DBH/2)^2)*0.0001*Count), #in m2/subpot
            ba_ABCO = sum(pi*((DBH[Species=="ABCO"]/2)^2)*0.0001 * 
                            Count[Species=="ABCO"]), #in m2
            ba_ABMA = sum(pi*((DBH[Species=="ABMA"]/2)^2)*0.0001 * 
                            Count[Species=="ABMA"]), #in m2
            ba_PICO = sum(pi*((DBH[Species=="PICO"]/2)^2)*0.0001 * 
                            Count[Species=="PICO"]), #in m2
            ba_PIJE = sum(pi*((DBH[Species=="PIJE"]/2)^2)*0.0001 * 
                            Count[Species=="PIJE"]) #in m2
  )
trees_summ_15.2 <- as.data.frame(trees_summ_15.2)

trees_summ_61.0 <- #Trees > 61.0 cm (24 in)
  trees %>%
  group_by(Year, Subplot) %>%
  summarise(n_trees = sum(Count[DBH>61]),
            ba = sum(pi*((DBH[DBH>61]/2)^2)*0.0001*Count[DBH>61]), #in m2/subpot
            ba_ABCO = sum(pi*((DBH[Species=="ABCO"&DBH>61]/2)^2)*0.0001 * 
                            Count[Species=="ABCO"&DBH>61]), #in m2
            ba_ABMA = sum(pi*((DBH[Species=="ABMA"&DBH>61]/2)^2)*0.0001 * 
                            Count[Species=="ABMA"&DBH>61]), #in m2
            ba_PICO = sum(pi*((DBH[Species=="PICO"&DBH>61]/2)^2)*0.0001 * 
                            Count[Species=="PICO"&DBH>61]), #in m2
            ba_PIJE = sum(pi*((DBH[Species=="PIJE"&DBH>61]/2)^2)*0.0001 * 
                            Count[Species=="PIJE"&DBH>61]) #in m2
  ) 
trees_summ_61.0 <- as.data.frame(trees_summ_61.0)

trees_summ_100 <- #Trees > 100 cm (39 in)
  trees %>%
  group_by(Year, Subplot) %>%
  summarise(n_trees = sum(Count[DBH>100]),
            ba = sum(pi*((DBH[DBH>100]/2)^2)*0.0001*Count[DBH>100]), #in m2/subpot
            ba_ABCO = sum(pi*((DBH[Species=="ABCO"&DBH>100]/2)^2)*0.0001 * 
                            Count[Species=="ABCO"&DBH>100]), #in m2
            ba_ABMA = sum(pi*((DBH[Species=="ABMA"&DBH>100]/2)^2)*0.0001 * 
                            Count[Species=="ABMA"&DBH>100]), #in m2
            ba_PICO = sum(pi*((DBH[Species=="PICO"&DBH>100]/2)^2)*0.0001 * 
                            Count[Species=="PICO"&DBH>100]), #in m2
            ba_PIJE = sum(pi*((DBH[Species=="PIJE"&DBH>100]/2)^2)*0.0001 * 
                            Count[Species=="PIJE"&DBH>100]) #in m2  ) 
  )
trees_summ_100 <- as.data.frame(trees_summ_100)

#Standardize by ha, and add to subplots df:
#Plots are 0.2 ac (0.0809371 ha) large
#Units will now be trees/ha and m2/ha
subplots$dens <- #density
  trees_summ$n_trees/0.0809371
subplots$dens_15.2 <- #density
  trees_summ_15.2$n_trees/0.0809371
subplots$dens_61.0 <- #density
  trees_summ_61.0$n_trees/0.0809371
subplots$dens_100 <- #density
  trees_summ_100$n_trees/0.0809371
subplots$ba <- #basal area
  trees_summ$ba/0.0809371
subplots$ba_15.2 <- #basal area
  trees_summ_15.2$ba/0.0809371
subplots$ba_61.0 <- #basal area
  trees_summ_61.0$ba/0.0809371
subplots$ba_100 <- #basal area
  trees_summ_100$ba/0.0809371
subplots$ABCO <- trees_summ$ba_ABCO/0.0809371
subplots$ABMA <- trees_summ$ba_ABMA/0.0809371
subplots$PICO <- trees_summ$ba_PICO/0.0809371
subplots$PIJE <- trees_summ$ba_PIJE/0.0809371
subplots$ABCO_15.2 <- trees_summ_15.2$ba_ABCO/0.0809371
subplots$ABMA_15.2 <- trees_summ_15.2$ba_ABMA/0.0809371
subplots$PICO_15.2 <- trees_summ_15.2$ba_PICO/0.0809371
subplots$PIJE_15.2 <- trees_summ_15.2$ba_PIJE/0.0809371
subplots$ABCO_61.0 <- trees_summ_61.0$ba_ABCO/0.0809371
subplots$ABMA_61.0 <- trees_summ_61.0$ba_ABMA/0.0809371
subplots$PICO_61.0 <- trees_summ_61.0$ba_PICO/0.0809371
subplots$PIJE_61.0 <- trees_summ_61.0$ba_PIJE/0.0809371
subplots$ABCO_100 <- trees_summ_100$ba_ABCO/0.0809371
subplots$ABMA_100 <- trees_summ_100$ba_ABMA/0.0809371
subplots$PICO_100 <- trees_summ_100$ba_PICO/0.0809371
subplots$PIJE_100 <- trees_summ_100$ba_PIJE/0.0809371

#Set up long form datasets for stacking BA by 4 main species:
subplots_ba_spp <-
  subplots %>%
  gather(key=Species,value=ba_spp,c(ABCO,ABMA,PICO,PIJE))
subplots_ba_spp_15.2 <-
  subplots %>%
  gather(key=Species,value=ba_spp,c(ABCO_15.2,ABMA_15.2,PICO_15.2,PIJE_15.2))
subplots_ba_spp_15.2$Species <- subplots_ba_spp$Species
subplots_ba_spp_61.0 <-
  subplots %>%
  gather(key=Species,value=ba_spp,c(ABCO_61.0,ABMA_61.0,PICO_61.0,PIJE_61.0))
subplots_ba_spp_61.0$Species <- subplots_ba_spp$Species
subplots_ba_spp_100 <-
  subplots %>%
  gather(key=Species,value=ba_spp,c(ABCO_100,ABMA_100,PICO_100,PIJE_100))
subplots_ba_spp_100$Species <- subplots_ba_spp$Species

#Create Change data frame
change <- subplots[subplots$Year==2017,c(2:11)]
change$dens_2017 <- subplots[subplots$Year==2017,"dens"]
change$dens_1970 <- subplots[subplots$Year==1970,"dens"]
change$dens_change <- change$dens_2017 - change$dens_1970
change$ba_2017 <- subplots[subplots$Year==2017,"ba"]
change$ba_1970 <- subplots[subplots$Year==1970,"ba"]
change$ba_change <- change$ba_2017 - change$ba_1970
change$shrub_gain <- ifelse(
  subplots[subplots$Year==1970,"Shrubs_Simple"]=="None" &
    subplots[subplots$Year==2017,"Shrubs_Simple"]!="None", 1,0
)
change$shrub_loss <- ifelse(
  subplots[subplots$Year==1970,"Shrubs_Simple"]!="None" &
    subplots[subplots$Year==2017,"Shrubs_Simple"]=="None", 1,0
)

#Other housekeeping:
subplots$Year <- factor(subplots$Year)
subplots$Shrubs_Any <- ifelse(subplots$Shrubs_Simple=="None",0,1)

####3.Statistical analyses and plots####
#hist(change$dens_change) #wide range of density change values (2017-1970), tendency for increase
#hist(change$ba_change) #same for BA, tendency for decrease
#hist(log(subplots$dens+0))

#Density, all sizes#:
m<-vector()
m[1] <- GetME_PVals(lmer(log(dens) ~ Year + (1|Subplot), data = subplots[subplots$Burned==0,]))$p.KR[2]
m[2] <- GetME_PVals(lmer(log(dens) ~ Year + (1|Subplot), data = subplots[subplots$Burned==1,]))$p.KR[2]
m[3] <- GetME_PVals(lmer(log(dens) ~ Year + (1|Subplot), data = subplots[subplots$Burned==2,]))$p.KR[2]
p.labs <- ifelse(m<0.05,"*","")
p3a<-
  ggplot(subplots,
       aes(x = Burned, y = dens, fill = factor(Year))) +
  geom_bar(stat="summary",fun.y = "mean",
           position = position_dodge(width=0.95)) +
  stat_summary(fun.data = "mean_se", geom = "errorbar",
               position = position_dodge(width=0.95), width = 0.5) +
  scale_fill_manual(values = c("1970" = "#fdc086", "2017" = "#386cb0"))+
  labs(x = "", y = expression(atop("trees > 7.6 cm","trees" ~ha ^ -1)), 
                      fill = "year", title = "tree density") +
  annotate(geom = "text", x = c(0,1,2), y = 360, label = p.labs, size = 6.5) +
  annotate(geom = "text", x = c(2.5), y = 360, label = "a", size = 6.5) +
  theme_bw() +
  theme(legend.position = "none")

#Basal area, all sizes#
#hist(log(subplots$ba))
m[1] <- GetME_PVals(lmer(log(ba) ~ Year + (1|Subplot), data = subplots[subplots$Burned==0,]))$p.KR[2]
m[2] <- GetME_PVals(lmer(log(ba) ~ Year + (1|Subplot), data = subplots[subplots$Burned==1,]))$p.KR[2]
m[3] <- GetME_PVals(lmer(log(ba) ~ Year + (1|Subplot), data = subplots[subplots$Burned==2,]))$p.KR[2]
p.labs <- ifelse(m<0.05,"*","")
p3b<-
ggplot(subplots,
       aes(x = Burned, y = ba, fill = factor(Year))) +
  geom_bar(stat="summary",fun.y = "mean",
           position = position_dodge(width=0.95)) +
  stat_summary(fun.data = "mean_se", geom = "errorbar",
               position = position_dodge(width=0.95), width = 0.5) +
  scale_fill_manual(values = c("1970" = "#fdc086", "2017" = "#386cb0"))+
  labs(x = "", y = expression(m^2*ha^-1), fill = "year", title = "basal area") +
  annotate(geom = "text", x = c(0,1,2), y = 65, label = p.labs, size = 6.5) +
  annotate(geom = "text", x = c(2.5), y = 65, label = "b", size = 6.5) +
  theme_bw() + 
  theme(legend.position = "none")

#Basal area by species, all sizes#
p3c<- 
  ggplot(subplots_ba_spp, aes(x = factor(Year), y = ba_spp, fill = Species)) +
  geom_bar(position = "stack", stat = "summary",fun.y = "mean") +
  scale_fill_manual(values = c("ABMA" = "#5e3c99", "ABCO" = "#b2abd2",
                               "PICO" = "#f4a582", "PIJE" = "#e66101"))+
  labs(x = "", y = expression(m^2*ha^-1), fill = "species:  ",
       subtitle = "basal area by species") +
  facet_wrap( ~ Burned, labeller = as_labeller(c(
    "0" = "0 times \nburned", "1" = "1 time \nburned", "2" = "2 times \nburned")
    ))+
  annotate(geom = "text", x = c(2.4), y = 40, label = c("", "", "c"), size = 6.5) +
  theme_bw() +
  guides(fill = guide_legend(nrow = 2)) +
  theme(plot.subtitle = element_text(hjust = 0.5), legend.position = "none")

#Density, trees > 15.2 cm dbh#
m[1] <- GetME_PVals(lmer(log(dens_15.2) ~ Year + (1|Subplot), 
                         data = subplots[subplots$Burned==0,]))$p.KR[2]
m[2] <- GetME_PVals(lmer(log(dens_15.2) ~ Year + (1|Subplot), 
                         data = subplots[subplots$Burned==1,]))$p.KR[2]
m[3] <- GetME_PVals(lmer(log(dens_15.2) ~ Year + (1|Subplot), 
                         data = subplots[subplots$Burned==2,]))$p.KR[2]
p.labs <- ifelse(m<0.05,"*","")
p3d<-
  ggplot(subplots,
         aes(x = Burned, y = dens_15.2, fill = factor(Year))) +
  geom_bar(stat="summary",fun.y = "mean",
           position = position_dodge(width=0.95)) +
  stat_summary(fun.data = "mean_se", geom = "errorbar",
               position = position_dodge(width=0.95), width = 0.5) +
  scale_fill_manual(values = c("1970" = "#fdc086", "2017" = "#386cb0"))+
  labs(x = "", y = expression(atop("trees > 15.2 cm","trees" ~ha ^ -1)), 
       fill = "year", main = "tree density") +
  annotate(geom = "text", x = c(0,1,2), y = 360, label = p.labs, size = 6.5) +
  annotate(geom = "text", x = c(2.5), y = 360, label = "d", size = 6.5) +
  theme_bw()+
  theme(legend.position = "none")

#Basal area, trees > 15.2 cm dbh#
m[1] <- GetME_PVals(lmer(log(ba_15.2) ~ Year + (1|Subplot), 
                         data = subplots[subplots$Burned==0,]))$p.KR[2]
m[2] <- GetME_PVals(lmer(log(ba_15.2) ~ Year + (1|Subplot), 
                         data = subplots[subplots$Burned==1,]))$p.KR[2]
m[3] <- GetME_PVals(lmer(log(ba_15.2) ~ Year + (1|Subplot), 
                         data = subplots[subplots$Burned==2,]))$p.KR[2]
p.labs <- ifelse(m<0.05,"*","")
p3e<-
  ggplot(subplots,
         aes(x = Burned, y = ba_15.2, fill = factor(Year))) +
  geom_bar(stat="summary",fun.y = "mean",
           position = position_dodge(width=0.95)) +
  stat_summary(fun.data = "mean_se", geom = "errorbar",
               position = position_dodge(width=0.95), width = 0.5) +
  scale_fill_manual(values = c("1970" = "#fdc086", "2017" = "#386cb0"))+
  labs(x = "", y = expression(m^2*ha^-1), fill = "year") +
  annotate(geom = "text", x = c(0,1,2), y = 65, label = p.labs, size = 6.5) +
  annotate(geom = "text", x = c(2.5), y = 65, label = "e", size = 6.5) +
  theme_bw()+
  theme(legend.position = "none")

#Basal area by species, trees > 15.2 cm dbh#
p3f<- 
  ggplot(subplots_ba_spp_15.2, aes(x = factor(Year), y = ba_spp, fill = Species)) +
  geom_bar(position = "stack", stat = "summary",fun.y = "mean") +
  scale_fill_manual(values = c("ABMA" = "#5e3c99", "ABCO" = "#b2abd2",
                               "PICO" = "#f4a582", "PIJE" = "#e66101"))+
  labs(x = "", y = expression(m^2*ha^-1), fill = "species:  ") +
  facet_wrap( ~ Burned)+
  annotate(geom = "text", x = c(2.4), y = 40, label = c("", "", "f"), size = 6.5) +
  theme_bw() +
  guides(fill = guide_legend(nrow = 2)) +
  theme(legend.position = "none")

#Density, trees > 61 cm dbh#
m[1] <- GetME_PVals(lmer(log(dens_61.0+1) ~ Year + (1|Subplot), 
                         data = subplots[subplots$Burned==0,]))$p.KR[2]
m[2] <- GetME_PVals(lmer(log(dens_61.0+1) ~ Year + (1|Subplot), 
                         data = subplots[subplots$Burned==1,]))$p.KR[2]
m[3] <- GetME_PVals(lmer(log(dens_61.0+1) ~ Year + (1|Subplot), 
                         data = subplots[subplots$Burned==2,]))$p.KR[2]
p.labs <- ifelse(m<0.05,"*","")
p3g<-
  ggplot(subplots,
         aes(x = Burned, y = dens_61.0, fill = factor(Year))) +
  geom_bar(stat="summary",fun.y = "mean",
           position = position_dodge(width=0.95)) +
  stat_summary(fun.data = "mean_se", geom = "errorbar",
               position = position_dodge(width=0.95), width = 0.5) +
  scale_fill_manual(values = c("1970" = "#fdc086", "2017" = "#386cb0"))+
  labs(x = "", y = expression(atop("trees > 61 cm","trees" ~ha ^ -1)), 
       fill = "year", main = "tree density") +
  annotate(geom = "text", x = c(0,1,2), y = 50, label = p.labs, size = 9.5) +
  annotate(geom = "text", x = c(2.5), y = 90, label = "g", size = 6.5) +
  theme_bw()+
  theme(legend.position = "none")

#Basal area, trees > 61 cm dbh#
m[1] <- GetME_PVals(lmer(log(ba_61.0+1) ~ Year + (1|Subplot), 
                         data = subplots[subplots$Burned==0,]))$p.KR[2]
m[2] <- GetME_PVals(lmer(log(ba_61.0+1) ~ Year + (1|Subplot), 
                         data = subplots[subplots$Burned==1,]))$p.KR[2]
m[3] <- GetME_PVals(lmer(log(ba_61.0+1) ~ Year + (1|Subplot), 
                         data = subplots[subplots$Burned==2,]))$p.KR[2]
p.labs <- ifelse(m<0.05,"*","")
p3h<-
  ggplot(subplots,
         aes(x = Burned, y = ba_61.0, fill = factor(Year))) +
  geom_bar(stat="summary",fun.y = "mean",
           position = position_dodge(width=0.95)) +
  stat_summary(fun.data = "mean_se", geom = "errorbar",
               position = position_dodge(width=0.95), width = 0.5) +
  scale_fill_manual(values = c("1970" = "#fdc086", "2017" = "#386cb0"))+
  labs(x = "", y = expression(m^2*ha^-1), fill = "year") +
  annotate(geom = "text", x = c(0,1,2), y = 30, label = p.labs, size = 9.5) +
  annotate(geom = "text", x = c(2.5), y = 65, label = "h", size = 6.5) +
  theme_bw()+
  theme(legend.position = "none")

#Basal area by species, trees > 61 cm dbh#
p3i<- 
  ggplot(subplots_ba_spp_61.0, aes(x = factor(Year), y = ba_spp, fill = Species)) +
  geom_bar(position = "stack", stat = "summary",fun.y = "mean") +
  scale_fill_manual(values = c("ABMA" = "#5e3c99", "ABCO" = "#b2abd2",
                               "PICO" = "#f4a582", "PIJE" = "#e66101"))+
  labs(x = "", y = expression(m^2*ha^-1), fill = "species:  ") +
  facet_wrap( ~ Burned)+
  annotate(geom = "text", x = c(2.4), y = 40, label = c("", "", "i"), size = 6.5) +
  theme_bw() +
  guides(fill = guide_legend(nrow = 2)) +
  theme(legend.position = "none")

#Density, trees > 100 cm dbh#
m[1] <- GetME_PVals(lmer(log(dens_100+1) ~ Year + (1|Subplot), 
                         data = subplots[subplots$Burned==0,]))$p.KR[2]
m[2] <- GetME_PVals(lmer(log(dens_100+1) ~ Year + (1|Subplot), 
                         data = subplots[subplots$Burned==1,]))$p.KR[2]
m[3] <- GetME_PVals(lmer(log(dens_100+1) ~ Year + (1|Subplot), 
                         data = subplots[subplots$Burned==2,]))$p.KR[2]
p.labs <- ifelse(m<0.05,"*","")
p3j<-
  ggplot(subplots,
         aes(x = Burned, y = dens_100, fill = factor(Year))) +
  geom_bar(stat="summary",fun.y = "mean",
           position = position_dodge(width=0.95)) +
  stat_summary(fun.data = "mean_se", geom = "errorbar",
               position = position_dodge(width=0.95), width = 0.5) +
  scale_fill_manual(values = c("1970" = "#fdc086", "2017" = "#386cb0"))+
  labs(x = "times burned", y = expression(atop("trees > 100 cm","trees" ~ha ^ -1)), 
       fill = "year", main = "tree density") +
  annotate(geom = "text", x = c(0,1,2), y = 30, label = p.labs, size = 9.5) +
  annotate(geom = "text", x = c(2.5), y = 45, label = "j", size = 6.5) +
  theme_bw()+
  theme(legend.position = "none")

#Basal area, trees > 100 cm dbh#
m[1] <- GetME_PVals(lmer(log(ba_100+1) ~ Year + (1|Subplot), 
                         data = subplots[subplots$Burned==0,]))$p.KR[2]
m[2] <- GetME_PVals(lmer(log(ba_100+1) ~ Year + (1|Subplot), 
                         data = subplots[subplots$Burned==1,]))$p.KR[2]
m[3] <- GetME_PVals(lmer(log(ba_100+1) ~ Year + (1|Subplot), 
                         data = subplots[subplots$Burned==2,]))$p.KR[2]
p.labs <- ifelse(m<0.05,"*","")
p3k<-
  ggplot(subplots,
         aes(x = Burned, y = ba_100, fill = factor(Year))) +
  geom_bar(stat="summary",fun.y = "mean",
           position = position_dodge(width=0.95)) +
  stat_summary(fun.data = "mean_se", geom = "errorbar",
               position = position_dodge(width=0.95), width = 0.5) +
  scale_fill_manual(values = c("1970" = "#fdc086", "2017" = "#386cb0"))+
  labs(x = "times burned", y = expression(m^2*ha^-1), fill = "year") +
  annotate(geom = "text", x = c(0,1,2), y = 30, label = p.labs, size = 9.5) +
  annotate(geom = "text", x = c(2.5), y = 65, label = "k", size = 6.5) +
  theme_bw()+
  theme(legend.position = c(0.2, 0.78))

#Basal area, trees > 100 cm dbh#
p3l<- 
  ggplot(subplots_ba_spp_100, aes(x = factor(Year), y = ba_spp, fill = Species)) +
  geom_bar(position = "stack", stat = "summary",fun.y = "mean") +
  scale_fill_manual(values = c("ABMA" = "#5e3c99", "ABCO" = "#b2abd2",
                               "PICO" = "#f4a582", "PIJE" = "#e66101"))+
  labs(x = "year", y = expression(m^2*ha^-1), fill = "species:  ") +
  facet_wrap( ~ Burned)+
  annotate(geom = "text", x = c(2.4), y = 40, label = c("", "", "l"), size = 6.5) +
  theme_bw() +
  guides(fill = guide_legend(nrow = 2)) +
  theme(legend.position = c(0.62,0.65), legend.title = element_text(size = 7),
    legend.text = element_text(size = 6), 
    legend.background = element_rect(fill = NA))

p4<-
  ggplot(subplots,
         aes(x = Burned, y = Shrubs_Any, fill = factor(Year))) +
  geom_bar(stat="summary",fun.y = function(x){sum(x)/length(x)},
           position = position_dodge(width=0.95)) +
  scale_fill_manual(values = c("1970" = "darkgreen", "2017" = "goldenrod"))+
  labs(x = "times burned", y = "proportion plots with shrubs", fill = "year") +
  theme_bw()+
  theme(legend.position = c(0.2, 0.8))


pdf("./Figures/Fig3.pdf", width = 9, height = 10)
grid.arrange(p3a,p3b,p3c,p3d,p3e,p3f,p3g,p3h,p3i,p3j,p3k,p3l, ncol = 3)
dev.off()

pdf("./Figures/Fig4.pdf", width = 3, height = 3)
p4
dev.off()

