

####1. Forestry analyses####
library(tidyr) #for gather(); v 0.7.2
library(readxl) #For readxl; v 1.0.0
library(dplyr) #For between, others; v 0.7.4
library(ggplot2) #for ggplot; v 2.2.1
library(gridExtra) #for gird.arrange(); v 2.3
library(lme4) #for lmer()
#library(Hmisc) #for stat_summary(); v 4.1-1; Don't need
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

####1a. Data processing####
trees <- #Read tree list data
  as.data.frame(read_xlsx("./Raw Data/SugarloafForestryPlotDataCombined.xlsx", 
            sheet = "Treelist"))
subplots <- #Read plot-level data
  as.data.frame(read_xlsx("./Raw Data/SugarloafForestryPlotDataCombined.xlsx", 
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

####1b. Structure calculations, data prep####
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

#Set up long form datasets for stacking BA by species:
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

####1c.Statistical analyses and plots####
#hist(change$dens_change) #wide range of density change values (2017-1970), tendency for increase
#hist(change$ba_change) #same for BA, tendency for decrease

#hist(log(subplots$dens+0))
m<-vector()
m[1] <- GetME_PVals(lmer(log(dens) ~ Year + (1|Subplot), data = subplots[subplots$Burned==0,]))$p.KR[2]
m[2] <- GetME_PVals(lmer(log(dens) ~ Year + (1|Subplot), data = subplots[subplots$Burned==1,]))$p.KR[2]
m[3] <- GetME_PVals(lmer(log(dens) ~ Year + (1|Subplot), data = subplots[subplots$Burned==2,]))$p.KR[2]
p.labs <- ifelse(m<0.05,"*","")
p4a<-
  ggplot(subplots,
       aes(x = Burned, y = dens, fill = factor(Year))) +
  geom_bar(stat="summary",fun.y = "mean",
           position = position_dodge(width=0.95)) +
  stat_summary(fun.data = "mean_se", geom = "errorbar",
               position = position_dodge(width=0.95), width = 0.5) +
  scale_fill_manual(values = c("1970" = "#fdc086", "2017" = "#386cb0"))+
  labs(x = "", y = expression(atop("trees > 7.6 cm","trees" ~ha ^ -1)), 
                      fill = "year", title = "tree density") +
  annotate(geom = "text", x = c(0,1,2), y = 370, label = p.labs, size = 6.5) +
  annotate(geom = "text", x = c(2.5), y = 370, label = "a", size = 6.5) +
  theme_bw() +
  theme(legend.position = "none")

#hist(log(subplots$ba))
m[1] <- GetME_PVals(lmer(log(ba) ~ Year + (1|Subplot), data = subplots[subplots$Burned==0,]))$p.KR[2]
m[2] <- GetME_PVals(lmer(log(ba) ~ Year + (1|Subplot), data = subplots[subplots$Burned==1,]))$p.KR[2]
m[3] <- GetME_PVals(lmer(log(ba) ~ Year + (1|Subplot), data = subplots[subplots$Burned==2,]))$p.KR[2]
p.labs <- ifelse(m<0.05,"*","")
p4b<-
ggplot(subplots,
       aes(x = Burned, y = ba, fill = factor(Year))) +
  geom_bar(stat="summary",fun.y = "mean",
           position = position_dodge(width=0.95)) +
  stat_summary(fun.data = "mean_se", geom = "errorbar",
               position = position_dodge(width=0.95), width = 0.5) +
  scale_fill_manual(values = c("1970" = "#fdc086", "2017" = "#386cb0"))+
  labs(x = "", y = expression(m^2*ha^-1), fill = "year", title = "basal area") +
  annotate(geom = "text", x = c(0,1,2), y = 70, label = p.labs, size = 6.5) +
  annotate(geom = "text", x = c(2.5), y = 70, label = "b", size = 6.5) +
  theme_bw() + 
  theme(legend.position = "none")

p4c<- 
  ggplot(subplots_ba_spp, aes(x = factor(Year), y = ba_spp, fill = Species)) +
  geom_bar(position = "stack", stat = "summary",fun.y = "mean") +
  scale_fill_manual(values = c("ABMA" = "#5e3c99", "ABCO" = "#b2abd2",
                               "PICO" = "#f4a582", "PIJE" = "#e66101"))+
  labs(x = "", y = expression(m^2*ha^-1), fill = "species:  ",
       subtitle = "basal area by species") +
  facet_wrap( ~ Burned, labeller = as_labeller(c(
    "0" = "0 times burned", "1" = "1 time burned", "2" = "2 times burned")
    ))+
  annotate(geom = "text", x = c(2.4), y = 43, label = c("", "", "c"), size = 6.5) +
  theme_bw() +
  guides(fill = guide_legend(nrow = 2)) +
  theme(plot.subtitle = element_text(hjust = 0.5), legend.position = "none")

#ggplot(change) + #deprecated
#  geom_bar(aes(x=Burned,y=dens_change),stat="summary",fun.y = "mean")
#ggplot(change) +
#  geom_bar(aes(x=Burned,y=ba_change),stat="summary",fun.y = "mean")
#ggplot(change) +
#  geom_bar(aes(x=Last_burned,y=dens_change),stat="summary",fun.y = "mean")

m[1] <- GetME_PVals(lmer(log(dens_15.2) ~ Year + (1|Subplot), 
                         data = subplots[subplots$Burned==0,]))$p.KR[2]
m[2] <- GetME_PVals(lmer(log(dens_15.2) ~ Year + (1|Subplot), 
                         data = subplots[subplots$Burned==1,]))$p.KR[2]
m[3] <- GetME_PVals(lmer(log(dens_15.2) ~ Year + (1|Subplot), 
                         data = subplots[subplots$Burned==2,]))$p.KR[2]
p.labs <- ifelse(m<0.05,"*","")
p4d<-
  ggplot(subplots,
         aes(x = Burned, y = dens_15.2, fill = factor(Year))) +
  geom_bar(stat="summary",fun.y = "mean",
           position = position_dodge(width=0.95)) +
  stat_summary(fun.data = "mean_se", geom = "errorbar",
               position = position_dodge(width=0.95), width = 0.5) +
  scale_fill_manual(values = c("1970" = "#fdc086", "2017" = "#386cb0"))+
  labs(x = "", y = expression(atop("trees > 15.2 cm","trees" ~ha ^ -1)), 
       fill = "year", main = "tree density") +
  annotate(geom = "text", x = c(0,1,2), y = 370, label = p.labs, size = 6.5) +
  annotate(geom = "text", x = c(2.5), y = 370, label = "d", size = 6.5) +
  theme_bw()+
  theme(legend.position = "none")

m[1] <- GetME_PVals(lmer(log(ba_15.2) ~ Year + (1|Subplot), 
                         data = subplots[subplots$Burned==0,]))$p.KR[2]
m[2] <- GetME_PVals(lmer(log(ba_15.2) ~ Year + (1|Subplot), 
                         data = subplots[subplots$Burned==1,]))$p.KR[2]
m[3] <- GetME_PVals(lmer(log(ba_15.2) ~ Year + (1|Subplot), 
                         data = subplots[subplots$Burned==2,]))$p.KR[2]
p.labs <- ifelse(m<0.05,"*","")
p4e<-
  ggplot(subplots,
         aes(x = Burned, y = ba_15.2, fill = factor(Year))) +
  geom_bar(stat="summary",fun.y = "mean",
           position = position_dodge(width=0.95)) +
  stat_summary(fun.data = "mean_se", geom = "errorbar",
               position = position_dodge(width=0.95), width = 0.5) +
  scale_fill_manual(values = c("1970" = "#fdc086", "2017" = "#386cb0"))+
  labs(x = "", y = expression(m^2*ha^-1), fill = "year") +
  annotate(geom = "text", x = c(0,1,2), y = 70, label = p.labs, size = 6.5) +
  annotate(geom = "text", x = c(2.5), y = 70, label = "e", size = 6.5) +
  theme_bw()+
  theme(legend.position = "none")

p4f<- 
  ggplot(subplots_ba_spp_15.2, aes(x = factor(Year), y = ba_spp, fill = Species)) +
  geom_bar(position = "stack", stat = "summary",fun.y = "mean") +
  scale_fill_manual(values = c("ABMA" = "#5e3c99", "ABCO" = "#b2abd2",
                               "PICO" = "#f4a582", "PIJE" = "#e66101"))+
  labs(x = "", y = expression(m^2*ha^-1), fill = "species:  ") +
  facet_wrap( ~ Burned)+
  annotate(geom = "text", x = c(2.4), y = 43, label = c("", "", "f"), size = 6.5) +
  theme_bw() +
  guides(fill = guide_legend(nrow = 2)) +
  theme(legend.position = "none")

m[1] <- GetME_PVals(lmer(log(dens_61.0+1) ~ Year + (1|Subplot), 
                         data = subplots[subplots$Burned==0,]))$p.KR[2]
m[2] <- GetME_PVals(lmer(log(dens_61.0+1) ~ Year + (1|Subplot), 
                         data = subplots[subplots$Burned==1,]))$p.KR[2]
m[3] <- GetME_PVals(lmer(log(dens_61.0+1) ~ Year + (1|Subplot), 
                         data = subplots[subplots$Burned==2,]))$p.KR[2]
p.labs <- ifelse(m<0.05,"*","")
p4g<-
  ggplot(subplots,
         aes(x = Burned, y = dens_61.0, fill = factor(Year))) +
  geom_bar(stat="summary",fun.y = "mean",
           position = position_dodge(width=0.95)) +
  stat_summary(fun.data = "mean_se", geom = "errorbar",
               position = position_dodge(width=0.95), width = 0.5) +
  scale_fill_manual(values = c("1970" = "#fdc086", "2017" = "#386cb0"))+
  labs(x = "", y = expression(atop("trees > 61 cm","trees" ~ha ^ -1)), 
       fill = "year", main = "tree density") +
  annotate(geom = "text", x = c(0,1,2), y = 100, label = p.labs, size = 6.5) +
  annotate(geom = "text", x = c(2.5), y = 100, label = "g", size = 6.5) +
  theme_bw()+
  theme(legend.position = "none")

m[1] <- GetME_PVals(lmer(log(ba_61.0+1) ~ Year + (1|Subplot), 
                         data = subplots[subplots$Burned==0,]))$p.KR[2]
m[2] <- GetME_PVals(lmer(log(ba_61.0+1) ~ Year + (1|Subplot), 
                         data = subplots[subplots$Burned==1,]))$p.KR[2]
m[3] <- GetME_PVals(lmer(log(ba_61.0+1) ~ Year + (1|Subplot), 
                         data = subplots[subplots$Burned==2,]))$p.KR[2]
p.labs <- ifelse(m<0.05,"*","")
p4h<-
  ggplot(subplots,
         aes(x = Burned, y = ba_61.0, fill = factor(Year))) +
  geom_bar(stat="summary",fun.y = "mean",
           position = position_dodge(width=0.95)) +
  stat_summary(fun.data = "mean_se", geom = "errorbar",
               position = position_dodge(width=0.95), width = 0.5) +
  scale_fill_manual(values = c("1970" = "#fdc086", "2017" = "#386cb0"))+
  labs(x = "", y = expression(m^2*ha^-1), fill = "year") +
  annotate(geom = "text", x = c(0,1,2), y = 70, label = p.labs, size = 6.5) +
  annotate(geom = "text", x = c(2.5), y = 70, label = "h", size = 6.5) +
  theme_bw()+
  theme(legend.position = "none")

p4i<- 
  ggplot(subplots_ba_spp_61.0, aes(x = factor(Year), y = ba_spp, fill = Species)) +
  geom_bar(position = "stack", stat = "summary",fun.y = "mean") +
  scale_fill_manual(values = c("ABMA" = "#5e3c99", "ABCO" = "#b2abd2",
                               "PICO" = "#f4a582", "PIJE" = "#e66101"))+
  labs(x = "", y = expression(m^2*ha^-1), fill = "species:  ") +
  facet_wrap( ~ Burned)+
  annotate(geom = "text", x = c(2.4), y = 43, label = c("", "", "i"), size = 6.5) +
  theme_bw() +
  guides(fill = guide_legend(nrow = 2)) +
  theme(legend.position = "none")

m[1] <- GetME_PVals(lmer(log(dens_100+1) ~ Year + (1|Subplot), 
                         data = subplots[subplots$Burned==0,]))$p.KR[2]
m[2] <- GetME_PVals(lmer(log(dens_100+1) ~ Year + (1|Subplot), 
                         data = subplots[subplots$Burned==1,]))$p.KR[2]
m[3] <- GetME_PVals(lmer(log(dens_100+1) ~ Year + (1|Subplot), 
                         data = subplots[subplots$Burned==2,]))$p.KR[2]
p.labs <- ifelse(m<0.05,"*","")
p4j<-
  ggplot(subplots,
         aes(x = Burned, y = dens_100, fill = factor(Year))) +
  geom_bar(stat="summary",fun.y = "mean",
           position = position_dodge(width=0.95)) +
  stat_summary(fun.data = "mean_se", geom = "errorbar",
               position = position_dodge(width=0.95), width = 0.5) +
  scale_fill_manual(values = c("1970" = "#fdc086", "2017" = "#386cb0"))+
  labs(x = "times burned", y = expression(atop("trees > 100 cm","trees" ~ha ^ -1)), 
       fill = "year", main = "tree density") +
  annotate(geom = "text", x = c(0,1,2), y = 50, label = p.labs, size = 6.5) +
  annotate(geom = "text", x = c(2.5), y = 50, label = "j", size = 6.5) +
  theme_bw()+
  theme(legend.position = "none")

m[1] <- GetME_PVals(lmer(log(ba_100+1) ~ Year + (1|Subplot), 
                         data = subplots[subplots$Burned==0,]))$p.KR[2]
m[2] <- GetME_PVals(lmer(log(ba_100+1) ~ Year + (1|Subplot), 
                         data = subplots[subplots$Burned==1,]))$p.KR[2]
m[3] <- GetME_PVals(lmer(log(ba_100+1) ~ Year + (1|Subplot), 
                         data = subplots[subplots$Burned==2,]))$p.KR[2]
p.labs <- ifelse(m<0.05,"*","")
p4k<-
  ggplot(subplots,
         aes(x = Burned, y = ba_100, fill = factor(Year))) +
  geom_bar(stat="summary",fun.y = "mean",
           position = position_dodge(width=0.95)) +
  stat_summary(fun.data = "mean_se", geom = "errorbar",
               position = position_dodge(width=0.95), width = 0.5) +
  scale_fill_manual(values = c("1970" = "#fdc086", "2017" = "#386cb0"))+
  labs(x = "times burned", y = expression(m^2*ha^-1), fill = "year") +
  annotate(geom = "text", x = c(0,1,2), y = 70, label = p.labs, size = 6.5) +
  annotate(geom = "text", x = c(2.5), y = 70, label = "k", size = 6.5) +
  theme_bw()+
  theme(legend.position = c(0.6, 0.6))

p4l<- 
  ggplot(subplots_ba_spp_100, aes(x = factor(Year), y = ba_spp, fill = Species)) +
  geom_bar(position = "stack", stat = "summary",fun.y = "mean") +
  scale_fill_manual(values = c("ABMA" = "#5e3c99", "ABCO" = "#b2abd2",
                               "PICO" = "#f4a582", "PIJE" = "#e66101"))+
  labs(x = "year", y = expression(m^2*ha^-1), fill = "species:  ") +
  facet_wrap( ~ Burned)+
  annotate(geom = "text", x = c(2.4), y = 43, label = c("", "", "l"), size = 6.5) +
  theme_bw() +
  guides(fill = guide_legend(nrow = 2)) +
  theme(legend.position = c(0.62,0.65), legend.title = element_text(size = 7),
    legend.text = element_text(size = 6), 
    legend.background = element_rect(fill = NA))

p5<-
  ggplot(subplots,
         aes(x = Burned, y = Shrubs_Any, fill = factor(Year))) +
  geom_bar(stat="summary",fun.y = function(x){sum(x)/length(x)},
           position = position_dodge(width=0.95)) +
  scale_fill_manual(values = c("1970" = "darkgreen", "2017" = "goldenrod"))+
  labs(x = "times burned", y = "proportion plots with shrubs", fill = "year") +
  theme_bw()+
  theme(legend.position = c(0.2, 0.8))


pdf("./Figures/MS/Fig4.pdf", width = 9, height = 10)
grid.arrange(p4a,p4b,p4c,p4d,p4e,p4f,p4g,p4h,p4i,p4j,p4k,p4l, ncol = 3)
dev.off()

pdf("./Figures/MS/Fig5.pdf", width = 3, height = 3)
p5
dev.off()

####2. Imagery####
library(raster)
library(rgdal)
library(igraph) #for clump(); version 1.2.2
#library(corrplot) #for corrplot(); version 0.84
#library(htmlTable)
library(ggcorrplot)
library(tableHTML)
library(gridExtra)
library(grid)
library(reshape2) #for melt(); version 1.4.3

####2a. Data processing- only run once, don't need to re-run####
#r73scb <- 
#  raster("./Processed Data/Classified Images/Final rasters/1973_raster_match_SCB.tif")

#Fix an error where shrubs (should be 1) and sparse meadow (should be 2) were reversed in 2014 (Only ran once, now fixed):
#r14scb <-
#  raster("./Processed Data/Classified Images/Final rasters/2014_raster_match_SCB.tif")
#tmp <- reclassify(r14scb,c(0.5,1.5,2,1.5,2.5,1))
#writeRaster(tmp,"./Processed Data/Classified Images/Final rasters/2014_raster_match_SCB_fixed.tif") 

#r14scb <-
#  raster("./Processed Data/Classified Images/Final rasters/2014_raster_match_SCB_fixed.tif")

#Take single pixels surrounded by other values and merge.
#r <- r73scb #Apply to eithe r14scb or r73scb
#active_cells <- which(!is.na(getValues(r)))
#a <- adjacent(r, active_cells, 4, pairs=TRUE) #CHECKME could eventually change 4 ->8
#for(c in unique(a[,1])){ #Change value of isolated cells.
#  print(c)
#  target_val <- #Value of a focal cell
#    getValues(r)[c]
#  comp_vals <- #Value of each of the 4 adjacent cells
#    getValues(r)[a[which(a[,1]==c),2]]
#  print(target_val)
#  print(comp_vals)
#  if(target_val != 3 & #ignore granite
#     #ignore isolated cells where all the border cells are granite:
#     #unique(length(unique(comp_vals))!=1 & unique(comp_vals)!=3) &
#     !any(is.na(comp_vals)) & #ignore border cells that have NA adjacent
#     !any(comp_vals%in%target_val) #if all border cells are different
#     ){
    #If all adjacent cells are of different values, and not NA:
    #adjust the cell c in question:
#    r[c] <- #Reassign focal cell to the most common of the adjacent cells.
#      ifelse(all(comp_vals==3), r[c], #If all neighbors are granite, don't change
#             ifelse(any(comp_vals==3), #Exclude granite from reclass options
#                    modal(comp_vals[-c(which(comp_vals==3))]),
#                    modal(comp_vals) )
#             )
#    print(paste("changing cell", c, "from", target_val, "to", modal(comp_vals)))
#  }
#}
#plot(r)

#writeRaster(r,"./Processed Data/Classified Images/Final rasters/1973_raster_match_SCB_analysis.tif", overwrite = TRUE)

#16:28:21
#3927
#16:28:37
#Took 16 seconds to do 1/63rd of the dataset. 16*63 = 1008 seconds, 17 minutes to do a full raster reclass.

####2b. Fire perimeter stuff####
###Read data
r73scb <- #Load the processed 1973 veg raster
  raster("./Processed Data/Classified Images/Final rasters/1973_raster_match_SCB_analysis.tif")
r14scb_new <- #Load the processed 2014 veg raster
  raster("./Processed Data/Classified Images/Final rasters/2014_raster_match_SCB_analysis_fixwetlandsGB.tif")
r14scb <- #Load the processed 2014 veg raster, before G fixed wetlands
  raster("./Processed Data/Classified Images/Final rasters/2014_raster_match_SCB_analysis.tif")
perims <- readOGR("./Processed Data/GIS/Sugarloaf Fires 1973-2003.shp")

###Process data
#Tip: http://r-sig-geo.2731867.n2.nabble.com/Efficient-way-to-obtain-gridded-count-of-overlapping-polygons-td6034590.html
perims <- spTransform(perims, CRSobj = crs(r73scb)) #reproject perimeters
#plot(perims, add=T) #Confirm alignment
r73_pts <- rasterToPoints(r73scb,spatial=TRUE)
r14_pts <- rasterToPoints(r14scb,spatial=TRUE)

##Below code identifies four missing values from revised layer and inserts them
##Only need to run once?
vals14<- na.exclude(getValues(r14scb))
vals14_new <- na.exclude(getValues(r14scb_new))
tmp <- vals14 - vals14_new
which(tmp!=0)
vals14[1715:1730]
vals14_new[1715:1730]
#missing new value 1 = #1717 (should be 4)
vals14_new <- c(vals14_new[1:1716],4,vals14_new[1717:length(vals14_new)])
tmp <- vals14 - vals14_new
which(tmp!=0)
#missing old value 1 = #1862 (should be 4)
vals14[1860:1880]
vals14_new[1860:1880]
vals14_new <- c(vals14_new[1:1861],4,vals14_new[1862:length(vals14_new)])
tmp <- vals14 - vals14_new
which(tmp!=0)
#missing old value 1 = #17301 (should be 4)
vals14[17300:17310]
vals14_new[17300:17310]
vals14_new <- c(vals14_new[1:17300],4,vals14_new[17301:length(vals14_new)])
tmp <- vals14 - vals14_new
which(tmp!=0)
#missing old value 1 = #17501 (should be 4)
vals14[17500:17520]
vals14_new[17500:17520]
vals14_new <- c(vals14_new[1:17500],4,vals14_new[17501:length(vals14_new)])
tmp <- vals14 - vals14_new
which(tmp!=0)
vals14[which(tmp!=0)]
vals14_new[which(tmp!=0)]
tmp <- getValues(r14scb)
blank <- rep(NA,length(tmp))
blank[which(!is.na(tmp))] <- vals14_new
tmp <- setValues(r14scb, blank)
#writeRaster(tmp,"./Processed Data/Classified Images/Final rasters/2014_raster_match_SCB_analysis.tif") 

which(tmp!=0)
getValues(r14scb)[which(tmp!=0)]
getValues(r14scb_old)[which(tmp!=0)]
#r73_tmp <- resample(r73scb,r14scb,method = "ngb")
#r73_pts_tmp <- rasterToPoints(r73_tmp,spatial=TRUE)

#tmp <- setValues(r14scb,getValues(r14scb_new))
r14_pts_new <- rasterToPoints(r14scb_new,spatial=TRUE)


r14_pts$n_fires <- r73_pts$n_fires <- #test alignment
  sapply(over(r73_pts, geometry(perims), returnList = TRUE), length)
r14_pts_old$n_fires <- r73_pts$n_fires <- #test alignment
  sapply(over(r73_pts, geometry(perims), returnList = TRUE), length)
r14_pts_tmp$n_fires <- r73_pts$n_fires <-  #test alignment
  sapply(over(r73_pts, geometry(perims), returnList = TRUE), length)


##

r14_pts$n_fires <- r73_pts$n_fires <- #Count number of overlapping fires
  sapply(over(r73_pts, geometry(perims), returnList = TRUE), length)
r73_pts$n_fires[r73_pts$n_fires > 2] <- #93 pixels had 4 fires, small n. 1360 pixels had 3 fires, chi-squared test was not converging. Converting to 3+4 burns to 2.
  2 #93 pixels at 0.16 ha/pixel = 14.88 ha
r14_pts$n_fires[r14_pts$n_fires > 2] <- #93 pixels had 4 fires, small n. Converting to 3.
  2

###Plot
#spplot(r73_pts["n_fires"], sp.layout=list("sp.polygons", perims, first=F)) #Option to plot w perims, don't use if gridded = TRUE
gridded(r73_pts) <- TRUE #For more efficient plotting, converts to "SpatialPixels"
pdf("../../GIS/Base Layers/TimesBurned.pdf", width = 3, height = 4)
spplot(r73_pts["n_fires"],
       col.regions = c("#0eb8f0", "#ffffbf", "#ff7f00","#ff7f00"),
       #col.regions = rainbow(100, start = 4/6, end = 1),
       #col.regions = brewer.pal(n = 5, name = "OrRd"),
       colorkey=list(at=c(0,1,2,4),
                     labels = list(at=c(0.5,1.5,3), labels = c("0","1","2-4")
                                   )
                     ),
       cuts = 3,
       main = "times burned 1973-2003"
       ) # "#0eb8f0", "#ffffbf", "#ff7f00"
#grid.text("times burned", x=unit(0.9, "npc"), y=unit(0.98, "npc"))
dev.off()

#r <- raster(r73_pts, "n_fires") #Only need to do once.
#p <- rasterToPolygons(r, dissolve = TRUE)
#writeOGR(p,"../../GIS/Raw Data/Sugarloaf Fire Perimeters","Nburned", driver = "ESRI Shapefile")
####2c. Analysis####

###Set up change analysis
sum_table <- 
  data.frame(veg_type = c(
    "shrub", "sparse meadow", "mixed conifer", "dense meadow"
  ))
scenarios <- c("watershed","no burns", "one burn", "two burns")
change_mat_total <- change_mat_exp <- gcp <- list()

for(s in scenarios){
  if(s == "watershed"){
    v73 <- getValues(r73scb)
    v14 <- getValues(r14scb)
    sum_table$pix_73 <- as.vector(table(v73[v73!=3 & v73!=5]))
    sum_table$pix_14 <- as.vector(table(v14[v14!=3 & v14!=5]))
  }
  if(s == "no burns"){
    v73 <- getValues(r73scb)[which(r73_pts$n_fires==0)]
    v14 <- getValues(r14scb)[which(r14_pts$n_fires==0)]
    sum_table$pix_73 <- as.vector(table(v73[v73!=3 & v73!=5]))
    sum_table$pix_14 <- as.vector(table(v14[v14!=3 & v14!=5]))
  }
  if(s == "one burn"){
    v73 <- getValues(r73scb)[which(r73_pts$n_fires==1)]
    v14 <- getValues(r14scb)[which(r14_pts$n_fires==1)]
    sum_table$pix_73 <- as.vector(table(v73[v73!=3 & v73!=5]))
    sum_table$pix_14 <- as.vector(table(v14[v14!=3 & v14!=5]))
  }
  if(s == "two burns"){
    v73 <- getValues(r73scb)[which(r73_pts$n_fires==2)]
    v14 <- getValues(r14scb)[which(r14_pts$n_fires==2)]
    sum_table$pix_73 <- as.vector(table(v73[v73!=3 & v73!=5]))
    sum_table$pix_14 <- as.vector(table(v14[v14!=3 & v14!=5]))
  }
  if(s == "three burns"){
    v73 <- getValues(r73scb)[which(r73_pts$n_fires==3)]
    v14 <- getValues(r14scb)[which(r14_pts$n_fires==3)]
    sum_table$pix_73 <- as.vector(table(v73[v73!=3 & v73!=5]))
    sum_table$pix_14 <- as.vector(table(v14[v14!=3 & v14!=5]))
  }
  
  change_mat_total[[s]] <- matrix(nrow = 6, ncol = 6)
  for(i in c(1,2,4,6)){ ##Observed matrix loop
    #1,2,4,6 = shrub, sparse meadow, mixed conifer, dense meadow
    change_mat_total[[s]][i,c(1,2,4,6)] <-
      c(length(which(v73==i & v14==1)),
        length(which(v73==i & v14==2)),
        length(which(v73==i & v14==4)),
        length(which(v73==i & v14==6))
      )
  }##End observed matrix loop
  
  change_mat_total[[s]] <- change_mat_total[[s]][-c(3,5),-c(3,5)]
  rownames(change_mat_total[[s]]) <- c("1973 shrub", "1973 sparse meadow", 
                                  "1973 mixed conifer", "1973 dense meadow")
  colnames(change_mat_total[[s]]) <- c("2014 shrub", "2014 sparse meadow", 
                                  "2014 mixed conifer", "2014 dense meadow")
  change_mat_exp[[s]] <- change_mat_total[[s]]
  
  for(i in 1:4){ ##Expected matrix loop
    if(i<4){
      change_mat_exp[[s]][i,i+1] <-
        change_mat_exp[[s]][i+1,i] <-
        round(mean(c(change_mat_exp[[s]][i,i+1],change_mat_exp[[s]][i+1,i])))
    }
    if(i<3){
      change_mat_exp[[s]][i,i+2] <-
        change_mat_exp[[s]][i+2,i] <-
        round(mean(c(change_mat_exp[[s]][i,i+2],change_mat_exp[[s]][i+2,i])))
    }
    if(i<2){
      change_mat_exp[[s]][i,i+3] <-
        change_mat_exp[[s]][i+3,i] <-
        round(mean(c(change_mat_exp[[s]][i,i+3],change_mat_exp[[s]][i+3,i])))
    }
  } ##End expected matrix loop
  
  
  #for(i in 1:ncol(change_mat_total)){
  #  change_mat_total[i,i] <- NA
  #} #Deprecated?
  
  X2 <- chisq.test(as.vector(change_mat_total[[s]]), 
                    p = as.vector(change_mat_exp[[s]]/sum(change_mat_exp[[s]])))
  print(X2)
  #deviance_ha <- round((change_mat_total[[s]]-change_mat_exp[[s]]) *0.16)  #0.16 ha/pixel
  deviance_prop <- (change_mat_total[[s]]-change_mat_exp[[s]]) / (change_mat_exp[[s]]+1)
  deviance_prop_melt <- list()
  deviance_prop_melt[[s]] <- melt(deviance_prop)
  names(deviance_prop_melt[[s]]) <- c("Y1973","Y2014","residual_prop")
  deviance_prop_melt[[s]]$Y1973 <- 
    factor(deviance_prop_melt[[s]]$Y1973, rev(levels(deviance_prop_melt[[s]]$Y1973) ) )
  deviance_prop_melt[[s]]$cells <- melt(change_mat_total[[s]])$value
  #Deprecated:
  #cp <- corrplot(deviance_ha, is.corr = FALSE)
  #gcp <- ggcorrplot(deviance_ha)
  #t <- tableHTML(change_mat_total)
  #t_grob <- tableGrob(change_mat_total)
  #grid.arrange(t_grob,ncol=1)
  #  ##Print Output
  #write_tableHTML(t,file = paste0("./Tables/Table1_",s,".html"))
  
  #Plot
  gcp[[s]] <- ggplot(deviance_prop_melt[[s]],aes(Y2014,Y1973)) +
    geom_tile(aes(fill = residual_prop),colour = "white") +
    geom_text(aes(label = cells)) +
    scale_fill_gradientn(colors = c("darkred", "white","darkblue"), limits = c(-1,1)) +
    theme_grey(base_size = 9) + 
    labs(title = s)+
    scale_x_discrete(position = "top") +
    theme(axis.ticks = element_blank(), axis.text.x = element_text(
      size = 9 *0.8, angle = 45, hjust = 0, colour = "grey50"))
  
  
}##End scenario loop

pdf("./Figures/MS/Fig3.pdf", width = 8, height = 6)
grid.arrange(gcp[[1]], gcp[[2]], gcp[[3]], gcp[[4]], 
             ncol = 2)
dev.off()



####3. Soil Moisture ####
library(readxl)
library(dplyr)
library(plotrix) #for std.err
library(ggplot2)
#library(raster) #for raster()
####3a. Data processing####
#r73scb <- #Load the processed 1973 veg raster
#  raster("./Processed Data/Classified Images/Final rasters/1973_raster_match_SCB_analysis.tif")
#r14scb <- #Load the processed 2014 veg raster
#  raster("./Processed Data/Classified Images/Final rasters/2014_raster_match_SCB_analysis.tif")
#perims <- readOGR("../../GIS/Raw Data/Sugarloaf Fire Perimeters/Sug_FRAP_FirePerims.shp")
sm1 <- 
  read_excel("./Raw Data/Soil Moisture/SoilMoistureForms_May2016_SEKI_complete_corrected.xlsx",
             sheet = "AllTDR")
sm1$Year <- "2016_05"
names(sm1)[names(sm1)=="Soil Sat. (%)"] = "sm"
sm2 <- 
  read_excel("./Raw Data/Soil Moisture/SoilMoistureForms_July2016_SEKI_complete_corrected.xlsx",
             sheet = "AllTDR")
sm2$Year <- "2016_07"
names(sm2)[names(sm2)=="Soil Sat. ( % )"] = "sm"
sm3 <- 
  read_excel("./Raw Data/Soil Moisture/SoilMoistureForms_SEKI_June2017_complete.xlsx",
             sheet = "AllTDR")
sm3$Year <- "2017_06"
names(sm3)[names(sm3)=="Soil Sat (%)"] = "sm"
names(sm3)[names(sm3)=="Veg (Conifer, Shrub, Sparse, Wetland)"] = "Veg"
sm4 <- 
  read_excel("./Raw Data/Soil Moisture/SoilMoistureForms_SEKI_July2017_complete.xlsx",
             sheet = "AllTDR")
sm4$Year <- "2017_07"
names(sm4)[names(sm4)=="Soil Sat (%)"] = "sm"
names(sm4)[names(sm4)=="Veg (Conifer, Shrub, Sparse, Wetland)"] = "Veg"
sm5 <- 
  read_excel("./Raw Data/Soil Moisture/SoilMoistureForms_SEKI_Jun2018_Final.xlsx",
             sheet = "AllTDR")
sm5$Year <- "2018_06"
names(sm5)[names(sm5)=="Soil Sat (%)"] = "sm"
names(sm5)[names(sm5)=="Veg (Conifer, Shrub, Sparse, Wetland)"] = "Veg"

sm <- rbind(
  sm1[,c("Year","LatPoint", "LonPoint", "sm" ,"Site","Veg")],
  sm2[,c("Year","LatPoint", "LonPoint", "sm" ,"Site","Veg")],
  sm3[,c("Year","LatPoint", "LonPoint", "sm" ,"Site","Veg")],
  sm4[,c("Year","LatPoint", "LonPoint", "sm" ,"Site","Veg")],
  sm5[,c("Year","LatPoint", "LonPoint", "sm" ,"Site","Veg")]
  
)

sm$sm[c(grep("water",sm$sm), grep("Water",sm$sm), grep("sat",sm$sm), grep("Sat",sm$sm),
        grep("999999900",sm$sm))] <- "100.00"
sm$sm <- as.numeric(sm$sm) #warning ok, we want to introduce NA's for text.
sm$sm[sm$sm>100 & !is.na(sm$sm)] <- 100 #Fix a few TDR glitches from wet soil or standng h2o

##Updated data read-in##
sm <- 
  read_excel("./Raw Data/Soil Moisture/SoilMoistureForms_SEKI_Combined_AllNums.xlsx",
             col_types = c("guess", rep("skip",2),"text", rep("skip",6),
                           rep("numeric",4), rep("guess",5), rep("skip",7)),
             col_names = c("MEASID","Site","LonPoint","LatPoint","sm","ID",
                           "Zone","Notes","Veg_Orig","Veg","Trip"),
             skip = 1, na = "NA",
             sheet = "AllTDR")
sm$Year <- NA
#sm$Veg <- NA #deprecate, added manually to spreadsheet
sm[grep("2016",sm$Trip),"Year"] <- "2016"
sm[grep("2017",sm$Trip),"Year"] <- "2017"
sm[grep("2018",sm$Trip),"Year"] <- "2018"
#sm$Year <- factor (sm$Year)

#Merge veg classes into target 4 classes:
sm$Veg[c(grep("Shrub", sm$Veg_Orig), grep("shrub",sm$Veg_Orig))] <- "shrub"
sm$Veg[c(grep("Wet meadow", sm$Veg_Orig),grep("Wet Meadow", sm$Veg_Orig),
         which((sm$Site=="GBW1" | sm$Site=="GBW3")
               & sm$Veg_Orig=="Meadow"))] <- "dense meadow"
sm$Veg[c(grep("grasses", sm$Veg_Orig), grep("open, meadow",sm$Veg_Orig),
         grep("meadow, open",sm$Veg_Orig),grep("Meadow, Open",sm$Veg_Orig),
         which(sm$Site=="GBM2"&sm$Veg_Orig=="Meadow"))] <- "sparse meadow"
sm$Veg[c(grep("Conifer", sm$Veg_Orig), grep("conifer",sm$Veg_Orig),
         grep("Forest",sm$Veg_Orig),grep("forest",sm$Veg_Orig))] <- "mixed conifer"


####3b. Analysis####
sum_table_sm <-
  sm[which(!is.na(sm$Veg)),] %>%
  group_by(Veg) %>%
  summarize(mean_sm = mean(sm, na.rm = T),
            sd_sm = sd(sm, na.rm = T),
            se_sm = std.error(sm, na.rm = T),
            n_plots = length(unique(Site))
  )

ggplot(sum_table_sm) +
  geom_col(aes(Veg, mean_sm))+
  geom_errorbar(aes(Veg, ymax = mean_sm + se_sm, ymin = mean_sm - se_sm),
                width = 0.5) +
  labs(y = "mean soil moisture (% VWC)") +
  theme_bw()

sum_table_sm_yr <-
  sm[which(!is.na(sm$Veg)),] %>%
  group_by(Veg, Year) %>%
  summarize(mean_sm = mean(sm, na.rm = T),
            sd_sm = sd(sm, na.rm = T),
            se_sm = std.error(sm, na.rm = T) ,
            sites = paste(unique(Site), collapse = ", ")
  )
dodge <- position_dodge(width=0.9)

pdf("./Figures/MS/Fig5_11_20.pdf", width = 8, height = 8)
ggplot(sum_table_sm_yr, aes(x = Veg, fill = Year)) +
  geom_col(aes(y=mean_sm), position = position_dodge())+
  geom_errorbar(aes(ymax = mean_sm + se_sm, ymin = mean_sm - se_sm),
                position = position_dodge()) +
  labs(y = "mean soil moisture (% VWC)") +
  theme_bw()
dev.off()

sum_table_sm_trip <-
  sm[which(!is.na(sm$Veg)),] %>%
  group_by(Veg, Trip) %>%
  summarize(mean_sm = mean(sm, na.rm = T),
            sd_sm = sd(sm, na.rm = T),
            se_sm = std.error(sm, na.rm = T) ,
            sites = paste(unique(Site), collapse = ", ")
  )

ggplot(sum_table_sm_trip, aes(x = Veg, fill = Trip)) +
  geom_col(aes(y=mean_sm), position = position_dodge())+
  geom_errorbar(aes(ymax = mean_sm + se_sm, ymin = mean_sm - se_sm),
                position = position_dodge()) +
  labs(y = "mean soil moisture (% VWC)")