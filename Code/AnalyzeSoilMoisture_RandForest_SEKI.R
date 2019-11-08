#Created by Gabrielle Boisrame
#Train models using measured soil moisture

#Be sure to set the working directory to this source file's location

AggData=1 #For grouping measurements that are the same veg type at the same site
Ghost=0   #For adding in "ghost" measurements to counteract effects of missing data
PICOsep=0  #For separating PICO from other conifers
Extrap=1 #If want to extrapolate to other non-measured points

library(nlme) # fit regression w/ spatially correlated errors
library(chron)
library(gstat) # classical geostatistics
library(MASS)
library(timeSeries)
library(sp)
library(matrixStats)
library(Hmisc)
library(plyr)
library(xts)
library(hydroGOF)
#library(tiff)

#Load data
SoilM <- read.csv('../Raw Data/Soil Moisture/SoilMoisture_SEKI_Combined_GISextract.csv',header=TRUE)    #('SoilMoistureMATLAB_All_10_02_15.csv')
#Choose early or late season
#SoilM<-SoilM[SoilM$DOY>180,] #Late Summer
#SoilM<-SoilM[(SoilM$DOY>153)&(SoilM$DOY<180),] #June

if(Extrap){
SoilMbig <- read.csv('../Processed Data/GIS/RasterToPointsValues.csv',header=TRUE)    #('SoilMoistureMATLAB_All_10_02_15.csv')
SoilMbig <- SoilMbig[SoilMbig$aspect>-2,] #Get rid of bad data points
SoilMbig <- SoilMbig[SoilMbig$SevNum<5,] #Get rid of bad data points
SoilMbig <- SoilMbig[(SoilMbig$Veg!=3)&(SoilMbig$Veg!=5)&(SoilMbig$Veg73!=3)&(SoilMbig$Veg73!=5)&SoilMbig$Veg>0,]

SoilMbig$Fire_Num[SoilMbig$Fire_Num<0] <- 0

SoilMbig$Veg14<-SoilMbig$Veg
SoilMbig$Veg[SoilMbig$Veg==1]<-"shrub"
#SoilMbig$Veg73[SoilMbig$Veg73==1]<-"shrub"
SoilMbig$Veg[SoilMbig$Veg==2]<-"sparse meadow"
#SoilMbig$Veg73[SoilMbig$Veg73==2]<-"sparse meadow"
SoilMbig$Veg[SoilMbig$Veg==4]<-"mixed conifer"
#SoilMbig$Veg73[SoilMbig$Veg73==4]<-"mixed conifer"
SoilMbig$Veg[SoilMbig$Veg==6]<-"dense meadow"
#SoilMbig$Veg73[SoilMbig$Veg73==6]<-"dense meadow"

SoilMbig$Veg<-as.factor(SoilMbig$Veg)
SoilMbig$Veg73<-as.factor(SoilMbig$Veg73)

  AspInd=.5*(1-cos(pi*(SoilMbig[,'aspect']-30)/180)) 
  AspInd[SoilMbig[,'aspect']<0]=0
  SoilMbig$AspectDeg=SoilMbig$aspect
  SoilMbig$Aspect=AspInd

    SoilMbig$Time_Since_Fire <- (2016-SoilMbig$Fire_Year)
    SoilMbig$Time_Since_Fire[SoilMbig$Time_Since_Fire>100]=100
    
      SoilMbig$TWI.10m<-log(SoilMbig$flow_acc_d/max(0.0001,tan(SoilMbig$slope_deg*pi/180)))
      SoilMbig$TWI.10m[SoilMbig$TWI.10m<0]=0

        SoilMbig$Upslope.Area<-SoilMbig$flow_acc_d
        SoilMbig$Upslope.Area[SoilMbig$flow_acc_d>5000]=5000
        
 hist(SoilMbig$slope_deg[SoilMbig$Veg14==4],freq=FALSE,col = rgb(0,1,0,.5))
 hist(SoilMbig$slope_deg[SoilMbig$Veg14==6],freq=FALSE,col = rgb(0,0,1,.5),add=TRUE)

 hist(SoilMbig$TWI.10m[SoilMbig$Veg14==4],freq=FALSE,col = rgb(0,1,0,.5))
 hist(SoilMbig$TWI.10m[SoilMbig$Veg14==6],freq=FALSE,col = rgb(0,0,1,.5),add=TRUE)
}

SoilM <-SoilM[SoilM$Site!='Calibration',] #Remove a few calibration points


thetaM <- as.numeric(SoilM$Soil_Sat)

#Remove NaNs (rocks)
SoilM=SoilM[!is.na(thetaM),]
thetaM=thetaM[!is.na(thetaM)]
SoilM=SoilM[thetaM>=0,]
thetaM=thetaM[thetaM>=0]

if (PICOsep){
  SoilM$Veg14[SoilM$VegNum==3.4]=1.5 #Separate out PICOs
}

#Transform aspect to be index from 0 to 1
AspInd=.5*(1-cos(pi*(SoilM[,'aspect']-30)/180)) 
AspInd[SoilM[,'aspect']<0]=0
SoilM$AspectDeg=SoilM$aspect
SoilM$Aspect=AspInd


#Get dates playing nice
SoilM$Date <- as.Date(SoilM$TimeStamp,format="%m/%d/%Y")
SoilM$Year <- 2015+as.numeric(years(SoilM$Date))
SoilM$Month <- months(SoilM$Date)
SoilM$DOY <- as.numeric(1+SoilM$Date-as.Date(sprintf('01/01/%i',SoilM$Year),format="%m/%d/%Y"))

#Make sure numbers are numbers
SoilM$Fire_Year<-as.integer(SoilM$Fire_Year)
SoilM$Fire_Num<-as.integer(SoilM$Fire_Num)
SoilM$Fire_Num[SoilM$Fire_Num<0]=0
#SoilM$SevNum<-as.integer(SoilM$SevNum)

SoilM$Veg73<-as.numeric(SoilM$X1973_veg)
SoilM$Veg14<-as.numeric(SoilM$X2014_veg)

#Get time since fire
SoilM$Time_Since_Fire <- (SoilM$Year-SoilM$Fire_Year)
SoilM$Time_Since_Fire[SoilM$Time_Since_Fire>100]=100 #Set max years since fire to 100



if (AggData){ 

  source("AggSubSites.R")
  SoilMa<-AggSubSites(SoilM)
  
  SoilM=SoilMa
  thetaM=SoilM$Soil_Sat
}

#VarName='Veg14'
#plot(SoilMa[,VarName],sqrt(SoilMVar[,VarName]),xlab='Mean',ylab='Std Dev',main=VarName)



#thetaM=log(thetaM)
#thetaM[thetaM<(-3)]=-3  #Get rid of neg infinity
#thetaM=sqrt(thetaM)
hist(thetaM)


#Get mean VWC for each site


#Add Ghost Sites 
if(Ghost){
  SoilM$VWC<-thetaM
SoilM_o = SoilM #Hold a copy of original data
source("AddGhosts.R")
SoilM<-AddGhosts(SoilM)

thetaM<-SoilM$VWC

thetaM<-SoilM$Soil_Sat

}

if(AggData==0){
  SiteMeansA<-aggregate(thetaM,list(SoilM$Site,SoilM$Veg),mean)
  SiteMeans<-aggregate(SiteMeansA$x,list(SiteMeansA$Group.1),max)
  SiteCounts<-aggregate(thetaM,list(SoilM$Site),length)
  #SiteMeans<-aggregate(thetaM,list(SoilM$SiteNum),mean)
  #SiteMeans<-aggregate(thetaM,list(SoilM$SiteNum),quantile,.6)
  #Divide sites into groups
}else{
  SiteMeansA<-SoilM[,c('SubSite','Soil_Sat')]
  SiteMeans<-aggregate(SiteMeansA$Soil_Sat,list(round(SiteMeansA$SubSite)),max)
  SiteCounts<-aggregate(thetaM,list(SoilM$Site),length)
}


GrpVWC<-c(10,25,50,60)
SoilM$VegChange<-(10*SoilM$Veg73+SoilM$Veg14)

#Calculate topographic wetness index
SoilM$flow_acc_d[SoilM$flow_acc_d<0]=0 #Replace noData values with zeros
SoilM$slope_deg[SoilM$slope_deg<0]=0 #Replace noData values with zeros

SoilM$TWI.10m<-log(SoilM$flow_acc_d/max(0.0001,tan(SoilM$slope_deg*pi/180)))
SoilM$TWI.10m[SoilM$TWI.10m<0]=0



#SiteVeg<-aggregate(SoilM$VegChange,list(SoilM$Site),median)
#GrpVeg<-sort(unique(SiteVeg$x))

SiteMeans$WetCat=1
SiteMeans$CountCat=0
SiteMeans$CountCat[SiteMeans$x<=GrpVWC[1]]<-sum(SiteMeans$x<=GrpVWC[1])

for(i in 1:(length(GrpVWC)-1)){
  SiteMeans$WetCat[SiteMeans$x>GrpVWC[i]]<-(i+1)
  SiteMeans$CountCat[SiteMeans$x>GrpVWC[i]]<-sum((SiteMeans$x>GrpVWC[i])&(SiteMeans$x<=GrpVWC[i+1]))
}

hist(SiteMeans$WetCat) #Look at how many sites are in each wetness category
plot(SiteMeans$WetCat,SiteMeans$x) #Look at mean site moisture within each wetness category
plot(SoilM$VegChange,thetaM,xlim=c(10,70))
#hist(SoilM$SiteNum)

#Correlation Matrix
if(AggData){
  CM=cor(SoilM[,c(5:15,26:34,37)])
}else{
CM=cor(SoilM[,c(12:14,32:43)])}
.01*round(CM*100)
sum((abs(CM)>.7)&(CM<1))/2
sum((abs(CM)>.4)&(CM<1))/2
((abs(CM)>.4)&(CM<1))


#SoilM1<-SoilM[SoilM$Year==2016,c(12:14,28:38)]

#SoilM1<-SoilM[SoilM$Year==2016,c(12:14,32:43)]

#cor(SoilM1)


#EarlyInd<-(SoilM$DOY<=160)
#MidInd<-(SoilM$DOY>170 & SoilM$DOY<190)
#LateInd<-(SoilM$DOY > 190)


SoilM$SevNum=SoilM$Fire_Max_Sev
#For now: set fire to minimum severity for fires where don't have severity data. 
SoilM$SevNum[(SoilM$Fire_Max_Sev<0)&SoilM$Fire_Num>0]=1 
SoilM$SevNum<-as.factor(SoilM$SevNum)

#Add some veg categories based on pre-fire veg.
#SoilM$Veg14b=SoilM$Veg14
#SoilM$Veg14b[SoilM$Veg14==4]=5
#SoilM$Veg14b[SoilM$VegChange==44]=6 #New wetlands =5, persistent wetlands = 6
#SoilM$Veg14b[SoilM$Veg14==5]=7
#SoilM$Veg14b[SoilM$VegChange==55]=8 #New aspen =7, persistent aspen =8
#SoilM$Veg14b[SoilM$Veg14==1.5]=2
#SoilM$Veg14b[SoilM$Veg14==2]=3
#SoilM$Veg14b[SoilM$Veg14==3]=4


thetaMs<-sqrt(thetaM)

library(randomForest)

SoilM$VWC<-thetaM

SoilM<-SoilM[SoilM$VegChange>0] #FOR NOW remove bad data points. Fix data later.

#SoilM<-SoilM[SoilM$VegChange>0] #FOR NOW remove bad data points. Fix data later.


#Set certain values as factors, since order doesn't matter.
SoilM$Veg14<-as.factor(SoilM$Veg14)
SoilM$Veg73<-as.factor(SoilM$Veg73)
SoilM$Year<-as.factor(SoilM$Year)


#load('ICB_RandomTree_26_04_17_AggSites_NoPICO_NoGhosts.rdata')
#load('ICB_RandomForest_16_03_17_SubSiteMeans.rdata')
load('ICB_RandomForest_2019_03_01_AggData.rdata')
Tfit_ICB<-Tfit

#load('SoilM_3_11_16.rdata')
#load('thetaM_3_11_16.rdata')

#thetaM<-SoilM$VWC

SoilM$Upslope.Area<-SoilM$flow_acc_d
SoilM$Upslope.Area[SoilM$flow_acc_d>5000]=5000



##Un-comment these if want to normalize
#SoilM$TimeSevInd<-(1-(SoilM$Time.Since.Fire/100))*(SoilM$SevNum/4)
#TempTPI<-(SoilM$TPI300m-min(SoilM$TPI300m))/(max(SoilM$TPI300m)-min(SoilM$TPI300m))
#TempSlope<-(SoilM$Slope-min(SoilM$Slope))/(max(SoilM$Slope)-min(SoilM$Slope))
#TempArea<-SoilM$Upslope.Area
#TempArea[TempArea>600]=625
#TempArea<-(TempArea-min(SoilM$Upslope.Area))/(625-min(SoilM$Upslope.Area))
#SoilM$TopInd<-.5*TempTPI+.3*TempSlope-0.2*TempArea
#SoilM$TopInd<-(.5*SoilM$TPI300m+.3*SoilM$Slope-0.2*SoilM$Upslope.Area)

sum(is.na(SoilM$Dist_From_River))

Tfit<-randomForest(VWC~Veg+Veg73+Year+DOY+Upslope.Area+slope_deg+Aspect+tpi_300m+TWI.10m+Time_Since_Fire+Fire_Num+SevNum+Elevation+Dist_From_River,data=SoilM,nodesize=5,ntree=500)


imp=Tfit$importance
barplot(xlab='Variable',ylab='Importance',names.arg=rownames(imp)[order(imp[, 1], decreasing=TRUE)],height=imp[order(imp[, 1], decreasing=TRUE),1],cex.names=.5)

a=partialPlot(x=Tfit,pred.data=SoilM,x.var='Year',ylab='VWC')
barplot(a$y,names=a$x,xlab='',ylab='VWC',ylim=c(0,20))
barplot(a$y[order(a$y, decreasing=FALSE)],names=a$x[order(a$y, decreasing=FALSE)],ylim=c(0,20),xlab='',ylab='VWC')
barplot(a$y[c(3,1,2,4)],names=a$x[c(3,1,2,4)],xlab='Year',ylab='VWC',ylim=c(0,20))


a=partialPlot(x=Tfit,pred.data=SoilM,x.var='Upslope.Area',ylab='VWC',ylim=c(0,20)) #TWI.10m, TPI300m, ='Dist.from.River')
#plot(a$x,a$y,ylim=c(0,20),lwd=2,xlab='Years Since Fire',ylab='Mean VWC',main='Modeled Effect of Variable on VWC',type='l')
points(a$x,a$y)

if(Extrap){
  SoilMbig$Year=factor(2018,levels=levels(SoilM$Year))
  SoilMbig$DOY=153 #205 #
  SoilMbig$SevNum<-as.factor(SoilMbig$SevNum)
  
  PredTreeBig<-predict(Tfit,SoilMbig,predict.all=TRUE)
  AllTrees<-PredTreeBig$individual
  StdTrees1<-apply(AllTrees,1,sd) #Get the standard deviation of the predicted values over all trees
  Q25Trees1<-apply(AllTrees,1,quantile,.25) #Get the standard deviation of the predicted values over all trees
  Q75Trees1<-apply(AllTrees,1,quantile,.75) #Get the standard deviation of the predicted values over all trees
  PredToday<-PredTreeBig$aggregate
  
  SoilMold=SoilMbig
  SoilMold$Veg14<-SoilMold$Veg73
  SoilMold$Time_Since_Fire<-100
  SoilMold$SevNum<-factor(0,levels=levels(SoilM$SevNum))
  SoilMold$Fire_Num<-0
  
  PredTreeOld<-predict(Tfit,SoilMold,predict.all=TRUE)
  AllTrees<-PredTreeOld$individual
  StdTrees2<-apply(AllTrees,1,sd) #Get the standard deviation of the predicted values over all trees
  Q25Trees2<-apply(AllTrees,1,quantile,.25) #Get the standard deviation of the predicted values over all trees
  Q75Trees2<-apply(AllTrees,1,quantile,.75) #Get the standard deviation of the predicted values over all trees
  PredUnburned<-PredTreeOld$aggregate
  
  plot(PredUnburned,PredToday)
  
  Fchange<-(SoilMbig$Veg73!=as.factor(SoilMbig$Veg14))
  plot(PredUnburned[Fchange],PredToday[Fchange],xlab="Modeled Unburned Soil Moisture (%)",ylab="Modeled Actual Soil Moisture (%)",main="June 2018",xlim=c(3,45),ylim=c(3,45))
  points(PredUnburned[SoilMbig$Veg73=='4'&SoilMbig$Veg14==6],PredToday[SoilMbig$Veg73=='4'&SoilMbig$Veg14==6],pch=15,col='blue') #Con-Wet
  points(PredUnburned[SoilMbig$Veg73=='4'&SoilMbig$Veg14==2],PredToday[SoilMbig$Veg73=='4'&SoilMbig$Veg14==2],pch=16,col='grey') #Con-Sparse
  points(PredUnburned[SoilMbig$Veg73=='6'&SoilMbig$Veg14==4],PredToday[SoilMbig$Veg73=='6'&SoilMbig$Veg14==4],pch=17,col=rgb(.2,.75,.2)) #Wet-Con
  points(PredUnburned[SoilMbig$Veg73=='4'&SoilMbig$Veg14==1],PredToday[SoilMbig$Veg73=='4'&SoilMbig$Veg14==1],pch=18,col='brown') #Con-Shrub
  lines(c(3,45),c(3,45),lty=2,lwd=2)
  legend(3,43,c("Conifer - Dense Mdw.","Conifer - Sparse Mdw.","Conifer-Shrub","Dense Mdw. - Conifer","Other"),pch=c(15,16,18,17,1),col=c("blue","grey","brown",rgb(.2,.75,.2),"black"))
  
  hist(PredToday[Fchange]-PredUnburned[Fchange],breaks=seq(-4.25,3.75,.5),main="June 2018",xlab="Actual-Unburned Soil Moisture (%)")
  hist(100*(PredToday[Fchange]-PredUnburned[Fchange])/PredUnburned[Fchange],breaks=seq(-62.5,62.5,5),main="July 2016",xlab="% Change in Actual-Unburned Soil Moisture")
  mean(PredToday[Fchange]-PredUnburned[Fchange])
}

PredTree<-predict(Tfit,SoilM,predict.all=TRUE)
AllTrees<-PredTree$individual
StdTrees<-apply(AllTrees,1,sd) #Get the standard deviation of the predicted values over all trees
PredTree<-PredTree$aggregate
plot(.01*SoilM$VWC,.01*PredTree)
errbar(.01*SoilM$VWC,.01*PredTree,.01*PredTree+.01*StdTrees,.01*PredTree-.01*StdTrees,xlab='Measured',ylab='Predicted',main='Random Forest Model Fit to Test Data',errbar.col='grey')
points(.01*SoilM$VWC,.01*PredTree)
lines(c(0,.65),c(0,.65),col='red')


#Separate by trip and veg
TripVegMat<-data.frame(cbind(Year=as.factor(c(2016,2016,2017,2017,2018)),Trip=c(1,2,1,2,1)),DenseMeadow=c(1:5),MixCon=c(1:5),Shrub=c(1:5),Sparse=c(1:5))
if(AggData){SoilMpred<-SoilM}else{SoilMpred<-SoilM[,5:53]}
SoilMpred<-SoilMpred[SoilMpred$Trip=='2017_Early',]
TempFrYrs<-(SoilMpred$Time_Since_Fire-1)
Yrs<-unique(SoilM$Year)
#SoilMpred$DOY<-180
#SoilMpred$Year<-2016
#SoilMpred<-unique(SoilMpred) #DOESNT WORK. I THINK BECAUSE TIME SINCE FIRE IS TIED TO YEAR
for (i in 1:5){
 Yr<-Yrs[ceiling(i/2)] 
 Trp<-TripVegMat$Trip[i]
SoilMpred$Year<-Yr
if(Trp==1){Dt<-178}else{Dt<-205}
SoilMpred$DOY<-Dt
SoilMpred$Time_Since_Fire<-(TempFrYrs+ceiling(-1+i/2))
SoilMpred$Time_Since_Fire[SoilMpred$Time_Since_Fire>90]<-100
a=partialPlot(x=Tfit,pred.data=SoilMpred,x.var='Veg',ylab='VWC')
# if(Trp==1){Dt<-(SoilM$DOY<190)}else{Dt<-SoilM$DOY>190}
#a=partialPlot(x=Tfit,pred.data=SoilM[SoilM$Year==Yr & Dt,],x.var='Veg',ylab='VWC')  #If want to look at modeled means for only sites measured in each trip
#barplot(0.01*a$y,names=a$x,ylim=c(0,.2),xlab='Dominant Veg',ylab='Mean VWC',main='Modeled Effect of Variable on VWC')
a$x
TripVegMat[i,3:6]<-.01*a$y
}
barplot(100*as.matrix(TripVegMat[,3:6]),beside=TRUE,ylab='Mean Volumetric Water Content (%)',legend.text=c("June 2016","July 2016","June 2017","July 2017","June 2018"),names.arg = c("Dense Meadow","Conifer","Shrub","Sparse"))


#Test alternate realities in terms of veg cover and fire history
SoilM_allburn<-SoilM
VegOps<-levels(SoilM$Veg)
SoilM_allburn$Veg[SoilM$Veg!='sparse meadow']<-VegOps[4]
#levels(SoilM_allburn$Veg)<-VegOps
SevOps<-levels(SoilM$SevNum)
SoilM_allburn$SevNum[SoilM$Veg!='sparse meadow']<-factor(SevOps[5])
#levels(SoilM_allburn$SevNum)<-SevOps
SoilM_allburn$Time_Since_Fire[(SoilM_allburn$Time_Since_Fire>90)&(SoilM$Veg!='sparse meadow')]<-13
SoilM_allburn$Fire_Num[(SoilM_allburn$Fire_Num==0)&(SoilM$Veg!='sparse meadow')]<-1
PredSparse<-predict(Tfit,SoilM_allburn)
PredMeas<-predict(Tfit,SoilM)
plot(PredMeas,PredSparse,type='p',xlim=c(0,60),ylim=c(0,60))
points(PredMeas[SoilM$Veg=='dense meadow'],PredSparse[SoilM$Veg=='dense meadow'])
points(PredMeas[SoilM$Veg=='dense meadow'],PredSparse[SoilM$Veg=='dense meadow'],col='blue')
points(PredMeas[SoilM$Veg=='sparse meadow'],PredSparse[SoilM$Veg=='sparse meadow'],col='red')
points(PredMeas[SoilM$Veg=='mixed conifer'],PredSparse[SoilM$Veg=='mixed conifer'],col='green')
lines(c(0,60),c(0,60))

SoilM_noburn<-SoilM
VegOps<-levels(SoilM$Veg)
SoilM_noburn$Veg[SoilM$Veg!='mixed conifer']<-VegOps[2]
#levels(SoilM_noburn$Veg)<-VegOps
SevOps<-levels(SoilM$SevNum)
SoilM_noburn$SevNum[SoilM$SevNum!=0]<-factor(0)
levels(SoilM_noburn$SevNum)<-SevOps
SoilM_noburn$Time_Since_Fire[(SoilM_noburn$Time_Since_Fire<90)]<-100
SoilM_noburn$Fire_Num[(SoilM_noburn$Fire_Num==0)]<-0
PredNoburn<-predict(Tfit,SoilM_noburn)
plot(PredMeas,PredNoburn,type='p',xlim=c(0,60),ylim=c(0,60))
points(PredMeas[SoilM$Veg=='dense meadow'],PredNoburn[SoilM$Veg=='dense meadow'])
points(PredMeas[SoilM$Veg=='dense meadow'],PredNoburn[SoilM$Veg=='dense meadow'],col='blue')
points(PredMeas[SoilM$Veg=='sparse meadow'],PredNoburn[SoilM$Veg=='sparse meadow'],col='red')
points(PredMeas[SoilM$Veg=='mixed conifer'],PredNoburn[SoilM$Veg=='mixed conifer'],col='green')
lines(c(0,60),c(0,60))


#Compare models from ICB and Sugarloaf

SoilM_match<-SoilM
SoilM_match$veg12<-SoilM$X2014_veg
 SoilM_match$veg12[SoilM$Veg=='dense meadow']=4
 SoilM_match$veg12[SoilM$Veg=='sparse meadow']=3
 SoilM_match$veg12[SoilM$Veg=='mixed conifer']=1
 SoilM_match$veg12[SoilM$Veg=='shrub']=2
#SoilM_match$veg12<-SoilM$X2014_veg
# SoilM_match$veg12[SoilM$X2014_veg==6]=4
# SoilM_match$veg12[SoilM$X2014_veg==1]=2
# SoilM_match$veg12[SoilM$X2014_veg==4]=1
# SoilM_match$veg12[SoilM$X2014_veg==2]=3
SoilM_match$veg69<-SoilM$X1973_veg
 SoilM_match$veg69[SoilM$X1973_veg==6]=4
 SoilM_match$veg69[SoilM$X1973_veg==1]=2
 SoilM_match$veg69[SoilM$X1973_veg==4]=1
 SoilM_match$veg69[SoilM$X1973_veg==2]=3
SoilM_match$Slope<-SoilM$slope_deg
SoilM_match$TPI300m<-SoilM$tpi_300m
SoilM_match$Dist.from.River<-SoilM_match$Dist_From_River
SoilM_match$Time.Since.Fire<-SoilM$Time_Since_Fire
#SoilM_match$Times.Burned<-SoilM$Fire_Num
SoilM_match$Times.Burned<-as.factor(SoilM_match$Fire_Num)
SoilM_match$Elev<-SoilM$Elevation
#SoilM_match<-SoilM_match[SoilM_match$Year==2016,]
SoilM_match<-SoilM_match[SoilM_match$veg12>0,]
SoilM_match$veg12<-factor(SoilM_match$veg12,levels=c(1:5))
SoilM_match$veg69<-factor(SoilM_match$veg69,levels=c(1:5))
SoilM_match$Year<-factor(SoilM_match$Year,levels=c(2014:2018))
#SoilM_match$Year<-as.numeric(SoilM_match$Year)
SoilM_match$SevNum<-as.numeric(SoilM_match$SevNum)
SoilM_match$SevNum<-(SoilM_match$SevNum-2)
SoilM_match$SevNum[SoilM_match$SevNum<0]<-0
SoilM_match$SevNum<-as.factor(SoilM_match$SevNum)

VWCpred_ICB<-predict(Tfit_ICB,SoilM_match)

VWCpred_SL<-predict(Tfit,SoilM)
plot(SoilM$VWC,VWCpred_SL,xlab='Measured SCB Moisture',ylab='Modeled')
lines(c(0,55),c(0,55),col='grey')
points(SoilM_match$VWC,VWCpred_ICB,col='red')
legend('bottomright',c('SCB Model','ICB Model'),col=c('black','red'),pch=1)

cor(SoilM$VWC,VWCpred_SL)
cor(SoilM_match$VWC,VWCpred_ICB)
  
hist(VWCpred_SL-SoilM$VWC,col=rgb(0,0,0,.5),main='Model Error',alpha=.5,xlim=c(-40,40),xlab='Modeled-Measured')
hist(VWCpred_ICB-SoilM_match$VWC,col=rgb(1,0,0,.5),add=TRUE)  
  

#Select random training points from all points:
#nt=round(length(thetaM)/4)
#Trn<-sample(1:length(thetaM),nt*3,replace=FALSE)



#Select random sites for training
ns=nrow(SiteMeans)
pctX=.70 #Proportion of sites to use
numX=min(SiteMeans$CountCat)-1 #Number of sites to use from each category
ByPct=1
ByVWC=0
RandSites=1 #Don't have both RandSites and ByVWC==1

#SoilMall=SoilM
#SoilM=SoilM[SoilM$Year==2016,]

#
NumSamps = 100 #Number of times to run with random training set
CorTrain = 0 #Hold correlation for each training set
CorTest = 0 #Hold correlation for each test set
RMSEtrain = 0
RMSEtest = 0
Top3=matrix(nrow=3,ncol=NumSamps)
YearMeanVWC = matrix(nrow=3,ncol=NumSamps) #Hold mean VWC for each year for each training set

SoilM$SiteNum<-floor(SoilM$SubSite)

for(g in 1:NumSamps){

if(RandSites){
  nt=round(pctX*ns)
  TrnS<-sample(1:ns,nt,replace=FALSE) #Select X% of the sites
  Trn<-which(SoilM$SiteNum %in% TrnS)   #Find which measurements are within those sites
  
}else{ if(ByVWC){
#Group by VWC
ng=length(GrpVWC)
#nt=round(ns*pctX)
Trn=c()
for(j in 1:ng){
  SubGp<-SiteMeans[SiteMeans$WetCat==j,]
  nsg=length(SubGp$x)
  if(ByPct){
  nt=round(nsg*pctX)}
  else{
    nt=numX
  }
  TrnS<-sample(1:nsg,nt,replace=FALSE)
  TrnT<-which(SoilM$SiteNum %in% SubGp$Group.1[sort(TrnS)])
  Trn<-c(Trn,TrnT)
}
}else{
#Select sites within veg groups instead of soil moisture
ng=length(GrpVeg)
#nt=round(ns*pctX)
Trn=c()
for(j in GrpVeg){
  SubGp<-SiteVeg[SiteVeg$x==j,]
  nsg=length(SubGp$x)
  nt=round(nsg*pctX)
  TrnS<-sample(1:nsg,nt,replace=FALSE)
  TrnT<-which(SoilM$SiteNum %in% SubGp$Group.1[sort(TrnS)])
  Trn<-c(Trn,TrnT)
}
}
}


Trn<-sort(Trn)
Tstx<-SoilM[-Trn,]
Tsty<-SoilM$VWC[-Trn]
Trnx<-SoilM[Trn,]
Trny<-SoilM$VWC[Trn]

length(Trny)/length(thetaM) #Proportion of msmts used in training 

#Tfit<-randomForest(thetaM~Aspect+Time.Since.Fire+SevNum+Slope+Elev+TWI.10m+Upslope.Area+Dist.from.River+VegNum+TPI300m+POINT_X+POINT_Y+TimeSevInd+TopInd,data=SoilM,strata="VegNames",nodesize=20)
#Tfit<-randomForest(thetaM~Aspect+Time.Since.Fire+SevNum+Slope+Elev+Upslope.Area+Dist.from.River+VegNum+TPI300m+POINT_X+POINT_Y+TimeSevInd+TopInd,data=SoilM,nodesize=20)
#Tfit<-randomForest(thetaM~DOY+VegNum+Upslope.Area+Slope+Aspect+TPI300m+TWI.10m+Dist.from.River+Time.Since.Fire+SevNum,data=SoilM,nodesize=5)
#Tfit<-randomForest(thetaM~VegNum+Upslope.Area+Slope+Aspect+TPI300m+TWI.10m+Dist.from.River+Time.Since.Fire+SevNum,data=SoilM,nodesize=5,ntree=100)
#Tfit<-randomForest(thetaM~VegNum+Year+DOY+Upslope.Area+Slope+Aspect+TPI300m+TWI.10m+Dist.from.River+Time.Since.Fire+SevNum+LatPoint+LonPoint,data=SoilM,subset=Trn,nodesize=5,ntree=500)
#Tfit<-randomForest(thetaM~VegNum+Year+DOY+Upslope.Area+Slope+Aspect+TPI300m+TWI.10m+Dist.from.River+Time.Since.Fire+Times.Burned+SevNum+LatPoint+LonPoint+veg69,data=SoilM,subset=Trn,nodesize=5,ntree=500)
#With no "observed" veg, latitude, or longitude:

#No Fire Data:
#Tfit<-randomForest(VWC~Veg14+Year+DOY+Upslope.Area+Slope+Aspect+TPI300m+TWI.10m+Dist.from.River+Elev+veg69,data=SoilM,subset=Trn,nodesize=5,ntree=500)
#All Data:
#NOTE: the best fit uses both Veg14 and observed veg, but that can't be upscaled as reliably
#No veg69, but more groups of Veg14:
#Tfit<-randomForest(VWC~Veg14b+Year+DOY+Upslope.Area+Slope+Aspect+TPI300m+TWI.10m+Dist.from.River+Time.Since.Fire+Times.Burned+SevNum+Elev,data=SoilM,subset=Trn,nodesize=5,ntree=500)

Tfit<-randomForest(VWC~Veg+Veg73+Year+DOY+Upslope.Area+slope_deg+Aspect+tpi_300m+TWI.10m+Time_Since_Fire+Fire_Num+SevNum+Elevation+Dist_From_River,data=SoilM,subset=Trn,nodesize=5,ntree=500)
#No Fire or Year: Tfit<-randomForest(VWC~Veg+Veg73+DOY+Upslope.Area+slope_deg+Aspect+tpi_300m+TWI.10m+Elevation+Dist_From_River,data=SoilM,subset=Trn,nodesize=5,ntree=500)


if(g==1 || g==NumSamps){
plot(Tfit)
a=partialPlot(x=Tfit,pred.data=SoilM,x.var='TWI.10m',ylab='VWC') #TWI.10m, TPI300m, ='Dist.from.River')
plot(a$x,.01*a$y,ylim=c(0,.2),lwd=2,xlab='TWI',ylab='Mean VWC',main='Modeled Effect of Variable on VWC',type='l')
points(a$x,.01*a$y)
lf<-lm(.01*a$y~a$x)
lines(c(0,14),lf$coefficients[1]+lf$coefficients[2]*c(0,14))

a=partialPlot(x=Tfit,pred.data=SoilM,x.var='Veg14',ylab='VWC')
barplot(0.01*a$y,names=a$x,ylim=c(0,.2),xlab='Dominant Veg',ylab='Mean VWC',main='Modeled Effect of Variable on VWC')
}

PredTree2<-predict(Tfit,Tstx,predict.all=TRUE)
AllTrees<-PredTree2$individual
StdTrees<-apply(AllTrees,1,sd) #Get the standard deviation of the predicted values over all trees
Q25Trees<-apply(AllTrees,1,quantile,.25) #Get the standard deviation of the predicted values over all trees
Q75Trees<-apply(AllTrees,1,quantile,.75) #Get the standard deviation of the predicted values over all trees

PredTree2<-PredTree2$aggregate
MeanTrees<-apply(AllTrees,1,mean)#Just a test: make sure get same as PredTree2$aggregate
max(abs(MeanTrees-PredTree2))

PredThetaT<-predict(Tfit,Trnx,predict.all=TRUE)
AllTtrees<-PredThetaT$individual
Q25Ttrees<-apply(AllTtrees,1,quantile,.25) #Get the standard deviation of the predicted values over all trees
Q75Ttrees<-apply(AllTtrees,1,quantile,.75) #Get the standard deviation of the predicted values over all trees
PredThetaT<-PredThetaT$aggregate

CorTrain[g] = cor(PredThetaT,Trny)
#plot(Trny,PredThetaT)

CorTest[g] = cor(Tsty,PredTree2)

RMSEtrain[g] = rmse(PredThetaT,Trny)
RMSEtest[g] = rmse(PredTree2,Tsty)

imp=Tfit$importance
Tlist=rownames(imp)[order(imp[, 1], decreasing=TRUE)]
Top3[,g]=Tlist[1:3]

a=partialPlot(x=Tfit,pred.data=SoilM,x.var='Year',ylab='VWC')
YearMeanVWC[,g]=.01*a$y

if(g==1){
  
errbar(.01*Trny,.01*PredThetaT,.01*Q75Ttrees,.01*Q25Ttrees,xlab='Measured',ylab='Fit to Training Data',main='Random Forest Model Fit to Training Data',errbar.col='grey')
lines(c(0,.55),c(0,.55),col='red')
  
#plot(.01*Tsty,.01*PredTree2,xlab='Measured',ylab='Predicted',main='Random Forest Model Fit to Test Data')
errbar(.01*Tsty,.01*PredTree2,.01*PredTree2+.01*StdTrees,.01*PredTree2-.01*StdTrees,xlab='Measured',ylab='Predicted',main='Random Forest Model Fit to Test Data',errbar.col='grey')
points(.01*Tsty,.01*PredTree2)
lines(c(0,.55),c(0,.55),col='red')

errbar(.01*Tsty,.01*PredTree2,.01*Q75Trees,.01*Q25Trees,xlab='Measured',ylab='Predicted',main='Random Forest Model Fit to Test Data',errbar.col='grey')
points(.01*Tsty,.01*PredTree2)
lines(c(0,.55),c(0,.55),lty=2,col='grey')

hist(.01*(PredTree2-Tsty),main='Random Forest Model Fit to Test Data',xlab='VWC Error (Predicted-Measured)')
sqrt(sum((PredTree2-Tsty)^2)/length(Tsty))
barplot(xlab='Variable',ylab='Importance',names.arg=rownames(imp)[order(imp[, 1], decreasing=TRUE)],height=imp[order(imp[, 1], decreasing=TRUE),1],cex.names=.5)

truehist(StdTrees-abs(PredTree2-Tsty)) #Look at std dev of predictions vs error
} 
}

truehist(CorTrain,xlim=c(0,1))
truehist(CorTest,xlim=c(0,1))
truehist(.01*RMSEtrain,xlim=c(0,.3))
truehist(.01*RMSEtest,xlim=c(0,.3))
Top3  #Show top 3 variables in terms of importance for each run

YM<-rowMeans(YearMeanVWC)
Yrs<-sort(unique(SoilM$Year))
Ystds<-rowSds(YearMeanVWC)
barplot(YM,names=Yrs,ylim=c(0,.2),xlab='Year',ylab='Mean VWC',main='Modeled Effect of Variable on VWC')
errbar(Yrs,YM,YM+Ystds,YM-Ystds,xlab='Year',ylab='VWC',main='Modeled Effect of Year on VWC, Cross Validated',errbar.col='grey')

#Look at ranges in differences between years:
hist(YearMeanVWC[3,]-YearMeanVWC[2,])
hist(YearMeanVWC[3,]-YearMeanVWC[1,])
hist(YearMeanVWC[2,]-YearMeanVWC[1,])
t.test(YearMeanVWC[3,],YearMeanVWC[1,])
#The cross validation shows that 2016 VWC was significantly different than 2014 (though 2015 and 2016 were not significantly different). 
#Should I add error bars to all these graphs, or just state differences? Error bars might get complicated with variables that have lots of possible values (unlike year, which only has 3).

#Model using full dataset
Tfit<-randomForest(VWC~Dist_From_River+LatPoint+LonPoint+Veg+Veg73+Year+DOY+Upslope.Area+slope_deg+Aspect+tpi_300m+TWI.10m+Time_Since_Fire+Fire_Num+SevNum+Elevation,data=SoilM,importance=TRUE,importanceSD=TRUE,nodesize=5,ntree=500)#save(Tfit,file='RandomTree_04_11_15_7veg_tst25.rdata')
imp<-Tfit$importance
impSE<-Tfit$importanceSD
barplot(xlab='Variable',ylab='Importance',names.arg=rownames(imp)[order(imp[, 1], decreasing=TRUE)],height=imp[order(imp[, 1], decreasing=TRUE),1],cex.names=.5)
errbar(rownames(imp),imp[,1],imp[,1]+impSE,imp[,1]-impSE,errbar.col='grey')




##map
palette(c('red','orange','green','blue'))
     nclrs=4
plot(SoilX$LonInd,SoilX$LatInd,col=round(thetaM*nclrs/10),type='p',pch=1)
  points(SoilX$LonInd,SoilX$LatInd,col=round(PredTree2*nclrs/10),type='p',pch=2)
  points(SoilX$LonInd,SoilX$LatInd,col=round(PredTree2*nclrs/10),type='p',pch=3)

palette(c('black','blue','purple','green','yellow','orange','red','magenta'))
     nclrs=8
     clrnums=abs(round((PredTree2-thetaM)*nclrs/40))
plot(jitter(SoilX$LonInd),jitter(SoilX$LatInd),col=clrnums,main='Map of Model Error',pch=21,bg=clrnums)
     legend('bottomleft',legend=c('5','10','15','20','25','30','35','40'),col=c(palette()),pch=21)

     


#Predict Soil Moisture over Whole Map
  #library(randomForest)
    # load('RandomForest_16_03_17_SubSiteMeans.rdata')
ICBmod<-Tfit

# library(tiff)   # Veg97map<-readTIFF('../GPSstuff/VegRasterClip_WHR_MaxArea1.tif',native=TRUE)
#library(adehabitat)    #image.plot

library(rgdal)
LoadSaved=1

Veg97map<-readGDAL('../GPSstuff/VegRaster_WHR_10m_MaxArea1.tif')
#image(Veg97map)
CoordVeg=coordinates(Veg97map)  #In meters.
#CoordVeg=coordinates(WGSmap)    #In degrees
CoordVeg[1:4,]


Veg97det<-readGDAL('../GPSstuff/cl_vegdetail_ProjectRaster_C1.tif') #Use this to change "Conifer Reproduction" from Shrub type
#WGSmap<-readGDAL('../GPSstuff/clip_10m_WGS') #Something wrong. Too big
VegEcogmap<-readGDAL('../GPSstuff/clip_Veg14_ICB.tif')   #veg2012_NAIP_rast10m.tif') 
Veg69map<-readGDAL('../GPSstuff/clip_veg69_ICB.tif') 

VegMorig=as.matrix(Veg97map)

 #coordinates(VegMorig)=CoordVeg #Doesnt work for matrices
     
     #image(VegMorig)
     Rlon=289923 #Rightmost limit of ICB 
     MxX=max(CoordVeg[,1])
     MnX=min(CoordVeg[,1])
     MxY=max(CoordVeg[,2])
     MnY=min(CoordVeg[,2])
     
     ICB_Ind<-CoordVeg[,1]<(.4*(MxX-MnX)+MnX)&CoordVeg[,2]>(.25*(MxY-MnY)+MnY)&CoordVeg[,2]<(.75*(MxY-MnY)+MnY)
     VegMorig=VegMorig[ICB_Ind]  #[round(.2*3484):round(.7*3484),0:round(.4*5295)]
     XYvals=CoordVeg[ICB_Ind,]
     names(VegMorig)='VegNum'
     names(XYvals)=c('x','y')
     VegMorig=as.data.frame(VegMorig)
     coordinates(VegMorig)<-XYvals

#ICB_IndR<-1:2000
#ICB_IndC<-1000:3000

     #VegMorig=VegMorig[ICB_IndR,ICB_IndC]
#spplot(VegMorig)

     VegDetOrig=as.matrix(Veg97det)    

if(LoadSaved<1){         
VegM=0*VegMorig$VegMorig
VegNom=VegM
  VegM[VegMorig$VegMorig>=98|VegMorig$VegMorig==2|VegMorig$VegMorig==20]=NaN  #Rocks and water are NaN for now
	#VegM[VegMorig$VegMorig==11]=3.1    #11=red fir, 
     #VegNom[VegMorig$VegMorig==11]='Conifer, mix fir'
	VegM[VegMorig$VegMorig==44|VegMorig$VegMorig==43|VegMorig$VegMorig==17|VegMorig$VegMorig==24]=2   #44=montane chaparall, 43=sagebrush  #17 and 24 are hardwood.  Put them here?
     VegNom[VegMorig$VegMorig==44|VegMorig$VegMorig==43|VegMorig$VegMorig==17|VegMorig$VegMorig==24]='Shrub'
     VegM[VegMorig$VegMorig==21]=5   #21=aspen
      VegNom[VegMorig$VegMorig==21]='aspen'
	VegM[VegMorig$VegMorig==12]=3.4   #12=lodgepole
     VegNom[VegMorig$VegMorig==12]='Conifer, pico'
      VegM[VegMorig$VegMorig==28]=7    #28=montane riparian
      VegNom[VegMorig$VegMorig==28]='Riparian'
	VegM[VegMorig$VegMorig==3|VegMorig$VegMorig==4]=6    #3=wet meadow
   VegNom[VegMorig$VegMorig==3|VegMorig$VegMorig==4]='Wet Meadow'
	VegM[VegMorig$VegMorig==6]=4    #6=perennial grassland (meadow or dry meadow)
   VegNom[VegMorig$VegMorig==6]='Meadow'
  VegM[VegMorig$VegMorig==45]=7    #45=willows (mostly wetland type, but may want to divide up later)
   VegNom[VegMorig$VegMorig==45]='Willow, rip.'
	VegM[VegMorig$VegMorig==11|VegMorig$VegMorig==14|VegMorig$VegMorig==16|VegMorig$VegMorig==10|VegMorig$VegMorig==13|VegMorig$VegMorig==15]=3.1  #14=white fir, #16=Jeffrey pine, #15=Doug fir (may fit better w/ pico?)
    VegNom[VegMorig$VegMorig==11|VegMorig$VegMorig==14|VegMorig$VegMorig==16|VegMorig$VegMorig==10|VegMorig$VegMorig==13|VegMorig$VegMorig==15]='Conifer, mix fir'
  
#Change conifer reproduction from shrub to "Conifer, mix con rec"
 VegDetOrig=VegDetOrig[ICB_Ind]
   VegM[VegDetOrig==7]=3.3   #7=Conifer reproduction
   VegNom[VegDetOrig==7]='Conifer, mix con rec'


#20=juniper (all growing in fairly barren areas)
		#10=subalpine conifer.  Didn't measure anything up here, but just set to driest conifer group
     
BigSoilX<-data.frame(cbind(VegNum=VegM,VegNames=VegNom)) #For normalizing
BigSoilM<-BigSoilX  #For keeping original values

veg69num<-as.matrix(Veg69map)
BigSoilM$veg69<-veg69num[ICB_Ind]
Veg14num<-as.matrix(VegEcogmap)
BigSoilM$Veg14<-Veg14num[ICB_Ind]
BigSoilM$veg69[BigSoilM$veg69>6]=NaN
BigSoilM$Veg14[BigSoilM$Veg14>6]=NaN
BigSoilM$Veg14[(BigSoilM$Veg14==1) & (VegMorig$VegMorig==12)]=1.5

DEM10m<-readGDAL('../GPSstuff/Clip_10mDEM_L')  #<-Lidar DEM  ('../GPSstuff/Clip_10mDEM')
  DEM10mM<-as.matrix(DEM10m)
  BigSoilX$Elev<-(DEM10mM[ICB_Ind]-min(SoilM[,'Elev']))/(max(SoilM[,'Elev'])-min(SoilM[,'Elev']))
  BigSoilM$Elev<-DEM10mM[ICB_Ind]
#CoordDEM<-coordinates(DEM10m)   
#image(DEM10m)

#TWImap<-readGDAL('../GPSstuff/Clip_TWI_L') # ('../GPSstuff/ClipTWI10shft') #May need to smooth this to get more sensible output
#TWIm<-as.matrix(TWImap)
#BigSoilX$TWI.10m<-(TWIm[ICB_Ind]-min(SoilM[,'TWI.10m']))/(max(SoilM[,'TWI.10m'])-min(SoilM[,'TWI.10m']))
# BigSoilM$TWI.10m<-TWIm[ICB_Ind]
# BigSoilM$TWI.10m[BigSoilM$TWI.10m<0]=0
#THERE ARE SOME ISSUES WITH THE TWI MAP. SWITCH TO DIRECT CALCULATION:
BigSoilM$TWI.10m<-log(BigSoilM$Upslope.Area/tan(BigSoilM$Slope*pi/180))
BigSoilM$TWI.10m[is.na(BigSoilM$TWI.10m)]=0 #Deal with NaN when both slope and upslope area are zero
BigSoilM$TWI.10m[BigSoilM$TWI.10m>20]=20 #Deal with infinity values
BigSoilM$TWI.10m[BigSoilM$TWI.10m<(-5)]=-5
BigSoilX$TWI.10m<-(BigSoilM$TWI.10m-min(SoilM[,'TWI.10m']))/(max(SoilM[,'TWI.10m'])-min(SoilM[,'TWI.10m']))

SlpMap<-readGDAL('../GPSstuff/Clip_Slope_L')     #Clip10mSlope')
 SlopeM<-as.matrix(SlpMap)
 BigSoilX$Slope<-(SlopeM[ICB_Ind]-min(SoilM[,'Slope']))/(max(SoilM[,'Slope'])-min(SoilM[,'Slope']))
 BigSoilM$Slope<-SlopeM[ICB_Ind]
Rdist<-readGDAL('../GPSstuff/Clip_RivDst_L')   #This version uses distance from area w/ at least 500m2 upslope area
 RdistM<-as.matrix(Rdist)
 BigSoilX$Dist.from.River<-(RdistM[ICB_Ind]-min(SoilM[,'Dist.from.River']))/(max(SoilM[,'Dist.from.River'])-min(SoilM[,'Dist.from.River']))
 BigSoilM$Dist.from.River<-RdistM[ICB_Ind]
SevNum<-readGDAL('../GPSstuff/clip_Mx_rdnbr')
SevNum<-as.matrix(SevNum)
SevNum[is.na(SevNum)]=0
 #BigSoilX$SevNum<-SevNum[ICB_Ind]/1000 #1000 is near the maximum of the bulk of severity values
 # BigSoilX$SevNum[BigSoilX$SevNum<0]=0
 # BigSoilX$SevNum[BigSoilX$SevNum>1]=1
TempSevNum<-SevNum[ICB_Ind]
BigSoilM$SevNum<-TempSevNum
BigSoilM$SevNum[TempSevNum>640]=3
BigSoilM$SevNum[TempSevNum<641]=2
BigSoilM$SevNum[TempSevNum<316]=1
BigSoilM$SevNum[TempSevNum<69]=0

 BigSoilM$RdNBR<-SevNum[ICB_Ind]
 BigSoilX$SevNum<-BigSoilM$SevNum/3
# BigSoilM$SevNum[BigSoilM$SevNum<0]=0
# BigSoilM$SevNum[BigSoilM$RdNBR>0&BigSoilM$RdNBR<69]=1 #Using categories from Miller 2007 for unch=1 to high=4
# BigSoilM$SevNum[BigSoilM$RdNBR>=69&BigSoilM$RdNBR<=315]=2
# BigSoilM$SevNum[BigSoilM$RdNBR>315&BigSoilM$RdNBR<=640]=3
# BigSoilM$SevNum[BigSoilM$RdNBR>640]=4
Aspct<-readGDAL('../GPSstuff/Clip_Aspct_L')      #Clip_Aspect10')
AspctM<-as.matrix(Aspct)
AspIndMap=.5*(1-cos(pi*(AspctM-30)/180)) 
AspIndMap[AspctM==-1]=0
#image(AspIndMap)
  BigSoilX$Aspect<-AspIndMap[ICB_Ind]
  BigSoilM$Aspect<-AspIndMap[ICB_Ind]
AreaAccum<-readGDAL('../GPSstuff/clip_flowac_L')
#Made a slight mistake, just need to remove first row and col of AreaAccum
A_upslope<-as.matrix(AreaAccum)
   A_upslope[A_upslope>5000]=5000
  BigSoilX$Upslope.Area<-A_upslope[ICB_Ind]/800
  BigSoilM$Upslope.Area<-A_upslope[ICB_Ind]
FireYr<-readGDAL('../GPSstuff/All_IllFires_19302009_Raster1.tif')
YearsSince<-as.matrix(FireYr)
YearsSince=2014-YearsSince
YearsSince[is.na(YearsSince)]=100
  BigSoilX$Time.Since.Fire<-YearsSince[ICB_Ind]/100
  BigSoilM$Time.Since.Fire<-YearsSince[ICB_Ind]
NumBurns<-readGDAL('../GPSstuff/clip_FireNum')
TimesBurned<-as.matrix(NumBurns)
TimesBurned[is.na(TimesBurned)]=0      #If no value, then it didn't burn
  BigSoilX$Times.Burned<-TimesBurned[ICB_Ind]/max(SoilM[,'Times.Burned'])
 BigSoilX$BurnInd<-BigSoilX$Times.Burned
  BigSoilX$BurnInd[BigSoilX$Times.Burned>1]=1
  BigSoilX$BurnInd[BigSoilX$Times.Burned==0]=1.25
  BigSoilX$BurnInd<-(BigSoilX$BurnInd-0.25)
 BigSoilM$Times.Burned<-TimesBurned[ICB_Ind]
 BigSoilM$BurnInd<-BigSoilX$BurnInd

TPI300<-readGDAL('../GPSstuff/clip_tpi300')
TPInum<-as.matrix(TPI300)
BigSoilM$TPI300m<-TPInum[ICB_Ind]

#Make sure to convert distances to be in degrees:
  #1 deg latitude is 111000 meters
  #1 deg longitude is 88200 meters at this latitude
BigXcoords<-XYvals
#BigXcoords[,1]=XYvals[,1]/88200
#BigXcoords[,2]=XYvals[,2]/111000

#MAY NEED TO DO THIS MORE PRECISELY!!!!!!
#BigSoilX$LonInd<-(BigXcoords[,1]-min(BigXcoords[,1]))/(max(BigXcoords[,1])-min(BigXcoords[,1])) #/(max(SoilM[,'LonPoint'])-min(SoilM[,'LonPoint']))
#BigSoilX$LatInd<-(BigXcoords[,2]-min(BigXcoords[,2]))/(max(BigXcoords[,2])-min(BigXcoords[,2])) #/(max(SoilM[,'LatPoint'])-min(SoilM[,'LatPoint']))
#BigSoilX$LonInd<-(XYvals[,1]-min(SoilM[,'POINT_X']))/(max(SoilM[,'POINT_X'])-min(SoilM[,'POINT_X']))
#BigSoilX$LatInd<-(XYvals[,2]-min(SoilM[,'POINT_Y']))/(max(SoilM[,'POINT_Y'])-min(SoilM[,'POINT_Y']))

BigSoilM$POINT_X<-XYvals[,1]
BigSoilM$POINT_Y<-XYvals[,2]



coordinates(BigSoilX)<-BigXcoords

#save(BigSoilM,file='BigSoilM_11_01_17.rdata')



#save(BigSoilX,file='BigSoilX_12_08_14.rdata') #LonInd and LatInd are normalized distances in meters. Included conifer reproduction from 97 veg map
#save(BigSoilX,file='BigSoilX_12_10_14.rdata') #Same but better distribution of fire severity
#load('BigSoilX_12_10_14.rdata') 
#load('BigSoilM97_23_03_15.rdata')

}else{
  load('Rfiles/BigSoilM_11_01_17.rdata')
  #load('BigSoilM_2_11_16.rdata')
  #load('BigSoilM_17_10_15.rdata')
  }
     
     
     

ModYear=1970
#ModYear=1997
#ModYear=2014
#ModYear=2019 #Code for all shrub
#ModYear=2020 #Code for all meadow

BigSoilM$VegNum<-as.numeric(BigSoilM$VegNum)
Veg97temp<-BigSoilM$VegNum

BigSoilM$VegChange=(10*BigSoilM$veg69+BigSoilM$Veg14)

BigSoilM$Veg14b=BigSoilM$Veg14
BigSoilM$Veg14b[BigSoilM$Veg14==4]=5
BigSoilM$Veg14b[BigSoilM$VegChange==44]=6 #New wetlands =5, persistent wetlands = 6
BigSoilM$Veg14b[BigSoilM$Veg14==5]=7
BigSoilM$Veg14b[BigSoilM$VegChange==55]=8 #New aspen =7, persistent aspen =8
BigSoilM$Veg14b[BigSoilM$Veg14==1.5]=2
BigSoilM$Veg14b[BigSoilM$Veg14==2]=3
BigSoilM$Veg14b[BigSoilM$Veg14==3]=4

if(PICOsep<1){BigSoilM$Veg14[BigSoilM$Veg14==1.5]=1}

if(ModYear<1974){
#For 1969 scenario: 

 BigSoilM$SevNum<-0*BigSoilM$Slope
 BigSoilM$Times.Burned<-0*BigSoilM$Slope
 BigSoilM$Time.Since.Fire=100
 
 BigSoilM$Veg14=BigSoilM$veg69 #Veg14 is really "current veg"
 if(PICOsep){BigSoilM$Veg14[VegMorig$VegMorig==12 & BigSoilM$veg69==1]=1.5}

 BigSoilM$DOY=160 #This should be altered for any year
 BigSoilM$TPI300m[BigSoilM$TPI300m>60]=60#This should be altered for any year
 BigSoilM$TPI300m[BigSoilM$TPI300m<(-60)]=-60
 BigSoilM$Dist.from.River[BigSoilM$Dist.from.River>500]=500#This should be altered for any year

 BigSoilM$Veg14b=BigSoilM$Veg14
 BigSoilM$Veg14b=BigSoilM$Veg14
 BigSoilM$Veg14b[BigSoilM$Veg14==4]=6
 #BigSoilM$Veg14b[BigSoilM$VegChange==44]=6 #New wetlands =5, persistent wetlands = 6
 BigSoilM$Veg14b[BigSoilM$Veg14==5]=7
 BigSoilM$Veg14b[BigSoilM$VegChange==55]=8 #New aspen =7, persistent aspen =8
 BigSoilM$Veg14b[BigSoilM$Veg14==1.5]=2
 BigSoilM$Veg14b[BigSoilM$Veg14==2]=3
 BigSoilM$Veg14b[BigSoilM$Veg14==3]=4
 
 }
if(ModYear==1997){
 Fire97<-readGDAL('../GPSstuff/clip_97fire') #Times burned before 1997
 Fire97M<-as.matrix(Fire97)
 Fire97M[is.na(Fire97M)]=0      #If no value, then it didn't burn
  BigSoilX$Times.Burned<-Fire97M[ICB_Ind]/max(SoilM[,'Times.Burned'])
 BigSoilX$BurnInd<-BigSoilX$Times.Burned
  BigSoilX$BurnInd[BigSoilX$Times.Burned>1]=1
  BigSoilX$BurnInd[BigSoilX$Times.Burned==0]=1.25
  BigSoilX$BurnInd<-(BigSoilX$BurnInd-0.25)
 BigSoilM$Times.Burned<-Fire97M[ICB_Ind]
 BigSoilM$BurnInd<-BigSoilX$BurnInd

 FireYr<-readGDAL('../GPSstuff/clip_p97fryr.tif')
 YearsSince<-as.matrix(FireYr)
	YearsSince<-YearsSince[,2:3485]
 YearsSince[is.na(YearsSince)]=1900
 YearsSince=(1997-YearsSince)
 BigSoilX$Time.Since.Fire<-YearsSince[ICB_Ind]/100
 BigSoilX$Time.Since.Fire[BigSoilX$Times.Burned==0]=0.97 
 BigSoilM$Time.Since.Fire<-YearsSince[ICB_Ind]
 BigSoilM$Time.Since.Fire[BigSoilM$Times.Burned==0]=97
 
    #FireSevPre97<-readGDAL('../GPSstuff/clip_Pre97nbr') NOT DONE
 BigSoilX$SevNum[BigSoilX$Times.Burned==0]=0  #Really should make a new severity map, since some more recent fires higher severity, but don't have nbr for older fires
  BigSoilM$SevNum[BigSoilM$Times.Burned==0]=0
  BigSoilM$SevNum[(BigSoilM$Times.Burned>0)&(BigSoilM$SevNum<1)]=1 #If burned but rdnbr<0, label as sev=1
#save(BigSoilM,file='BigSoilM97_23_03_15.rdata')
}
if((ModYear>2011)&(ModYear<2017))
{
BigSoilM$Year=ModYear
#BigSoilM$Veg14[is.na(Veg97temp)]=NaN #Block out cells outside of ICB
BigSoilM$DOY=220 #This should be altered for any year
BigSoilM$Time.Since.Fire=BigSoilM$Time.Since.Fire+(ModYear-2014)
#BigSoilM$TPI300m<-as.numeric(BigSoilM$TPI300m)
BigSoilM$TPI300m[BigSoilM$TPI300m>60]=60#This should be altered for any year
BigSoilM$TPI300m[BigSoilM$TPI300m<(-60)]=-60
BigSoilM$Dist.from.River[BigSoilM$Dist.from.River>500]=500#This should be altered for any year
}
if(ModYear==2019){
  BigSoilX$VegNames[BigSoilX$VegNum<7] = 'Shrub'
  #BigSoilX$SevNum[BigSoilX$SevNum==0]=0.5
  #BigSoilX$Time.Since.Fire[BigSoilX$Time.Since.Fire>.7]=0.5

}
if(ModYear==2020){
  BigSoilX$VegNames[BigSoilX$VegNum<7] = 'Wet Meadow'
  BigSoilX$BurnInd[BigSoilX$Times.Burned==0]=0.5
  BigSoilX$SevNum[BigSoilX$SevNum==0]=0.5
  BigSoilX$Time.Since.Fire[BigSoilX$Time.Since.Fire>.7]=0.5
}

#Calculate Indices
#BigSoilM$TimeSevInd<-(1-(BigSoilM$Time.Since.Fire/100))*(BigSoilM$SevNum/4)
# TempTPI<-(BigSoilM$TPI300m-min(SoilM$TPI300m))/(max(SoilM$TPI300m)-min(SoilM$TPI300m))
# TempSlope<-(BigSoilM$Slope-min(SoilM$Slope))/(max(SoilM$Slope)-min(SoilM$Slope))
# TempArea<-BigSoilM$Upslope.Area
#  TempArea[TempArea>600]=625
#  TempArea<-(TempArea-min(SoilM$Upslope.Area))/(625-min(SoilM$Upslope.Area))
#BigSoilM$TopInd<-.5*TempTPI+.3*TempSlope-0.2*TempArea
# BigSoilM$TopInd<-(.5*BigSoilM$TPI300m+.3*BigSoilM$Slope-0.2*BigSoilM$Upslope.Area)

#Xgood=BigSoilX[!is.na(BigSoilX$VegNum)&as.numeric(BigSoilX$VegNum)>0,]
##XgoodCoords<-BigXcoords[!is.na(BigSoilX$VegNum)&as.numeric(BigSoilX$VegNum)>0,]
#XgoodCoords<-XYvals[!is.na(BigSoilX$VegNum)&as.numeric(BigSoilX$VegNum)>0,]


AllwaysVeg<-(!is.na(BigSoilM$Veg14)&(BigSoilM$Veg14>0)&!is.na(BigSoilM$veg69)&(BigSoilM$veg69>0))

Xgood=BigSoilM[AllwaysVeg,]
#Xgood$VegNum<-as.numeric(Xgood$VegNum)
#Xgood$TopInd[Xgood$TopInd>0.6]=0.6
XgoodCoords<-XYvals[AllwaysVeg,]

ToPredict<-seq(from=1,to=nrow(Xgood),by=4)  #<-runif(n=50000,min=1,max=nrow(Xgood)) #c(199000:200000)  #(Xgood$VegNum==4) #c(5000:7000)
thetaMod<-predict(ICBmod,Xgood[ToPredict,])
  #thetaMod<-thetaMod^2

#thetaMod[thetaMod<0]=0   #Make sure nothing <0
#thetaMod[thetaMod>60]=60

     #thetaMod[VegMorig==2]=100  #2 is the code for water

#Add rocks and water bodies
#thetaMod2=rbind(thetaMod,rep(100,times=length(WaterInd),rep(0,times=RockInd))
	
palette(c('red','orange','green','blue'))
nclrs=4
plot(XgoodCoords[ToPredict,1],XgoodCoords[ToPredict,2],col=round(thetaMod*nclrs/10),type='p',pch=1)
#plot(Xgood$VegNames[ToPredict],thetaMod)

#boxplot(thetaMod~names(thetaMod))


#thetaMod=as.data.frame(thetaMod)
#coordinates(thetaMod)<-XgoodCoords[ToPredict,]
#     spplot(thetaMod)   #May need to be in a grid to make this work

#Export thetaMod along with corresponding coordinates from XYvals to create points on a GIS map.
XYpred<-XgoodCoords[ToPredict,]
LonPred=XYpred[,1]
LatPred=XYpred[,2]
Veg14Pred=Xgood$Veg14[ToPredict]
Veg69Pred=Xgood$veg69[ToPredict]
TWIPred=Xgood$TWI.10m[ToPredict]
SevPred=Xgood$SevNum[ToPredict]
SlpPred=Xgood$Slope[ToPredict]
ElevPred=Xgood$Elev[ToPredict]
#save(Veg14Pred,file='ModeledVWC/SubSiteMeans/Veg14pred_Y2014D220_clim2014_RandForest_04_26_17_NoPICO_NoGhost',ascii=TRUE)
#save(Veg69Pred,file='ModeledVWC/SubSiteMeans/Veg69pred_Y2014D220_clim2014_RandForest_04_26_17_NoPICO_NoGhost',ascii=TRUE)
#save(LonPred,file='ModeledVWC/SubSiteMeans/Lonpred_Y2014D220_clim2014_RandForest_04_26_17_NoPICO_NoGhost',ascii=TRUE)
#save(LatPred,file='ModeledVWC/SubSiteMeans/Latpred_Y2014D220_clim2014_RandForest_04_26_17_NoPICO_NoGhost',ascii=TRUE)
#save(TWIPred,file='ModeledVWC/SubSiteMeans/TWIpred_Y2014D220_clim2014_RandForest_04_26_17_NoPICO_NoGhost',ascii=TRUE)
#save(thetaMod,file='ModeledVWC/SubSiteMeans/VWCpred_Y1970D160_clim2014_RandForest_04_26_17_NoPICO_NoGhost',ascii=TRUE)



#Look at relationship between different variables

#Correlation Matrix to check for changed collinearity
CMB=cor(Xgood[,c(3:6,8:12,14,17:18)])   
.01*round(CMB*100)
sum((abs(CMB)>.7)&CMB<1)/2
sum((abs(CMB)>.4)&CMB<1)/2


SN<-SevNum[(TimesBurned>0)&~is.na(TimesBurned)]
TB<-TimesBurned[(TimesBurned>0)&~is.na(TimesBurned)]
cor(SN,TB,use="p")
plot(TB,SN,ylim=c(0,3000))
hist(SN[TB>3&SN>0])
hist(SN[TB<3&SN>0&SN<=1100])


plot(Xgood[ToPredict,'VegNum'],thetaMod)


