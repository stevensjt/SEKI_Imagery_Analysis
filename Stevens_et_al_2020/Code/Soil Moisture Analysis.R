##This is a component of analysis code for Stevens et al. 2020: 
##Forest vegetation change and its impacts on soil water following 47 years of managed wildfire. 
##Ecosystems, In Press.
##Lead Author: Gabrielle Boisrame
##Contact: Jens Stevens, stevensjt@gmail.com
##Soil moisture analyses: This code Train models using measured soil moisture
##And produces figures 7-8 in the manuscript and figures D1-D5 in the supporting material.

####0. Load libraries and initialize parameters####

library(chron) #for years(); version 2.3-54
library(randomForest) #for randomForest(); version 4.6-14
library(Hmisc) #for errbar(); version 4.2-0
library(matrixStats) #for rowSds(); version 0.55.0
library(hydroGOF) #for rmse(); version 0.3-10
library(MASS) #for truehist(); version 7.3-51.4


AggData <- 1 #For grouping measurements that are the same veg type at the same site
Extrap <- #Set to "1" if want to extrapolate to other non-measured points. 
  #To only work with observation locations, set this to zero.
  1 
Validate <- #Set to "1" if want to run validation of the random forest model (takes a long time to run). 
  #To skip this step, set to zero.
  1 

#Be sure to set the working directory to this source file's location
setwd("./Code")

####1. Load and process data####
SoilM <- read.csv('./Data/SugarloafSoilMoisture.csv',header=TRUE) #Observed Data


if(Extrap){ #Set up gridded predictor values for predicting soil moisture in unmeasured areas of the watershed 
  SoilMbig <- #gridded dataframe
    read.csv('./Data/Geospatial/RasterToPointsValues.csv',header=TRUE)  #Gridded data across full watershed.
  SoilMbig <- SoilMbig[SoilMbig$aspect>-2,] #Get rid of bad data points
  SoilMbig <- SoilMbig[SoilMbig$SevNum<5,] #Get rid of bad data points
  SoilMbig <- SoilMbig[(SoilMbig$Veg!=3) & 
                         (SoilMbig$Veg!=5) & 
                         (SoilMbig$Veg73!=3) & 
                         (SoilMbig$Veg73!=5) & 
                         SoilMbig$Veg>0,
                       ]
  
  #If the number of fires in the map is <0 ("no data" code), set the number of fires to zero.
  SoilMbig$Fire_Num[SoilMbig$Fire_Num<0] <- 0 
  
  #Change variable names to match those in the fitted model.
  SoilMbig$Veg14 <- SoilMbig$Veg
  SoilMbig$Veg[SoilMbig$Veg==1] <- "shrub"
  SoilMbig$Veg[SoilMbig$Veg==2] <- "sparse meadow"
  SoilMbig$Veg[SoilMbig$Veg==4] <- "mixed conifer"
  SoilMbig$Veg[SoilMbig$Veg==6] <- "dense meadow"
  SoilMbig$Veg <- as.factor(SoilMbig$Veg)
  SoilMbig$Veg73 <- as.factor(SoilMbig$Veg73)

  AspInd <- #Convert aspect to an index between 0 and 1
    .5*(1-cos(pi*(SoilMbig[,'aspect']-30)/180)) 
  AspInd[SoilMbig[,'aspect']<0] <- 0
  SoilMbig$AspectDeg <- SoilMbig$aspect
  SoilMbig$Aspect <- AspInd

  #Calculate time since fire for the gridded dataset 
  #For sampling points where no fires had been observed since 1970, 
    #set time since fire to 100 yrs
  SoilMbig$Time_Since_Fire <- (2016-SoilMbig$Fire_Year)
  SoilMbig$Time_Since_Fire[SoilMbig$Time_Since_Fire>100] <- 100
  
  #Calculate Topographic Wetness Index (TWI)
  #(flow accumulation area)/tan(Slope)
  SoilMbig$TWI.10m <- #Set minimum of tan(slope) to 0.0001 to avoid division by zero.
    log(SoilMbig$flow_acc_d/max(0.0001,tan(SoilMbig$slope_deg*pi/180))) 
  SoilMbig$TWI.10m[SoilMbig$TWI.10m<0] <- 0

  #Change variable names to match fitted random forest model.
  SoilMbig$Upslope.Area <- SoilMbig$flow_acc_d
  SoilMbig$Upslope.Area[SoilMbig$flow_acc_d>5000] <- 5000 #Set max flow accumulation area to 5000 
  #(to avoid being too far outside or range used to fit the model) 
  
  #Option to plot: 
  #1) Distribution of slope of mixed conifer plots
  #hist(SoilMbig$slope_deg[SoilMbig$Veg14==4],freq=FALSE,col = rgb(0,1,0,.5))
  #2) Distribution of slope of dense meadow plots; blue color indicates overlapping bars
  #hist(SoilMbig$slope_deg[SoilMbig$Veg14==6],freq=FALSE,col = rgb(0,0,1,.5),add=TRUE)
  #3) Distribution of topographic wetness index (TWI) of mixed conifer plots
  #hist(SoilMbig$TWI.10m[SoilMbig$Veg14==4],freq=FALSE,col = rgb(0,1,0,.5))
  #4) Distribution of TWI of dense meadow plots; blue color indicates overlapping bars
  #hist(SoilMbig$TWI.10m[SoilMbig$Veg14==6],freq=FALSE,col = rgb(0,0,1,.5),add=TRUE)
}

##Additional data processing

#Remove a few calibration points
SoilM <-
  SoilM[SoilM$Site!='Calibration',] 

#Remove NaNs (rocks)
thetaM <- #Make sure values are numeric so can use ">=" equations
  #thetaM is key response variable of soil moisture (= volumetric water content)
  as.numeric(SoilM$Soil_Sat)
SoilM <- SoilM[!is.na(thetaM),]
thetaM <-thetaM[!is.na(thetaM)]
SoilM <- SoilM[thetaM>=0,]
thetaM <- thetaM[thetaM>=0]

#Transform aspect to be index from 0 to 1
AspInd <- .5*(1-cos(pi*(SoilM[,'aspect']-30)/180)) 
AspInd[SoilM[,'aspect']<0] <- 0
SoilM$AspectDeg <- SoilM$aspect
SoilM$Aspect <- AspInd

#Get dates playing nice
SoilM$Date <- as.Date(SoilM$TimeStamp,format="%m/%d/%Y")
SoilM$Year <- 2015+as.numeric(years(SoilM$Date))
SoilM$Month <- months(SoilM$Date)
SoilM$DOY <- as.numeric(1+SoilM$Date-as.Date(sprintf('01/01/%i',SoilM$Year),format="%m/%d/%Y"))

#Make sure numbers are numbers
SoilM$Fire_Year <- as.integer(SoilM$Fire_Year)
SoilM$Fire_Num <- as.integer(SoilM$Fire_Num)
SoilM$Fire_Num[SoilM$Fire_Num<0] <- 0

#Ensure that R is reading values as numbers rather than text
SoilM$Veg73 <- as.numeric(SoilM$X1973_veg)
SoilM$Veg14 <- as.numeric(SoilM$X2014_veg)

#Calculate time since fire for the gridded dataset 
#For sampling points where no fires had been observed since 1970, 
  #set time since fire to 100 yrs
SoilM$Time_Since_Fire <- (SoilM$Year-SoilM$Fire_Year)
SoilM$Time_Since_Fire[SoilM$Time_Since_Fire>100] <- 100 

if (AggData){ 
  #This function takes the full dataset, which has multiple sample points per site 
  #it then calculates both the mean and maximum values for soil moisture (volumetric water content) 
  #and other variables of interest within each dominant vegetation type of each site.
  source("./Code/Functions/AggSubSites.R")
  SoilM <- AggSubSites(SoilM)
  thetaM <- SoilM$Soil_Sat
}


#Get mean VWC for each site, into three new data frames
if(AggData==0){
  SiteMeansA <- aggregate(thetaM,list(SoilM$Site,SoilM$Veg),mean)
  SiteMeans <- aggregate(SiteMeansA$x,list(SiteMeansA$Group.1),max)
  SiteCounts <- aggregate(thetaM,list(SoilM$Site),length)
  #Divide sites into groups
}else{
  SiteMeansA <- SoilM[,c('SubSite','Soil_Sat')]
  SiteMeans <- aggregate(SiteMeansA$Soil_Sat,list(round(SiteMeansA$SubSite)),max)
  SiteCounts <- aggregate(thetaM,list(SoilM$Site),length)
}


#Replace no-data values for flow accumulation and slope with zeros
SoilM$flow_acc_d[SoilM$flow_acc_d<0] <- 0 #Replace noData values with zeros
SoilM$slope_deg[SoilM$slope_deg<0] <- 0 #Replace noData values with zeros

SoilM$TWI.10m <- #Calculate topographic wetness index
  log(SoilM$flow_acc_d/max(0.0001,tan(SoilM$slope_deg*pi/180)))
SoilM$TWI.10m[SoilM$TWI.10m<0] <- 0


#Change names of fire severity variables to match model
SoilM$SevNum <- SoilM$Fire_Max_Sev
#Set fire to minimum severity for fires where don't have severity data. 
SoilM$SevNum[(SoilM$Fire_Max_Sev<0)&SoilM$Fire_Num>0] <- 1 
SoilM$SevNum <- as.factor(SoilM$SevNum)

SoilM$VWC <- thetaM #VWC = volumetric water content

#Set certain values as factors, since order doesn't matter.
SoilM$Veg14<-as.factor(SoilM$Veg14)
SoilM$Veg73<-as.factor(SoilM$Veg73)
SoilM$Year<-as.factor(SoilM$Year)

####2. Random Forest modeling####

#Load random forest model trained to Illilouette Creek Basin data
load('./Code/Functions/ICB_RandomForestModel.rdata')
Tfit_ICB<-Tfit #Save the Tfit model just loaded under a separate name.

#Change variable names to match model.
SoilM$Upslope.Area<-SoilM$flow_acc_d
SoilM$Upslope.Area[SoilM$flow_acc_d>5000]=5000 #Set max flow accumulation area to 5000 
#(to avoid being too far outside or range used to fit the model)

#Fit random forest model to Sugarloaf data
Tfit <- randomForest(VWC ~ 
                       Veg + Veg73 + Year + DOY + Upslope.Area + 
                       slope_deg + Aspect + tpi_300m + TWI.10m + 
                       Time_Since_Fire + Fire_Num + SevNum + Elevation + 
                       Dist_From_River,
                     data=SoilM, nodesize=5, ntree=500)

####3. Review output from Random Forest modeling####

#Create importance plot showing which variables affect the model results the most (Figure D1). 
#This may appear slightly different from the manuscript due to stochasticity in the model training.
imp <- Tfit$importance
barplot(xlab='Variable', ylab='Importance', names.arg=rownames(imp)[order(imp[, 1], decreasing=TRUE)], 
        height=imp[order(imp[, 1], decreasing=TRUE),1], cex.names=.5)

#Create plots showing how individual variables (in this case, year) 
  #affect the soil moisture independently of other variables (used to create barplot panels in Figure D3). 
a <- partialPlot(x=Tfit,pred.data=SoilM,x.var='Year',ylab='VWC')
barplot(a$y,names=a$x,xlab='',ylab='VWC',ylim=c(0,20))

#Create plots showing how individual variables (in this case, upslope area) 
  #affect the soil moisture independently of other variables (used to create components of Figures D2 and D3).
a <- #Variables can include tpi_300m, Elevation, Aspect, TWI.10m, slope_deg, Dist_From_River, Upslope.Area
  partialPlot(x=Tfit, pred.data=SoilM, x.var='Upslope.Area', ylab='VWC', ylim=c(0,20)) 
points(a$x,a$y)

####4. Extrapolate soil moisture model across the entire watershed####

#Model soil moisture across the entire watershed, 
  #both under current conditions and under an unburned scenario with 1973 vegetation
if(Extrap){
  SoilMbig$Year <- factor(2018,levels=levels(SoilM$Year)) #Set the "year" in the model to 2018
  SoilMbig$DOY <- 153 #205 # Set the date (in terms of days since December 31) you wish to model.
  SoilMbig$SevNum <- as.factor(SoilMbig$SevNum)
  
  PredTreeBig <- #Estimate soil moisture at the points given in SoilMbig using the random forest model
    predict(Tfit,SoilMbig,predict.all=TRUE)
  AllTrees <- PredTreeBig$individual
  StdTrees1 <- #Get the standard deviation of the predicted values over all trees
    apply(AllTrees,1,sd)
  Q25Trees1 <- #Get the 25th percentile of the predicted values over all trees 
    apply(AllTrees,1,quantile,.25,na.rm=TRUE)
  Q75Trees1 <- #Get the 75th percentile of the predicted values over all trees  
    apply(AllTrees,1,quantile,.75,na.rm=TRUE) 
  PredToday <- PredTreeBig$aggregate
  
  #Set up a new matrix of data that will represent unburned conditions with 1973 vegetation.
  SoilMold <- SoilMbig
  SoilMold$Veg14 <- SoilMold$Veg73
  SoilMold$Time_Since_Fire <- 100
  SoilMold$SevNum <- factor(0,levels=levels(SoilM$SevNum))
  SoilMold$Fire_Num <- 0
  
  PredTreeOld <- #Estimate soil moisture at the points given in SoilMbig using the random forest model
    predict(Tfit,SoilMold,predict.all=TRUE) 
  AllTrees <- PredTreeOld$individual
  StdTrees2 <- #Get the standard deviation of the predicted values over all trees
    apply(AllTrees,1,sd) 
  Q25Trees2 <- #Get the 25th percentil of the predicted values over all trees 
    apply(AllTrees,1,quantile,.25,na.rm=TRUE) 
  Q75Trees2 <- #Get the 75th percentile of the predicted values over all trees
    apply(AllTrees,1,quantile,.75,na.rm=TRUE) 
  PredUnburned <- PredTreeOld$aggregate
  
  #Plot modeled soil moisture under current conditions, 
    #versus what soil mositure would be under unburned conditions.
    #This is use to create Figure 8 in the manuscript.
  #The next two plots may appear slightly different from Figure 8 in the manuscript due to stochasticity in the model training.
  Fchange <- (SoilMbig$Veg73!=as.factor(SoilMbig$Veg14))
  plot(PredUnburned[Fchange], PredToday[Fchange], xlab="Modeled Unburned Soil Moisture (%)", 
       ylab="Modeled Actual Soil Moisture (%)", main="June 2018", xlim=c(3,45), ylim=c(3,45))
  points(PredUnburned[SoilMbig$Veg73=='4'&SoilMbig$Veg14==6], 
         PredToday[SoilMbig$Veg73=='4'&SoilMbig$Veg14==6],
         pch=15, col='blue') #Conifer-Wet transition
  points(PredUnburned[SoilMbig$Veg73=='4'&SoilMbig$Veg14==2], 
         PredToday[SoilMbig$Veg73=='4'&SoilMbig$Veg14==2],
         pch=16, col='grey') #Conifer-Sparse transition
  points(PredUnburned[SoilMbig$Veg73=='6'&SoilMbig$Veg14==4],
         PredToday[SoilMbig$Veg73=='6'&SoilMbig$Veg14==4],
         pch=17, col=rgb(.2,.75,.2)) #Wet-Conifer transition
  points(PredUnburned[SoilMbig$Veg73=='4'&SoilMbig$Veg14==1],
         PredToday[SoilMbig$Veg73=='4'&SoilMbig$Veg14==1],
         pch=18, col='brown') #Conifer-Shrub transition
  lines(c(3,45),c(3,45),lty=2,lwd=2)
  legend(3,43,c("Conifer - Dense Mdw.", "Conifer - Sparse Mdw.", "Conifer-Shrub", 
                "Dense Mdw. - Conifer", "Other"),
         pch=c(15,16,18,17,1), col=c("blue","grey","brown", rgb(.2,.75,.2), "black"))
  
  #Create histogram of change, included in Figure 8
  hist(PredToday[Fchange] - PredUnburned[Fchange], 
       breaks=seq(-4.25,4.25,.5), main="June 2018", xlab="Actual-Unburned Soil Moisture (%)") 
  hist(100*(PredToday[Fchange]-PredUnburned[Fchange])/PredUnburned[Fchange], 
       main="July 2016",xlab="% Change in Actual-Unburned Soil Moisture")
  print("Mean difference in soil moisture under observed fire minus simulated no fire scenario:")
  print(paste(round(mean(PredToday[Fchange]-PredUnburned[Fchange]),2), "%")) 
}

####5. Model validation and plotting soil moisture across SCB subsites####

#Create model validation plot, showing predicted versus measured soil moisture values.
PredTree <- predict(Tfit,SoilM,predict.all=TRUE)
AllTrees <- PredTree$individual
StdTrees <- apply(AllTrees,1,sd) #Get the standard deviation of the predicted values over all trees
PredTree <- PredTree$aggregate
errbar(.01*SoilM$VWC, .01*PredTree, .01*PredTree+.01*StdTrees, .01*PredTree-.01*StdTrees, 
       xlab='Measured', ylab='Predicted', main='Random Forest Model Fit to Test Data', errbar.col='grey')
points(.01*SoilM$VWC,.01*PredTree)
lines(c(0,.65),c(0,.65), col='red')


#Separate by trip and veg, for creating Figure 7
TripVegMat <- matrix(nrow=5,ncol=6)
TripVegMat[,1] <- c(2016,2016,2017,2017,2018)
TripVegMat[,2] <- c(1,2,1,2,1)
TripVegMatSD <- TripVegMat
TripVegMatMx <- TripVegMat
TripVegMatMn <- TripVegMat
TripVegMatObs <- TripVegMat


if(AggData){
  SoilMpred<-SoilM
  } else { 
    SoilMpred<-SoilM[,5:53]
  }


SoilMpred <- #Calculate across all sites that were measured in June 2016
  SoilMpred[SoilMpred$Trip=='2016_Early' | 
              (SoilMpred$Site=='GBW7') | 
              ((SoilMpred$Trip=='2017_Early') & 
                 ((SoilMpred$Site=='HighSevDryMET') | 
                    (SoilMpred$Site=='HighSevWetMET') | 
                    (SoilMpred$Site=='LowSevMET')
                  )
               ),
            ] 

TempFrYrs <- (SoilMpred$Time_Since_Fire) 
Yrs <- unique(SoilM$Year)

nsites <- length(SoilMpred$VWC)
AllModMat <- #Create matrix to hold model output
  data.frame(Year=rep(0,5*nsites), DOY=rep(0,5*nsites), VWC=rep(0,5*nsites), Veg=rep(0,5*nsites))

#Loop through the 5 site visits and calculate VWC at each site for each date
for (i in 1:5){
  Yr<-Yrs[ceiling(i/2)] #Set the yera
  Trp<-TripVegMat[i,2]
  SoilMpred$Year<-Yr
  if(Trp==1){ #Set DOY to 178 (June) or 205 (July)
    Dt<-178   
    } else {
      Dt<-205
      }
  SoilMpred$DOY <- Dt
  SoilMpred$Time_Since_Fire <- (TempFrYrs+ceiling(-1+i/2))
  SoilMpred$Time_Since_Fire[SoilMpred$Time_Since_Fire>90] <- 100
  
  a <- predict(Tfit,SoilMpred)
  
  #Calculate means within each vegetation class
  
  TripVegMat[i,3]<-.01*mean(a[SoilMpred$Veg=='dense meadow'])
  TripVegMatSD[i,3]<-.01*sd(a[SoilMpred$Veg=='dense meadow'])
  TripVegMatMx[i,3]<-.01*max(a[SoilMpred$Veg=='dense meadow'])
  TripVegMatMn[i,3]<-.01*min(a[SoilMpred$Veg=='dense meadow'])
  
  TripVegMat[i,4]<-.01*mean(a[SoilMpred$Veg=='mixed conifer'])
  TripVegMatSD[i,4]<-.01*sd(a[SoilMpred$Veg=='mixed conifer'])
  TripVegMatMx[i,4]<-.01*max(a[SoilMpred$Veg=='mixed conifer'])
  TripVegMatMn[i,4]<-.01*min(a[SoilMpred$Veg=='mixed conifer'])
  
  TripVegMat[i,5]<-.01*mean(a[SoilMpred$Veg=='shrub'])
  TripVegMatSD[i,5]<-.01*sd(a[SoilMpred$Veg=='shrub'])
  TripVegMatMx[i,5]<-.01*max(a[SoilMpred$Veg=='shrub'])
  TripVegMatMn[i,5]<-.01*min(a[SoilMpred$Veg=='shrub'])
  
  TripVegMat[i,6]<-.01*mean(a[SoilMpred$Veg=='sparse meadow'])
  TripVegMatSD[i,6]<-.01*sd(a[SoilMpred$Veg=='sparse meadow'])
  TripVegMatMx[i,6]<-.01*max(a[SoilMpred$Veg=='sparse meadow'])
  TripVegMatMn[i,6]<-.01*min(a[SoilMpred$Veg=='sparse meadow'])
  
  if(Trp==1){
    
    TripVegMatObs[i,3]<-mean(SoilM$VWC[SoilM$Veg=='dense meadow' & SoilM$Year==Yr & SoilM$DOY<190])
    TripVegMatObs[i,4]<-mean(SoilM$VWC[SoilM$Veg=='mixed conifer' & SoilM$Year==Yr & SoilM$DOY<190])
    TripVegMatObs[i,5]<-mean(SoilM$VWC[SoilM$Veg=='shrub' & SoilM$Year==Yr & SoilM$DOY<190])
    TripVegMatObs[i,6]<-mean(SoilM$VWC[SoilM$Veg=='sparse meadow' & SoilM$Year==Yr & SoilM$DOY<190])
    }
  else{
    TripVegMatObs[i,3]<-mean(SoilM$VWC[SoilM$Veg=='dense meadow' & SoilM$Year==Yr & SoilM$DOY>190]) 
    TripVegMatObs[i,4]<-mean(SoilM$VWC[SoilM$Veg=='mixed conifer' & SoilM$Year==Yr & SoilM$DOY>190])
    TripVegMatObs[i,5]<-mean(SoilM$VWC[SoilM$Veg=='shrub' & SoilM$Year==Yr & SoilM$DOY>190])
    TripVegMatObs[i,6]<-mean(SoilM$VWC[SoilM$Veg=='sparse meadow' & SoilM$Year==Yr & SoilM$DOY>190])
    }
  
  AllModMat$Year[(1+nsites*(i-1)):(nsites*i)] <- Yr
  AllModMat$DOY[(1+nsites*(i-1)):(nsites*i)] <- Dt
  AllModMat$VWCmod[(1+nsites*(i-1)):(nsites*i)] <- a
  AllModMat$Veg[(1+nsites*(i-1)):(nsites*i)] <-SoilMpred$Veg
} #End For loop

AllModMat$Trip <- (AllModMat$Year+2015+AllModMat$DOY/365)

#Figure 7 #
pdf("./Figures/Fig7.pdf", width = 6.57, height = 3.92)
par(mar=c(2.6,4,1,1))
boxplot(VWCmod~Trip+Veg, data=AllModMat, ann=FALSE, col=gray.colors(5,start=.3, end=.95), 
        at=c(1:5,7:11,13:17,19:23), xaxt="n",
        names=c("","","Dense Meadow","","","","","Conifer","","","","",
                "Shrub","","","","","Sparse","","")
        )
mtext("Volumetric Water Content (%)", side=2, line=2.5, cex=1.4)
mtext(at=c(3,9,15,21),c("Dense Meadow","Conifer","Shrub","Sparse"), 
      side=1, line=0.6, cex=1.2)
points(x=c(1:5,7:11,13:17,19:23),y=100*TripVegMat[,3:6],pch=15)
legend(18,54,c("June 2016","July 2016","June 2017","July 2017","June 2018"),
       fill=gray.colors(5,start=.3, end=.95))
dev.off()
#dev.copy2pdf(file="./Figures/tmp.pdf") #other option to print.

####6. Compare models from ICB and Sugarloaf ####
#This creates Figures D4 and D5

SoilM_match <- SoilM

#Change variable names and codes to match other model.
SoilM_match$veg12 <- SoilM$X2014_veg
 SoilM_match$veg12[SoilM$Veg=='dense meadow'] = 4
 SoilM_match$veg12[SoilM$Veg=='sparse meadow'] = 3
 SoilM_match$veg12[SoilM$Veg=='mixed conifer'] = 1
 SoilM_match$veg12[SoilM$Veg=='shrub'] = 2
SoilM_match$veg69 <- SoilM$X1973_veg
 SoilM_match$veg69[SoilM$X1973_veg==6] = 4
 SoilM_match$veg69[SoilM$X1973_veg==1] = 2
 SoilM_match$veg69[SoilM$X1973_veg==4] = 1
 SoilM_match$veg69[SoilM$X1973_veg==2] = 3
SoilM_match$Slope <- SoilM$slope_deg
SoilM_match$TPI300m <- SoilM$tpi_300m
SoilM_match$Dist.from.River <- SoilM_match$Dist_From_River
SoilM_match$Time.Since.Fire <- SoilM$Time_Since_Fire
SoilM_match$Times.Burned <- 
  factor(SoilM_match$Fire_Num,levels=Tfit_ICB$forest$xlevels$Times.Burned)
SoilM_match$Elev <- SoilM$Elevation
SoilM_match <- SoilM_match[SoilM_match$veg12>0,]
SoilM_match$veg12 <- factor(SoilM_match$veg12,levels=c(1:5))
SoilM_match$veg69 <- factor(SoilM_match$veg69,levels=c(1:5))
SoilM_match$Year <- factor(SoilM_match$Year,levels=c(2014:2018))
SoilM_match$SevNum <- as.numeric(SoilM_match$SevNum)
SoilM_match$SevNum <- (SoilM_match$SevNum-2)
SoilM_match$SevNum[SoilM_match$SevNum<0] <- 0
SoilM_match$SevNum <- as.factor(SoilM_match$SevNum)

#Create Figure D4
VWCpred_ICB <- predict(Tfit_ICB,SoilM_match)
VWCpred_SL <- predict(Tfit,SoilM)

pdf("./Figures/FigD4.pdf", width = 6.8, height = 6.0)
par(mar=c(4.4,4.1,1,1))
plot(SoilM$VWC,VWCpred_SL,xlab='Measured SCB Moisture',ylab='Modeled',
     cex.lab = 1.1)
lines(c(0,55),c(0,55),col='grey')
points(SoilM_match$VWC,VWCpred_ICB,col='red')
legend('bottomright',c('SCB Model','ICB Model'),col=c('black','red'),pch=1)
dev.off()

cor(SoilM$VWC,VWCpred_SL)
cor(SoilM_match$VWC,VWCpred_ICB)

#Create Figure D5  
pdf("./Figures/FigD5.pdf", width = 5, height = 4)
hist(VWCpred_SL-SoilM$VWC, col=rgb(0,0,0,.5), main='Model Error', 
     xlim=c(-40,40), xlab='Modeled-Measured Volumetric Water Content',
     cex.lab = 0.8, cex.axis = 0.8, cex.main = 1)
hist(VWCpred_ICB-SoilM_match$VWC, col=rgb(1,0,0,.5), add=TRUE, breaks=c(-8:8)*5)  
dev.off()


#Validate the random forest model by training on subsets of data then testing on remaining data
#Takes ~30 seconds
if(Validate){
  

  #Select random sites for training
  ns <- nrow(SiteMeans)
  pctX <- 0.70 #Proportion of sites to use
  NumSamps <- 100 #Number of times to run with random training set
  CorTrain <- 0 #Hold correlation for each training set
  CorTest <- 0 #Hold correlation for each test set
  RMSEtrain <- 0
  RMSEtest <- 0
  Top3 <- matrix(nrow=3,ncol=NumSamps)
  YearMeanVWC <- matrix(nrow=3,ncol=NumSamps) #Hold mean VWC for each year for each training set
  
  SoilM$SiteNum <- #Each subsite is in the format X.Y, 
    #where X is the site number and Y is the subsite identifier, 
    #so rounding down gives the site number.
    floor(SoilM$SubSite) 
  
  for(g in 1:NumSamps){ 
    #Train the model numsamps times on different subsets of the data
    nt <- round(pctX*ns)
    TrnS <- sample(1:ns,nt,replace=FALSE) #Select X% of the sites
    Trn <- which(SoilM$SiteNum %in% TrnS)   #Find which measurements are within those sites
    
    Trn <- sort(Trn)
    Tstx <- SoilM[-Trn,]
    Tsty <- SoilM$VWC[-Trn]
    Trnx <- SoilM[Trn,]
    Trny <- SoilM$VWC[Trn]
    
    length(Trny)/length(thetaM) #Proportion of msmts used in training 
    
    #Fit random forest model to subset of observations.
    Tfit <- randomForest(VWC ~
                           Veg + Veg73 + Year + DOY + Upslope.Area + 
                           slope_deg + Aspect + tpi_300m + TWI.10m + 
                           Time_Since_Fire + Fire_Num + SevNum + Elevation + 
                           Dist_From_River,
                         data=SoilM, subset=Trn, nodesize=5, ntree=500)
    
    PredTree2 <- predict(Tfit,Tstx,predict.all=TRUE)
    AllTrees <- PredTree2$individual
    StdTrees <- apply(AllTrees,1,sd) #Get the standard deviation of the predicted values over all trees
    Q25Trees <- #Get the 25th percentile of the predicted values over all trees 
      apply(AllTrees,1,quantile,.25) 
    Q75Trees <- #Get the 75th percentile of the predicted values over all trees 
      apply(AllTrees,1,quantile,.75) 
    
    PredTree2 <- PredTree2$aggregate
    MeanTrees <- apply(AllTrees,1,mean) #Just a test: make sure get same as PredTree2$aggregate
    max(abs(MeanTrees-PredTree2))
    
    PredThetaT <- predict(Tfit,Trnx,predict.all=TRUE)
    AllTtrees <- PredThetaT$individual
    Q25Ttrees <- #Get the 25th percentile of the predicted values over all trees
      apply(AllTtrees,1,quantile,.25) 
    Q75Ttrees <- #Get the 75th percentile of the predicted values over all trees
      apply(AllTtrees,1,quantile,.75) 
    PredThetaT <- PredThetaT$aggregate
    
    #Save the correlation coefficient of modeled vs. observed values for this 
      #round of training (CorTrain) and test (CorTest) datasets
    CorTrain[g] = cor(PredThetaT,Trny)
    CorTest[g] = cor(Tsty,PredTree2)
    
    #Save the root mean square error of modeled vs. observed values for this 
      #round of training (CorTrain) and test (CorTest) datasets
    RMSEtrain[g] = rmse(PredThetaT,Trny)
    RMSEtest[g] = rmse(PredTree2,Tsty)
    
    #Record which variables were selected as the top 3 most important predictors in this 
      #round of training
    imp <- Tfit$importance
    Tlist <- rownames(imp)[order(imp[, 1], decreasing=TRUE)]
    Top3[,g] <- Tlist[1:3]
    
    #Keep track of how year affects soil moisture in each training run
    a <- partialPlot(x=Tfit,pred.data=SoilM,x.var='Year',ylab='VWC', plot = FALSE)  
    YearMeanVWC[,g] <- .01*a$y
    
  }
  
  truehist(CorTrain,xlim=c(0,1))
  truehist(CorTest,xlim=c(0,1))
  truehist(.01*RMSEtrain,xlim=c(0,.3))
  truehist(.01*RMSEtest,xlim=c(0,.3))
  Top3  #Show top 3 variables in terms of importance for each run
  
  YM <- rowMeans(YearMeanVWC)
  Yrs <- sort(unique(SoilM$Year))
  Ystds <- rowSds(YearMeanVWC)
  barplot(YM,names=Yrs,ylim=c(0,.2), xlab='Year', ylab='Mean VWC', 
          main='Modeled Effect of Variable on VWC')
  errbar(Yrs, YM, YM+Ystds, YM-Ystds, 
         xlab='Year', ylab='VWC', main='Modeled Effect of Year on VWC, Cross Validated', 
         errbar.col='grey')
}
