AggSubSites <- function(SoilM){
  AggList=list(SoilM$Site,SoilM$SubSite,SoilM$Trip,SoilM$Veg)
  #AggList=list(SoilM$Site,SoilM$Year,SoilM$DOY,SoilM$veg12,SoilM$veg69)
  SoilMa=aggregate(SoilM[,c(12:14,31:35,40:41)],by=AggList,FUN="mean")
  SoilMa2=aggregate(SoilM[,c(29:30,36:38,43:48)],by=AggList,FUN="median") #Get values that must be integers
  CountIn=aggregate(SoilM[,'Soil_Sat'],by=AggList,length)
  SoilMa$Count=CountIn$x
  SoilMa$Site=SoilMa$Group.1
  SoilMa$Veg=SoilMa$Group.4
  SoilMa$SubSite=SoilMa$Group.2
  #SoilMa$Times.Burned=SoilMa2$Times.Burned
  #SoilMa$Time.Since.Fire=SoilMa2$Time.Since.Fire
  
  
  SoilMa=cbind(SoilMa,SoilMa2)
  SoilMVar=aggregate(SoilM[,c(12:14,31:38)],by=AggList,FUN="var")
  #SoilMse=sqrt(SoilMVar$VWC)/sqrt(CountIn$x)     #Standard Error
  #plot(SoilMa$VWC,SoilMse)
  
  #SiteNamesCheck=aggregate(SoilM[,c("Site","SubSite","Veg")],by=AggList,FUN="last")
  
  return(SoilMa)
}

