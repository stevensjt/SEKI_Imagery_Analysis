AggSubSites <- function(SoilM){
  AggList=list(SoilM$Site,SoilM$SubSite,SoilM$Trip,SoilM$Veg)
  #AggList=list(SoilM$Site,SoilM$Year,SoilM$DOY,SoilM$veg12,SoilM$veg69)

  SoilMa=aggregate(SoilM[,c(12:14,30:34,40:42)],by=AggList,FUN="mean")
  SoilMa2=aggregate(SoilM[,c(28:29,35:38,44,46:49)],by=AggList,FUN="median") #Get values that must be integers
  SoilMa2$Trip<-SoilMa2$Group.3

  CountIn=aggregate(SoilM[,'Soil_Sat'],by=AggList,length)
  SoilMa$Count=CountIn$x
  SoilMa$Site=SoilMa$Group.1
  SoilMa$Veg=SoilMa$Group.4
  SoilMa$SubSite=SoilMa$Group.2
  #SoilMa$Times.Burned=SoilMa2$Times.Burned
  #SoilMa$Time.Since.Fire=SoilMa2$Time.Since.Fire
  
  
  SoilMa=cbind(SoilMa,SoilMa2)
  SoilMVar=aggregate(SoilM[,c(12:14,30:34,40:42)],by=AggList,FUN="var")
  #SoilMse=sqrt(SoilMVar$VWC)/sqrt(CountIn$x)     #Standard Error
  #plot(SoilMa$VWC,SoilMse)
  
  #SiteNamesCheck=aggregate(SoilM[,c("Site","SubSite","Veg")],by=AggList,FUN="last")
  
  return(SoilMa)
}

