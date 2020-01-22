##This is a component of analysis code for Stevens et al. 2020: 
##Forest vegetation change and its impacts on soil water following 47 years of managed wildfire. 
##Ecosystems, In Press.
##Lead Author: Gabrielle Boisrame
##Contact: Jens Stevens, stevensjt@gmail.com
##AggSubSites: This function takes the full dataset, which has multiple sample points per site 
#it then calculates either
#the mean (for continuous variables such as for soil moisture; Soil_Sat = volumetric water content) and
#or the maximum (for discrete variables) value of a suite of
#variables of interest within a subsite (area of fairly uniform vegetation type within a site)
#In this dataframe, "subsites" are essentially equivalent to "sites" except for a few sites where sampling
#Was conducted linearly across ecotones to sample different vegetation types.

AggSubSites <- function(SoilM){
  AggList <- #Consolidate relevant information into a list
    list(SoilM$Site,SoilM$SubSite,SoilM$Trip,SoilM$Veg)

  SoilMa <- #Aggregate continuous variables by mean function
    aggregate(SoilM[,c("LonPoint", "LatPoint", "Soil_Sat", "aspect", 
                       "tpi_300m", "flow_acc_d", "slope_deg", "Elevation", 
                       "Dist_From_River", "AspectDeg", "Aspect")],
              by = AggList, FUN = "mean")
  SoilMa2 <- #Organize categorical variables by max function
    aggregate(SoilM[,c("X1973_veg", "X2014_veg", "Fire_Year", "Fire_Num", 
                       "Fire_Max_Sev", "Fire_Max_Sev_Obs", "Year", "DOY", 
                       "Veg73", "Veg14", "Time_Since_Fire")],
              by = AggList, FUN = "max")
     
  SoilMa2$Trip <- SoilMa2$Group.3

  CountIn <- #Get the number of measurements included in each aggregated data point
    aggregate(SoilM[,'Soil_Sat'],by=AggList,length)
  SoilMa$Count <- CountIn$x
  SoilMa$Site <- SoilMa$Group.1
  SoilMa$Veg <- SoilMa$Group.4
  SoilMa$SubSite <- SoilMa$Group.2

  SoilMa <- #Remove first four columns which are now duplicates
    cbind(SoilMa[,5:ncol(SoilMa)],SoilMa2[,5:ncol(SoilMa2)])
  
  return(SoilMa)
}