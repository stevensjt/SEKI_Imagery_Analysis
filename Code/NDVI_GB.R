library(readxl)
library(tidyverse)
library(raster)
library(rgdal)

##Process Spatial Data
#1. Read
r_dir <- ("../Processed Data/GIS/NDVI/1984/")
v_dir <- ("../Processed Data/GIS/NDVI")
r_mask <- readOGR(v_dir,layer = "SCB_ICB_WGS")
rs <- stack()
for(i in 1:length(list.files(r_dir))){
  if(i == 1){
    rs <- stack(
      raster(paste0(r_dir,list.files(r_dir)[i]),band = 1)
    )
  } else {
    rs[[i]] <- raster(paste0(r_dir,list.files(r_dir)[i]),band = 1)
  }
  #Values=255 (nodata) and <15 (likely granite) set to NA
  rs[[i]] [!between(getValues(rs[[i]]),15,100)]<-NA
}

#2. Create layer of NA values 
#'sum' produces NA for any cell that has NA at any date
r_NA <- overlay(rs, fun = sum)
rs[is.na(r_NA)]<-NA

#3. Calculate statistics within two watersheds:
# Extract raster values to list object
r.vals <- extract(rs, r_mask)

# Use list apply to calculate mean for each polygon
#r.mean <- lapply(r.vals, FUN=mean) #deprecate, can't handle NAs

#Heavy handed approach:
ICB_vals<-as.data.frame(r.vals[[1]])
length(which(is.na(ICB_vals[,1])))/nrow(ICB_vals) #proportion of dataset that is NA
ICB_vals<-ICB_vals[-which(is.na(ICB_vals[,1])),]
ICB_means <- colMeans(ICB_vals)

SCB_vals<-as.data.frame(r.vals[[2]])
length(which(is.na(SCB_vals[,1])))/nrow(SCB_vals) #proportion of dataset that is NA
SCB_vals<-SCB_vals[-which(is.na(SCB_vals[,1])),]
SCB_means <- colMeans(SCB_vals)

##Visualize
#d <- read_excel("./Processed Data/WatershedProductivity.xlsx", sheet = "ggplot") #original 1985 data
d<-data.frame(Date=rep(as.Date(c(
  "1985-05-09 UTC", "1985-05-25 UTC", "1985-06-10 UTC", 
  "1985-06-26 UTC", "1985-07-28 UTC")),2),
  Watershed = c(rep("ICB",5),rep("SCB",5))
)
d$NDVI <- c(ICB_means/100,SCB_means/100) #Put back in normal NDVI units.

ggplot(d) +
  geom_line(aes(x = Date, y = NDVI, col = Watershed))+
  labs(title = "1984")

hist(ICB_vals/100,main='',xlab='NDVI',freq=FALSE,col=rgb(1,.1,.1,.5),ylim=c(0,4))
hist(SCB_vals/100,main='',xlab='NDVI',freq=FALSE,col=rgb(0,.2,1,.3),add=TRUE)
legend(.60,4,legend=c("ICB","SCB"),col=c(rgb(1,.1,.1,.5),rgb(0,.2,1,.3)),pch=15)