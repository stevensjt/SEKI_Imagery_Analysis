library(readxl)
library(tidyverse)
library(raster)
library(rgdal)

##Process Spatial Data
#1. Read
r_dir <- ("./Processed Data/GIS/NDVI/1984/")
v_dir <- ("./Processed Data/GIS/NDVI/")
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
#r.vals_1985 <- r.vals #if working with 85 data
#r.vals_1984 <- r.vals #if working with 84 data

# Use list apply to calculate mean for each polygon
#r.mean <- lapply(r.vals, FUN=mean) #deprecate, can't handle NAs

#Heavy handed approach:
ICB_vals<-as.data.frame(r.vals_1985[[1]])
length(which(is.na(ICB_vals[,1])))/nrow(ICB_vals) #proportion of dataset that is NA
ICB_vals<-ICB_vals[-which(is.na(ICB_vals[,1])),]
ICB_means <- colMeans(ICB_vals)

SCB_vals<-as.data.frame(r.vals_1985[[2]])
length(which(is.na(SCB_vals[,1])))/nrow(SCB_vals) #proportion of dataset that is NA
SCB_vals<-SCB_vals[-which(is.na(SCB_vals[,1])),]
SCB_means <- colMeans(SCB_vals)

##Visualize
#d <- read_excel("./Processed Data/WatershedProductivity.xlsx", sheet = "ggplot") #original 1985 data

#new 1984 data:
d.1984<-data.frame(Date=rep(as.Date(c(
  "1984-05-09 UTC", "1984-05-25 UTC", "1984-06-10 UTC", 
  "1984-06-26 UTC", "1984-07-28 UTC")),2),
  Watershed = c(rep("ICB",5),rep("SCB",5))
)
d.1984$NDVI <- c(ICB_means/100,SCB_means/100) #Put back in normal NDVI units.
d.1984$Year="1984"
#d.1984$Date <- format(d.1984$Date,format = "%m-%d")
#new 1985 data:
d.1985<-data.frame(Date=rep(as.Date(c(
  "1984-05-09 UTC", "1984-05-25 UTC", "1984-06-10 UTC", 
  "1984-06-26 UTC", "1984-07-12 UTC")),2),
  Watershed = c(rep("ICB",5),rep("SCB",5))
)
d.1985$NDVI <- c(ICB_means/100,SCB_means/100) #Put back in normal NDVI units.
d.1985$Year="1985"
#d.1985$Date <- format(d.1985$Date,format = "%m-%d")

d<- full_join(d.1984,d.1985)

pdf("./Figures/MS/NDVI.pdf",height=4,width=4)
ggplot(d) +
  geom_point(aes(x = Date, y = NDVI, col = Watershed, pch = Year))+
  geom_smooth(aes(x = Date, y = NDVI, col = Watershed)) +
  scale_color_manual(values=c("deepskyblue3","darkgoldenrod2"))+
  theme_bw()+
  theme(legend.position = c(0.8,0.3),
        axis.text.x  = element_text(hjust=0.85))
dev.off()

