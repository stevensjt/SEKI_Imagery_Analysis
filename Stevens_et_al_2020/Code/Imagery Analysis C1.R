##This is a component of analysis code for Stevens et al. 2020: Forest vegetation change and its impacts on soil water following 47 years of managed wildfire. Ecosystems, In Press.
##Contact: Jens Stevens, stevensjt@gmail.com
##Imagery analyses: This code produces Figure C1

####0. Load libraries####
require(devtools)
#install_version("ggplot2", version = "2.2.1", repos = "http://cran.us.r-project.org") #skip updates
library(ggplot2) #for ggplot; v 2.2.1
library(raster) #for raster(); version 2.8-4
library(rgdal) #for readOGR(); version 1.3-6
library(reshape2) #for melt(); version 1.4.3
library(gridExtra) #for grid.arrange(); version 2.3
#library(igraph) #for clump(); version 1.2.2
#library(ggcorrplot)
#library(tableHTML)

#library(grid) #for grid.text(); version 3.6.1; deprecated?

####1. Data processing- only run once, don't need to re-run####
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

####2. Read and process data####
#Read data#
r73scb <- #Load the processed 1973 veg raster
  raster("./Data/Geospatial/1973_veg_classes.tif")
r14scb <- #Load the processed 2014 veg raster
  raster("./Data/Geospatial/2014_veg_classes.tif")
#Metadata: EPSG Coordinate Ref System: 26911
#Metadata: pixel values for veg classes are: 
#1 = shrub, 2 = sparse meadow, 3 = granite
#4 = mixed conifer, 5 = water, 6 = dense meadow
perims <- readOGR("./Data/Geospatial/Sugarloaf Fires 1973-2003.shp")
#Metadata: EPSG Coordinate Ref System: 26911

#Process data#
perims <- #reproject perimeters (already the same but just to be safe)
  spTransform(perims, CRSobj = crs(r73scb)) 
r73_pts <- rasterToPoints(r73scb,spatial=TRUE)
r14_pts <- rasterToPoints(r14scb,spatial=TRUE)
r14_pts$n_fires <- r73_pts$n_fires <- #Count number of overlapping fires
  sapply(over(r73_pts, geometry(perims), returnList = TRUE), length)
#93 pixels had 4 fires, small n (93 pixels at 0.16 ha/pixel = 14.88 ha). 1360 pixels had 3 fires, chi-squared test was not converging. Converting to 3+4 burns to 2. 
r73_pts$n_fires[r73_pts$n_fires > 2] <- 2 
r14_pts$n_fires[r14_pts$n_fires > 2] <- 2 
length(r73_pts$n_fires[r73_pts$n_fires==0])* 0.16 #Area unburned = 5724 ha
#slight discrepancy from published paper (5707) 
#is due to realignment of rasters to UTM11 from unknown projection, for metadata purposes.
length(r73_pts$n_fires[r73_pts$n_fires==1])* 0.16 #Area burning once
length(r73_pts$n_fires[r73_pts$n_fires==2])* 0.16 #Area burning twice or more

#Plot data#
gridded(r73_pts) <- TRUE #For more efficient plotting, converts to "SpatialPixels"
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


####3. Statistical analyses and plots####

###Set up change analysis
sum_table <- 
  data.frame(
    veg_type = c("shrub", "sparse meadow", "mixed conifer", "dense meadow")
    )
scenarios <- c("watershed","no burns", "one burn", "two burns")
change_mat_total <- change_mat_exp <- gcp <- list()

for(s in scenarios){ #Begin scenario loop (different numbers of times burned)
  if(s == "watershed"){
    v73 <- getValues(r73scb)
    v14 <- getValues(r14scb)
    sum_table$pix_73 <- as.vector(table(v73[v73!=3 & v73!=5]))
    sum_table$pix_14 <- as.vector(table(v14[v14!=3 & v14!=5]))
    l = "a"
  }
  if(s == "no burns"){
    v73 <- getValues(r73scb)[which(r73_pts$n_fires==0)]
    v14 <- getValues(r14scb)[which(r14_pts$n_fires==0)]
    sum_table$pix_73 <- as.vector(table(v73[v73!=3 & v73!=5]))
    sum_table$pix_14 <- as.vector(table(v14[v14!=3 & v14!=5]))
    l = "b"
  }
  if(s == "one burn"){
    v73 <- getValues(r73scb)[which(r73_pts$n_fires==1)]
    v14 <- getValues(r14scb)[which(r14_pts$n_fires==1)]
    sum_table$pix_73 <- as.vector(table(v73[v73!=3 & v73!=5]))
    sum_table$pix_14 <- as.vector(table(v14[v14!=3 & v14!=5]))
    l = "c"
  }
  if(s == "two burns"){
    v73 <- getValues(r73scb)[which(r73_pts$n_fires==2)]
    v14 <- getValues(r14scb)[which(r14_pts$n_fires==2)]
    sum_table$pix_73 <- as.vector(table(v73[v73!=3 & v73!=5]))
    sum_table$pix_14 <- as.vector(table(v14[v14!=3 & v14!=5]))
    l = "d"
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
  
  change_mat_total[[s]] <- change_mat_total[[s]][c(4,1,2,6),c(4,1,2,6)]
  rownames(change_mat_total[[s]]) <- 
    colnames(change_mat_total[[s]]) <- c("mixed \nconifer", "shrub", 
                                         "sparse \nmeadow", "dense \nmeadow")
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
  
  X2 <- chisq.test(as.vector(change_mat_total[[s]])+1, #+1 to deal with convergence issues at 0's
                    p = as.vector((change_mat_exp[[s]]+1)/sum(change_mat_exp[[s]]+1)))
  print(X2)
  deviance_prop <- (change_mat_total[[s]]-change_mat_exp[[s]]) / (change_mat_exp[[s]]+1)
  deviance_prop_melt <- list()
  deviance_prop_melt[[s]] <- melt(deviance_prop)
  names(deviance_prop_melt[[s]]) <- c("Y1973","Y2014","residual_prop")
  deviance_prop_melt[[s]]$Y1973 <- 
    factor(deviance_prop_melt[[s]]$Y1973, rev(levels(deviance_prop_melt[[s]]$Y1973) ) )
  deviance_prop_melt[[s]]$cells <- melt(change_mat_total[[s]])$value

  #Plot within loop, for later combination
  gcp[[s]] <- 
    ggplot(deviance_prop_melt[[s]],aes(Y2014,Y1973)) +
    geom_tile(aes(fill = residual_prop),colour = "white") +
    geom_text(aes(label = cells)) +
    scale_fill_gradientn(name = "residual \nproportion", 
                         colors = c("darkred", "white","cornflowerblue"), 
                         limits = c(-1,1)) +
    theme_grey(base_size = 11) + 
    labs(title = paste(l,"          ",s), x = "2014", y = 1973) +
    scale_x_discrete(position = "top") +
    theme(axis.ticks = element_blank(), axis.text.x = element_text(
      size = 11 *0.9, angle = 45, hjust = 0, colour = "grey50"))
  
  
}##End scenario loop

pdf("./Figures/FigC1.pdf", width = 8, height = 6)
grid.arrange(gcp[[1]], gcp[[2]], gcp[[3]], gcp[[4]], 
             ncol = 2)
dev.off()
#Note: Slight changes in pixel counts relative to Fig. C1 in Manuscript
#is due to realignment of rasters to UTM11 from unknown projection, for metadata purposes.
#Conclusions are unchanged.


