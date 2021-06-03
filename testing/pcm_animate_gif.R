setwd("/home/jason/R/runoptGPP/")


# Load Packages and Data ####################################################################
library(runoptGPP)
library(raster)
library(rgdal)
library(sp)
library(Rsagacmd)
library(RColorBrewer)

# Functions

runoutGeom <- function(runout_plys, elev, ID = NULL) {
  #Create a dataframe to store the ID, length, width and (landslide) area
  if(is.null(ID)){
    bbDf <- data.frame(fid = 0:(length(runout_plys)-1), id = 1:length(runout_plys),
                       width = NA, length = NA, area = NA, surfacearea = NA,
                       maxelev = NA, minelev = NA, reachangle = NA)
  } else {
    bbDf <- data.frame(fid = 0:(length(runout_plys)-1), id = NA,
                       width = NA, length = NA, area = NA, surfacearea = NA,
                       maxelev = NA, minelev = NA, reachangle = NA)
  }


  n_features <- length(runout_plys@polygons)

  for (i in 1:n_features){
    runout_ply <- runout_plys[i,]

  #par(mfrow = c(2,2))
  #for (i in 1:100){
  #  runout_ply <- runout_polygons[i,]

    # Get vertices of

    if(!is.null(ID)){bbDf[i,]$id <- runout_ply@data[,ID]}
    #Get the vertices coordinates of from SpatialPolygons
    #pnts <- runout_ply@polygons[[1]]@Polygons[[1]]@coords
    pnts <- getVertices(runout_ply)
    #Calculate minimum area rectangle
    bb <- minbb(pnts)

    #Calcluate the difference in elevation between points
    #elevPly <- raster::extract(elev, pnts)
    #maxElev <- max(elevPly)


    #plot(bb$box)
    #points(pnts)
    #text(bb$box, labels=1:4 , cex=3, font=2)
  #}

    elevExt <- raster::extract(elev, bb$box[1:4,])
    dElev12 <- sqrt( (elevExt[1] - elevExt[2])^2)
    dElev14 <- sqrt( (elevExt[1] - elevExt[4])^2)
    #Determine (and calculate) the length and width based on delta elevation

    if(dElev12 > dElev14) {
      length <- EuclDist(bb$box[1,1], bb$box[1,2], bb$box[2,1], bb$box[2,2])
      width <- EuclDist(bb$box[1,1], bb$box[1,2], bb$box[4,1], bb$box[4,2])
    } else {
      length <- EuclDist(bb$box[1,1], bb$box[1,2], bb$box[4,1], bb$box[4,2])
      width <- EuclDist(bb$box[1,1], bb$box[1,2], bb$box[2,1], bb$box[2,2])
    }

    bbDf[i,]$length <- length #planar
    bbDf[i,]$width <- width #planar
    bbDf[i,]$area <- rgeos::gArea(runout_ply) #area of landslide (*not area of bbox) {Rgeos}


    #Calculate the 'true' surface area
    elevCrop <- raster::crop(elev, runout_ply) #crop and mask used to speed up calculation
    elevMask <- raster::mask(elevCrop, runout_ply)
    elevMaskSGDF <- as(elevMask , "SpatialGridDataFrame")
    bbDf[i,]$surfacearea <- sp::surfaceArea(elevMaskSGDF) #surfacearea of landslide {sp}
    #Calculate the max. and min. elevation
    elevPnts <- raster::extract(elev, pnts)
    bbDf[i,]$maxelev <- max(elevPnts)
    bbDf[i,]$minelev <- min(elevPnts)
    bbDf[i,]$reachangle <- 90 - (atan( bbDf[i,]$length / ( bbDf[i,]$maxelev- bbDf[i,]$minelev))*180/pi)

  }
  return(bbDf)
}


sp.pred <- raster::rasterToPolygons(gpp$gpp.parea, n = 4, dissolve = TRUE, na.rm=TRUE)
geom_pred <- runoutGeom(sp.pred, elev = dem)
geom_act <- runoutGeom(obs_poly, elev = dem)


plot.PCMBBox <- function(dem, slide_ply, slide_src, src_thresh=0.5, MD, MU, hillshade,
                         map_ext = 3000){

  gpp <- pcmPerformance(dem, slide_plys = slide_ply, slide_src = slide_src,
                        rw_slp = 30, rw_ex = 3, rw_per = 2,
                        pcm_mu = MU, pcm_md = MD,
                        gpp_iter = 1000, buffer_ext = map_ext, buffer_source = 50,
                        plot_eval = FALSE, return_features = TRUE, saga_lib = saga)

  # Get bounding boxes for obs. and simulated runout
  pred_thres <- rasterCdf(gpp$gpp.parea)
  pred_thres[pred_thres >= src_thresh] = 1
  pred_thres[pred_thres < src_thresh] = NA
  pred_poly <- raster::rasterToPolygons(pred_thres, n = 4, dissolve = TRUE, na.rm=TRUE)

  bb_est <- minBBoxSpatialPolygons(pred_poly)
  bb_obs <- minBBoxSpatialPolygons(slide_ply)

  colors <- brewer.pal(n = 9, name = "Blues")

  hs <- crop(hillshade, gpp$dem)

  plot(hs,
       col=grey.colors(100, start=1, end=0),
       legend=F,
       #cex.axis = 0.7, cex.main = 0.8,
       axes = FALSE, box = FALSE)

  mtext(paste("µ:", MU, "\nRelErr:", round(gpp$length.relerr, digits = 3)),
        side=3, adj=0, line=0, cex=1)


  #plot(gpp$dem,
  #     col=terrain.colors(100),
  #     alpha=0.3,
  #     add=T,
  #     legend=F)

  plot( rasterCdf(gpp$gpp.parea), col = viridis::viridis(10, direction = -1), add = T, alpha = 0.7)
  plot(bb_est, add = TRUE, lty = 2)
  plot(bb_obs, add = TRUE, border = "red", lty = 2)
  plot(slide_ply, add = TRUE, lty = 3)

  # add the DSM on top of the hillshade


}


#plot.PCMBBox(dem, slide_ply = runout_polygon, slide_src = source_point, src_thresh = 0,
#             MD = 40, MU = 0.1, hillshade, map_ext = 3000)
# LOAD OBJECTS #############################################

# Initiate a SAGA-GIS geoprocessing object
saga <- saga_gis(opt_lib = "sim_geomorphology")

# Set workspace
setwd("/home/jason/Data/Chile/")

# Load digital elevation model (DEM)
dem <- raster("elev_alos_12_5m.tif")

# Load runout source points
source_points <- readOGR(".", "debris_flow_source_points")

# Load runout track polygons and assign object ID based on row number
runout_polygons <- readOGR(".", "debris_flow_polys_sample")
runout_polygons$objectid <- 1:length(runout_polygons)


setwd("saga_data")
hillshade <- raster("hs_filtered_filled.sdat")

setwd("/home/jason/Scratch/Figures")




# EXPLORE LANDSLIDE SHAPES ######################################
#par(mfrow = c(3, 6))
#for(i in 1:length(runout_polygons)){
#  runout_poly <- runout_polygons[i,]
#  plot(runout_poly, main = i)
#}


# Examples of slides that change direction
#slides <- c(30, 6)


# SELECT SLIDE FOR EXAMPLE ######################################

#runout_polygon <- runout_polygons[77,] # Oringinal comparison
#runout_polygon <- runout_polygons[86,] # Good example of topographic control of Perla model - not box.

id = 77

runout_polygon <- runout_polygons[id,] # Good example of topographic control of Perla model - not box.

plot(runout_polygon)

sel_source_point  <- over(source_points, runout_polygon)
source_point <- source_points[!is.na(sel_source_point$objectid),]
plot(source_point, add = TRUE)


# TEST DIFFERENT PCM PARAMETERS ##################################
numFun <- function(x){

  # for numbers < 10000
  num_digi <- nchar(trunc(abs(x)))

  if( num_digi == 1){
    num <- paste0("000", x)
  } else if(num_digi == 2){
    num <- paste0("00", x)
  } else if(num_digi == 3) {
    num <- paste0("0", x)
  } else {
    num <- paste0(x)
  }

  return(num)

}

setwd("/home/jason/Scratch/Figures/pcm_animate")
freq_thresh <- 0
massdrag <- 40
mu_vec <- c(0.04, 0.1, 0.3, 0.6)
mu_vec <- seq(0.04, 0.6, by=0.04)


options(digits = 4)


#par(mfrow=c(2,2))
for(i in 1:length(mu_vec)){

  img_name <- paste0("pcm_", id, "_sim_lbl_it_", numFun(length(mu_vec) - i + 1), ".png")

  #png(filename=img_name, res = 300, width = 5.2, height = 4,
  #    units = "in", pointsize = 11)

  par(family = "Arial", mfrow = c(1,1), mar = c(0.5, 1, 4, 0.5))

  plot.PCMBBox(dem, slide_ply = runout_polygon, slide_src = source_point, src_thresh = 0.5,
               MD = massdrag, MU = mu_vec[i], hillshade, map_ext = 1500)

  #dev.off()
}


#https://gifmaker.me/
# PLOT SINGLE


# Explore friction coefficient sensitivity #################################
# In model, tan(theta) > mu : acceleration
#           tan(theta) < mu : deceleration
# https://www.cambridge.org/core/journals/journal-of-glaciology/article/twoparameter-model-of-snowavalanche-motion/B87923FFC6ADAF61B0079EEBCBD96F19



mu <- seq(0.04, 0.6, by=0.04)
slps <- atan(mu)*180/pi

df_control <- data.frame(mu = mu, theta = slps)
df_control

atan(0.11)*180/pi
atan(0.18)*180/pi


# Look at velocity #######################################################


runout_polygon <- runout_polygons[77,] # example in paper
#runout_polygon <- runout_polygons[86,]
plot(runout_polygon)

sel_source_point  <- over(source_points, runout_polygon)
source_point <- source_points[!is.na(sel_source_point$objectid),]
mds <- seq(20, 150, by=5)

max_vel <- rep(NA,length(mds))
min_vel <- rep(NA,length(mds))
med_vel <- rep(NA, length(mds))
iqr_vel <- rep(NA, length(mds))
avg_vel <- rep(NA, length(mds))
sd_vel <- rep(NA, length(mds))


for(i in 1:length(mds)){
  md <- mds[i]
  gpp <- pcmPerformance(dem, slide_plys = runout_polygon, slide_src = source_point,
                        rw_slp = 30, rw_ex = 3, rw_per = 2,
                        pcm_mu = 0.11, pcm_md = md,
                        gpp_iter = 1000, buffer_ext = 500, buffer_source = 50,
                        plot_eval = FALSE, return_features = TRUE, saga_lib = saga)

  #plot(gpp$gpp.maxvel)
  max_vel[i] <- max(getValues(gpp$gpp.maxvel), na.rm = TRUE)
  avg_vel[i] <- mean(getValues(gpp$gpp.maxvel), na.rm = TRUE)
  min_vel[i] <- min(getValues(gpp$gpp.maxvel), na.rm = TRUE)
  med_vel[i] <- median(getValues(gpp$gpp.maxvel), na.rm = TRUE)
  iqr_vel[i] <- sd(getValues(gpp$gpp.maxvel), na.rm = TRUE)
  sd_vel[i] <- IQR(getValues(gpp$gpp.maxvel), na.rm = TRUE)

}

df_vel <- data.frame(MD = mds,
                     max_vel = max_vel,
                     min_vel = min_vel,
                     avg_vel = avg_vel,
                     sd_vel = sd_vel,
                     med_vel = med_vel,
                     iqr_vel = iqr_vel)
df_vel

#save(df_vel, file = "md_velocity_df.Rd")

(load("md_velocity_df.Rd"))

plot(df_vel$MD, df_vel$max_vel, ylab = "Max. Velocity (m/s)", xlab = "Mass-to-drag ratio (m)")
plot(df_vel$MD, df_vel$avg_vel, ylab = "Mean Velocity (m/s)", xlab = "Mass-to-drag ratio (m)",
     col = "Red", type = "l",
     ylim = c(0, max(df_vel$avg_vel + df_vel$sd_vel*1.96 )))
#plot(mds, sd_vel, ylab = "Max. Velocity (m/s)", xlab = "Mass-to-drag ratio (m)")
#plot(mds, iqr_vel, ylab = "Mean Velocity (m/s)", xlab = "Mass-to-drag ratio (m)")
lines(df_vel$MD, df_vel$avg_vel + df_vel$sd_vel*1.96)
lines(df_vel$MD, df_vel$avg_vel - df_vel$sd_vel*1.96)

plot(df_vel$MD, df_vel$max_vel, ylab = "Max. Velocity (m/s)", xlab = "Mass-to-drag ratio (m)",
     ylim = c(0, max(max_vel)))
lines(df_vel$MD, df_vel$avg_vel)
lines(df_vel$MD, df_vel$min_vel)

# Test All ###################################################################


numFun <- function(x){

  # for numbers < 10000
  num_digi <- nchar(trunc(abs(x)))

  if( num_digi == 1){
    num <- paste0("000", x)
  } else if(num_digi == 2){
    num <- paste0("00", x)
  } else if(num_digi == 3) {
    num <- paste0("0", x)
  } else {
    num <- paste0(x)
  }

  return(num)

}

setwd("/home/jason/Scratch/Figures/pcm_animate")



# Run loop

ids <- 1:100
freq_thresh <- 0
massdrag <- 40
mu_vec <- c(0.04, 0.1, 0.3, 0.6)
mu_vec <- seq(0.04, 0.6, by=0.04)
options(digits = 4)


for(ID in 1:length(ids)){

  runout_polygon <- runout_polygons[ID,] # Good example of topographic control of Perla model - not box.

  plot(runout_polygon)

  sel_source_point  <- over(source_points, runout_polygon)
  source_point <- source_points[!is.na(sel_source_point$objectid),]
  plot(source_point, add = TRUE)


  #par(mfrow=c(2,2))
  for(i in 1:length(mu_vec)){

    img_name <- paste0("pcm_", ID, "_sim_lbl_it_", numFun(length(mu_vec) - i + 1), ".png")

    png(filename=img_name, res = 300, width = 5.2, height = 4,
        units = "in", pointsize = 11)

    par(family = "Arial", mfrow = c(1,1), mar = c(0.5, 1, 4, 0.5))

    plot.PCMBBox(dem, slide_ply = runout_polygon, slide_src = source_point, src_thresh = 0.5,
                 MD = massdrag, MU = mu_vec[i], hillshade, map_ext = 1500)

    dev.off()
  }


  #https://gifmaker.me/

}




