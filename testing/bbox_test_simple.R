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


    plot(bb$box)
    points(pnts)
    text(bb$box, labels=1:4 , cex=3, font=2)
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
                        rw_slp = 40, rw_ex = 3, rw_per = 1.5,
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
       main = paste("frqcut:", src_thresh,
                    "MD:", MD,
                    "MU:", MU,
                    "RelErr:", round(gpp$length.relerr, digits = 4)),
       cex.axis = 0.7, cex.main = 0.8)

  #plot(gpp$dem,
  #     col=terrain.colors(100),
  #     alpha=0.3,
  #     add=T,
  #     legend=F)

  plot( rasterCdf(gpp$gpp.parea), col = colors, add = T, alpha = 0.7)
  plot(bb_est, add = TRUE, lty = 2)
  plot(bb_obs, add = TRUE, border = "red", lty = 2)
  plot(slide_ply, add = TRUE, lty = 3)

  # add the DSM on top of the hillshade


}


plot.MedianDist <- function(dem, slide_ply, slide_src, src_thresh=0.5, MD, MU, hillshade,
                         map_ext = 3000){

  gpp <- pcmPerformance(dem, slide_plys = slide_ply, slide_src = slide_src,
                        rw_slp = 40, rw_ex = 3, rw_per = 1.5,
                        pcm_mu = MU, pcm_md = MD,
                        gpp_iter = 1000, buffer_ext = map_ext, buffer_source = 50,
                        plot_eval = FALSE, return_features = TRUE, saga_lib = saga)

  stop_pnts <- rasterToPoints(gpp$gpp.stop)
  stop_dists <- sp::spDistsN1(stop_pnts[,1:2], slide_src, longlat = FALSE)
  stop_pnts <- as.data.frame(stop_pnts)
  stop_pnts$dists <- stop_dists
  names(stop_pnts) <- c("x", "y", "freq", "dist")

  est_dist <- median(rep(stop_pnts$dist, times=stop_pnts$freq))

  obs_dist <- runoutGeom(slide_ply, gpp$dem)$length

  stop_relerr <- abs(obs_dist - est_dist) / obs_dist

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
       main = paste("frqcut:", src_thresh,
                    "MD:", MD,
                    "MU:", MU,
                    "\nBBre:", round(gpp$length.relerr, digits = 4),
                    "SPre:", round(stop_relerr, digits = 4)),
       cex.axis = 0.7, cex.main = 0.8)

  #plot(gpp$dem,
  #     col=terrain.colors(100),
  #     alpha=0.3,
  #     add=T,
  #     legend=F)

  plot( rasterCdf(gpp$gpp.stop), col = colors, add = T, alpha = 0.7)
  #plot(bb_est, add = TRUE, lty = 2)
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

runout_polygon <- runout_polygons[77,]
#runout_polygon <- runout_polygons[6,]
plot(runout_polygon)

sel_source_point  <- over(source_points, runout_polygon)
source_point <- source_points[!is.na(sel_source_point$objectid),]
plot(source_point, add = TRUE)


# TEST DIFFERENT PCM PARAMETERS ##################################

freq_thresh <- 0
massdrag <- 150
mu_vec <- c(0.04, 0.1, 0.3, 0.6)


par(mfrow=c(2,2))
for(i in 1:4){
  plot.PCMBBox(dem, slide_ply = runout_polygon, slide_src = source_point, src_thresh = 0.5,
               MD = massdrag, MU = mu_vec[i], hillshade, map_ext = 2500)
}




par(mfrow=c(2,2))
for(i in 1:4){
  plot.MedianDist(dem, slide_ply = runout_polygon, slide_src = source_point, src_thresh = 0.5,
               MD = massdrag, MU = mu_vec[i], hillshade, map_ext = 2500)
}



