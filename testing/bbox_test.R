setwd("/home/jason/R/runoptGPP/")


# Load Packages and Data ####################################################################
library(runoptGPP)
library(raster)
library(rgdal)
library(sp)
library(Rsagacmd)

# Functions

pretty_PCM_example <- function(dem, slide_ply, slide_src, src_thresh=0.5, MD, MU, hillshade){

  require(broom)
  require(tibble)
  require(dplyr)


  gpp <- pcmPerformance(dem, slide_plys = slide_ply, slide_src = slide_src,
                        rw_slp = 40, rw_ex = 3, rw_per = 1.5,
                        pcm_mu = MU, pcm_md = MD,
                        gpp_iter = 1000, buffer_ext = 600, buffer_source = 50,
                        plot_eval = TRUE, return_features = TRUE, saga_lib = saga)

  # Get bounding boxes for obs. and simulated runout
  pred_thres <- rasterCdf(gpp$gpp.parea)
  pred_thres[pred_thres >= src_thresh] = 1
  pred_thres[pred_thres < src_thresh] = NA
  pred_poly <- raster::rasterToPolygons(pred_thres, n = 4, dissolve = TRUE, na.rm=TRUE)

  bb_est <- minBBoxSpatialPolygons(pred_poly)
  bb_obs <- minBBoxSpatialPolygons(slide_ply)

  obs_df <- tidy(gpp$actual.poly, region = "objectid" )
  bbEst_df <- tidy(bb_est, region = "value" )
  bbObs_df <- tidy(bb_obs, region = "value" )

  obs_df <- obs_df %>% add_column(label = "Obs. debris flow runout")
  bbObs_df <- bbObs_df %>% add_column(label = "Obs. min. area bounding box")
  bbEst_df <- bbEst_df %>% add_column(label = "Est. min. area bounding box")

  obs_df <- obs_df %>% add_column(fill_ply = 1)
  bbObs_df <- bbObs_df %>% add_column(fill_ply = 2)
  bbEst_df <- bbEst_df %>% add_column(fill_ply = 3)

  layers <- bind_rows(list(obs_df, bbObs_df, bbEst_df))

  parea_df <- as.data.frame(rasterCdf(gpp$gpp.parea), xy = TRUE)
  names(parea_df)[3] <- "value"
  parea_df <- parea_df[!is.na(parea_df$value),]

  hs <- crop(hillshade, gpp$dem)

  hs_df <-  as.data.frame(hs, xy = TRUE)

  return(list(
    parea_df = parea_df,
    hs_df = hs_df,
    layers = layers,
    mu = MU,
    md = MD,
    roc = gpp$roc,
    length.relerr = gpp$length.relerr
  ))
}


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

# Select a debris flow and source point for this example


par(mfrow = c(3, 6))
for(i in 1:length(runout_polygons)){
  runout_poly <- runout_polygons[i,]
  plot(runout_poly, main = i)
}


# Examples of slides that change direction
slides <- c(30, 6)

runout_polygon <- runout_polygons[77,]
#runout_polygon <- runout_polygons[86,]
plot(runout_polygon)

sel_source_point  <- over(source_points, runout_polygon)
source_point <- source_points[!is.na(sel_source_point$objectid),]
plot(source_point, add = TRUE)


par(mfrow = c(2, 2))
pcm <- pcmPerformance(dem, slide_plys = runout_polygon, slide_src = source_point,
                      rw_slp = 40, rw_ex = 3, rw_per = 1.5,
                      pcm_mu = 0.15, pcm_md = 120,
                      gpp_iter = 1000, buffer_ext = 500, buffer_source = 50,
                      plot_eval = TRUE, return_features = TRUE, saga_lib = saga)


pcm <- pcmPerformance(dem, slide_plys = runout_polygon, slide_src = source_point,
                      rw_slp = 40, rw_ex = 3, rw_per = 1.5,
                      pcm_mu = 0.04, pcm_md = 120,
                      gpp_iter = 1000, buffer_ext = 1000, buffer_source = 50,
                      plot_eval = TRUE, return_features = TRUE, saga_lib = saga)

plot(pcm$)
