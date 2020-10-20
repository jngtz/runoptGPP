# FUNCTIONS ####################################################################

RasterCdfRescale <- function(x){
  x_ecdf <- ecdf(getValues(x))
  prob_x <- setValues(x, x_ecdf(getValues(x)))
  prob_x

}


RasterThreshold <- function(x, prob_range = c(0.9, Inf)){
  if(prob_range[2] != Inf){

    m <- c(-Inf, prob_range[1], NA, prob_range, 1, prob_range[2], Inf, NA)
    rclmat <- matrix(m, ncol=3, byrow=TRUE)

  } else {

    m <- c(-Inf, prob_range[1], NA, prob_range, 1)
    rclmat <- matrix(m, ncol=3, byrow=TRUE)

  }

  rc <- reclassify(x, rclmat)
  rc

}


gppSourceThreshold <- function(sub_catch.id, sub_catch_vec, dem, source_pred, cutoffs,
                               rw_slp, rw_exp, rw_per, pcm_mu, pcm_md,
                               gpp_iter = 1000){

  dem_grid <- crop(dem, sub_catch_vec[sub_catch.id,])
  #dem_grid <- mask(dem_grid, sub_catch_vec[sub_catch.id,] ) # Don't mask border to

  source_grid <- crop(source_pred, sub_catch_vec[sub_catch.id,])
  source_grid <- mask(source_grid, sub_catch_vec[sub_catch.id,])

  for(i in 1:length(cutoffs)){

    source_threshold <- RasterThreshold(source_grid, prob_range = c(cutoffs[i], Inf))

    #just in case no source cells are present...
    if(freq(source_threshold, value = 1) != 0) {
      gpp <- saga$sim_geomorphology$gravitational_process_path_model(dem = dem_grid,
                                                                     release_areas = source_threshold,
                                                                     process_path_model = 1,
                                                                     rw_slope_thres = rw_slp,
                                                                     rw_exponent = rw_exp,
                                                                     rw_persistence = rw_per,
                                                                     gpp_iterations = gpp_iter,
                                                                     friction_model = 5,
                                                                     friction_mu = pcm_mu,
                                                                     friction_mass_to_drag = pcm_md)

      write_name <- paste("gpp_subcatch_", sub_catch.id, "_cutoff_", cutoffs[i], sep="")

      writeRaster(gpp$process_area, paste(write_name, "_parea.tif", sep = ""),
                  format = "GTiff", overwrite = TRUE)
      writeRaster(gpp$stop_positions, paste(write_name, "_stopp.tif", sep = ""),
                  format = "GTiff", overwrite = TRUE)

      rm(gpp)

    }

  }
}


# LOAD PACKAGES AND DATA #######################################################

library(runout.opt)
library(raster)
library(rgdal)
library(foreach)

library(Rsagacmd)
saga <- saga_gis()

# Initiate a SAGA-GIS geoprocessing object
saga <- saga_gis(opt_lib = "sim_geomorphology")

# Set workspace
setwd("/home/jason/Data/Chile/")

# Load digital elevation model (DEM)
dem_all <- raster("dem_alos_12_5m _no sinks.tif")

#load shapefile of subcatchment polygons
sub_catch.vec = readOGR("sub_catchments.shp")

# Load runout source points
source_points <- readOGR(".", "debris_flow_source_points")

# Load runout track polygons
runout_polygons <- readOGR(".", "debris_flow_polys_sample")

# Load sampling mask
mask <- raster("mask_dflow_repo.tif")

# source area prediction [GAM]
source_pred <- raster("gam_dflow_v2repo_elv_slp_carea_plcrv_dstflt_meshdenoise.tif")

# resample source_pred to dem - seems to be a problem in extent when loading to
# SAGA
#source_pred <- resample(source_pred, dem_all, method = "ngb")
#writeRaster(source_pred, "rsmp_gam_dflow_v2repo_elv_slp_carea_plcrv_dstflt_meshdenoise.tif", format = "GTiff", overwrite = TRUE)

source_pred_all <- raster("rsmp_gam_dflow_v2repo_elv_slp_carea_plcrv_dstflt_meshdenoise.tif")


# LOAD OPTIMAL PARAMETER SET FOR RW AND PCM ####################################

setwd("/home/jason/Scratch/GPP_PCM_Paper")
(load("pcm_opt_params.Rd"))

setwd("/home/jason/Scratch/GPP_RW_Paper")
(load("rw_opt_params.Rd"))

# RUN RW PCM MODEL FOR DIFFERENT SOURCE THRESHOLDS #############################

# By catchment area for parallel processing
setwd("/home/jason/Scratch/GPP_Subcatch_Paper")
sub_catch.ids <- as.numeric(sub_catch.vec$OBJECTID)

cl <- parallel::makeCluster(32)
doParallel::registerDoParallel(cl)

results.list <-
  foreach(sub_catch_id=sub_catch.ids, .packages=c('rgdal','raster', 'Rsagacmd')) %dopar% {
    gppSourceThreshold(sub_catch.id = sub_catch_id,
                       sub_catch_vec = sub_catch.vec,
                       dem = dem_all,
                       source_pred = source_pred_all,
                       cutoffs = seq(.5,.95, by=0.05),
                       rw_slp = rw_opt$rw_slp_opt,
                       rw_exp = rw_opt$rw_exp_opt,
                       rw_per = rw_opt$rw_per_opt,
                       pcm_mu = pcm_opt$pcm_mu,
                       pcm_md = pcm_opt$pcm_md,
                       gpp_iter = 1000)
  }
parallel::stopCluster(cl)


# MERGE PROCESS AREAS ##########################################################
sub_catch.ids <- as.numeric(sub_catch.vec$OBJECTID)

cutoffs = seq(.5,.9, by=0.05)
gpp_raster <- list()

# Determine if runout raster exits...

for(i in 1:length(cutoffs)){

  raster_exists <- rep(NA, length(sub_catch.ids))

  for(k in 1:length(sub_catch.ids)){

    raster_name <- paste("gpp_subcatch_", sub_catch.ids[k], "_cutoff_", cutoffs[i], "_parea.tif", sep="")

    if(file.exists(raster_name)){
      raster_exists[k] <- TRUE
    } else {
      raster_exists[k] <- FALSE
    }

  }

  #select only sub catchments that exists
  sub_catch_exists <- sub_catch.ids[raster_exists]

  raster_list <- list()
  for(j in 1:length(sub_catch_exists)){

    raster_name <- paste("gpp_subcatch_", sub_catch_exists[j], "_cutoff_", cutoffs[i], "_parea.tif", sep="")
    raster_list[[j]] <- raster(raster_name)

  }

  merge_rasters <- do.call(merge,raster_list)
  gpp_raster[[i]] <- merge_rasters
  write_name <- paste("gpp_all_cutoff", cutoffs[i], ".tif", sep="")
  writeRaster(merge_rasters, write_name, format = "GTiff", overwrite = TRUE)

}

gpp_raster


# LOAD RUNOUT FROM DIFFERENT SOURCE THRESHOLDS #################################

setwd("/home/jason/Scratch/GPP_Subcatch_Paper")
sub_catch.ids <- as.numeric(sub_catch.vec$OBJECTID)

# Define thresholds (cutoffs) that were used for runout simulation
cutoffs = seq(.5,.95, by=0.05)
gpp_raster <- list()

for(i in 1:length(cutoffs)){

  raster_name <- paste("gpp_all_cutoff", cutoffs[i], ".tif", sep="")
  gpp_raster[[i]] <- raster(raster_name)

}

gpp_raster


# COMPUTE SOURCE THRESHOLD PERFORMANCE #########################################

auroc_cutoffs <- rep(NA, length(cutoffs))
area_cutoffs <- rep(NA, length(cutoffs))

# Convert runout polygons to raster for sampling
runout_area <- rasterize(runout_polygons, dem, field=1, background = NA)

for(i in 1:length(cutoffs)){
  print(cutoffs[i])
  auroc_cutoffs[i] <- rocParea(gpp_raster[[i]], runout_area, mask, smp_size = 1000)
  area_cutoffs[i] <- areaPerParea(gpp_raster[[i]], dem)
}


cutoffs_df <- data.frame(cutoffs = cutoffs, auroc = auroc_cutoffs,
                         per_area = area_cutoffs)

cutoffs_df

save(cutoffs_df, file="performance_sourcearea_thresholds.Rd")


# PLOT RESULTS #################################################################

setwd("/home/jason/Scratch/Figures")

(load("performance_sourcearea_thresholds.Rd"))

png(filename="src_area_auroc_threshold.png", res = 300, width = 7.5, height = 3,
    units = "in", pointsize = 11)

par(family = "Arial", mfrow = c(1,2), mar = c(4, 3.5, 1, 1),
    mgp = c(2.1, 0.75, 0))

plot(cutoffs_df$cutoffs, cutoffs_df$auroc,
     xlab = "Source area threshold",
     ylab = "AUROC",
     cex.axis = 1, cex.lab = 1, type = "l", main = "" )
points(cutoffs_df$cutoffs, cutoffs_df$auroc, pch = 20)

plot(cutoffs_df$cutoffs, cutoffs_df$per_area,
     xlab = "Source area threshold",
     ylab = "Runout area (%)",
     cex.axis = 1, cex.lab = 1, type = "l", main = "" )
points(cutoffs_df$cutoffs, cutoffs_df$per_area, pch = 20)


