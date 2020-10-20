setwd("/home/jason/R/runout.opt/")

library(devtools)
build()
install()


# Load packages and data #######################################################
library(runout.opt)
library(raster)
library(rgdal)
library(Rsagacmd)

# Initiate a SAGA-GIS geoprocessing object
saga <- saga_gis(opt_lib = "sim_geomorphology")

# Set workspace
setwd("/home/jason/Data/Chile/")

# Load digital elevation model (DEM)
dem <- raster("elev_alos_12_5m.tif")

# Load runout source points
source_points <- readOGR(".", "debris_flow_source_points")

# Load runout track polygons
runout_polygons <- readOGR(".", "debris_flow_polys_sample")

#load shapefile of subcatchment polygons
sub_catchments <- readOGR("sub_catchments.shp")

# Load sampling mask
mask <- raster("mask_dflow_repo.tif")


# Load predicted runout rasters calc. from different source area thresholds  ###

setwd("/home/jason/Scratch/GPP_Subcatch_Paper")
sub_catch.ids <- as.numeric(sub_catchments$OBJECTID)

# Define thresholds (cutoffs) that were used for runout simulation
cutoffs = seq(.5,.95, by=0.05)
gpp_raster <- list()

for(i in 1:length(cutoffs)){

  raster_name <- paste("gpp_all_cutoff", cutoffs[i], ".tif", sep="")
  gpp_raster[[i]] <- raster(raster_name)

}

gpp_raster


# Compute performance and coverage area for  thresholds ########################

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

#(load("source_pred_roc_cutoffs.Rd"))

# Plot results #################################################################

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


