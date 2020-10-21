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

# Load shapefile of sub catchment polygons
sub_catchments <- readOGR("sub_catchments.shp")

# Load and crop DEM for a single sub catchment
dem <- raster("elev_alos_12_5m.tif")
dem <- crop(dem, sub_catchments[sub_catchments$ID == 56,])
dem <- mask(dem_crop, sub_catchments[sub_catchments$ID == 56,])

# Load and crop source area prediction raster
source_pred <- raster("source_pred_gam.tif")
source_pred <- crop(source_pred, sub_catchments[sub_catchments$ID == 56,])
source_pred <- mask(source_pred, sub_catchments[sub_catchments$ID == 56,])

# For visualization make a hillshade model from the DEM
slope <- terrain(dem, opt='slope')
aspect <- terrain(dem, opt='aspect')
hillshade <- hillShade(slope, aspect, angle=40, direction=270)


# Plot source prediction area
library(ggplot2)
library(patchwork) # Tile multiple ggplots
library(ggnewscale) # Allow multiple scales per ggplot

# Aggregate data for faster mapping
hillshade_agg <- aggregate(hillshade, fact = 5, fun = mean)
source_agg <- aggregate(source_pred, fact = 5, fun = mean)

hillshade_df <- as.data.frame(hillshade_agg, xy = TRUE)
hillshade_df <- hillshade_df[!is.na(hillshade_df[,3]),]

source_df <- as.data.frame(source_agg, xy = TRUE)
source_df <- source_df[!is.na(source_df[,3]),]
names(source_df) <- c("x", "y", "pred")

map.srcpred <- ggplot() +
  geom_tile(data=hillshade_df, aes(x=x, y=y, fill = layer),
              show.legend = FALSE) +
  scale_fill_gradient(high = "white", low = "black", na.value = "#FFFFFF") +

  new_scale("fill") +
  geom_tile(data=source_df, aes(x=x, y=y, fill = pred) ) +
  scale_fill_viridis_c(name = "Source area\nprediction",  alpha = 0.6, direction = -1) +

  xlab("Easting (m)") +
  ylab("Northing (m)") +
  coord_fixed() +
  theme_light() +
  theme(text = element_text(size = 9), axis.title = element_text(size = 9),
        axis.text = element_text(size = 6), axis.text.y = element_text(angle = 90))
map.srcpred

# Get optimal parameters #######################################################
setwd("/home/jason/Scratch/GPP_PCM_Paper")
(load("pcm_gridsearch_multi.Rd"))

pcm_opt <- pcmGetOpt(pcm_gridsearch_multi)

setwd("/home/jason/Scratch/GPP_RW_Paper")
(load("rw_gridsearch_multi.Rd"))

rw_opt <- rwGetOpt(rw_gridsearch_multi)

# Classify source area #########################################################

# Will identify source areas with a prediction probability > 0.7
source_area <- rasterThreshold(source_pred, c(0.7, Inf))

# Plot source area
sourcearea_agg <- aggregate(source_area, fact = 5, fun = mean)
sourcearea_df <- as.data.frame(sourcearea_agg, xy = TRUE)
sourcearea_df <- sourcearea_df[!is.na(sourcearea_df[,3]),]
names(sourcearea_df) <- c("x", "y", "pred")

map.srcarea <- ggplot() +
  geom_tile(data=hillshade_df, aes(x=x, y=y, fill = layer),
              show.legend = FALSE) +
  scale_fill_gradient(high = "white", low = "black", na.value = "#FFFFFF") +

  new_scale("fill") +
  geom_tile(data=sourcearea_df, alpha = 0.6, aes(x=x, y=y, fill = "") ) +
  scale_fill_manual(name = "Source area", values = "#e74c3c" ) +

  xlab("Easting (m)") +
  ylab("Northing (m)") +
  coord_fixed() +
  theme_light() +
  theme(text = element_text(size = 9), axis.title = element_text(size = 9),
        axis.text = element_text(size = 6), axis.text.y = element_text(angle = 90))
map.srcarea


# Apply GPP RW PCM runout model ################################################

gpp <- saga$sim_geomorphology$gravitational_process_path_model(dem = dem, release_areas = source_area,
                                                               process_path_model = 1,
                                                               rw_slope_thres = rw_opt$rw_slp_opt,
                                                               rw_exponent = rw_opt$rw_exp_opt,
                                                               rw_persistence = rw_opt$rw_per_opt,
                                                               gpp_iterations = 1000,
                                                               friction_model = 5,
                                                               friction_mu = pcm_opt$pcm_mu,
                                                               friction_mass_to_drag = pcm_opt$pcm_md)

# Save/Export results
setwd("/home/jason/Scratch/GPP_PCM_Paper")
writeRaster(gpp$process_area, filename="gpp_process_area.tif", format = "GTiff")
writeRaster(gpp$max_velocity, filename="gpp_max_velocity.tif", format = "GTiff")
writeRaster(gpp$stop_positions, filename="gpp_stop_positions.tif", format = "GTiff")

parea <- raster("gpp_process_area.tif")
maxvel <- raster("gpp_max_velocity.tif")
stopp <- raster("gpp_stop_positions.tif")

# Quantile classification of runout and stop position frequencies
parea_cdf <- rasterCdf(parea)
stopp_cdf <- rasterCdf(stopp)

# Aggregate for quicker mapping
parea_cdf <- aggregate(parea_cdf, fact = 5, fun = mean)
stopp_cdf <- aggregate(stopp_cdf, fact = 5, fun = mean)
maxvel <- aggregate(maxvel, fact = 5, fun = mean)

# Format to dataframe for ggplot2
parea_df <- as.data.frame(parea_cdf, xy = TRUE)
parea_df <- parea_df[!is.na(parea_df[,3]),]

stopp_df <- as.data.frame(stopp_cdf, xy = TRUE)
stopp_df <- stopp_df[!is.na(stopp_df[,3]),]

maxvel_df <- as.data.frame(maxvel, xy = TRUE)
maxvel_df <- maxvel_df[!is.na(maxvel_df[,3]),]

# Process area map

map.parea <- ggplot() +
  geom_tile(data=hillshade_df, aes(x=x, y=y, fill = layer),
              show.legend = FALSE) +
  scale_fill_gradient(high = "white", low = "black", na.value = "#FFFFFF") +

  new_scale("fill") +
  geom_tile(data=parea_df, aes(x=x, y=y, fill = gpp_process_area) ) +
  scale_fill_viridis_c(name = "Runout impact\nindex",  alpha = 0.6, direction = -1) +

  xlab("Easting (m)") +
  ylab("Northing (m)") +
  coord_fixed() +
  theme_light() +
  theme(text = element_text(size = 9), axis.title = element_text(size = 9),
        axis.text = element_text(size = 6), axis.text.y = element_text(angle = 90))
map.parea
setwd("/home/jason/Scratch/Figures")
ggsave("regional_runout_impact_map.png", dpi = 300, width = 4, height = 6, units = "in")


# Stop positions map

map.stopp <- ggplot() +
  geom_tile(data=hillshade_df, aes(x=x, y=y, fill = layer),
            show.legend = FALSE) +
  scale_fill_gradient(high = "white", low = "black", na.value = "#FFFFFF") +

  new_scale("fill") +
  geom_tile(data=stopp_df, aes(x=x, y=y, fill = gpp_stop_positions) ) +
  scale_fill_viridis_c(name = "Runout stop\npositions\n(Quantiles)",  alpha = 0.6, direction = -1) +

  xlab("Easting (m)") +
  ylab("Northing (m)") +
  coord_fixed() +
  theme_light() +
  theme(text = element_text(size = 9), axis.title = element_text(size = 9),
        axis.text = element_text(size = 6), axis.text.y = element_text(angle = 90))
map.stopp

# Max. velocities

map.maxvel <- ggplot() +
  geom_tile(data=hillshade_df, aes(x=x, y=y, fill = layer),
            show.legend = FALSE) +
  scale_fill_gradient(high = "white", low = "black", na.value = "#FFFFFF") +

  new_scale("fill") +
  geom_tile(data=maxvel_df, aes(x=x, y=y, fill = gpp_max_velocity) ) +
  scale_fill_viridis_c(name = "Maximum\nvelocities\n(m/s)",  alpha = 0.6, direction = -1) +

  xlab("Easting (m)") +
  ylab("Northing (m)") +
  coord_fixed() +
  theme_light() +
  theme(text = element_text(size = 9), axis.title = element_text(size = 9),
        axis.text = element_text(size = 6), axis.text.y = element_text(angle = 90))
map.maxvel

map.parea + map.maxvel

map.srcarea + map.parea
