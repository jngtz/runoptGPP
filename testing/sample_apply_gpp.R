# LOAD PACKAGES ################################################################

library(runoptGPP)
library(rgdal)
library(raster)
library(rgeos)
library(sf)

library(Rsagacmd)

# Initiate a SAGA-GIS geoprocessing object
saga <- saga_gis(opt_lib = "sim_geomorphology")

# LOAD DATA ####################################################################

setwd("/home/jason/Data/Chile/")
# elevation model
# Load digital elevation model (DEM)
dem <- raster("elev_alos_12_5m_no_sinks.tif")

# Load runout source points
source_points <- readOGR(".", "debris_flow_source_points")

# Load runout track polygons and assign object ID based on row number
runout_polygons <- readOGR(".", "debris_flow_polys_sample")
runout_polygons$objectid <- 1:length(runout_polygons)

#crs(slide_point_vec) <- crs(slide_poly_vec)

buffer_source = 50

source_buffer <- rgeos::gBuffer(source_points, width = buffer_source)
source_grid <- raster::rasterize(source_buffer, dem, field=1 )
source_grid <- raster::mask(source_grid, runout_polygons)
source_plot <- raster::rasterToPolygons(source_grid)

gpp_rw <- saga$sim_geomorphology$gravitational_process_path_model(dem = dem_grid, release_areas = source_grid,
                                                                   process_path_model = 1,
                                                                   rw_slope_thres = 40,
                                                                   rw_exponent = 3.0,
                                                                   rw_persistence = 1.9,
                                                                   gpp_iterations = gpp_iter)

cdf_rw_parea <- rasterCdf(gpp_rw$process_area)

gpp_pcm <- saga$sim_geomorphology$gravitational_process_path_model(dem = dem, release_areas = source_grid,
                                                                   #friction_mu_grid = mu.grid,
                                                                   process_path_model = 1,
                                                                   rw_slope_thres = 40,
                                                                   rw_exponent = 3.0,
                                                                   rw_persistence = 1.9,
                                                                   gpp_iterations = 1000,
                                                                   friction_model = 5,
                                                                   friction_mu = 0.11,
                                                                   friction_mass_to_drag = 40)

# Rescale to values from 0 to 1
cdf_pcm_parea <- rasterCdf(gpp_pcm$process_area)
