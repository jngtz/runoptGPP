setwd("/home/jason/R/runout.opt/")

library(devtools)
document()
build()
install()


# Load Packages and Data ####################################################################
library(runout.opt)
library(raster)
library(rgdal)
library(sp)
library(Rsagacmd)

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
runout_polygon <- runout_polygons[77,]

sel_source_point  <- over(source_points, runout_polygon)
source_point <- source_points[!is.na(sel_source_point$objectid),]


dem_grid <- raster::crop(dem, raster::extent(runout_polygon) + 1000)
source_grid <- raster::rasterize(matrix(sp::coordinates(source_point)[1:2], ncol = 2), dem_grid, field = 1)

s_time <- Sys.time()
gpp <- saga$sim_geomorphology$gravitational_process_path_model(dem = dem_grid, release_areas = source_grid,
                                                               #friction_mu_grid = mu.grid,
                                                               process_path_model = 1,
                                                               rw_slope_thres = 40,
                                                               rw_exponent = 1.9,
                                                               rw_persistence = 2,
                                                               gpp_iterations = 1000,
                                                               friction_model = 5,
                                                               friction_mu = 0.11,
                                                               friction_mass_to_drag = 40,
                                                               process_area = "test_parea.sgrd",
                                                               .all_outputs = FALSE)
                                                               # max_velocity = "test_max_vel")
Sys.time() - s_time


river <- readOGR("Custom_DrainageLine.shp")
river_buff <- rgeos::gBuffer(river, width = 30)

# Try to extract only the parea freq values at border of river polygon...
river_buff <- as(river_buff, 'SpatialLines')
plot(gpp$process_area)
plot(river_buff, add = TRUE)
river_inter <- unlist(extract(gpp$process_area, river_buff))
river_inter/1000
sum(river_inter/1000, na.rm=TRUE)

# try with sf
library(sf)

# Convert to sf
sf_river <- st_as_sf(river)
sf_riverbf <- st_as_sf(river_buff)

sf_crop <- st_crop(sf_riverbf, extent(dem_grid))

# Splite
sf_riversplit <- st_difference(sf_riverbf, st_buffer(st_intersection(sf_riverbf, sf_river), dist=1e-12))

river_split <- as_Spatial(sf_riversplit)

plot(gpp$process_area)
plot(as_Spatial(sf_riversplit), add=TRUE)

