# LOAD PACKAGES ################################################################

library(rgdal)
library(raster)
library(Rsagacmd)

# Initiate a SAGA-GIS geoprocessing object
saga <- saga_gis()

# LOAD DATA ####################################################################

setwd("/home/jason/Data/Chile/")
# elevation model
dem <- raster("elev_alos_12_5m_no_sinks.tif")

# slide start/source points
slide_point_vec <- readOGR(".", "debris_flow_source_points")

# actual/mapped debris flow polygons
slide_poly_vec <- readOGR(".", "debris_flow_polys_sample")
slide_poly_vec$objectid <- 1:100

crs(slide_point_vec) <- crs(slide_poly_vec)

# FIND MODULE FOR LOADING INTO SAGA ###############

saga$io_gdal$import_raster(files = "elev_alos_12_5m_no_sinks.tif")

saga$
