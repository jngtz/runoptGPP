library(raster)
library(rgdal)

setwd("/home/jason/Data/Chile/")

# Load digital elevation model (DEM)
dem <- raster("elev_alos_12_5m.tif")
slp <- terrain(dem, opt="slope", unit="degrees")

# Load runout track polygons and assign object ID based on row number
runout_polygons <- readOGR(".", "debris_flow_polys_sample")
runout_polygons$objectid <- 1:length(runout_polygons)

# Select a debris flow and source point for this example
ply_deposit <- runout_polygons[77,]

slp_deposit <- raster::extract(slp, ply_deposit)[[1]]
dem_deposit <- raster::extract(dem, ply_deposit)[[1]]

# Compute friction angle from point of start of deposition
theta <- slp_deposit[dem_deposit == max(dem_deposit)]
mu = tan(theta*pi/180)
mu

# Compute friction angle from average slope in deposition area
theta_avg <- mean(slp_deposit)

# Convert to friction coefficient mu
mu_avg = tan(theta_avg*pi/180)
mu_avg
