setwd("/home/jason/R/runoptGPP/")

#library(devtools)
#document()
#build()
#install()

#check(vignettes = FALSE)


# Load Packages and Data ####################################################################
library(runoptGPP)
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


# GPP random walk simulation ################################

rwPerformance(dem, slide_plys = runout_polygon, slide_src = source_point,
                   slp = 30, ex = 3, per = 2,
                   gpp_iter = 1000, buffer_ext = 500, buffer_source = 50,
                   plot_eval = TRUE, saga_lib = saga)


crop_dem <- crop(dem, extent(runout_polygon)*2)
plot(crop_dem)
plot(runout_polygon, add = TRUE)
plot(source_point, add = TRUE)

# Custom random walk ########################################
#https://www.r-bloggers.com/2018/05/simulating-animal-movements-and-habitat-use/

# Functions #############




slp_th <- 30 # slope threshold

# get cell res
res_cell <- res(dem)[1]
# get start cell
str_cell <- cellFromXY(dem, xy = c(source_point$POINT_X, source_point$POINT_Y))
ngh_cells <- adjacent(dem, str_cell, directions = 8, pairs = FALSE)

elv_ngh <- extract(dem, ngh_cells)
elv_str <- extract(dem, str_cell)

slp_cells <-  tan( (elv_ngh - elv_str)  / res_cell) * 180/pi


df <- data.frame(cell = ngh_cells,
                     elv = elv_ngh,
                     slp = slp_cells)

df$gamma <- tan(df$slp) / tan(slp_th)

# pick neighbors that are below slope threshold


