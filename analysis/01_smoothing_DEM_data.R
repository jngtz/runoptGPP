setwd("/home/jason/Data/Chile/")

# FUNCTIONS ###################################################################
geoTiffToSGRD_custom <- function(file_name, write_name){
  #Convert a tif to a SAGA format
  require(rgdal)
  setwd("../")
  spgdf <- readGDAL(file_name)
  setwd("saga_data")
  writeGDAL(spgdf, write_name, drivername = "SAGA")
}


# RSAGA SETTINGS ###############################################################

library(RSAGA)

myenv <- rsaga.env()

# IMPORT DEM ################################################################
setwd("saga_data")

#Convert gtiff to sgrd for processing in SAGA
geoTiffToSGRD_custom("elev_alos_12_5m.tif",  "dem.sdat")



# FILTER (SAR) NOISE ###########################################################

# https://www.sciencedirect.com/science/article/pii/S0169555X09002967?via%3Dihub

#rsaga.get.modules("grid_filter")

rsaga.get.usage("grid_filter", "Mesh Denoise") # Sun et al 2007
rsaga.geoprocessor("grid_filter", "Mesh Denoise", env=myenv,
                    param=list(
                      INPUT = "dem.sgrd",
                      OUTPUT = "dem_filtered.sgrd",
                      SIGMA = 0.9, # Threshold
                      ITER = 70, # Number of iterations for Normal Updating
                      VITER = 100, # Number of iterations for Vertex Updating
                      NB_CV = 0)) # Common vertex



# FILL SINKS ###################################################################
#rsaga.get.modules("ta_preprocessor", env=myenv)
#rsaga.get.usage("ta_preprocessor", "Fill Sinks (Planchon/Darboux, 2001)", env=myenv)

rsaga.geoprocessor("ta_preprocessor", "Fill Sinks (Planchon/Darboux, 2001)",
                   param=list(
                     DEM = "dem_filtered.sgrd",
                     RESULT = "dem_filtered_filled.sgrd"),
                   env=myenv)


out_dem <- raster("dem_filtered_filled.sgrd")
setwd("../")
writeRaster(out_dem, filename = "dem_filtered_filled.tiff", format="GTiff", overwrite=TRUE)


setwd("saga_data")


# COMPUTE TERRAIN ATTRIBUTES ###################################################

# slope
# aspect
# general curvature
# plan curvature
# profile curvature
# catchment slope
# catchment area (flow accumulation) (modified - saga wetness index)
# saga wetness index


# slope, aspect, curvature (general, plan and profile)

#rsaga.get.modules("ta_morphometry", env=myenv)
#rsaga.get.usage("ta_morphometry", "Slope, Aspect, Curvature", env=myenv)

rsaga.geoprocessor("ta_morphometry", "Slope, Aspect, Curvature", env=myenv,
                   param = list(
                     ELEVATION = "dem_filtered_filled.sgrd",
                     SLOPE = "slope_filtered_filled.sgrd",
                     ASPECT = "aspect_filtered_filled.sgrd",
                     C_GENE = "gencrv_filtered_filled.sgrd",
                     C_PROF = "profcrv_filtered_filled.sgrd",
                     C_PLAN = "plancrv_filtered_filled.sgrd",
                     UNIT_SLOPE = 1, #degrees
                     UNIT_ASPECT = 1 #degrees
                   ))


# catchment area/contributing upslope area/flow accumulation [m^2]

#rsaga.get.modules("ta_hydrology", env=myenv)
#rsaga.get.usage("ta_hydrology", "SAGA Wetness Index", env=myenv)

rsaga.geoprocessor("ta_hydrology", "SAGA Wetness Index", env=myenv,
                   param = list(
                     DEM = "dem_filtered_filled.sgrd",
                     AREA = "carea_filtered_filled.sgrd",
                     SLOPE = "cslope_filtered_filled.sgrd",
                     AREA_MOD = "modcarea_filtered_filled.sgrd",
                     TWI = "swi_filtered_filled.sgrd"
                   ))


# hillshade

#rsaga.get.modules("ta_lighting", env=myenv)
#rsaga.get.usage("ta_lighting", "Analytical Hillshading", env=myenv)

rsaga.geoprocessor("ta_lighting", "Analytical Hillshading", env=myenv,
                   param = list(
                     ELEVATION = "dem_filtered_filled.sgrd",
                     SHADE = "hs_filtered_filled.sgrd"
                   ))


# terrain surface convexity (tsconv)
#rsaga.get.modules("ta_morphometry", env=myenv)
#rsaga.get.usage("ta_morphometry", "Terrain Surface Convexity", env=myenv)

rsaga.geoprocessor("ta_morphometry", "Terrain Surface Convexity", env=myenv,
                   param = list(
                     DEM = "dem_filtered_filled.sgrd",
                     CONVEXITY = "tsconv_filtered_filled.sgrd"
                   ))

# terrain rugedness index (tri)
#rsaga.get.modules("ta_morphometry", env=myenv)
#rsaga.get.usage("ta_morphometry", "Terrain Ruggedness Index (TRI)", env=myenv)

rsaga.geoprocessor("ta_morphometry", "Terrain Ruggedness Index (TRI)", env=myenv,
                   param = list(
                     DEM = "dem_filtered_filled.sgrd",
                     TRI = "tri_filtered_filled.sgrd"
                   ))

