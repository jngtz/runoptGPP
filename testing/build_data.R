setwd("inst/extdata")
dflow_runout_ply <- rgdal::readOGR("dflow_runout_ply.shp")
dflow_source_pnt <- rgdal::readOGR("dflow_source_pnt.shp")
elev_model <- raster::raster("elev_12_5m.tif")

usethis::use_data(dflow_runout_ply)
usethis::use_data(dflow_source_pnt)
usethis::use_data(elev_model)
