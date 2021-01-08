
# LOAD PACKAGES AND DATA #######################################################

library(rgdal)
library(raster)
library(sp)
library(runoptGPP)
library(ggplot2)

setwd("/home/jason/Data/Chile/")

#name of landslide point shapefile
nm_slidepoint <- "debris_flow_source_points"

#name of landslide polygon shapefile
nm_slidepoly <- "debris_flow_runout_polys"

#load slide points and polygons
slide_points <- readOGR(".", nm_slidepoint)

#Use only debris flows
slide_polys <- readOGR(".", nm_slidepoly)
slide_polys$value <- 1

# SUMMARIZE DEBRIS FLOWS ##########################

geom_slide <- runoutGeom(slide_polys)
#save(geom_slide, file = "geom_slide_polys.Rd")
(load("geom_slide_polys.Rd"))

summary(geom_slide)
boxplot(geom_slide$surfacearea)

sapply(geom_slide, FUN=median)
sapply(geom_slide, FUN=IQR)
sapply(geom_slide, FUN=min)
sapply(geom_slide, FUN=max)
