# LOAD PACKAGES AND DATA ####################################################################

library(rgdal)
library(raster)
library(sp)
library(mgcv)
library(sperrorest)

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

#Note Filtered result provides a much more homogenous source zones, but
# does not work well for runout model.. -

setwd("saga_data/meshdenoise")
#load raster files
ta_list <- list("dem_meshdenoise_filled.sdat",
                "slope_meshdenoise_filled.sdat",
                "aspect_meshdenoise_filled.sdat",
                "gencrv_meshdenoise_filled.sdat",
                "profcrv_meshdenoise_filled.sdat",
                "plancrv_meshdenoise_filled.sdat",
                "carea_meshdenoise_filled.sdat",
                "cslope_meshdenoise_filled.sdat",
                "modcarea_meshdenoise_filled.sdat",
                "swi_meshdenoise_filled.sdat",
                "tsconv_meshdenoise_filled.sdat",
                "tri_meshdenoise_filled.sdat")

ta_names <- c("elev", "slope", "aspect", "gen_curv", "prof_curv", "plan_curv",
              "carea", "cslope", "mod_carea", "swi", "convexity", "tri")

layers <- stack(ta_list)
names(layers) <- ta_names

setwd("../../")

layers[["cslope"]] <- layers[["cslope"]] * 180/pi
layers[["logcarea"]] <- log10(layers[["carea"]])
layers[["dist_faults"]] <- raster("dist_faults.tif")



# MASK AREAS OUTSIDE OF LANDSLIDE MAPPING AREAS ################################

# create mask layer for sampling
area_mask <- mask(layers[['elev']], layers[['elev']], inverse=TRUE, updatevalue=1)
slide_mask <- mask(area_mask, slide_polys, inverse=TRUE)

# rasterize selected subcatchments
sel_subcatchments <- rasterize(sub_catchment_poly, dem, field = "Seleccion")
sub_catch_vec = readOGR("sub_catchments.shp")

# create mask for sampling non-landslide areas
mask <- mask(slide_mask, sel_subcatchments, maskvalue=1)
writeRaster(mask, "mask_dflow_repo.tif", format="GTiff", overwrite=TRUE)


# SPATIAL SAMPLE OF LANDSLIDE AND NON-LANDSLIDE AREAS ##########################

# to provide data for our spatially predictive models we need to create
# sample of the response variable: landslide and non-landslides

# sample of landslide points
smp.slides <- slide_points
smp.slides$x <- coordinates(smp.slides)[,1]
smp.slides$y <- coordinates(smp.slides)[,2]
xy.slides <- data.frame(x = coordinates(smp.slides)[,1], y=coordinates(smp.slides)[,2])

# randomly sample cell locations of non-landslide cells
set.seed(1234)
smp.size <- nrow(xy.slides)

smp.noslides <- as.data.frame(sampleRandom(mask, size = smp.size, xy = TRUE))
smp.noslides[,3] = NULL
points(smp.noslides)

# extract predictor variable values from grids
# for the landslide samples
df.slides <- extract(layers, xy.slides, sp = TRUE)
df.slides <- as.data.frame(df.slides)
# add xy values
df.slides$x <- smp.slides$x
df.slides$y <- smp.slides$y
# add variable indicating that these are landslides
df.slides$landslide <- "TRUE"
# remove landslides outside test area (area of interest)
df.slides <- df.slides[!is.na(df.slides$elev),]

# for the non-landslide sample
df.noslides <- as.data.frame(extract(layers, smp.noslides, sp = TRUE))
df.noslides$x <- smp.noslides$x
df.noslides$y <- smp.noslides$y
df.noslides$landslide <- "FALSE"

# combine the slides and no slides data into one dataframe
df <- rbind(df.slides, df.noslides)

# encode the landslide variable as a factor
df$landslide <- as.factor(df$landslide)

# Remove NA related to land-use
df <- df[!is.na(df$landuse),]
df$landuse <- as.factor(df$landuse)


# BUILD SOURCE AREA PREDICTION MODELS ##########################################

# Transform predictor variables
df$logcarea <- log10(df$carea)
df$plan_curv[df$plan_curv > .1] <- .1
df$plan_curv[df$plan_curv < -.1] <- -.1
df$logcarea <- log10(df$carea)

# Create formula
fml <- landslide ~ elev + slope + logcarea + plan_curv + dist_faults

# fit a logistic regression model
model.lr <- glm(formula = fml, family = 'binomial', data = df)
summary(model.lr)

# to automatically fit a GAM with smoothing functions based on the formula used
# for the GLM ("fml"), we use the following function
my.gam <- function(formula, data, family = binomial, k = 5) {
  response.name <- as.character(formula)[2]
  predictor.names <- labels(terms(formula))
  categorical <- sapply(data[,predictor.names], is.logical) |
    sapply(data[,predictor.names], is.factor)
  formula <- paste(response.name, "~",
                   paste(predictor.names[categorical], "+"),
                   paste("s(", predictor.names[!categorical], ", k=", k, ")", collapse = "+"))
  formula <- as.formula(formula)
  fit <- gam(formula, data, family = family, select=TRUE)
  return(fit)
}


# fit a generalized additive model (GAM)
model.gam <- my.gam(fml, df)
# display gam model formula
formula(model.gam)
# visualize smoothing functions in the GAM
par(mfrow = c(2,3), mar = c(4, 4, 1, 1))
plot(model.gam, all.terms = TRUE)



# MODEL VALIDATION AND COMPARISON ##############################################

# to make certain that performance of our model is not by chance,
# we can perform repeated cross-validation, or in this example,
# repated cross-validation based on splitting the data into
# spatially defined subsets (i.e., spatial cross-validation)

# perform 1000 repeated, 5-fold spatial cv for the generalized additive model

gam.sp.results <- sperrorest(formula = fml, data = df, coords = c("x","y"),
                          model_fun = my.gam,
                          pred_args = list(type="response"),
                          smp_fun = partition_kmeans,
                          smp_args = list(repetition = 1:1000, nfold = 5))

gam.sp.auroc.test <- unlist(summary(gam.sp.results$error_rep,level=1)[,"test_auroc"])
mean(gam.sp.auroc.test)
sd(gam.sp.auroc.test)

median(gam.sp.auroc.test)
IQR(gam.sp.auroc.test)

# and resampling for spatial cross-validation using k-means
resamp.spcv = partition_kmeans(df, nfold=5, repetition=1)
plot(resamp.spcv, df)



# APPLY SPATIAL PREDICTIONS TO GRIDDED DATA ##################################

# now that we have validated our model peformance, we can create
# a spatial model prediction by applying the fitted model predictions
# to a gridded dataset (our RasterStack object)

#Since there is a bias in sample: high elevaiton areas were not sampled,
#a constant  value for elevation will be used in the prediction to grid
# the mean  value from the sample mean(df$elev) was used (2821 masl)

# Transformations
layers[['plan_curv']][layers[['plan_curv']] > 0.1] <- 0.1
layers[['plan_curv']][layers[['plan_curv']] < -0.1] <- -0.1


#layers[['elev']] <- mask(layers[['elev']], layers[['elev']], inverse=TRUE, updatevalue=2821)

#NOTE takes a while to process... look at spatial.tools package for improving speed
#pred.gam<- raster::predict(layers, model.gam, type = "response", progress = TRUE)

#Make it quicker (easier on memory) by removing unecessary layers from layers
predictor.names <- labels(terms(fml))
newdata.layers <- layers[[predictor.names]]

library(parallel)
beginCluster(15) # 5 for my notebook b/c only enough RAM for 4
clusterGAM.pred <- clusterR(newdata.layers, predict, args=list(model=model.gam, type="response"))
endCluster()

# export to a raster format
writeRaster(clusterGAM.pred, "src_area_pred.tif", format = "GTiff", overwrite = TRUE)

