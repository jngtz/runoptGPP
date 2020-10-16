setwd("/home/jason/R/runout.opt/")

library(devtools)
build()
install()

# Load Packages and Data ####################################################################
library(runout.opt)
library(raster)
library(rgdal)
library(Rsagacmd)

# Initiate a SAGA-GIS geoprocessing object
saga <- saga_gis(opt_lib = "sim_geomorphology")

# Data
setwd("/home/jason/Data/Chile/")
dem <- raster("elev_alos_12_5m.tif")

# Runout source points
source_points <- readOGR(".", "debris_flow_source_points")

# Runout track polygons
runout_polygons <- readOGR(".", "debris_flow_polys_sample")
# Assign an object ID to each row of the SpatialPolygonsDataFrame
runout_polygons$objectid <- 1:100


# Example of GPP random walk simulation using R ################################

rwPerformance(dem, slide_plys = runout_polygons, slide_src = source_points,
                   slide_id = 77, slp = 30, ex = 3, per = 2,
                   gpp_iter = 1000, buffer_ext = 500, buffer_source = 50,
                   plot_eval = TRUE)

# Define grid search values
steps <- 3
rwexp_vec <- seq(1.3, 3, len=steps)
rwper_vec <- seq(1.5, 2, len=steps)
rwslp_vec <- seq(20, 40, len=steps)

rw_gridsearch <- rwGridsearch(dem, slide_plys = runout_polygons, slide_src = source_points,
               slide_id = 77,
               #Input random walk grid search space
               slp_v = rwslp_vec, ex_v = rwexp_vec, per_v = rwper_vec,
               #Set number of simulation iterations
               gpp_iter = 1000,
               #Define processing extent size (m)
               buffer_ext = 500,
               #(Optional) Define size of buffer to make source area from point
               buffer_source = 50)

rw_gridsearch

rw_opt_single <- rwGetOpt_single(rw_gridsearch)
rw_opt_single

# Run RW optimization for multiple runout tracks  ##############################

## DO NOT RUN
polyid_vec <- 1:100

# Use parallel processing for faster computations
library(foreach)
cl <- parallel::makeCluster(4)
doParallel::registerDoParallel(cl)

rw_grisearch_multi <-
  foreach(poly_id=polyid_vec, .packages=c('rgdal','raster', 'rgeos', 'ROCR', 'Rsagacmd', 'sf', 'runout.opt')) %dopar% {

    .GlobalEnv$saga <- saga

    rwGridsearch(dem, slide_plys = runout_polygons, slide_src = source_points,
                   slide_id = poly_id, slp_v = rwslp_vec, ex_v = rwexp_vec, per_v = rwper_vec,
                   gpp_iter = 1000, buffer_ext = 500, buffer_source = 50, save_res = FALSE,
                   plot_eval = FALSE)

  }

parallel::stopCluster(cl)


# Get optimal random walk parameter set ########################################

## HIDE
setwd("/home/jason/Scratch/GPP_RW_Paper")
(load("rw_gridsearch_multi.Rd"))
##

# from object
rwGetOpt(rw_gridsearch_multi, measure = median)


# Validate transferability using spatial cross validation ######################

rw_spcv <- rwSPCV(x = rw_gridsearch_multi, slide_plys = runout_polygons,
                  n_folds = 5, repetitions = 10)

freq_rw <- rwPoolSPCV(rw_spcv, plot_freq = TRUE)


# Example of GPP PCM simulation #######################

pcm <- pcmPerformance(dem, slide_plys = runout_polygons, slide_src = source_points,
               slide_id = 77, rw_slp = 40, rw_ex = 3, rw_per = 1.5,
               pcm_mu = 0.15, pcm_md = 120,
               gpp_iter = 1000, buffer_ext = 500, buffer_source = 50,
               plot_eval = TRUE, return_features = TRUE)

# Runout distance relative error
pcm$length.relerr

# Plot GPP PCM runout modelling ouputs
gpp_output <- stack(pcm$gpp.parea, pcm$gpp.stop, pcm$gpp.maxvel)
names(gpp_output) <- c("Process_area", "Stop_positions", "Max_velocity")
plot(gpp_output)

# Grid search for optimal PCM parameter set ###############
pcmmd_vec <- seq(20, 120, by=20)
pcmmu_vec <- seq(0.05, 0.3, by=0.1)

pcm_gridsearch <- pcmGridsearch(dem,
                       slide_plys = runout_polygons, slide_src = source_points, slide_id = 77,
                       #Plug-in random walk optimal parameters
                       rw_slp = rw_opt_single$rw_slp_opt,
                       rw_ex = rw_opt_single$rw_exp_opt,
                       rw_per = rw_opt_single$rw_per_opt,
                       #Input PCM grid search space
                       pcm_mu_v = pcmmu_vec,
                       pcm_md_v = pcmmd_vec,
                       #Set number of simulation iterations
                       gpp_iter = 1000,
                       #Define processing extent size (m)
                       buffer_ext = 500,
                       #(Optional) Define size of buffer to make source area from point
                       buffer_source = 50)

pcmGetOpt_single(pcm_gridsearch)


# Apply grid search to multiple events #########

polyid_vec <- 1:4

library(foreach)
cl <- parallel::makeCluster(4)
doParallel::registerDoParallel(cl)

pcm_gridsearch_multi <-
  foreach(poly_id=polyid_vec, .packages=c('rgdal','raster', 'rgeos', 'ROCR', 'Rsagacmd', 'sf', 'runout.opt')) %dopar% {

    .GlobalEnv$saga <- saga

    pcmGridsearch(dem,
                  slide_plys = runout_polygons, slide_src = source_points, slide_id = poly_id,
                  rw_slp = rw_opt$rw_slp_opt, rw_ex = rw_opt$rw_exp_opt, rw_per = rw_opt$rw_per_opt,
                  pcm_mu_v = pcmmu_vec, pcm_md_v = pcmmd_vec,
                  gpp_iter = 1000,
                  buffer_ext = 500, buffer_source = NULL,
                  predict_threshold = 0.5,
                  plot_eval = FALSE)

  }

parallel::stopCluster(cl)

# GET PCM OPTIMAL PARAMETERS #######################################

pcmGetOpt(pcm_gridsearch_multi, performance = "relerr", measure = "median", plot_opt = TRUE)

library(ggplot2)
library(metR)
library(reshape2)

pcm_grid <- pcmGetGrid(pcm_gridsearch_multi, performance = "auroc", measure = "IQR")
pcm_grid_df <- melt(pcm_grid)

ggplot(data = pcm_grid_df, aes(x=Var2, y=Var1, z=value)) +
  geom_tile(aes(fill = value)) +
  ylab(expression(paste("Sliding friction coefficient"))) +
  xlab("Mass-to-drag ratio (m)") +
  labs(fill="Median relative\nrunout distance\nerror") +
  scale_fill_viridis_c(direction = 1) +
  theme_light() +
  theme(text = element_text(family = "Arial", size = 8), axis.title = element_text(size = 9),
        axis.text = element_text(size = 8),legend.position = "right")


# PCM spatial cross validation #######################################

pcm_spcv <- pcmSPCV(pcm_gridsearch_multi, slide_plys = runout_polygons,
                    n_folds = 5, repetitions = 100, from_save = FALSE)

freq_pcm <- pcmPoolSPCV(pcm_spcv, plot_freq = TRUE)




# STOP Vignette here ... Determine source area prediction threshold #############################
setwd("/home/jason/Data/Chile/")
source_pred <- raster::raster("source_area_prediction.tif")

cutoffs <- seq(.5,.95, by=0.05)

library(foreach)
cl <- parallel::makeCluster(4)
doParallel::registerDoParallel(cl)

gpp_pareas <-
  foreach(src_thrsh=cutoffs, .packages=c('rgdal','raster', 'rgeos', 'ROCR', 'Rsagacmd', 'sf', 'runout.opt')) %dopar% {

    .GlobalEnv$saga <- saga

    runoutPareaPredict(source_pred, dem, source_threshold = src_thrsh,
                     rw_slp, rw_exp, rw_per, pcm_mu, pcm_md,
                     gpp_iter = 1000)

  }

parallel::stopCluster(cl)


# Compute AUROC for each process area
parea_aurocs<- rep(NA, length(cutoffs))

for(i in 1:length(cutoffs)){
  parea_aurocs[i] <- rocParea(gpp_pareas[[i]], runout_polygons)
}

# Sample size analysis #################################################

