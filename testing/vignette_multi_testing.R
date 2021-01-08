setwd("/home/jason/R/runoutGPP/")

library(devtools)
document()
build()
install()


# Load Packages and Data ####################################################################
library(runout.opt)
library(raster)
library(rgdal)
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


# Run RW optimization for multiple runout tracks  ##############################

polyid.vec <- 1:100
steps <- 11
rwexp_vec <- seq(1.3, 3, len=steps)
rwper_vec <- seq(1.5, 2, len=steps)
rwslp_vec <- seq(20, 40, len=steps)

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


# Apply grid search to multiple events #########

pcmmd_vec <- seq(20, 150, by=5)
pcmmu_vec <- seq(0.04, 0.6, by=0.01)

polyid_vec <- 1:100

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

## HIDE
setwd("/home/jason/Scratch/GPP_PCM_Paper")
(load("pcm_gridsearch_multi.Rd"))
##


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


