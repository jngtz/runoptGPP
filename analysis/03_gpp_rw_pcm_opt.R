
# LOAD PACKAGES AND DATA #######################################################

library(runoptGPP)
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


# DEFINE RW AND PCM GRID SEARCH SPACE ##########################################

steps <- 11
rwexp_vec <- seq(1.3, 3, len=steps)
rwper_vec <- seq(1.5, 2, len=steps)
rwslp_vec <- seq(20, 40, len=steps)

pcmmd_vec <- seq(20, 150, by=5)
pcmmu_vec <- seq(0.04, 0.6, by=0.01)

polyid_vec <- 1:100


# RW GRIDSEARCH OPTIMIZATION W PARALLELIZATION  ################################
setwd("/home/jason/Scratch/GPP_RW_Paper")

library(foreach)
cl <- parallel::makeCluster(4)
doParallel::registerDoParallel(cl)

rw_grisearch_multi <-
  foreach(poly_id=polyid_vec, .packages=c('rgdal','raster', 'rgeos', 'ROCR', 'Rsagacmd', 'sf', 'runout.opt')) %dopar% {

    .GlobalEnv$saga <- saga

    rwGridsearch(dem, slide_plys = runout_polygons, slide_src = source_points,
                 slide_id = poly_id, slp_v = rwslp_vec, ex_v = rwexp_vec, per_v = rwper_vec,
                 gpp_iter = 1000, buffer_ext = 500, buffer_source = 50, save_res = TRUE,
                 plot_eval = FALSE)

  }

parallel::stopCluster(cl)

# Get RW optimal parameter set
rw_opt <- rwGetOpt(rw_gridsearch_multi, measure = median)
save(rw_opt, file = "rw_opt_params.Rd")
save(rw_gridsearch_multi, file = "rw_gridsearch_multi.Rd")


# RW PARAM VALIDATION W SPATIAL CV #############################################

rw_spcv <- rwSPCV(x = rw_gridsearch_multi, slide_plys = runout_polygons,
                  n_folds = 5, repetitions = 1000)

freq_rw <- rwPoolSPCV(rw_spcv, plot_freq = TRUE)
freq_rw


# PCM GRIDSEARCH OPTIMIZATION W PARALLELIZATION  ################################
setwd("/home/jason/Scratch/GPP_PCM_Paper")

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
                  predict_threshold = 0.5, save_res = TRUE,
                  plot_eval = FALSE)

  }

parallel::stopCluster(cl)

# Get PCM optimal parameter set
pcm_opt <- pcmGetOpt(pcm_gridsearch_multi, performance = "relerr", measure = "median", plot_opt = TRUE)
save(pcm_opt, file = "pcm_opt_params.Rd")
save(pcm_gridsearch_multi, file = "pcm_gridsearch_multi.Rd")


# PCM PARAM VALIDATION W SPATIAL CV #############################################

pcm_spcv <- pcmSPCV(pcm_gridsearch_multi, slide_plys = runout_polygons,
                    n_folds = 5, repetitions = 1000, from_save = FALSE)

freq_pcm <- pcmPoolSPCV(pcm_spcv, plot_freq = TRUE)
freq_pcm


