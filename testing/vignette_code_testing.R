setwd("/home/jason/R/runout.opt/")


# Load Data ####################################################################
setwd("/home/jason/Data/Chile/")
dem <- raster::raster("dem_alos_12_5m _no sinks.tif")

# slide start/source points
source_points <- rgdal::readOGR(".", "dflow_points_v1_reposition")

# actual/mapped debris flow polygons
runout_polygons<- rgdal::readOGR(".", "dflow_polygons_v1_reposition_sample_100")
runout_polygons$objectid <- 1:100

raster::crs(source_points) <- raster::crs(runout_polygons)

# initalize SAGA GIS geoprocessing
saga <- Rsagacmd::saga_gis(opt_lib = "sim_geomorphology")


# Example of GPP random walk simulation using R ################################

library(runout.opt)

rwPerformance(dem, slide_plys = runout_polygons, source_pnts = source_points,
                   slide_id = 10, slp = 33, ex = 3, per = 2,
                   gpp_iter = 1000, buffer_ext = 500, buffer_source = NULL,
                   plot_eval = TRUE)

steps <- 3
rwexp_vec <- seq(1.3, 3, len=steps)
rwper_vec <- seq(1.5, 2, len=steps)
rwslp_vec <- seq(20, 40, len=steps)

rw_gridsearch <- rwGridsearch(dem, slide_plys = runout_polygons, source_pnts = source_points,
               slide_id = 2, slp_v = rwslp_vec, ex_v = rwexp_vec, per_v = rwper_vec,
               gpp_iter = 1000, buffer_ext = 500, buffer_source = 50,
               save_res = FALSE, plot_eval = TRUE)

rw_gridsearch


rwGetOpt_single(rw_gridsearch)


# Run RW optimization for multiple runout tracks  ##############################

# Use parallel processing for faster computations

library(foreach)
cl <- parallel::makeCluster(4)
doParallel::registerDoParallel(cl)

Sys.time()
rw_grisearch_multi <-
  foreach(poly_id=polyid_vec, .packages=c('rgdal','raster', 'rgeos', 'ROCR', 'Rsagacmd', 'sf', 'runout.opt')) %dopar% {

    .GlobalEnv$saga <- saga

    rwGridsearch(dem, slide_plys = runout_polygons, source_pnts = source_points,
                   slide_id = poly_id, slp_v = rwslp_vec, ex_v = rwexp_vec, per_v = rwper_vec,
                   gpp_iter = 1000, buffer_ext = 500, buffer_source = 50, save_res = FALSE,
                   plot_eval = FALSE)

  }

parallel::stopCluster(cl)
Sys.time()


# Get optimal random walk parameter set ########################################

# from object
rwGetOpt(rw_gridsearch_multi, measure = median)

# from saved files
setwd("/home/jason/Scratch/GPP_RW_Paper")
rwGetOpt(measure = median, from_save = TRUE)



# Validate transferability using spatial cross validation ######################


rw_spcv <- rwSPCV(x = rw_gridsearch_multi, slide_plys = runout_polygons[1:4,],
       n_folds = 3, repetitions = 10)


setwd("/home/jason/Scratch/GPP_RW_Paper")

rw_spcv <- rwSPCV(slide_plys = runout_polygons,
            n_folds = 5,
            repetitions = 10,
            from_save = TRUE)

freq_rw <- rwPoolSPCV(rw_spcv, plot.freq = TRUE)


# Visualize freq of SPCV optimal parameter sets ################################
setwd("/home/jason/Scratch/GPP_RW_Paper")
# Load saga gpp random walk settings
load("gridsearch_rw_settings.Rd")

rwexp_vec <- rw_settings$vec_rwexp
rwper_vec <- rw_settings$vec_rwper
rwslp_vec <- rw_settings$vec_rwslp


library(ggplot2)

#Pick a slope threshold (slice) of grid search space
slope_thresh <- 40
freq_rw <- freq_rw[order(freq_rw$rel_freq, decreasing = TRUE),]

library(ggplot2)
breaks_bubble <- c(10, 30, 60)

gg_freq_rw <- freq_rw[freq_rw$slp == slope_thresh,]

ggplot(gg_freq_rw, aes(x=per, y=exp)) +
  ggtitle(paste("Slope threshold:", slope_thresh ))+
  # Can improve by making different colors for different slope thresholds...

  geom_point(alpha=0.7, aes(colour = median_auroc, size = rel_freq)) +
  scale_size(name="Relative\nfrequency (%)",
            breaks = breaks_bubble) +
  scale_colour_gradient(low = "#1B4F72", high = "#85C1E9",
                        name = "Median AUROC") +
  scale_x_continuous(expression(paste("Persistence factor")),
                     limits = c(min(rwper_vec), max = max(rwper_vec))) +
  scale_y_continuous(expression(paste("Exponent of divergence")),
                     limits = c(min(rwexp_vec), max = max(rwexp_vec)+.1)) +
  theme_light() +
  theme(text = element_text(family = "Arial", size = 8), axis.title = element_text(size = 9),
        axis.text = element_text(size = 8), title = element_text(size = 8))


# PCM MODEL OPTIMIZATION #######################
# Adapt pcm results similar to RW to work with object or saved files


workspace_dir <- "/home/jason/Scratch/F_Testing"
setwd(workspace_dir)

result <- pcmPerformance(dem, slide_plys = runout_polygons, source_pnts = source_points,
               slide_id = 8, rw_slp = 33, rw_ex = 3, rw_per = 2,
               pcm_mu = 0.11, pcm_md = 20,
               gpp_iter = 1000, buffer_ext = 500, buffer_source = 50,
               plot_eval = TRUE, return_features = TRUE)

pcmmd_vec <- seq(50, 150, by=50)
pcmmu_vec <- seq(0.04, 0.6, by=0.25)

pcm_result <- pcmGridsearch(dem, workspace = workspace_dir,
                       slide_plys = runout_polygons, source_pnts = source_points, slide_id = 5,
                       rw_slp = 40, rw_ex = 1.8, rw_per = 2,
                       pcm_mu_v = pcmmu_vec, pcm_md_v = pcmmd_vec,
                       gpp_iter = 1000,
                       buffer_ext = 500, buffer_source = NULL,
                       predict_threshold = 0.5,
                       plot_eval = FALSE)

for(i in 1:10){
  pcmGridsearch(dem, workspace = workspace_dir,
             slide_plys = runout_polygons, source_pnts = source_points, slide_id = i,
             rw_slp = 40, rw_ex = 1.8, rw_per = 2,
             pcm_mu_v = pcmmu_vec, pcm_md_v = pcmmd_vec,
             gpp_iter = 1000,
             buffer_ext = 500, buffer_source = NULL,
             predict_threshold = 0.5,
             plot_eval = FALSE)
}

# GET PCM OPTIMAL PARAMETERS #######################################

setwd(workspace_dir)

pcmGetOpt(pcm_md_vec = pcmmd_vec, pcm_mu_vec = pcmmu_vec, n_train = 10,
          performance = "relerr", measure = "median")

pcm_spcv <- pcmSPCV(slide_plys = runout_polygons[1:10,],
                    n_folds = 3, repetitions = 100, pcm_mu_v = pcmmu.vec,
                    pcm_md_v = pcmmd.vec)

# Test pooling
setwd("/home/jason/Scratch/GPP_PCM_Paper")

(load("gridsearch_pcm_settings.Rd"))

polyid.vec <- 1:pcm_settings$n_train
pcmmu.vec <- pcm_settings$vec_pcmmu
pcmmd.vec <- pcm_settings$vec_pcmmd



pcm_spcv <- pcmSPCV(slide_plys = runout_polygons[1:100,],
                    n_folds = 3, repetitions = 10, pcm_mu_v = pcmmu.vec,
                    pcm_md_v = pcmmd.vec)

(load("repeated_spcv_PCM.Rd"))

freq_pcm <- pcmPoolSPCV(pcm_spcv)


