setwd("/home/jason/R/runout.opt/")

setwd("/home/jason/Data/Chile/")
# elevation model
dem <- raster::raster("dem_alos_12_5m _no sinks.tif")

# slide start/source points
slide_point_vec <- rgdal::readOGR(".", "dflow_points_v1_reposition")

# actual/mapped debris flow polygons
slide_poly_vec <- rgdal::readOGR(".", "dflow_polygons_v1_reposition_sample_100")
slide_poly_vec$objectid <- 1:100

raster::crs(slide_point_vec) <- raster::crs(slide_poly_vec)

saga <- Rsagacmd::saga_gis()


rwPerformance(dem, slide_plys = slide_poly_vec, source_pnts = slide_point_vec,
                   slide_id = 2, slp = 33, ex = 3, per = 2,
                   gpp_iter = 1000, buffer_ext = 500, buffer_source = 50,
                   plot_eval = TRUE)


steps <- 3
rwexp_vec <- seq(1.3, 3, len=steps)
rwper_vec <- seq(1.5, 2, len=steps)
rwslp_vec <- seq(20, 40, len=steps)

workspace_dir <- "/home/jason/Scratch/F_Testing"

rwGridsearch(dem, workspace = workspace_dir, slide_plys = slide_poly_vec, source_pnts = slide_point_vec,
               slide_id = 2, slp_v = rwslp_vec, ex_v = rwexp_vec, per_v = rwper_vec,
               gpp_iter = 1000, buffer_ext = 500, buffer_source = 50,
               plot_eval = TRUE)

(load("result_rw_roc_2.Rd"))
roc_result


# RUN OPTIMIZATION #############################################################

# Define which runout polygons to run grid search for

polyid_vec <- 1:length(slide_poly_vec)
polyid_vec <- 1:4

# Save RW GridSearch Settings

setwd(workspace_dir)

rw_settings <- list(
  n_train = length(polyid_vec),
  vec_rwexp = rwexp_vec,
  vec_rwper = rwper_vec,
  vec_rwslp = rwslp_vec
)

save(rw_settings, file = "gridsearch_rw_settings.Rd")


# Run in parellel

library(foreach)
cl <- parallel::makeCluster(4)
doParallel::registerDoParallel(cl)

Sys.time()
results.list <-
  foreach(poly_id=polyid_vec, .combine='c', .packages=c('rgdal','raster', 'rgeos', 'ROCR', 'Rsagacmd', 'sf')) %dopar% {
    rwGridsearch(dem, workspace = "/home/jason/Scratch/F_Testing", slide_plys = slide_poly_vec, source_pnts = slide_point_vec,
                   slide_id = poly_id, slp_v = rwslp_vec, ex_v = rwexp_vec, per_v = rwper_vec,
                   gpp_iter = 1000, buffer_ext = 500, buffer_source = 50,
                   plot_eval = FALSE)
  }
parallel::stopCluster(cl)
Sys.time()


# GET RW OPTIMAL VALUES #############################################################

(load("gridsearch_rw_settings.Rd"))

rwGetOpt(workspace = getwd(),
                    rwslp_vec = rw_settings$vec_rwslp,
                    rwexp_vec = rw_settings$vec_rwexp,
                    rwper_vec = rw_settings$vec_rwper,
                    n_train = rw_settings$n_train,
                    measure = median)


# SPCV RW ###########################################################################

setwd("/home/jason/Scratch/GPP_RW_Paper")
# Load saga gpp random walk settings
load("gridsearch_rw_settings.Rd")

rwexp_vec <- rw_settings$vec_rwexp
rwper_vec <- rw_settings$vec_rwper
rwslp_vec <- rw_settings$vec_rwslp

#polyid.vec <- 1:rw_settings$n_train

rw_spcv <- rwSPCV(slide_plys = slide_poly_vec,
            n_folds = 3,
            repetitions = 5,
            rwslp_vec = rw_settings$vec_rwslp,
            rwexp_vec = rw_settings$vec_rwexp,
            rwper_vec = rw_settings$vec_rwper)

freq_rw <- rwPoolSPCV(rw_spcv)


# VISUALIZE RW SPCV OPT PARAM SET FREQ #######################

#Pick a slope threshold (slice) of grid search space
slope_thresh <- 40
freq_rw <- freq_rw[order(freq_rw$rel_freq, decreasing = TRUE),]

require(ggplot2)
breaks_bubble <- c(round(min(freq_rw$rel_freq)),
                   round(median(freq_rw$rel_freq)),
                   round(max(freq_rw$rel_freq)))

gg_freqRW <- freqRW[freq_rw$slp == slope_thresh,]

ggplot(freq_rw, aes(x=per, y=exp)) +
  ggtitle(paste("Slope threshold:", slope_thresh ))+
  # Can improve by making different colors for different slope thresholds...

  geom_point(alpha=0.7, aes(colour = median_auroc, size = rel_freq)) +
  scale_size(range = c(2, 10), name="Relative\nfrequency (%)",
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
workspace_dir <- "/home/jason/Scratch/F_Testing"

result <- pcmPerformance(dem, slide_plys = slide_poly_vec, source_pnts = slide_point_vec,
               slide_id = 8, rw_slp = 33, rw_ex = 3, rw_per = 2,
               pcm_mu = 0.11, pcm_md = 20,
               gpp_iter = 1000, buffer_ext = 500, buffer_source = 50,
               plot_eval = TRUE, return_features = TRUE)

pcmmd_vec <- seq(50, 150, by=50)
pcmmu_vec <- seq(0.04, 0.6, by=0.25)

pcm_result <- pcmGridsearch(dem, workspace = workspace_dir,
                       slide_plys = slide_poly_vec, source_pnts = slide_point_vec, slide_id = 5,
                       rw_slp = 40, rw_ex = 1.8, rw_per = 2,
                       pcm_mu_v = pcmmu_vec, pcm_md_v = pcmmd_vec,
                       gpp_iter = 1000,
                       buffer_ext = 500, buffer_source = NULL,
                       predict_threshold = 0.5,
                       plot_eval = FALSE)

for(i in 1:10){
  pcmGridsearch(dem, workspace = workspace_dir,
             slide_plys = slide_poly_vec, source_pnts = slide_point_vec, slide_id = i,
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

pcm_spcv <- pcmSPCV(slide_plys = slide_poly_vec[1:10,],
                    n_folds = 3, repetitions = 100, pcm_mu_v = pcmmu.vec,
                    pcm_md_v = pcmmd.vec)

# Test pooling
setwd("/home/jason/Scratch/GPP_PCM_Paper")

(load("gridsearch_pcm_settings.Rd"))

polyid.vec <- 1:pcm_settings$n_train
pcmmu.vec <- pcm_settings$vec_pcmmu
pcmmd.vec <- pcm_settings$vec_pcmmd

(load("repeated_spcv_PCM.Rd"))

freq_pcm <- pcmPoolSPCV(rep_spcv_pcm)


