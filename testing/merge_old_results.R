
# Merge PCM relerr and auroc into one object ################

setwd("/home/jason/Scratch/GPP_PCM_Paper")
# Load saga gpp random walk settings
load("gridsearch_pcm_settings.Rd")

n_train <- 100

for(i in 1:n_train){
  roc_res_nm <- paste("result_roc_", i, ".Rd", sep="")
  load(roc_res_nm) #res
  roc_res <- roc_result

  relerr_res_nm <- paste("result_relerr_length_", i, ".Rd", sep="")
  load(relerr_res_nm) #res
  relerr_res <- relerr_length_result

  result_pcm <- list(
    auroc = roc_result,
    relerr = relerr_res
  )

  save(result_pcm, file = paste("result_pcm_gridsearch_", i, ".Rd", sep=""))
}

# Create a list of all PCM results ###########

setwd("/home/jason/Scratch/GPP_PCM_Paper")
# Load saga gpp random walk settings
load("gridsearch_pcm_settings.Rd")

n_train <- 100

pcm_gridsearch_multi <- list()

for(i in 1:n_train){
  roc_res_nm <- paste("result_roc_", i, ".Rd", sep="")
  load(roc_res_nm) #res
  roc_res <- roc_result

  relerr_res_nm <- paste("result_relerr_length_", i, ".Rd", sep="")
  (load(relerr_res_nm)) #res
  relerr_res <- relerr_length_result

  result_pcm <- list(
    auroc = roc_result,
    relerr = relerr_res
  )

  pcm_gridsearch_multi[[i]] <- result_pcm
}

save(pcm_gridsearch_multi, file="pcm_gridsearch_multi.Rd")

# Create a list of all RW results ############

setwd("/home/jason/Scratch/GPP_RW_Paper")

n_train <- 100

rw_gridsearch_multi <- list()

for(i in 1:n_train){
  roc_res_nm <- paste("result_rw_roc_", i, ".Rd", sep="")
  load(roc_res_nm) #res
  rw_gridsearch_multi[[i]] <- roc_result
}

# Create source area as polygon ##############
source_list <- list()

for(i in 1:10){
  slide_poly_single <- runout_polygons[i,]

  # Crop dem to slide polygon
  dem_grid <- raster::crop(dem, raster::extent(slide_poly_single) + 500)

  # Subset corresponding source/start point of runout
  sel_over_start_point  <- sp::over(source_points, slide_poly_single)
  sel_start_point <- source_points[!is.na(sel_over_start_point$objectid),]


  # Create a buffer around source point to create a source/release area
  source_buffer <- rgeos::gBuffer(sel_start_point, width = 20)
  # Clip polygon using border of intersecting runout polygon
  source_list[[i]] <- raster::intersect(source_buffer, slide_poly_single)
}

source_areas <- do.call(raster::bind, source_list)




source_grid <- raster::rasterize(source_buffer, dem_grid, field=1 )
source_grid <- raster::mask(source_grid, slide_poly_single )



# For visualize rw freq testing ##########

test_rw <- data.frame(
  slp = c(20, 20, 30, 30, 30, 10),
  per = c(1.95, 1.95, 1.80, 1.80, 2.0, 2.0),
  exp = c(3, 2.3, 2.66, 1.6, 1.9, 1.7),
  freq = c(10, 15, 1, 16, 20, 4),
  rel_freq = rep(NA, 6),
  median_auroc = c(.92, .8, .71, .82, .95, .90),
  iqr_auroc = rep(0, 6)
)
test_freq <- rbind(test_rw, freq_rw)
test_freq$rel_freq = test_freq$freq/sum(test_freq$freq) * 100
freq_rw <- test_freq


# Visualize freq of SPCV optimal parameter sets ################################

library(ggplot2)

#Pick a slope threshold (slice) of grid search space
slope_thresh <- 40
freq_rw <- freq_rw[order(freq_rw$rel_freq, decreasing = TRUE),]

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
                     limits = c(min(rw_spcv$settings$rwper_vec), max = max(rw_spcv$settings$rwper_vec))) +
  scale_y_continuous(expression(paste("Exponent of divergence")),
                     limits = c(min(rw_spcv$settings$rwexp_vec), max = max(rw_spcv$settings$rwexp_vec)+.1)) +
  theme_light() +
  theme(text = element_text(family = "Arial", size = 8), axis.title = element_text(size = 9),
        axis.text = element_text(size = 8), title = element_text(size = 8))

