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
  source_list[[i]]<- raster::intersect(source_buffer, slide_poly_single)
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
