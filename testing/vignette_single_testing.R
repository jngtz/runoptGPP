setwd("/home/jason/R/runoptGPP/")

#library(devtools)
#document()
#build()
#install()

#check(vignettes = FALSE)


# Load Packages and Data ####################################################################
library(runoptGPP)
library(raster)
library(rgdal)
library(sp)
library(Rsagacmd)

# Initiate a SAGA-GIS geoprocessing object
saga <- saga_gis(opt_lib = "sim_geomorphology")

# Set workspace
setwd("/home/jason/Data/Chile/")

# Load digital elevation model (DEM)
dem <- raster("elev_alos_12_5m.tif")
fd <- terrain(dem, opt = "flowdir")

# Load runout source points
source_points <- readOGR(".", "debris_flow_source_points")

# Load runout track polygons and assign object ID based on row number
runout_polygons <- readOGR(".", "debris_flow_polys_sample")
runout_polygons$objectid <- 1:length(runout_polygons)

# Select a debris flow and source point for this example
runout_polygon <- runout_polygons[77,]
#runout_polygon <- runout_polygons[86,]
plot(runout_polygon)


par(mfrow = c(3, 6))
for(i in 1:length(runout_polygons)){
  runout_poly <- runout_polygons[i,]
  plot(runout_poly, main = i)
}


# Select a debris flow and source point for this example

#runout_polygon <- runout_polygons[92,]
runout_polygon <- runout_polygons[77,] # example in paper
#runout_polygon <- runout_polygons[86,]
plot(runout_polygon)

sel_source_point  <- over(source_points, runout_polygon)
source_point <- source_points[!is.na(sel_source_point$objectid),]


# TESTING flow path #############################
setwd("/home/jason/Data/Chile/saga_data")
slp <- raster("slope_alos_12_5m_filled.sdat")
crop_slp <- crop(slp, runout_polygon)
mask_slp <- mask(crop_slp, runout_polygon)

crop_dem <- crop(dem, runout_polygon)
mask_dem <- mask(crop_dem, runout_polygon)

crop_fd <- crop(fd, runout_polygon)

plot(mask_dem %in% minValue(mask_dem))

plot(mask_slp)
plot(source_point, add = TRUE)


flow_pth <- flowPath(fd, coordinates(source_point))

xy <- xyFromCell(fd, flow_pth)
plot(crop_dem)
plot(runout_polygon, add = TRUE)
lines(xy)

LineLength(xy) / 1000 # in km



points <- sp::SpatialPoints(xy)

# use as to convert to line
sp_line <- as(points,"SpatialLines")
crs(sp_line) <- crs(runout_polygon)

ll_sp_line <- spTransform(sp_line, CRS("+proj=longlat +datum=WGS84"))
SpatialLinesLengths(ll_sp_line, longlat = TRUE)

# GPP random walk simulation ################################


dem
slide_plys = runout_polygon
slide_src = source_point
slp = 33
ex = 3
per = 2
gpp_iter = 1000
buffer_ext = 500
buffer_source = 50
measure = "auroc"
plot_eval = TRUE
saga_lib = saga

rwPerformance_test <- function(dem, slide_plys, slide_src,
                          slide_id = 2, slp = 33, ex = 3, per = 2,
                          gpp_iter = 1000, buffer_ext = 500, buffer_source = 30,
                          measure = "auroc", plot_eval = TRUE, saga_lib)
{

  # If single runout polygon as input, assign slide_id value of 1
  if(length(slide_plys) == 1){
    slide_id <- 1
  }

  slide_plys$objectid <- 1:length(slide_plys)
  # Subset a single slide polygon
  slide_poly_single <- slide_plys[slide_id,]

  # Crop dem to slide polygon
  dem_grid <- raster::crop(dem, raster::extent(slide_poly_single) + buffer_ext)


  if(class(slide_src) == "SpatialPointsDataFrame"){
    # Subset corresponding source/start point of runout
    sel_over_start_point  <- sp::over(slide_src, slide_poly_single)
    sel_start_point <- slide_src[!is.na(sel_over_start_point$objectid),]

    if(!is.null(buffer_source)){
      # Create a buffer around source point to create a source/release area
      source_buffer <- rgeos::gBuffer(sel_start_point, width = buffer_source)
      source_grid <- raster::rasterize(source_buffer, dem_grid, field=1 )
      source_grid <- raster::mask(source_grid, slide_poly_single )
      source_plot <- raster::rasterToPolygons(source_grid)
    } else {
      # Just use source point
      source_plot <- sel_start_point
      source_grid <- raster::rasterize(matrix(sp::coordinates(sel_start_point)[1:2], ncol = 2), dem_grid, field = 1)
    }
  }

  if(class(slide_src) == "SpatialPolygonsDataFrame" ){
    sel_over_start_poly <- sp::over(slide_src, slide_poly_single)
    sel_start_poly <- slide_src[!is.na(sel_over_start_poly$objectid),]
    source_plot <- sel_start_poly
    source_grid <- raster::rasterize(sel_start_poly, dem_grid, field=1 )

  }





  # Run runout model using Rsagacmd (faster than RSAGA)
  gpp <- saga_lib$sim_geomorphology$gravitational_process_path_model(dem = dem_grid, release_areas = source_grid,
                                                                     process_path_model = 1,
                                                                     rw_slope_thres = slp,
                                                                     rw_exponent = ex,
                                                                     rw_persistence = per,
                                                                     gpp_iterations = gpp_iter)



  rescale_process_area <- rasterCdf(gpp$process_area)

  # AUROC
  #NEW
  #rescale_process_area[is.na(rescale_process_area)] <- 0
  #buff_slide <- raster::buffer(slide_poly_single, width = 200)
  #rescale_process_area <- raster::mask(rescale_process_area, buff_slide)

  # Down sampling negatives

  plot(dem_grid)
  max_elv <- max(getValues(raster::mask(dem_grid, source_grid)), na.rm = TRUE)
  mask_elv <- dem_grid
  mask_elv[mask_elv > max_elv] <- NA


  pred_values <- raster::getValues(rescale_process_area)
  pred_values[is.na(pred_values)] <- 0

  slide_area <- raster::rasterize(slide_poly_single, rescale_process_area, field=1, background = 0)
  #slide_area <- raster::mask(slide_area, buff_slide)

  obs_values <- raster::getValues(slide_area)
  pred_area <- ROCR::prediction(predictions = pred_values[!is.na(pred_values)], labels = obs_values[!is.na(obs_values)])

  auroc_area <- ROCR::performance(pred_area, "auc")
  roc <- auroc_area@y.values[[1]]

  aucpr_area <- ROCR::performance(pred_area, "aucpr")
  aucpr_area@y.values[[1]]

  aucpr_area <- ROCR::performance(pred_area, "acc")
  aucpr_area@y.values[[1]]

  # area under precision recall curve
  #roc<-PRROC::roc.curve(scores.class0 = pred_values[!is.na(pred_values)], weights.class0 = obs_values[!is.na(obs_values)])
  pr<-PRROC::pr.curve(scores.class0 = pred_values[!is.na(pred_values)], weights.class0 = obs_values[!is.na(obs_values)], curve = TRUE)

  # aucpr curve is relative to the % of positives....
  per_pos <- sum(obs_values[!is.na(obs_values)]) / length(obs_values[!is.na(obs_values)]) * 100
  #https://glassboxmedicine.com/2019/03/02/measuring-performance-auprc/



  #iou
  pred_runout <- gpp$process_area
  pred_runout[pred_runout > 0] = 1
  pred_runout[is.na(pred_runout)] = 0

  obs_runout <- raster::rasterize(slide_poly_single, dem_grid)
  obs_runout[is.na(obs_runout)] = 0

  intersect_runout <- pred_runout + obs_runout

  sum_inter <- raster::freq(intersect_runout, value = 2)
  sum_union <- sum_inter + raster::freq(intersect_runout, value = 1)

  iou <- sum_inter / sum_union


  # Plot evaluation results
  if(plot_eval){
    raster::plot(rescale_process_area,
                 main = paste("id", slide_id, "iou", round(iou, digits = 3),
                              "auroc", round(roc, digits=3), "\n",
                              "Ex", ex, "Pr", per, "Slp", slp),
                 cex.main = 0.7, cex.axis = 0.7, cex=0.7)
    sp::plot(slide_poly_single, add=TRUE)
    sp::plot(source_plot, add=TRUE)
  }


  if(measure == "iou"){
    rw_perf = iou
  }

  if(measure == "pr"){
    rw_perf = as.numeric(pr[2])
  }

  if(measure == "auroc"){
    rw_perf = roc
  }


  return(rw_perf)

}




rwGridsearch_test <- function(dem, slide_plys, slide_src,
                         slide_id = NULL, slp_v, ex_v, per_v,
                         gpp_iter = 1000, buffer_ext = 500, buffer_source = NULL,
                         measure = "iou", save_res = FALSE, plot_eval = FALSE, saga_lib)

{

  if(is.null(slide_id)){
    slide_id = 1
  }

  roc_result.nm <- paste("result_rw_roc_", slide_id, ".Rd", sep="")

  column.names <- ex_v
  row.names <- slp_v
  matrix.names <- per_v

  roc_result <- array(NA ,dim = c(length(slp_v), length(ex_v), length(per_v)),dimnames = list(row.names, column.names, matrix.names))

  #roc[row, col, matrix]
  #roc[rwslp, rwexp, rwper]

  for(j in 1:length(per_v)){
    for(i in 1:length(ex_v)){
      for(k in 1:length(slp_v)){
        #res[k, i, j] <- paste(pcmmd[k], rwexp[i], rwper[j] )
        roc <- rwPerformance_test(dem,
                             slide_plys = slide_plys,
                             slide_src = slide_src,
                             slide_id = slide_id,
                             slp = slp_v[k],
                             ex = ex_v[i],
                             per = per_v[j],
                             gpp_iter = gpp_iter,
                             buffer_ext = buffer_ext,
                             buffer_source = buffer_source,
                             measure = measure,
                             plot_eval = plot_eval,
                             saga_lib = saga_lib)

        roc_result[k, i, j] <- roc

      }}}


  if(save_res){
    save(roc_result, file=roc_result.nm)
  }

  return(roc_result)
}



setwd("/home/jason/Scratch/NullExt_GPP_Compare")



for(i in 1:length(runout_polygons)){
  runout_poly <- runout_polygons[i,]
  sel_source_point  <- over(source_points, runout_poly)
  source_point <- source_points[!is.na(sel_source_point$objectid),]

  png_nm = paste0("rw_perf_compare", i, ".png")
  png(png_nm, width = 480, height = 800)
  par(mfrow = c(2,1))

  rwPerformance(dem, slide_plys = runout_poly, slide_src = source_point,
                slp = 40, ex =  3 , per = 1.9,
                gpp_iter = 1000, buffer_ext = 500, buffer_source = 50,
                plot_eval = TRUE, saga_lib = saga)

  # Opt using poly extent
  rwPerformance(dem, slide_plys = runout_poly, slide_src = source_point,
                slp = 38, ex =  2.32 , per = 1.7,
                gpp_iter = 1000, buffer_ext = 500, buffer_source = 50,
                plot_eval = TRUE, saga_lib = saga)

  dev.off()
}

# Opt with 500 m buffer


# Test measures ################################
# Define grid search values
steps <- 3
rwexp_vec <- seq(1.3, 3, len=2)
rwper_vec <- seq(1.5, 2, len=2)
rwslp_vec <- seq(30, 60, len=5)

gs_auroc <- rwGridsearch_test(dem, slide_plys = runout_polygon, slide_src = source_point,
                           #Input random walk grid search space
                           slp_v = rwslp_vec, ex_v = rwexp_vec, per_v = rwper_vec,
                           #Set number of simulation iterations
                           gpp_iter = 1000,
                           #Define processing extent size (m)
                           buffer_ext = 0,
                           #(Optional) Define size of buffer to make source area from point
                           buffer_source = 50,
                           measure = "auroc", saga_lib = saga)

opt_auroc <- rwGetOpt_single(gs_auroc)
opt_auroc

gs_iou <- rwGridsearch_test(dem, slide_plys = runout_polygon, slide_src = source_point,
               #Input random walk grid search space
               slp_v = rwslp_vec, ex_v = rwexp_vec, per_v = rwper_vec,
               #Set number of simulation iterations
               gpp_iter = 1000,
               #Define processing extent size (m)
               buffer_ext = 500,
               #(Optional) Define size of buffer to make source area from point
               buffer_source = 50,
               measure = "iou", saga_lib = saga)

gs_pr <- rwGridsearch_test(dem, slide_plys = runout_polygon, slide_src = source_point,
                           #Input random walk grid search space
                           slp_v = rwslp_vec, ex_v = rwexp_vec, per_v = rwper_vec,
                           #Set number of simulation iterations
                           gpp_iter = 1000,
                           #Define processing extent size (m)
                           buffer_ext = 500,
                           #(Optional) Define size of buffer to make source area from point
                           buffer_source = 50,
                           measure = "pr", saga_lib = saga)


opt_par <- rbind(rwGetOpt_single(gs_auroc), rwGetOpt_single(gs_iou), rwGetOpt_single(gs_pr))
opt_par$metric <- c("auroc", "iou", "pr")
opt_par

opt_auroc <- rwGetOpt_single(gs_auroc)
opt_iou <- rwGetOpt_single(gs_iou)
opt_pr <- rwGetOpt_single(gs_pr)

par(mfrow = c(3,1))
#opt auroc
rwPerformance(dem, slide_plys = runout_polygon, slide_src = source_point,
              slp = opt_auroc$rw_slp_opt, ex =  opt_auroc$rw_exp_opt , per = opt_auroc$rw_per_opt,
              gpp_iter = 1000, buffer_ext = 500, buffer_source = 50,
              plot_eval = TRUE, saga_lib = saga)

#opt iou
rwPerformance(dem, slide_plys = runout_polygon, slide_src = source_point,
              slp = opt_iou$rw_slp_opt, ex =  opt_iou$rw_exp_opt , per = opt_iou$rw_per_opt,
              gpp_iter = 1000, buffer_ext = 500, buffer_source = 50,
              plot_eval = TRUE, saga_lib = saga)


#opt pr
rwPerformance(dem, slide_plys = runout_polygon, slide_src = source_point,
              slp = opt_pr$rw_slp_opt, ex =  opt_pr$rw_exp_opt , per = opt_pr$rw_per_opt,
              gpp_iter = 1000, buffer_ext = 500, buffer_source = 50,
              plot_eval = TRUE, saga_lib = saga)


# GPP PCM simulation #######################

pcm <- pcmPerformance(dem, slide_plys = runout_polygon, slide_src = source_point,
               rw_slp = 40, rw_ex = 3, rw_per = 1.5,
               pcm_mu = 0.15, pcm_md = 120,
               gpp_iter = 1000, buffer_ext = 500, buffer_source = 50,
               plot_eval = TRUE, return_features = TRUE, saga_lib = saga)

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
                       slide_plys = runout_polygon, slide_src = source_point,
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
                       buffer_source = 50, saga_lib = saga)

# Get optimal parameters
pcmGetOpt_single(pcm_gridsearch)


