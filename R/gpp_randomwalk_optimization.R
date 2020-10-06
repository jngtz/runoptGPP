
#' Cumulative Distribution Function for Rasters
#'
#' Compute an empirical cumulative distribution function for a raster.
#' @param x A RasterLayer object
#' @return RasterLayer
#' @examples
#' ## Not run:
#' # Initalize a saga object
#' saga <- saga_gis()
#'

SagaGppRndWlk <- function(dem, slide_plys = slide_poly_vec, release_pnts = slide_point_vec,
                          slide_id = 2, rw_slp = 33, rw_ex = 3, rw_per = 2,
                          buffer_ext = 500, buffer_source = 50, gpp_iter = 1000,
                          plot_eval = FALSE)
{

  poly_id = slide_id
  slide_poly_vec = slide_plys
  slide_point_vec = release_pnts

  # Subset a single slide polygon
  slide_poly_single <- slide_poly_vec[poly_id,]

  # Crop dem to slide polygon
  dem_grid <- raster::crop(dem, extent(slide_poly_single) + buffer_ext)

  # Subset corresponding source/start point of runout
  sel_over_start_point  <- sp::over(slide_point_vec, slide_poly_single)
  sel_start_point <- slide_point_vec[!is.na(sel_over_start_point$objectid),]

  # Create a buffer around source point to create a source/release area
  source_buffer <- gBuffer(sel_start_point, width = buffer_source)
  source_grid <- raster::rasterize(source_buffer, dem_grid, field=1 )
  source_grid <- raster::mask(source_grid, slide_poly_single )

  # Run runout model using Rsagacmd (faster than RSAGA)
  gpp <- saga$sim_geomorphology$gravitational_process_path_model(dem = dem_grid, release_areas = source_grid,
                                                                 #friction_mu_grid = mu.grid,
                                                                 process_path_model = 1,
                                                                 rw_slope_thres = rw_slp,
                                                                 rw_exponent = rw_ex,
                                                                 rw_persistence = rw_per,
                                                                 gpp_iterations = gpp_iter)


  # Rescale to values from 0 to 1
  #rescale_process_area <- RasterMinMaxRescale(gpp$process_area)
  rescale_process_area <- RasterCdf(gpp$process_area)

  # AUROC
  pred_values <- getValues(rescale_process_area)
  pred_values[is.na(pred_values)] <- 0

  slide_area <- rasterize(slide_poly_single, rescale_process_area, field=1, background = 0)

  pred_area <- prediction(predictions = pred_values, labels = getValues(slide_area))
  perf_area <- performance(pred_area, "tpr", "fpr")

  auroc_area <- performance(pred_area, "auc")
  roc <- auroc_area@y.values[[1]]

  # Plot evaluation results
  if(plot_eval){
    plot(rescale_process_area,
         main = paste("id", poly_id,
                      "roc", round(roc, digits=2), "\n",
                      "Ex", rw_ex, "Pr", rw_per, "Slp", rw_slp),
         cex.main = 0.7, cex.axis = 0.7, cex=0.7)
    plot(slide_poly_single, add=TRUE)
  }

  return(roc)

}

