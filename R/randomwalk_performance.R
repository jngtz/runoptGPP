


#' Random walk runout path performance
#'
#' Computes the area under the receiver operating characteristic curve (AUROC) for
#'      runout paths simuluated using the random walk model component of the
#'      GPP tool in SAGA-GIS. The AUROC compares a runout polygon to the simulated
#'      path.
#' @param dem A DEM as a RasterLayer object
#' @param slide_plys Runout tracks as a SpatialPolygonsDataFrame
#' @param source_pnts Source points as a SpatialPointsDataFrame
#' @param slide_id Selects a single runout polygon from slide_plys by row
#' @param slp Random walk slope threshold - below lasteral spreading is modelled
#' @param ex Random walk exponent controlling lateral spread
#' @param per Random walk persistence factor to weight flow direction consistency
#' @param gpp_iter Number of random walk model iterations
#' @param buffer_ext Defines buffer distance (in meters) around runout polygon
#'      to crop source DEM. This helps to reduce computational time
#' @param buffer_source Can define a buffer distance (in meters) to extend source
#'      point to a source area
#' @param plot_eval If TRUE, will plot random walk path and runout polygon
#' @return the area under the receiver operating characteristic
#' @details Runout is either simulated from a single source point or a buffered
#'      area round the source point.
#' @examples
#' ## Not run:
#' # Initalize a saga object
#' saga <- Rsagacmd::saga_gis()
#'

performanceRndWalk <- function(dem, slide_plys, source_pnts,
                          slide_id = 2, slp = 33, ex = 3, per = 2,
                          gpp_iter = 1000, buffer_ext = 500, buffer_source = NULL,
                          plot_eval = FALSE)
{

  # Subset a single slide polygon
  slide_poly_single <- slide_plys[slide_id,]

  # Crop dem to slide polygon
  dem_grid <- raster::crop(dem, raster::extent(slide_poly_single) + buffer_ext)

  # Subset corresponding source/start point of runout
  sel_over_start_point  <- sp::over(source_pnts, slide_poly_single)
  sel_start_point <- source_pnts[!is.na(sel_over_start_point$objectid),]

  if(!is.null(buffer_source)){
    # Create a buffer around source point to create a source/release area
    source_buffer <- rgeos::gBuffer(sel_start_point, width = buffer_source)
    source_grid <- raster::rasterize(source_buffer, dem_grid, field=1 )
    source_grid <- raster::mask(source_grid, slide_poly_single )
  } else {
    # Just use source point
    source_grid <- raster::rasterize(matrix(sp::coordinates(sel_start_point)[1:2], ncol = 2), dem_grid, field = 1)
  }



  # Run runout model using Rsagacmd (faster than RSAGA)
  gpp <- saga$sim_geomorphology$gravitational_process_path_model(dem = dem_grid, release_areas = source_grid,
                                                                 #friction_mu_grid = mu.grid,
                                                                 process_path_model = 1,
                                                                 rw_slope_thres = slp,
                                                                 rw_exponent = ex,
                                                                 rw_persistence = per,
                                                                 gpp_iterations = gpp_iter)



  rescale_process_area <- rasterCdf(gpp$process_area)

  # AUROC
  pred_values <- raster::getValues(rescale_process_area)
  pred_values[is.na(pred_values)] <- 0

  slide_area <- raster::rasterize(slide_poly_single, rescale_process_area, field=1, background = 0)

  pred_area <- ROCR::prediction(predictions = pred_values, labels = raster::getValues(slide_area))

  auroc_area <- ROCR::performance(pred_area, "auc")
  roc <- auroc_area@y.values[[1]]

  # Plot evaluation results
  if(plot_eval){
    raster::plot(rescale_process_area,
         main = paste("id", slide_id,
                      "roc", round(roc, digits=3), "\n",
                      "Ex", rw_ex, "Pr", rw_per, "Slp", rw_slp),
         cex.main = 0.7, cex.axis = 0.7, cex=0.7)
    sp::plot(slide_poly_single, add=TRUE)
  }

  return(roc)

}


