

stops <- rasterToPoints(rasterCdf(gpp$gpp.stop))

#' PCM Runout Distance Error
#'
#' Using min. area bounding boxes to calculate predicted length, it measures
#'     the relative error, relative difference and error.
#' @param obs_poly The observed/mapped runout track as a SpatialPolygonsDataFrame
#' @param pred_raster A RasterLayer object predicting runout
#' @param dem A DEM as a RasterLayer object
#' @return A list of error measures

errMinBboxLength <- function(obs_poly, pred_raster, dem){
  #Calculates the relative error b/w bbox length estimates of slides
  sp.pred <- raster::rasterToPolygons(pred_raster, n = 4, dissolve = TRUE, na.rm=TRUE)
  geom_pred <- runoutGeom(sp.pred, elev = dem)
  geom_act <- runoutGeom(obs_poly, elev = dem)

  err <- geom_pred$length - geom_act$length
  rel_err <- abs(geom_act$length - geom_pred$length) / geom_act$length
  rel_diff<- (geom_pred$length - geom_act$length) / geom_act$length

  return(
    list(rel_error = rel_err,
         rel_difference = rel_diff,
         error = err
    ))
}



#' PCM runout distance performance
#'
#' Computes the error for runout distances simuluated using the random walk
#'      and PCM model components of the GPP tool in SAGA-GIS.
#' @param dem A DEM as a RasterLayer object
#' @param slide_plys Runout tracks as a SpatialPolygonsDataFrame
#' @param slide_src Source points as a SpatialPointsDataFrame or source areas
#'      as a SpatialPolygonsDataFrame
#' @param slide_id Selects a single runout polygon from slide_plys by row
#' @param rw_slp Random walk slope threshold - below lateral spreading is modelled
#' @param rw_ex Random walk exponent controlling lateral spread
#' @param rw_per Random walk persistence factor to weight flow direction consistency
#' @param pcm_mu PCM model sliding friction coefficient
#' @param pcm_md PCM model mass-to-drag ratio (m)
#' @param gpp_iter Number of model iterations
#' @param buffer_ext (Optional) Defines buffer distance (in meters) around runout polygon
#'      to crop source DEM. This helps to reduce computational time
#' @param buffer_source (Optional) Can define a buffer distance (in meters) to extend source
#'      point to a source area
#' @param predict_threshold A cutoff value to define what quantile of simulated runout
#'      frequencies is the predicted runout.
#' @param plot_eval logical. If TRUE will plot simulated runout and runout polygon
#' @param return_features logical. If TRUE, returned list will include GPP input and output
#'      data, in addition to a list of error measures.
#' @param saga_lib The initiated SAGA-GIS geoprocessor object
#' @return A list of runout distance performance measures.
#' @examples
#' \dontrun{
#' # Initialize a saga object
#' saga <- Rsagacmd::saga_gis()
#'
#' # Load elevation model (DEM)
#' dem <- raster(system.file("extdata/elev_12_5m.tif", package="runout.opt"))
#'
#' # Load runout polygons and source points
#' runout_plys <- rgdal::readOGR(system.file("extdata/dflow_runout_ply.shp", package="runout.opt"))
#' source_pnts <- rgdal::readOGR(system.file("extdata/dflow_source_pnt.shp", package="runout.opt"))
#'
#' # Run GPP PCM model for a rounout polygon
#' pcm <- pcmPerformance(dem, slide_plys = runout_plys[1,], slide_src = source_pnts,
#'   rw_slp = 40, rw_ex = 3, rw_per = 1.5,
#'   pcm_mu = 0.15, pcm_md = 120,
#'   gpp_iter = 1000, buffer_ext = 500, buffer_source = 50,
#'   plot_eval = TRUE, return_features = TRUE)
#'
#' # Runout distance relative error
#' pcm$length.relerr
#'
#' # Plot GPP PCM runout modelling ouputs
#' gpp_output <- stack(pcm$gpp.parea, pcm$gpp.stop, pcm$gpp.maxvel)
#' names(gpp_output) <- c("Process_area", "Stop_positions", "Max_velocity")
#' plot(gpp_output)
#'
#' }

pcmPerformance <- function(dem, slide_plys, slide_src, slide_id = 1,
                             rw_slp = 33, rw_ex = 3, rw_per = 2,
                             pcm_mu = 0.3, pcm_md = 75,
                             buffer_ext = 500, buffer_source = 50, gpp_iter = 1000,
                             predict_threshold = 0.5, plot_eval = FALSE,
                             return_features = FALSE, saga_lib)
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
                                                                 #friction_mu_grid = mu.grid,
                                                                 process_path_model = 1,
                                                                 rw_slope_thres = rw_slp,
                                                                 rw_exponent = rw_ex,
                                                                 rw_persistence = rw_per,
                                                                 gpp_iterations = gpp_iter,
                                                                 friction_model = 5,
                                                                 friction_mu = pcm_mu,
                                                                 friction_mass_to_drag = pcm_md)


  # Rescale to values from 0 to 1
  rescale_process_area <- rasterCdf(gpp$process_area)

  # AUROC
  pred_values <- raster::getValues(rescale_process_area)
  pred_values[is.na(pred_values)] <- 0

  slide_area <- raster::rasterize(slide_poly_single, gpp$process_area, field=1, background = 0)

  pred_thres <- rescale_process_area
  pred_thres[pred_thres >= predict_threshold] = 1
  pred_thres[pred_thres < predict_threshold] = NA

  pred_area <- ROCR::prediction(predictions = pred_values, labels = raster::getValues(slide_area))
  perf_area <- ROCR::performance(pred_area, "tpr", "fpr")
  auroc_area <- ROCR::performance(pred_area, "auc")
  roc <- auroc_area@y.values[[1]]

  # Length loss
  errMinBox <- errMinBboxLength(obs_poly = slide_poly_single,
                                pred_raster = pred_thres,
                                dem = dem)

  length_relerr <- errMinBox[["rel_error"]]
  length_reldiff <- errMinBox[["rel_difference"]]
  length_error <- errMinBox[["error"]]

  # Plot evaluation results
  if(plot_eval){
    sp::plot(rescale_process_area,
         main = paste("id", slide_id,
                      "roc", round(roc, digits=2), "\n",
                      "err.L", round(length_relerr, digits = 2),
                      "Ex", rw_ex, "Pr", rw_per, "Mu", pcm_mu, "Md", pcm_md),
         cex.main = 0.7, cex.axis = 0.7, cex=0.7)
    sp::plot(slide_poly_single, add=TRUE)
    sp::plot(source_plot, add = TRUE)
  }


  if(return_features){
    return(
      list(id = slide_id,
           roc = roc,
           length.relerr = length_relerr,
           length.reldiff = length_reldiff,
           length.error = length_error,
           dem = dem_grid,
           actual.poly = slide_poly_single,
           source.pnt = sel_start_point,
           source.area = source_grid,
           gpp.parea = gpp$process_area,
           gpp.stop = gpp$stop_positions,
           gpp.maxvel = gpp$max_velocity))
  } else {
    return(
      list(id =  slide_id,
           roc = roc,
           length.relerr = length_relerr,
           length.reldiff = length_reldiff,
           length.error = length_error))
  }

}




