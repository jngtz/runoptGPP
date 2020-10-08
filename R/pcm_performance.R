
#' PCM Runout Distance Error
#'
#' Using min. area bouding boxes to calculate predicted length, it measures
#'     the relative error, relative difference and error.
#' @param obs_poly A RasterLayer object
#' @return RasterLayer

errMinBboxLength <- function(obs_poly, pred_raster, elev_raster){
  #Calculates the relative error b/w bbox length estimates of slides
  sp.pred <- rasterToPolygons(pred_raster, n = 4, dissolve = TRUE, na.rm=TRUE)
  geom_pred <- runoutGeom(sp.pred, elev = elev_raster)
  geom_act <- runoutGeom(obs_poly, elev = elev_raster)

  err <- geom_pred$length - geom_act$length
  rel_err <- abs(geom_act$length - geom_pred$length) / geom_act$length
  rel_diff<- (geom_pred$length - geom_act$length) / geom_act$length

  return(
    list(rel_error = rel_err,
         rel_difference = rel_diff,
         error = error
    ))
}



