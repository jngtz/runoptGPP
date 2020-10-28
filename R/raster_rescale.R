#' Cumulative Distribution Function for Rasters
#'
#' Compute an empirical cumulative distribution function for a raster.
#' @param x A RasterLayer object
#' @return RasterLayer
#' @examples
#' r <- raster(system.file("external/rlogo.grd", package="raster"))
#' r_ecdf <- rasterCdf(logo)
#' plot(r_ecdf)

rasterCdf <- function(x){

  x_ecdf <- ecdf(raster::getValues(x))
  prob_x <- raster::setValues(x, x_ecdf(raster::getValues(x)))
  prob_x

}



#' Min-max normalization (Rescaling)
#'
#' Rescale a raster values in the range from 0 to 1 using min-max normalization.
#' @param x A RasterLayer object
#' @return RasterLayer
#' @examples
#' r <- raster(system.file("external/rlogo.grd", package="raster"))
#' r_minmax <- rasterMinMaxRescale(logo)
#' plot(r_minmax)

rasterMinMaxRescale <- function(x){

  ( x - min(raster::getValues(x), na.rm=TRUE)) / (max(raster::getValues(x), na.rm=TRUE) -  min(raster::getValues(x), na.rm=TRUE) )

}
