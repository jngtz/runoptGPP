#' Grid search optimization for random walk runout path
#'
#' Computes the area under the receiver operating characteristic curve (AUROC) for
#'      runout paths simuluated using the random walk model component of the
#'      GPP tool in SAGA-GIS. The AUROC compares a runout polygon to the simulated
#'      path.
#' @param dem A DEM as a RasterLayer object
#' @param workspace The file path where to save performance results for each runout polygon
#' @param slide_plys Runout tracks as a SpatialPolygonsDataFrame
#' @param source_pnts Source points as a SpatialPointsDataFrame
#' @param slide_id Selects a single runout polygon from slide_plys by row
#' @param slp_v A vector of random walk slope thresholds - below lasteral spreading is modelled
#' @param ex_v A vector or random walk exponents controlling lateral spread
#' @param per_v A vector or random walk persistence factors to weight flow direction consistency
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


gridOptRndWalk <- function(dem, workspace = getwd(), steps = 1, slide_plys, source_pnts,
                      slide_id, slp_v, ex_v, per_v,
                      gpp_iter = 1000, buffer_ext = 500, buffer_source = NULL,
                      plot_eval = FALSE)

{

  revert_wd <- getwd()
  setwd(workspace)

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
        roc <- performanceRndWalk(dem,
                             slide_plys = slide_plys,
                             source_pnts = source_pnts,
                             slide_id = slide_id,
                             slp = slp_v[k],
                             ex = ex_v[i],
                             per = per_v[j],
                             gpp_iter = gpp_iter,
                             buffer_ext = buffer_ext,
                             buffer_source = buffer_source,
                             plot_eval = plot_eval)

        roc_result[k, i, j] <- roc

      }}}

  save(roc_result, file=roc_result.nm)
  setwd(revert_wd)
}





#' Get random walk grid search optimal parameters
#'
#' @param workspace The file path where the performance results were saved
#' @param slide_id Selects a single runout polygon from slide_plys by row
#' @param rwslp_vec The vector of random walk slope thresholds used in the grid search
#' @param rwexp_vec The vector or random walk exponents used in the grid search
#' @param rwper_vec The vector or random walk persistence factors used in the grid search
#' @param n_train Number of runout polygons processed with grid search
#' @param measure A measure (e.g. median, mean) to find optimal parameter across grid search space
#' @return A dataframe  with the optimal parameter set and AUROC performance

getRndWalkOptParams <- function(workspace = getwd(), rwslp_vec, rwexp_vec, rwper_vec, n_train,
                                measure = median){

  roc_list <- list()

  for(PER in 1:length(rwper_vec)){

    per_list <- list()

    for(i in 1:n_train){
      res_nm <- paste("result_rw_roc_", i, ".Rd", sep="")
      load(res_nm) #res
      res <- roc_result
      per_list[[i]] <- as.matrix(res[,,PER])

    }

    Y <- do.call(cbind, per_list)
    Y <- array(Y, dim=c(dim(per_list[[1]]), length(per_list)))

    Y.measure <- apply(Y, c(1, 2), measure, na.rm = TRUE)

    roc_list[[PER]] <- Y.measure

  }

  roc <- do.call(cbind, roc_list)
  roc <- array(roc, dim=c(dim(roc_list[[1]]), length(roc_list)))

  # Find which combination of RW parameters had the highest auroc
  ROCwh <- which(roc==max(roc), arr.ind=T)
  ROCwh
  roc[ROCwh]

  opt_params <- data.frame(
    rw_slp_opt = rwslp_vec[ROCwh[1]],
    rw_exp_opt = rwexp_vec[ROCwh[2]],
    rw_per_opt = rwper_vec[ROCwh[3]],
    rw_auroc = roc[ROCwh]
  )

  opt_params

}




ggplot(sets, aes(x=per, y=exp)) +
  geom_point(alpha=0.7, aes(colour = median_auroc, size = rel_freq)) +
  scale_size(range = c(2, 10), name="Relative\nfrequency (%)",
             breaks = c(5, 15, 60)) +
  scale_colour_gradient(low = "#1B4F72", high = "#85C1E9",
                        name = "Median AUROC") +
  scale_x_continuous(expression(paste("Persistence factor")),
                     limits = c(min(rwper_vec), max = max(rwper_vec))) +
  scale_y_continuous(expression(paste("Exponent of divergence")),
                     limits = c(min(rwexp_vec), max = max(rwexp_vec)+.1)) +
  theme_light() +
  theme(text = element_text(family = "Arial", size = 8), axis.title = element_text(size = 9),
        axis.text = element_text(size = 8))

