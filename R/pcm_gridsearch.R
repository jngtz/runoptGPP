#' Grid search optimization for PCM runout distance
#'
#' Computes performance measures for
#'      runout paths simuluated using the random walk and PCM model components of the
#'      GPP tool in SAGA-GIS.
#' @param dem A DEM as a RasterLayer object
#' @param workspace The file path where to save performance results for each runout polygon
#' @param slide_plys Runout tracks as a SpatialPolygonsDataFrame
#' @param source_pnts Source points as a SpatialPointsDataFrame
#' @param slide_id Selects a single runout polygon from slide_plys by row
#' @param rw_slp Random walk slope threshold - below lasteral spreading is modelled
#' @param rw_ex Random walk exponent controlling lateral spread
#' @param rw_per Random walk persistence factor to weight flow direction consistency
#' @param pcm_mu_v A vector of PCM model sliding friction coefficients
#' @param pcm_md_v A vector of PCM model mass-to-drag ratios (m)
#' @param gpp_iter Number of random walk model iterations
#' @param buffer_ext Defines buffer distance (in meters) around runout polygon
#'      to crop source DEM. This helps to reduce computational time
#' @param buffer_source Can define a buffer distance (in meters) to extend source
#'      point to a source area
#' @param predict_threshold A cutoff value to define what quantile of simulated runout
#'      frequencies is the predicted runout.
#' @param plot_eval If TRUE, will plot simulated runout and runout polygon
#' @return A list with performances for each parameter combinations across grid search space.
#'      The performances measures include relative error, relative difference, error and AUROC.
#' @details Runout is either simulated from a single source point or a buffered
#'      area round the source point.

pcmGridsearch <- function(dem, workspace = getwd(),
                       slide_plys, source_pnts, slide_id,
                       rw_slp, rw_ex, rw_per,
                       pcm_mu_v, pcm_md_v,
                       gpp_iter = 1000,
                       buffer_ext = 500, buffer_source = NULL,
                       predict_threshold = 0.5,
                       plot_eval = FALSE){

  revert_wd <- getwd()
  setwd(workspace)

  column.names <- pcm_md_v
  row.names <- pcm_mu_v

  roc_result <- array(NA ,dim = c(length(pcm_mu_v),length(pcm_md_v)),dimnames = list(row.names, column.names))
  relerr_length_result <- roc_result
  reldiff_length_result <- roc_result
  err_length_result <- roc_result

  #roc[row, col]
  #roc[pcmmu, pcmmd]

  for(i in 1:length(pcm_md_v)){
    for(k in 1:length(pcm_mu_v)){

      result <- pcmPerformance(dem = dem,
                                 slide_plys = slide_plys,
                                 source_pnts = source_pnts,
                                 slide_id = slide_id,
                                 rw_slp = rw_slp,
                                 rw_ex = rw_ex,
                                 rw_per = rw_per,
                                 pcm_mu = pcm_mu_v[k],
                                 pcm_md = pcm_md_v[i],
                                 buffer_ext = buffer_ext,
                                 buffer_source = buffer_source,
                                 gpp_iter = gpp_iter,
                                 predict_threshold = predict_threshold,
                                 plot_eval = FALSE)

      roc_result[k, i] <- result$roc
      relerr_length_result[k, i] <- result$length.relerr
      reldiff_length_result[k, i] <- result$length.reldiff
      err_length_result[k, i] <- result$length.error

      #roc_result[k, i] <- paste(pcmmd[i], pcmmu[k])
    }}

  result_pcm <- list(
    auroc = roc_result,
    relerr = relerr_length_result,
    reldiff = reldiff_length_result,
    error = err_length_result)

  save(result_pcm, file = paste("result_pcm_gridsearch_", slide_id, ".Rd", sep=""))

  setwd(revert_wd)
  return(result_pcm_gridsearch)


}




#' Get PCM runout distance grid search optimal parameters
#'
#' @param workspace The file path where the performance results were saved
#' @param slide_id Selects a single runout polygon from slide_plys by row
#' @param pcm_mu_vec The vector of random walk slope thresholds used in the grid search
#' @param pcm_mu_vec The vector or random walk exponents used in the grid search
#' @param n_train Number of runout polygons processed with grid search
#' @param performace Performance measure "relerr" relative error
#' @param measure A measure (e.g. "median", "mean") to find optimal parameter across grid search space
#' @return A dataframe  with the optimal parameter set and performance

pcmGetOpt <- function(workspace = getwd(), pcm_md_vec, pcm_mu_vec, n_train,
                      performance = "relerr", measure = "median"){
  revert_wd <- getwd()
  setwd(workspace)

  err_list <- list()
  roc_list <- list()

  for(i in 1:n_train){

    res_nm <- paste("result_pcm_gridsearch_", i, ".Rd", sep="")
    load(res_nm) #res

    err_list[[i]] <- result_pcm_gridsearch[[performance]]
    roc_list[[i]] <- result_pcm_gridsearch[['auroc']]

  }

  err <- apply(simplify2array(err_list), 1:2, get(measure))
  roc <- apply(simplify2array(roc_list), 1:2, get(measure))

  err_wh <- which(err==min(err), arr.ind=T)
  err_wh
  err[err_wh]

  opt_md <- pcm_md_vec[err_wh[2]]
  opt_mu <- pcm_mu_vec[err_wh[1]]

  opt_gpp_par <- data.frame(
    pcm_mu = opt_mu,
    pcm_md = opt_md
    )

  opt_gpp_par[paste0(measure, "_", performance)] <- err[err_wh]
  opt_gpp_par[paste0(measure, "_", "auroc")] <- roc[err_wh]

  setwd(revert_wd)

  return(opt_gpp_par)

}






