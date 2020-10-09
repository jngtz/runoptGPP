

gridOptPCM <- function(dem, workspace = getwd(),
                       slide_plys, source_pnts, slide_id,
                       rw_slp, rw_ex, rw_per,
                       pcm_mu_v, pcm_md_v,
                       gpp_iter = 1000,
                       buffer_ext = 500, buffer_source = NULL,
                       predict_threshold = 0.5,
                       plot_eval = FALSE){

  roc_result.nm <- paste("result_roc_", slide_id, ".Rd", sep="")
  relerr_length_result.nm <- paste("result_relerr_length_", slide_id, ".Rd", sep="")
  reldiff_length_result.nm <- paste("result_reldiff_length_", slide_id, ".Rd", sep="")
  err_length_result.nm <- paste("result_err_length_", slide_id, ".Rd", sep="")

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

      result <- performancePCM(dem = dem,
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

  save(roc_result, file=roc_result.nm)
  save(relerr_length_result, file=relerr_length_result.nm)
  save(reldiff_length_result, file=reldiff_length_result.nm)
  save(err_length_result, file=err_length_result.nm)

  result_list <- list(
    auroc = roc_result,
    relerr = relerr_length_result,
    reldiff = reldiff_length_result,
    error = err_length_result)

  return(result_list)
}

