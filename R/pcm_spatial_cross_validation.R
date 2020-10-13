

#' Perform (repeated) spatial cross-validation of PCM optimization
#'
#' @param x A SpatialPolygonsDataframe
#' @param n_folds The number of cross-validated folds (i.e. partitions)
#' @param repetitions Number of cross-validation repetitions
#' @param pcm_mu_v A vector of PCM model sliding friction coefficients
#' @param pcm_md_v A vector of PCM model mass-to-drag ratios (m)
#' @return A vector containing numeric labels defining each fold

pcmSPCV <- function(slide_plys, n_folds, repetitions, pcm_mu_v, pcm_md_v){

  rep_spcv_pcm <- list()

  for(n in 1:repetitions){

    rep_seed <- 11082017 + n

    fold_labels <- spcvFoldsPoly(slide_plys, n_folds = n_folds, seed_cv = rep_seed)

    spcv_results <- list()

    spcv_pcm <- list()
    # Folds ###############
    for(k in 1:n_folds){

      train_relerr_list <- list()

      polyid_train.vec <- which(fold_labels != k)
      polyid_test.vec <- which(fold_labels == k)

      for(i in 1:length(polyid_train.vec)){
        obj.id <- polyid_train.vec[i]
        res_nm <- paste("result_pcm_gridsearch_", obj.id, ".Rd", sep="")
        load(res_nm) #res
        res <- result_pcm[['relerr']]
        train_relerr_list[[i]] <- res
        # calc median for these using apply

      }


      rel_err <- apply(simplify2array(train_relerr_list), 1:2, median)
      iqr_relerr<- apply(simplify2array(train_relerr_list), 1:2, IQR)

      rel_err_wh <- which(rel_err==min(rel_err), arr.ind=T)
      rel_err_wh # an arrayInd() ...
      rel_err[rel_err_wh]

      #In case two optimal params found, take first
      if(length(rel_err_wh) > 2){
        rel_err_wh <- rel_err_wh[1,]
      }

      opt_md <- pcm_md_v[rel_err_wh[2]]
      opt_mu <- pcm_mu_v[rel_err_wh[1]]

      # Test
      test_relerr_list <- list()

      for(i in 1:length(polyid_test.vec)){
        obj.id <- polyid_test.vec[i]
        res_nm <- paste("result_pcm_gridsearch_", obj.id, ".Rd", sep="")
        load(res_nm) #res
        res <- result_pcm[['relerr']]
        test_relerr_list[[i]] <- res
        # calc median for these using apply

      }


      test_rel_err <- apply(simplify2array(test_relerr_list), 1:2, median)
      test_rel_err[rel_err_wh[1], rel_err_wh[2]]
      test_iqr_relerr<- apply(simplify2array(test_relerr_list), 1:2, IQR)


      opt_gpp_par <- list(
        n_train = length(polyid_train.vec),
        n_test = length(polyid_test.vec),
        pcm_mu = opt_mu,
        pcm_md = opt_md,
        train_relerr = rel_err[rel_err_wh[1], rel_err_wh[2]],
        test_relerr = test_rel_err[rel_err_wh[1], rel_err_wh[2]]
      )

      spcv_pcm[[k]] <- opt_gpp_par

    }

    spcv_pcm_folds <- data.frame(matrix(unlist(spcv_pcm), nrow = length(spcv_pcm), byrow = T))
    names(spcv_pcm_folds) <- names(spcv_pcm[[1]])

    rep_spcv_pcm[[n]] <- spcv_pcm_folds

  }

  rep_spcv_pcm

}





#' Pool PCM spatial cross-validation results
#'
#' Creates a summary of the frequency and performance (relative error) of optimal parameter
#'     sets across spatial cross-validation repetitions.
#' @param x The (list) result of pcmSPCV()
#' @return A data frame summarizing the frequency and performance of each
#'     optimal parameter sets


pcmPoolSPCV<- function(x){

  pool_pcm <- do.call(rbind, x)

  # summarize by frequency and median score for optimal parameter sets
  opt_sets <- data.frame(mu = pool_pcm$pcm_mu, md = pool_pcm$pcm_md)
  freq_sets <- table(opt_sets)

  md_nm <- as.numeric(attributes(freq_sets)$dimnames$md)
  mu_nm <- as.numeric(attributes(freq_sets)$dimnames$mu)

  # Get array index of pairs
  set_ind <- which(freq_sets !=0, arr.ind = TRUE)
  sets <- data.frame(mu = mu_nm[set_ind[,1]],
                     md = md_nm[set_ind[,2]],
                     freq = freq_sets[set_ind])

  sets$rel_freq <- sets$freq/nrow(pool_pcm) * 100

  #Find median AUROC values
  set_id <- paste(sets$mu, sets$md)
  pool_pcm$set_id <- paste(pool_pcm$pcm_mu, pool_pcm$pcm_md)

  for(i in 1:length(set_id)){

    relerr_i <- pool_pcm$test_relerr[which(pool_pcm$set_id == set_id[i])]
    sets$rel_err[i] <- median(relerr_i)
    sets$iqr_relerr[i] <- IQR(relerr_i, na.rm = TRUE)

  }

  return(sets)

}




