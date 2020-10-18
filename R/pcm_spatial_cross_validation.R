

#' Perform (repeated) spatial cross-validation of PCM optimization
#'
#' @param x A SpatialPolygonsDataframe
#' @param n_folds The number of folds (i.e. partitions)
#' @param repetitions Number of repetitions
#' @param from_save (Logical) if TRUE, will load save files from current working directory
#' @return A vector containing numeric labels defining each fold

pcmSPCV <- function(x, slide_plys, n_folds, repetitions, from_save = FALSE){


  if(from_save){
    x <- list()

    files <- list.files(pattern = "result_pcm_gridsearch_")

    for(i in 1:length(files)){
      res_nm <- paste("result_pcm_gridsearch_", i, ".Rd", sep="")
      load(res_nm)
      x[[i]] <- result_pcm
    }

  }

  pcm_md_vec <- as.numeric(colnames(x[[1]][[1]]))
  pcm_mu_vec <- as.numeric(rownames(x[[1]][[1]]))

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
        res <- x[[obj.id]][['relerr']]
        train_relerr_list[[i]] <- res
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

      opt_md <- pcm_md_vec[rel_err_wh[2]]
      opt_mu <- pcm_mu_vec[rel_err_wh[1]]

      # Test
      test_relerr_list <- list()

      for(i in 1:length(polyid_test.vec)){
        obj.id <- polyid_test.vec[i]
        res <- x[[obj.id]][['relerr']]
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

  rep_spcv_pcm$settings <- list(pcm_md_vec = pcm_md_vec, pcm_mu_vec = pcm_mu_vec)

  rep_spcv_pcm


}





#' Pool PCM spatial cross-validation results
#'
#' Creates a summary of the frequency and performance (relative error) of optimal parameter
#'     sets across spatial cross-validation repetitions.
#' @param x The (list) result of pcmSPCV()
#' @param plot_freq (Logical) if TRUE, will a create bubble plot of optimal parameter
#'      set frequencies across grid search space
#' @return A data frame summarizing the frequency and performance of each
#'     optimal parameter sets


pcmPoolSPCV<- function(x, plot_freq = FALSE){

  pool_pcm <- do.call(rbind, x[1:(length(x) - 1)])

  pcm_md_vec <- x$settings$pcm_md_vec
  pcm_mu_vec <- x$settings$pcm_mu_vec

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

  if(plot_freq){
    gg <-  ggplot2::ggplot(sets, ggplot2::aes(x=md, y=mu)) +
      ggplot2::geom_point(alpha=0.7, ggplot2::aes(colour = rel_err, size = rel_freq)) +

      ggplot2::scale_colour_gradient(high = "#1B4F72", low = "#85C1E9",
                            name = "Median relative\nerror") +

      ggplot2::labs(size = "Relative\nfrequency (%)") +

      ggplot2::scale_x_continuous(expression(paste("Mass-to-drag ratio (m)")),
                         limits = c(min(pcm_md_vec), max = max(pcm_md_vec))) +
      ggplot2::scale_y_continuous(expression(paste("Sliding friction coefficient")),
                         limits = c(min(pcm_mu_vec), max = max(pcm_mu_vec))) +
      ggplot2::theme_light() +
      ggplot2::theme(text = ggplot2::element_text(family = "Arial", size = 8),
                     axis.title = ggplot2::element_text(size = 9),
                     axis.text = ggplot2::element_text(size = 8))

    print(gg)
  }

  return(sets)

}


