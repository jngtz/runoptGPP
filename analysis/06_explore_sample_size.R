
# LOAD PACKAGES AND DATA #######################################################

library(runoptGPP)
library(rgeos)
library(sperrorest)
library(raster)
library(rgdal)

setwd("/home/jason/Data/Chile/")
# elevation model
dem <- raster("elev_alos_12_5m_no_sinks.tif")

# slide start/source points
slide_point_vec <- readOGR(".", "debris_flow_source_points")

# actual/mapped debris flow polygons
slide_poly_vec <- readOGR(".", "debris_flow_polys_sample")
slide_poly_vec$objectid <- 1:100
crs(slide_point_vec) <- crs(slide_poly_vec)


# SPATIAL CROSS VALIDATION SETTINGS ###########################

repetitions <- 100
n_folds <- 5
n_smps <- seq(from = 10, to = 80, by = 10)


# RW SAMPLE SIZE ###############################################################

n_smps <- seq(from = 10, to = 80, by = 10)

setwd("/home/jason/Scratch/GPP_RW_Paper")
(load("rw_gridsearch_multi.Rd"))

rwslp_vec <- as.numeric(dimnames(rw_gridsearch_multi[[1]])[[1]])
rwexp_vec <- as.numeric(dimnames(rw_gridsearch_multi[[1]])[[2]])
rwper_vec <- as.numeric(dimnames(rw_gridsearch_multi[[1]])[[3]])


smp_sizes_rw <- vector("list", length = length(n_smps))


for(SMP in 1:length(n_smps)){

  n_smp <- n_smps[SMP]

  rep_spcv_rw <- vector("list", length = repetitions)

  for(n in 1:repetitions){

    print(paste("Sample size:", n_smp, "Rep:", n))

    rep_seed <- 11082017 + n

    # Random sample of slides
    #rndm_smp <- sample(length(polyid.vec), n_smp)
    #smp_slide_poly_vec <- slide_poly_vec[rndm_smp,]
    #plot(smp_slide_poly_vec)

    # Define folds for sample
    fold_labels <- spcvFoldsPoly(slide_poly_vec, n_folds = n_folds, seed = rep_seed)

    spcv_results <- list()



    # Folds ###############
    for(k in 1:n_folds){

      roc_median <- vector("list", length = n_folds)
      roc_iqr <- vector("list", length = n_folds)
      #roc_mean <- vector("list", length = n_folds)
      #roc_sd <- vector("list", length = n_folds)

      polyid_train.vec <- which(fold_labels != k)


      if(n_smp < length(polyid_train.vec)){

        polyid_train.vec <- sample(polyid_train.vec, n_smp, replace = FALSE)

      }

      #print(length(polyid_train.vec))




      polyid_test.vec <- which(fold_labels == k)

      for(PER in 1:length(rwper_vec)){

        per_list <- list()

        for(i in 1:length(polyid_train.vec)){
          obj.id <- polyid_train.vec[i]
          res <- rw_gridsearch_multi[[obj.id]]
          per_list[[i]] <- as.matrix(res[,,PER])

        }

        Y <- do.call(cbind, per_list)
        Y <- array(Y, dim=c(dim(per_list[[1]]), length(per_list)))

        Y.median <- apply(Y, c(1, 2), median, na.rm = TRUE)
        Y.iqr<- apply(Y, c(1, 2), IQR, na.rm = TRUE)
        #Y.mean <- apply(Y, c(1, 2), mean, na.rm = TRUE)
        #Y.sd<- apply(Y, c(1, 2), sd, na.rm = TRUE)

        roc_median[[PER]] <- Y.median
        roc_iqr[[PER]] <- Y.iqr
        #roc_mean[[PER]] <- Y.mean
        #roc_sd[[PER]] <- Y.sd

      }

      rocMedian <- do.call(cbind, roc_median)
      rocMedian <- array(rocMedian, dim=c(dim(roc_median[[1]]), length(roc_median)))

      rocIQR <- do.call(cbind, roc_iqr)
      rocIQR <- array(rocIQR, dim=c(dim(roc_iqr[[1]]), length(roc_iqr)))

      #rocMean <- do.call(cbind, roc_mean)
      #ocMean <- array(rocMean, dim=c(dim(roc_mean[[1]]), length(roc_mean)))

      #rocSD <- do.call(cbind, roc_sd)
      #rocSD <- array(rocSD, dim=c(dim(roc_sd[[1]]), length(roc_sd)))


      # Find which combination of RW parameters had the highest auroc
      ROCwh <- which(rocMedian==max(rocMedian), arr.ind=T)
      ROCwh

      #train_auc <- rocMedian[ROCwh]

      train_auc_median <- rocMedian[ROCwh]
      train_auc_iqr <- rocIQR[ROCwh]
      #train_auc_mean <- rocMean[ROCwh]
      #train_auc_sd <- rocSD[ROCwh]

      slp_opt <- rwslp_vec[ROCwh[1]]
      exp_opt <- rwexp_vec[ROCwh[2]]
      per_opt <- rwper_vec[ROCwh[3]]

      # to get test results
      #roc_median <- list()
      #roc_iqr <- list()
      #roc_mean <- list()
      #roc_sd <- list()

      for(PER in 1:length(rwper_vec)){

        per_list <- list()

        for(i in 1:length(polyid_test.vec)){
          obj.id <- polyid_test.vec[i]
          res <- rw_gridsearch_multi[[obj.id]]
          per_list[[i]] <- as.matrix(res[,,PER])

        }

        Y <- do.call(cbind, per_list)
        Y <- array(Y, dim=c(dim(per_list[[1]]), length(per_list)))

        Y.median <- apply(Y, c(1, 2), median, na.rm = TRUE)
        Y.iqr<- apply(Y, c(1, 2), IQR, na.rm = TRUE)
        #Y.mean <- apply(Y, c(1, 2), mean, na.rm = TRUE)
        #Y.sd<- apply(Y, c(1, 2), sd, na.rm = TRUE)

        roc_median[[PER]] <- Y.median
        roc_iqr[[PER]] <- Y.iqr
        #roc_mean[[PER]] <- Y.mean
        #roc_sd[[PER]] <- Y.sd

      }

      test_rocMedian <- do.call(cbind, roc_median)
      test_rocMedian <- array(test_rocMedian, dim=c(dim(roc_median[[1]]), length(roc_median)))

      test_rocIQR <- do.call(cbind, roc_iqr)
      test_rocIQR <- array(test_rocIQR, dim=c(dim(roc_iqr[[1]]), length(roc_iqr)))

      #test_rocMean <- do.call(cbind, roc_mean)
      #test_rocMean <- array(test_rocMean, dim=c(dim(roc_mean[[1]]), length(roc_mean)))

      #test_rocSD <- do.call(cbind, roc_sd)
      #test_rocSD <- array(test_rocSD, dim=c(dim(roc_sd[[1]]), length(roc_sd)))

      #test_roc <- do.call(cbind, roc_median)
      #test_roc<- array(test_roc, dim=c(dim(roc_median[[1]]), length(roc_median)))

      test_auc_median <- test_rocMedian[ROCwh]
      test_auc_iqr <- test_rocIQR[ROCwh]
      #test_auc_mean <- test_rocMean[ROCwh]
      #test_auc_sd <- test_rocSD[ROCwh]

      # result for test dataset
      result_list <- list(
        n_train = length(polyid_train.vec),
        n_test = length(polyid_test.vec),
        rw_slp_opt = slp_opt,
        rw_exp_opt = exp_opt,
        rw_per_opt = per_opt,
        test_auroc_median = test_auc_median,
        train_auroc_median = train_auc_median
      )

      spcv_results[[k]] <- result_list

    }


    spcv_rw_folds <- data.frame(matrix(unlist(spcv_results), nrow = length(spcv_results), byrow = T))
    names(spcv_rw_folds) <- names(spcv_results[[1]])

    rep_spcv_rw[[n]] <- spcv_rw_folds

  }

  smp_sizes_rw[[SMP]] <- rep_spcv_rw

}

#save(smp_sizes_rw, file = "smp_size_repeated_spcv_RW.Rd")

# PCM SAMPLE SIZE ##############################################################

setwd("/home/jason/Scratch/GPP_PCM_Paper")
(load("pcm_gridsearch_multi.Rd"))

pcmmd.vec <- as.numeric(colnames(pcm_gridsearch_multi[[1]][[1]]))
pcmmu.vec <- as.numeric(rownames(pcm_gridsearch_multi[[1]][[1]]))

spcv_pcm <- list()
smp_sizes_pcm <- vector("list", length = length(n_smps))


for(SMP in 1:length(n_smps)){

  n_smp <- n_smps[SMP]

  rep_spcv_pcm <- vector("list", length = length(n_smps))

  for(n in 1:repetitions){

    print(paste("Sample size:", n_smp, "Rep:", n))

    rep_seed <- 11082017 + n

    # Define folds for sample
    fold_labels <- spcvFoldsPoly(slide_poly_vec, n_folds = n_folds, seed = rep_seed)


    # Folds ###############
    for(k in 1:n_folds){

      train_relerr_list <- vector("list", length = n_folds)

      polyid_train.vec <- which(fold_labels != k)

      if(n_smp < length(polyid_train.vec)){

        polyid_train.vec <- sample(polyid_train.vec, n_smp, replace = FALSE)

      }

      polyid_test.vec <- which(fold_labels == k)

      for(i in 1:length(polyid_train.vec)){
        obj.id <- polyid_train.vec[i]
        res <- pcm_gridsearch_multi[[obj.id]]$relerr
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

      opt_md <- pcmmd.vec[rel_err_wh[2]]
      opt_mu <- pcmmu.vec[rel_err_wh[1]]

      # Test
      test_relerr_list <- vector("list", length = length(n_folds))

      for(i in 1:length(polyid_test.vec)){
        obj.id <- polyid_test.vec[i]
        res <- pcm_gridsearch_multi[[obj.id]]$relerr
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

  smp_sizes_pcm[[SMP]] <- rep_spcv_pcm
}


#save(smp_sizes_pcm, file = "smp_size_repeated_spcv_PCM.Rd")

# SUMMARY RW RESULTS ######################

setwd("/home/jason/Scratch/GPP_RW_Paper")
(load("smp_size_repeated_spcv_RW.Rd"))

smp_size_rw_l <- list()

for(SMP in 1:length(n_smps)){

  mean_rep <- rep(NA, repetitions)

  for(k in 1:repetitions){
    # summarize for each repetition
    mean_rep[k] <- mean(smp_sizes_rw[[SMP]][[k]]$test_auroc_median)
  }

  smp_size_rw_l[[SMP]] <- mean_rep

}

mean_auc <- sapply(smp_size_rw_l, FUN = mean)
sd_auc <- sapply(smp_size_rw_l, FUN = sd)
plot(n_smps, mean_auc,
     xlab = "Sample size", ylab = "AUROC",
     ylim = c(0.91, 0.94), type = "l")


epsilon = 0.3
for(i in 1:length(n_smps)) {
  up = mean_auc[i] + sd_auc[i]
  low = mean_auc[i] - sd_auc[i]
  segments(n_smps[i],low , n_smps[i], up)
  segments(n_smps[i]-epsilon, up , n_smps[i]+epsilon, up)
  segments(n_smps[i]-epsilon, low , n_smps[i]+epsilon, low)
}

points(n_smps, mean_auc, pch = 16)



# SUMMARY PCM RESULTS ####################

setwd("/home/jason/Scratch/GPP_PCM_Paper")

(load("smp_size_repeated_spcv_PCM.Rd"))

n_smps <- seq(from = 10, to = 80, by = 10)

smp_size_pcm_l <- list()
smp_size_pcm_sd_l <- list()

# test summarized mean k and mean rep.

for(SMP in 1:length(smp_sizes_pcm)){

  mean_rep <- rep(NA, length(smp_sizes_pcm[[1]]))

  for(k in 1:length(smp_sizes_pcm[[1]])){
    # summarize for each repetition
    mean_rep[k] <- mean(smp_sizes_pcm[[SMP]][[k]]$test_relerr)
  }

  smp_size_pcm_l[[SMP]] <- mean_rep

}

mean_relerr <- sapply(smp_size_pcm_l, FUN = mean)
sd_relerr <- sapply(smp_size_pcm_l, FUN = sd)
plot(n_smps, mean_relerr,
     xlab = "Sample size", ylab = "Relative error",
     ylim = c(0.1, 0.3), type = "l")


epsilon = 0.3
for(i in 1:length(n_smps)) {
  up = mean_relerr[i] + sd_relerr[i]
  low = mean_relerr[i] - sd_relerr[i]
  segments(n_smps[i],low , n_smps[i], up)
  segments(n_smps[i]-epsilon, up , n_smps[i]+epsilon, up)
  segments(n_smps[i]-epsilon, low , n_smps[i]+epsilon, low)
}

points(n_smps, mean_relerr, pch = 16)



smp_size_pcm_df <- data.frame(smp_size = n_smps, test_median_relerr = NA,
                              test_iqr_relerr = NA, opt_mu = NA, opt_md = NA, opt_relfreq = NA,
                              opt_median_mu = NA, opt_iqr_mu = NA,
                              opt_median_md = NA, opt_iqr_md = NA)

pairs_pcm <- list()


for(SMP in 1:length(n_smps)){
  pool_pcm <- do.call(rbind, smp_sizes_pcm[[SMP]])


  smp_size_pcm_df$test_median_relerr[SMP] <- median(pool_pcm$test_relerr)
  smp_size_pcm_df$test_iqr_relerr[SMP] <- IQR(pool_pcm$test_relerr)

  # Relative frequency of parameter combinations
  opt_pairs <- data.frame(mu = pool_pcm$pcm_mu, md = pool_pcm$pcm_md)
  freq_pairs <- table(opt_pairs)
  freq_pairs != 0

  mu_nm <- as.numeric(rownames(freq_pairs))
  md_nm <- as.numeric(colnames(freq_pairs))

  # Get array index of pairs
  pair_ind <- which(freq_pairs !=0, arr.ind = TRUE)
  pairs <- data.frame(mu = mu_nm[pair_ind[,1]],
                      md = md_nm[pair_ind[,2]],
                      freq = freq_pairs[pair_ind])
  pairs$rel_freq <- pairs$freq/nrow(pool_pcm) * 100


  #Find relative error values
  pair_id <- paste(pairs$mu, pairs$md)
  pool_pcm$pair_id <- paste(pool_pcm$pcm_mu, pool_pcm$pcm_md)

  for(i in 1:length(pair_id)){

    relerr_i <- pool_pcm$test_relerr[which(pool_pcm$pair_id == pair_id[i])]
    pairs$rel_err[i] <- median(relerr_i)
    pairs$iqr_relerr[i] <- IQR(relerr_i, na.rm = TRUE)

  }

  pairs_pcm[[SMP]] <- pairs

  opt_pool <- pairs[pairs$rel_freq == max(pairs$rel_freq),]

  smp_size_pcm_df$opt_mu[SMP] <- opt_pool$mu
  smp_size_pcm_df$opt_md[SMP] <- opt_pool$md
  smp_size_pcm_df$opt_relfreq[SMP] <- opt_pool$rel_freq

  smp_size_pcm_df$opt_median_mu[SMP] <- median(pool_pcm$pcm_mu)
  smp_size_pcm_df$opt_median_md[SMP] <- median(pool_pcm$pcm_md)

  smp_size_pcm_df$opt_iqr_mu[SMP] <- IQR(pool_pcm$pcm_mu)
  smp_size_pcm_df$opt_iqr_md[SMP] <- IQR(pool_pcm$pcm_md)

}



# PLOT RW AND PCM RESULTS ######################################################

setwd("/home/jason/Scratch/Figures")

# Sample size performance and variation ###

png(filename="sampe_size_performance.png", res = 300, width = 7.5, height = 3,
    units = "in", pointsize = 11)

par(family = "Arial", mfrow = c(1,2), mar = c(4, 4, 1, 0.5),
    mgp = c(2, 0.75, 0))


#AUROC

mean_auc <- sapply(smp_size_rw_l, FUN = mean)
sd_auc <- sapply(smp_size_rw_l, FUN = sd)
plot(n_smps, mean_auc,
     xlab = "Sample size", ylab = "Runout path AUROC",
     ylim = c(0.91, 0.94), type = "l")


epsilon = 0.3
for(i in 1:length(n_smps)) {
  up = mean_auc[i] + sd_auc[i]
  low = mean_auc[i] - sd_auc[i]
  segments(n_smps[i],low , n_smps[i], up)
  segments(n_smps[i]-epsilon, up , n_smps[i]+epsilon, up)
  segments(n_smps[i]-epsilon, low , n_smps[i]+epsilon, low)
}

points(n_smps, mean_auc, pch = 16)


#Relative error
mean_relerr <- sapply(smp_size_pcm_l, FUN = mean)
sd_relerr <- sapply(smp_size_pcm_l, FUN = sd)
plot(n_smps, mean_relerr,
     xlab = "Sample size", ylab = "Runout distance relative error",
     ylim = c(0.1, 0.25), type = "l")


epsilon = 0.3
for(i in 1:length(n_smps)) {
  up = mean_relerr[i] + sd_relerr[i]
  low = mean_relerr[i] - sd_relerr[i]
  segments(n_smps[i],low , n_smps[i], up)
  segments(n_smps[i]-epsilon, up , n_smps[i]+epsilon, up)
  segments(n_smps[i]-epsilon, low , n_smps[i]+epsilon, low)
}

points(n_smps, mean_relerr, pch = 16)


dev.off()



# REL FREQ PCM PARAMS / SAMPLE SIZE SCATTER ####################################

setwd("/home/jason/Scratch/Figures")

pairs_pcm_df <- pairs_pcm[[1]]
pairs_pcm_df$smp_size <- n_smps[1]

for(i in 2:length(pairs_pcm)){
  tmp_pcm_df <- pairs_pcm[[i]]
  tmp_pcm_df$smp_size <- n_smps[i]

  pairs_pcm_df <- rbind(pairs_pcm_df, tmp_pcm_df)

}

summary(as.factor(pairs_pcm_df$smp_size))

pairs_pcm_df$labels <- as.factor(pairs_pcm_df$smp_size)

mybreaks <- c(1, 10, 50 )

ggplot(pairs_pcm_df, aes(x=md, y=mu)) +
  geom_point(alpha=0.7, aes(colour = smp_size, size = rel_freq)) +


  scale_size(name="Relative\nfrequency (%)", range = c(2, 6), breaks = mybreaks) +

  #scale_colour_gradient(high = "#1B4F72", low = "#85C1E9", name = "Sample size") +
  scale_color_viridis_c(name = "Sample size")+

  xlab("Mass-to-drag ratio (m)") +
  ylab("Sliding friction coefficient") +


  theme_light() +
  theme(text = element_text(family = "Arial", size = 8), axis.title = element_text(size = 9),
        axis.text = element_text(size = 8)) +
  guides(
    colour = guide_colourbar(order = 1),
    size = guide_legend(order = 2))

ggsave("sample_size_scatter_opt_pcm.png", dpi = 300, width = 5.5, height = 3.25, units = "in")
ggsave("sample_size_scatter_opt_pcmw7h4.png", dpi = 300, width = 7.5, height =4.5, units = "in")

