library(rgeos)
library(sperrorest)
library(raster)

# LABEL SP CROSS VALIDATION #############################################
# Create a cross validation of the results

setwd("/home/jason/Data/Chile/")
# elevation model
dem <- raster("dem_alos_12_5m _no sinks.tif")

# slide start/source points
slide_point_vec <- readOGR(".", "dflow_points_v1_reposition")

# actual/mapped debris flow polygons
slide_poly_vec <- readOGR(".", "dflow_polygons_v1_reposition_sample_100")
slide_poly_vec$objectid <- 1:100
crs(slide_point_vec) <- crs(slide_poly_vec)

# convert polygon to single point
centroid_poly <- gCentroid(slide_poly_vec, byid=TRUE, id = slide_poly_vec$objectid)
plot(slide_poly_vec)
plot(centroid_poly, add = TRUE)

# Use kmeans clustering to spatially cluster data
set.seed(11082017)

poly_crds <- coordinates(centroid_poly)
n_folds <- 5
sp_resamp <- partition_kmeans(poly_crds, nfold = n_folds) # spatial cv
nsp_resamp <- partition_cv(poly_crds, nfold = n_folds) # non-spatial cv

plot(sp_resamp, poly_crds)

cv_labels <- data.frame(objectid = slide_poly_vec$objectid, spcv_lab = NA,
                        cv_lab = NA)
# Spatial cv labels
for(i in 1:n_folds){
  test_label <- sp_resamp[[1]][[i]]$test
  cv_labels$spcv_lab[test_label] <- i
}

summary(as.factor(cv_labels$spcv_lab))

# (Non-spatial) cv labels
for(i in 1:n_folds){
  test_label <- nsp_resamp[[1]][[i]]$test
  cv_labels$cv_lab[test_label] <- i
}

summary(as.factor(cv_labels$cv_lab))

setwd("/home/jason/Scratch/GPP_RW_Paper")
slide_poly_vec$spcv_label <- cv_labels$spcv_lab
slide_poly_vec$cv_label <- cv_labels$cv_lab

writeOGR(slide_poly_vec, dsn = ".",
         layer = "dflow_polygons_v1_reposition_sample_100_CV",
         driver = "ESRI Shapefile", overwrite_layer = TRUE
)


# GET RW OPTIMAL PARAMETERS ############################################
setwd("/home/jason/Scratch/GPP_RW_Paper")
# Load saga gpp random walk settings
load("gridsearch_rw_settings.Rd")

rwexp_vec <- rw_settings$vec_rwexp
rwper_vec <- rw_settings$vec_rwper
rwslp_vec <- rw_settings$vec_rwslp

polyid.vec <- 1:rw_settings$n_train
#polyid.vec <- 1:10

# GLOBAL PARAMETERS ###########################################

roc_median <- list()
roc_iqr <- list()
roc_mean <- list()
roc_sd <- list()

for(PER in 1:length(rwper_vec)){

  per_list <- list()

  for(i in 1:length(polyid.vec)){

    res_nm <- paste("result_rw_roc_", i, ".Rd", sep="")
    load(res_nm) #res

    # res is the result of different parameters tested for a single
    # landslide

    res <- roc_result
    per_list[[i]] <- as.matrix(res[,,PER])

  }

  Y <- do.call(cbind, per_list)
  Y <- array(Y, dim=c(dim(per_list[[1]]), length(per_list)))

  Y.median <- apply(Y, c(1, 2), median, na.rm = TRUE)
  Y.iqr<- apply(Y, c(1, 2), IQR, na.rm = TRUE)
  Y.mean <- apply(Y, c(1, 2), mean, na.rm = TRUE)
  Y.sd<- apply(Y, c(1, 2), sd, na.rm = TRUE)

  roc_median[[PER]] <- Y.median
  roc_iqr[[PER]] <- Y.iqr
  roc_mean[[PER]] <- Y.mean
  roc_sd[[PER]] <- Y.sd

}

rocMedian <- do.call(cbind, roc_median)
rocMedian <- array(rocMedian, dim=c(dim(roc_median[[1]]), length(roc_median)))

rocIQR <- do.call(cbind, roc_iqr)
rocIQR <- array(rocIQR, dim=c(dim(roc_iqr[[1]]), length(roc_iqr)))

rocMean <- do.call(cbind, roc_mean)
rocMean <- array(rocMean, dim=c(dim(roc_mean[[1]]), length(roc_mean)))

rocSD <- do.call(cbind, roc_sd)
rocSD <- array(rocSD, dim=c(dim(roc_sd[[1]]), length(roc_sd)))


# Find which combination of RW parameters had the highest auroc
ROCwh <- which(rocMedian==max(rocMedian), arr.ind=T)
ROCwh

#train_auc <- rocMedian[ROCwh]

global_rw <- data.frame(
  median_auc = rocMedian[ROCwh],
  iqr_auc = rocIQR[ROCwh],
  mean_auc = rocMean[ROCwh],
  sd_auc = rocSD[ROCwh],
  slp_opt_median = rwslp_vec[ROCwh[1]],
  exp_opt_median = rwexp_vec[ROCwh[2]],
  per_opt_median = rwper_vec[ROCwh[3]]
)

global_w

# CROSS VALIDATED PARAMETERS ###########################################
spcv_results <- list()

for(k in 1:n_folds){
  # Aggregate results for each landslide by median
  roc_median <- list()
  roc_iqr <- list()
  roc_mean <- list()
  roc_sd <- list()

  polyid_train.vec <- which(cv_labels$spcv_lab != k)
  polyid_test.vec <- which(cv_labels$spcv_lab == k)

  for(PER in 1:length(rwper_vec)){

    per_list <- list()

    for(i in 1:length(polyid_train.vec)){
      obj.id <- polyid_train.vec[i]
      res_nm <- paste("result_rw_roc_", obj.id, ".Rd", sep="")
      load(res_nm) #res

      # res is the result of different parameters tested for a single
      # landslide

      res <- roc_result
      per_list[[i]] <- as.matrix(res[,,PER])

    }

    Y <- do.call(cbind, per_list)
    Y <- array(Y, dim=c(dim(per_list[[1]]), length(per_list)))

    Y.median <- apply(Y, c(1, 2), median, na.rm = TRUE)
    Y.iqr<- apply(Y, c(1, 2), IQR, na.rm = TRUE)
    Y.mean <- apply(Y, c(1, 2), mean, na.rm = TRUE)
    Y.sd<- apply(Y, c(1, 2), sd, na.rm = TRUE)

    roc_median[[PER]] <- Y.median
    roc_iqr[[PER]] <- Y.iqr
    roc_mean[[PER]] <- Y.mean
    roc_sd[[PER]] <- Y.sd

  }

  rocMedian <- do.call(cbind, roc_median)
  rocMedian <- array(rocMedian, dim=c(dim(roc_median[[1]]), length(roc_median)))

  rocIQR <- do.call(cbind, roc_iqr)
  rocIQR <- array(rocIQR, dim=c(dim(roc_iqr[[1]]), length(roc_iqr)))

  rocMean <- do.call(cbind, roc_mean)
  rocMean <- array(rocMean, dim=c(dim(roc_mean[[1]]), length(roc_mean)))

  rocSD <- do.call(cbind, roc_sd)
  rocSD <- array(rocSD, dim=c(dim(roc_sd[[1]]), length(roc_sd)))


  # Find which combination of RW parameters had the highest auroc
  ROCwh <- which(rocMedian==max(rocMedian), arr.ind=T)
  ROCwh

  #train_auc <- rocMedian[ROCwh]

  train_auc_median <- rocMedian[ROCwh]
  train_auc_iqr <- rocIQR[ROCwh]
  train_auc_mean <- rocMean[ROCwh]
  train_auc_sd <- rocSD[ROCwh]

  slp_opt <- rwslp_vec[ROCwh[1]]
  exp_opt <- rwexp_vec[ROCwh[2]]
  per_opt <- rwper_vec[ROCwh[3]]

  # to get test results
  roc_median <- list()
  roc_iqr <- list()
  roc_mean <- list()
  roc_sd <- list()

  for(PER in 1:length(rwper_vec)){

    per_list <- list()

    for(i in 1:length(polyid_test.vec)){
      obj.id <- polyid_test.vec[i]
      res_nm <- paste("result_rw_roc_", obj.id, ".Rd", sep="")
      load(res_nm) #res
      res <- roc_result
      per_list[[i]] <- as.matrix(res[,,PER])

    }

    Y <- do.call(cbind, per_list)
    Y <- array(Y, dim=c(dim(per_list[[1]]), length(per_list)))

    Y.median <- apply(Y, c(1, 2), median, na.rm = TRUE)
    Y.iqr<- apply(Y, c(1, 2), IQR, na.rm = TRUE)
    Y.mean <- apply(Y, c(1, 2), mean, na.rm = TRUE)
    Y.sd<- apply(Y, c(1, 2), sd, na.rm = TRUE)

    roc_median[[PER]] <- Y.median
    roc_iqr[[PER]] <- Y.iqr
    roc_mean[[PER]] <- Y.mean
    roc_sd[[PER]] <- Y.sd

  }

  test_rocMedian <- do.call(cbind, roc_median)
  test_rocMedian <- array(test_rocMedian, dim=c(dim(roc_median[[1]]), length(roc_median)))

  test_rocIQR <- do.call(cbind, roc_iqr)
  test_rocIQR <- array(test_rocIQR, dim=c(dim(roc_iqr[[1]]), length(roc_iqr)))

  test_rocMean <- do.call(cbind, roc_mean)
  test_rocMean <- array(test_rocMean, dim=c(dim(roc_mean[[1]]), length(roc_mean)))

  test_rocSD <- do.call(cbind, roc_sd)
  test_rocSD <- array(test_rocSD, dim=c(dim(roc_sd[[1]]), length(roc_sd)))

  #test_roc <- do.call(cbind, roc_median)
  #test_roc<- array(test_roc, dim=c(dim(roc_median[[1]]), length(roc_median)))

  test_auc_median <- test_rocMedian[ROCwh]
  test_auc_iqr <- test_rocIQR[ROCwh]
  test_auc_mean <- test_rocMean[ROCwh]
  test_auc_sd <- test_rocSD[ROCwh]

  # result for test dataset
  result_list <- list(
    n_train = length(polyid_train.vec),
    n_test = length(polyid_test.vec),
    rw_slp_opt = slp_opt,
    rw_exp_opt = exp_opt,
    rw_per_opt = per_opt,
    test_auroc_median = test_auc_median,
    test_auroc_iqr = test_auc_iqr,
    test_auroc_mean = test_auc_mean,
    test_auroc_sd = test_auc_sd,
    train_auroc_median = train_auc_median,
    train_auroc_iqr = train_auc_iqr,
    train_auroc_mean = train_auc_mean,
    train_auroc_sd = train_auc_sd
  )

  spcv_results[[k]] <- result_list

}

spcv_results

spcv_fold_results <- data.frame(matrix(unlist(spcv_results), nrow = length(spcv_results), byrow = T))
names(spcv_fold_results) <- names(spcv_results[[1]])
spcv_fold_results

spcv_auroc_summary <- data.frame(train_median = median(spcv_fold_results$train_auroc_median),
                                 train_IQR = IQR(spcv_fold_results$train_auroc_median),
                                 train_mean = mean(spcv_fold_results$train_auroc_median),
                                 train_sd = sd(spcv_fold_results$train_auroc_median),

                                 test_median = median(spcv_fold_results$test_auroc_median),
                                 test_IQR = IQR(spcv_fold_results$test_auroc_median),
                                 test_mean = mean(spcv_fold_results$test_auroc_median),
                                 test_sd = sd(spcv_fold_results$test_auroc_median))

spcv_auroc_summary

spcv_rw_folds <- spcv_fold_results
# GET PCM OPTIMAL PARAMETERS #######################################

# This is results of 2 parameter optimization of PCM with MD and Mu

# GLOBAL PCM OPT PARAM #############################################
setwd("/home/jason/Scratch/GPP_PCM_Paper")
setwd("D:/Scratch/GPP_PCM")

(load("gridsearch_pcm_settings.Rd"))

polyid.vec <- 1:pcm_settings$n_train
pcmmu.vec <- pcm_settings$vec_pcmmu
pcmmd.vec <- pcm_settings$vec_pcmmd

relerr_list <- list()

for(i in 1:pcm_settings$n_train){

  res_nm <- paste("result_relerr_length_", i, ".Rd", sep="")
  load(res_nm) #res
  res <- relerr_length_result
  relerr_list[[i]] <- res
  # calc median for these using apply

}


rel_err <- apply(simplify2array(relerr_list), 1:2, median)
iqr_relerr<- apply(simplify2array(train_relerr_list), 1:2, IQR)


rel_err_wh <- which(rel_err==min(rel_err), arr.ind=T)
rel_err_wh
rel_err[rel_err_wh]

opt_md <- pcmmd.vec[rel_err_wh[2]]
opt_mu <- pcmmu.vec[rel_err_wh[1]]

global_pcm <- list(
  pcm_mu = opt_mu,
  pcm_md = opt_md,
  rel_err = rel_err[rel_err_wh],
  iqr_rel_err = iqr_relerr[rel_err_wh]
)

#save(opt_gpp_par, file = "optimal_gpp_parameters.Rd")


# SPCV PCM OPT PARAM #############################################
setwd("/home/jason/Scratch/GPP_PCM_Paper")
setwd("D:/Scratch/GPP_PCM")

(load("gridsearch_pcm_settings.Rd"))

polyid.vec <- 1:pcm_settings$n_train
pcmmu.vec <- pcm_settings$vec_pcmmu
pcmmd.vec <- pcm_settings$vec_pcmmd

spcv_pcm <- list()

for(k in 1:n_folds){

  train_relerr_list <- list()

  polyid_train.vec <- which(cv_labels$spcv_lab != k)
  polyid_test.vec <- which(cv_labels$spcv_lab == k)


  for(i in 1:length(polyid_train.vec)){
    obj.id <- polyid_train.vec[i]
    res_nm <- paste("result_relerr_length_", obj.id, ".Rd", sep="")
    load(res_nm) #res
    res <- relerr_length_result
    train_relerr_list[[i]] <- res
    # calc median for these using apply

  }


  rel_err <- apply(simplify2array(train_relerr_list), 1:2, median)
  iqr_relerr<- apply(simplify2array(train_relerr_list), 1:2, IQR)

  rel_err_wh <- which(rel_err==min(rel_err), arr.ind=T)
  rel_err_wh
  rel_err[rel_err_wh]

  opt_md <- pcmmd.vec[rel_err_wh[2]]
  opt_mu <- pcmmu.vec[rel_err_wh[1]]

  # Test
  test_relerr_list <- list()

  for(i in 1:length(polyid_test.vec)){
    obj.id <- polyid_test.vec[i]
    res_nm <- paste("result_relerr_length_", obj.id, ".Rd", sep="")
    load(res_nm) #res
    res <- relerr_length_result
    test_relerr_list[[i]] <- res
    # calc median for these using apply

  }


  test_rel_err <- apply(simplify2array(test_relerr_list), 1:2, median)
  test_rel_err[rel_err_wh]
  test_iqr_relerr<- apply(simplify2array(test_relerr_list), 1:2, IQR)


  opt_gpp_par <- list(
    pcm_mu = opt_mu,
    pcm_md = opt_md,
    train_relerr = rel_err[rel_err_wh],
    train_relerr_iqr = iqr_relerr[rel_err_wh],
    test_relerr = test_rel_err[rel_err_wh],
    test_relerr_iqr = test_iqr_relerr[rel_err_wh]
  )

  spcv_pcm[[k]] <- opt_gpp_par

}

spcv_pcm

spcv_pcm_folds <- data.frame(matrix(unlist(spcv_pcm), nrow = length(spcv_pcm), byrow = T))
names(spcv_pcm_folds) <- names(spcv_pcm[[1]])

spcv_rw_folds
spcv_pcm_folds

global_rw
global_pcm
