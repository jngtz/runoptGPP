
setwd("/home/jason/Data/Chile/")
# elevation model
dem <- raster::raster("dem_alos_12_5m _no sinks.tif")

# slide start/source points
slide_point_vec <- rgdal::readOGR(".", "dflow_points_v1_reposition")

# actual/mapped debris flow polygons
slide_poly_vec <- rgdal::readOGR(".", "dflow_polygons_v1_reposition_sample_100")
slide_poly_vec$objectid <- 1:100

raster::crs(slide_point_vec) <- raster::crs(slide_poly_vec)

saga <- Rsagacmd::saga_gis()


performanceRndWalk(dem, slide_plys = slide_poly_vec, source_pnts = slide_point_vec,
                   slide_id = 2, slp = 33, ex = 3, per = 2,
                   gpp_iter = 1000, buffer_ext = 500, buffer_source = 50,
                   plot_eval = TRUE)



steps <- 3
rwexp_vec <- seq(1.3, 3, len=steps)
rwper_vec <- seq(1.5, 2, len=steps)
rwslp_vec <- seq(20, 40, len=steps)

workspace_dir <- "/home/jason/Scratch/F_Testing"

gridOptRndWalk(dem, workspace = workspace_dir, slide_plys = slide_poly_vec, source_pnts = slide_point_vec,
               slide_id = 2, slp_v = rwslp_vec, ex_v = rwexp_vec, per_v = rwper_vec,
               gpp_iter = 1000, buffer_ext = 500, buffer_source = 50,
               plot_eval = TRUE)

(load("result_rw_roc_2.Rd"))
roc_result


# RUN OPTIMIZATION #############################################################

# Define which runout polygons to run grid search for

polyid_vec <- 1:length(slide_poly_vec)
polyid_vec <- 1:4

# Save RW GridSearch Settings

setwd(workspace_dir)

rw_settings <- list(
  n_train = length(polyid_vec),
  vec_rwexp = rwexp_vec,
  vec_rwper = rwper_vec,
  vec_rwslp = rwslp_vec
)

save(rw_settings, file = "gridsearch_rw_settings.Rd")


# Run in parellel

library(foreach)
cl <- parallel::makeCluster(4)
doParallel::registerDoParallel(cl)

Sys.time()
results.list <-
  foreach(poly_id=polyid_vec, .combine='c', .packages=c('rgdal','raster', 'rgeos', 'ROCR', 'Rsagacmd', 'sf')) %dopar% {
    gridOptRndWalk(dem, workspace = "/home/jason/Scratch/F_Testing", slide_plys = slide_poly_vec, source_pnts = slide_point_vec,
                   slide_id = poly_id, slp_v = rwslp_vec, ex_v = rwexp_vec, per_v = rwper_vec,
                   gpp_iter = 1000, buffer_ext = 500, buffer_source = 50,
                   plot_eval = FALSE)
  }
parallel::stopCluster(cl)
Sys.time()


# GET RW OPTIMAL VALUES #############################################################
(load("gridsearch_rw_settings.Rd"))

getRndWalkOptParams(workspace = getwd(),
                    rwslp_vec = rw_settings$vec_rwslp,
                    rwexp_vec = rw_settings$vec_rwexp,
                    rwper_vec = rw_settings$vec_rwper,
                    n_train = rw_settings$n_train,
                    measure = median)

