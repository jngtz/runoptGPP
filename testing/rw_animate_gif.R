

# Load Packages and Data ####################################################################
library(runoptGPP)
library(raster)
library(rgdal)
library(sp)
library(Rsagacmd)



# Functions #########################################################



# Load data ########################################################

# Initiate a SAGA-GIS geoprocessing object
saga <- saga_gis(opt_lib = "sim_geomorphology")
#saga <- saga_gis(saga_bin = "D:/Software/saga-6.4.0_x64/saga_cmd.exe", temp_path = 'D:/Temp/SAGAtmp', opt_lib = "sim_geomorphology")

# Set workspace
setwd("/home/jason/Data/Chile/")
#setwd("D:\\JasonGoetz\\Research\\Chile\\R-Project\\Data")

# Load digital elevation model (DEM)
dem <- raster("elev_alos_12_5m.tif")
#dem <- raster("dem_alos_12_5m _no sinks.tif")


# Load runout source points
source_points <- readOGR(".", "debris_flow_source_points")

# Load runout track polygons and assign object ID based on row number
runout_polygons <- readOGR(".", "debris_flow_polys_sample")
runout_polygons$objectid <- 1:length(runout_polygons)

# Select a debris flow and source point for this example
runout_polygon <- runout_polygons[77,]

sel_source_point  <- over(source_points, runout_polygon)
source_point <- source_points[!is.na(sel_source_point$objectid),]

setwd("saga_data")
hillshade <- raster("hs_filtered_filled.sdat")

setwd("/home/jason/Scratch/Figures")


# GPP random walk simulation ################################

test <- rwPerformance(dem, slide_plys = runout_polygon, slide_src = source_point,
              slp = 30, ex = 3, per = 2,
              gpp_iter = 1000, buffer_ext = 500, buffer_source = NULL,
              plot_eval = TRUE, saga_lib = saga)




crop_dem <- crop(dem, extent(runout_polygon)+500)
crop_hs <- crop(hillshade, crop_dem)
plot(crop_dem)

plot(crop_hs,
     col=grey.colors(100, start=1, end=0),
     legend=F,
     cex.axis = 0.7, cex.main = 0.8)
plot(runout_polygon, add = TRUE)
plot(source_point, add = TRUE)

# Custom random walk ########################################
#https://www.r-bloggers.com/2018/05/simulating-animal-movements-and-habitat-use/

start_time <- Sys.time()
dem = crop_dem
#slope = terrain(dem, opt="slope", unit = "degrees")
reps = 1000

beta_thres <- 30 # slope threshold - ! could be something weird going on here.. check rad to deg conversions...
alpha <- 3 # exponent for divergent flow # alpha >= 1
p <- 2 # persistence factor

# get initial cell center
cntr_cell <- cellFromXY(dem, xy = coordinates(source_point))
i = 0
n = 0
sim_paths <- vector(mode = "list", length = reps)
cell_pos <- 9999

# Get distance to each cell
ngh_cells <- adjacent(dem, cntr_cell, directions = 8, pairs = FALSE, id = TRUE)
ngh_points <- xyFromCell(dem, ngh_cells)
cntr_point <- xyFromCell(dem, cntr_cell)

euclideanDistance <- function(p1, p2){
  sqrt( (p1[1] - p2[1])^2 + (p1[2]- p2[2])^2 )
}

cell_dist <- apply(ngh_points, 1, euclideanDistance, p2 = cntr_point)

# Start repeated simulations

for(k in 1:reps){

  print(k)
  path_cells <- list()
  cntr_cell <- cellFromXY(dem, xy = coordinates(source_point))
  i = 0
  n = 0
  prv_pos <- 9999


  while(n == 0){

    i = i + 1

    ngh_cells <- adjacent(dem, cntr_cell, directions = 8, pairs = FALSE, id = TRUE)

    if(length(ngh_cells) < 8){
      break
    }


    elv_values <- getValues(dem)[c(ngh_cells, cntr_cell)]
    elv_ngh <- elv_values[1:8]
    elv_cntr <- elv_values[9]

    lower_elv <- elv_ngh <= elv_cntr
    # get elevations
    #elv_ngh <- extract(dem, ngh_cells)
    #elv_cntr <- extract(dem, cntr_cell)


    if(!any(lower_elv)){
      break # stop if no lower elevations
    }


    # compute slope/beta (in degrees)
    beta_ngh <-  atan( (elv_cntr - elv_ngh)  / cell_dist) * 180/pi

    #beta_ngh <- extract(slope, ngh_cells)

    if(anyNA(beta_ngh)){
      break
    }


    f <- rep(1, 8)

    if(prv_pos < 9){
      f[prv_pos] <- p
    }

    f <- f[lower_elv]
    cells <- ngh_cells[lower_elv]
    beta_ngh <- beta_ngh[lower_elv]

    gamma_i <- tan(beta_ngh*pi/180) / tan(beta_thres*pi/180)

    fj <- f * tan(beta_ngh*pi/180)

    prob <-  f*tan(beta_ngh*pi/180) / sum(fj)
    prob <- prob/sum(prob)

    gamma_max <- max(gamma_i)


    if(gamma_max > 1){
      # if gamma > select only steepest neighbor - can have ties
      N <- cells[gamma_max == gamma_i]

      if(length(N) > 1){
        nxt_cell <- sample(N, size = 1)
        prv_pos <- which(nxt_cell == ngh_cells)
      } else {
        nxt_cell <- N
        prv_pos <- which(nxt_cell == ngh_cells)
      }

    } else {
      # otherwise mfdf criterion
      N <- cells[gamma_i >= gamma_max^alpha]
      trans_prob <- prob[gamma_i >= gamma_max^alpha]

      if(length(N) > 1){
        nxt_cell <- sample(N, size = 1, prob = trans_prob)
        prv_pos <- which(nxt_cell == ngh_cells)
      } else {
        nxt_cell <- N
        prv_pos <- which(nxt_cell == ngh_cells)
      }

    }

    path_cells[[i]] <- nxt_cell
    cntr_cell <- nxt_cell


  }

  sim_paths[[k]] <- unlist(path_cells)

}

# Merge paths
r_blank <- setValues(dem, value = 0)
r_sims <- r_blank

for(k in 1:reps){
  r_path <- r_blank
  r_path[sim_paths[[k]]] <- 1

    r_sims <- r_path + r_sims
}

print(Sys.time() - start_time )

r_sims[r_sims == 0] <- NA
plot(r_sims)
plot(source_point, add = TRUE)
plot(runout_polygon, add = TRUE)

plot(rasterCdf(r_sims))
plot(source_point, add = TRUE)
plot(runout_polygon, add = TRUE)


setwd("/home/jason/Scratch/Figures")
#save(sim_paths, file = "rw_sim_paths_it1000.Rd")
(load("rw_sim_paths_it1000.Rd"))


# Create images #####################

r_blank <- setValues(dem, value = 0)
r_prev <- r_blank
r_first <- r_blank
r_first[r_first == 0] <- NA

plot(r_first, add = FALSE,
     legend = FALSE,
     axes = FALSE, box = FALSE, alpha = 0.7)


setwd("/home/jason/Scratch/Figures/rw_animate")


numFun <- function(x){

  # for numbers < 10000
  num_digi <- nchar(trunc(abs(x)))

  if( num_digi == 1){
    num <- paste0("000", x)
  } else if(num_digi == 2){
    num <- paste0("00", x)
  } else if(num_digi == 3) {
    num <- paste0("0", x)
  } else {
    num <- paste0(x)
  }

  return(num)

}


for(i in 1:100){


  img_name <- paste0("rwsim_it_", numFun(i), ".png")

  r_sim <- r_blank
  r_sim[sim_paths[[i]]] <- 1
  r_sim <- r_sim + r_prev
  r_sim[r_sim == 0] <- NA

  png(filename=img_name, res = 300, width = 7.5, height = 4,
      units = "in", pointsize = 11)

  par(family = "Arial", mfrow = c(1,1), mar = c(0.5, 1, 0.5, 0.5))
      #mgp = c(1, 0.75, 0))

  plot(r_sim, add = FALSE, col = viridis::viridis(10, direction = -1),
       legend = TRUE,
       axes = FALSE, box = FALSE, alpha = 0.7)

  plot(run_ply, add = TRUE, lty = 3)
  plot(src_pnt, add = TRUE)

  dev.off()

  r_prev <- r_sim
  r_prev[is.na(r_prev)] <- 0

}


# https://gifmaker.me/


# N labelled ##############
setwd("/home/jason/Scratch/Figures/rw_animate_label")



r_blank <- setValues(dem, value = 0)
r_prev <- r_blank
r_first <- r_blank
r_first[r_first == 0] <- NA

plot(r_first, add = FALSE,
     legend = FALSE,
     axes = FALSE, box = FALSE, alpha = 0.7)



numFun <- function(x){

  # for numbers < 10000
  num_digi <- nchar(trunc(abs(x)))

  if( num_digi == 1){
    num <- paste0("000", x)
  } else if(num_digi == 2){
    num <- paste0("00", x)
  } else if(num_digi == 3) {
    num <- paste0("0", x)
  } else {
    num <- paste0(x)
  }

  return(num)

}

for(i in 1:100){


  img_name <- paste0("rwsim_lbl_it_", numFun(i), ".png")

  r_sim <- r_blank
  r_sim[sim_paths[[i]]] <- 1
  r_sim <- r_sim + r_prev
  r_sim[r_sim == 0] <- NA

  png(filename=img_name, res = 300, width = 7.5, height = 4,
      units = "in", pointsize = 11)

  par(family = "Arial", mfrow = c(1,1), mar = c(0.5, 1, 2, 0.5))
  #mgp = c(1, 0.75, 0))

  plot(r_sim, add = FALSE, col = viridis::viridis(10, direction = -1),
       legend = TRUE,
       axes = FALSE, box = FALSE, alpha = 0.7)

  mtext(paste("n =", i), side=3, adj=0, line=-3, cex=1.7)

  plot(run_ply, add = TRUE, lty = 3)
  plot(src_pnt, add = TRUE)

  dev.off()

  r_prev <- r_sim
  r_prev[is.na(r_prev)] <- 0

}

# Animate Particle #####################


setwd("/home/jason/Scratch/Figures")
#save(sim_paths, file = "rw_sim_paths_it1000.Rd")
(load("rw_sim_paths_it1000.Rd"))

setwd("rw_particle_animate")
#pick one instance
path <- sim_paths[[1]]
r_blank <- setValues(crop_dem, NA)
r_path <- r_blank

for(i in 1:length(path)){

  img_name <- paste0("rw_particle_", numFun(i), ".png")


  png(filename=img_name, res = 300, width = 7.5, height = 4,
      units = "in", pointsize = 11)

  par(family = "Arial", mfrow = c(1,1), mar = c(0.5, 1, 2, 0.5))

  plot(r_blank, col=grey.colors(100, start=1, end=0), box = FALSE, axes = FALSE,
       legend = FALSE)
  r_path[path[i]] <- 1

  plot(r_path, add = TRUE, col = viridis::viridis(10, direction = -1), legend = FALSE)
  plot(run_ply, add = TRUE, lty = 3)
  plot(src_pnt, add = TRUE)

  dev.off()
}


# Animate optimization ######################


rwVisOpt <- function(dem, slide_plys = runout_polygon, slide_src = source_point,
                     src_thresh=0.5, SLP = 30, PER = 2, EXP = 3, hillshade,
                     map_ext = 3000, buffer_ext = 1500, buffer_source = 50){

  auroc <- rwPerformance(dem, slide_plys, slide_src,
                slp = SLP, ex = EXP, per = PER,
                gpp_iter = 1000, buffer_ext = map_ext, buffer_source = NULL,
                plot_eval = FALSE, saga_lib = saga)


  # If single runout polygon as input, assign slide_id value of 1
  if(length(slide_plys) == 1){
    slide_id <- 1
  }

  slide_plys$objectid <- 1:length(slide_plys)
  # Subset a single slide polygon
  slide_poly_single <- slide_plys[slide_id,]

  # Crop dem to slide polygon
  dem_grid <- raster::crop(dem, raster::extent(slide_poly_single) + buffer_ext)


  if(class(slide_src) == "SpatialPointsDataFrame"){
    # Subset corresponding source/start point of runout
    sel_over_start_point  <- sp::over(slide_src, slide_poly_single)
    sel_start_point <- slide_src[!is.na(sel_over_start_point$objectid),]

    if(!is.null(buffer_source)){
      # Create a buffer around source point to create a source/release area
      source_buffer <- rgeos::gBuffer(sel_start_point, width = buffer_source)
      source_grid <- raster::rasterize(source_buffer, dem_grid, field=1 )
      source_grid <- raster::mask(source_grid, slide_poly_single )
      source_plot <- raster::rasterToPolygons(source_grid)
    } else {
      # Just use source point
      source_plot <- sel_start_point
      source_grid <- raster::rasterize(matrix(sp::coordinates(sel_start_point)[1:2], ncol = 2), dem_grid, field = 1)
    }
  }

  if(class(slide_src) == "SpatialPolygonsDataFrame" ){
    sel_over_start_poly <- sp::over(slide_src, slide_poly_single)
    sel_start_poly <- slide_src[!is.na(sel_over_start_poly$objectid),]
    source_plot <- sel_start_poly
    source_grid <- raster::rasterize(sel_start_poly, dem_grid, field=1 )

  }


  # Run runout model using Rsagacmd (faster than RSAGA)
  gpp <- saga$sim_geomorphology$gravitational_process_path_model(dem = dem_grid, release_areas = source_grid,
                                                                     process_path_model = 1,
                                                                     rw_slope_thres = SLP,
                                                                     rw_exponent = EXP,
                                                                     rw_persistence = PER,
                                                                     gpp_iterations = 1000)

  r_cdf <- rasterCdf(gpp$process_area)

  hs <- crop(hillshade, dem_grid)

  plot(hs,
       col=grey.colors(100, start=1, end=0),
       legend=F,
       #cex.axis = 0.7, cex.main = 0.8,
       axes = FALSE, box = FALSE)

  mtext(paste("Slp. thresh:", SLP, "\nAUROC:", round(auroc, digits = 3)),
        side=3, adj=0, line=0, cex=1)


  plot( r_cdf, col = viridis::viridis(10, direction = -1), add = T, alpha = 0.7)
  plot(slide_plys, add = TRUE, lty = 3)


}



steps <- 11
rwexp_vec <- seq(1.3, 3, len=steps)
rwper_vec <- seq(1.5, 2, len=steps)
rwslp_vec <- seq(20, 40, len=steps)

for(i in 1:length(rwslp_vec)){

  par(family = "Arial", mfrow = c(1,1), mar = c(0.5, 1, 4, 0.5))

  rwVisOpt(dem, slide_plys = runout_polygon, slide_src = source_point,
           src_thresh=0.5, SLP = rwslp_vec[i], PER = 2, EXP = 2, hillshade,
           map_ext = 0, buffer_ext = 1500, buffer_source = 30)

}

