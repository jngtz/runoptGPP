

# Load Packages and Data ####################################################################
library(raster)
library(rgdal)

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

crop_dem <- crop(dem, extent(runout_polygon)+500)

# random walk implementation in R (it's slow...) ########################################

euclideanDistance <- function(p1, p2){
  sqrt( (p1[1] - p2[1])^2 + (p1[2]- p2[2])^2 )
}

randomWalk <- function(dem, source_point, slp_thresh = 30, exp_div = 3, per_fct = 2, reps = 100){

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

  cell_dist <- apply(ngh_points, 1, euclideanDistance, p2 = cntr_point)

  # Start repeated simulations

  for(k in 1:reps){

    #print(k)
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

      if(anyNA(beta_ngh)){
        break
      }

      f <- rep(1, 8)

      if(prv_pos < 9){
        f[prv_pos] <- per_fct
      }

      f <- f[lower_elv]
      cells <- ngh_cells[lower_elv]
      beta_ngh <- beta_ngh[lower_elv]

      gamma_i <- tan(beta_ngh*pi/180) / tan(slp_thresh*pi/180)

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
        N <- cells[gamma_i >= gamma_max^exp_div]
        trans_prob <- prob[gamma_i >= gamma_max^exp_div]

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

  return(r_sims)

}

start_time <- Sys.time()

r_sim <- randomWalk(crop_dem, source_point, reps = 100)

print(Sys.time() - start_time )

r_sims[r_sims == 0] <- NA
plot(r_sims)
plot(source_point, add = TRUE)
plot(runout_polygon, add = TRUE)

plot(rasterCdf(r_sims))
plot(source_point, add = TRUE)
plot(runout_polygon, add = TRUE)


