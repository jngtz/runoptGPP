# Notes
# - using custom adjcent cells function makes it much faster
#   but cellFromXY now slowing things down...
# - this code may already be 'fast' enough to check for cell connectivity



# Load Packages and Data ####################################################################
library(raster)
library(rgdal)
library(profvis)

# Set workspace
setwd("/home/jason/Data/Chile/")
setwd("C:\\Projects\\cetaqua\\data")
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

plot(crop_dem)
plot(runout_polygon, add = TRUE)
plot(source_point, add = TRUE)

# random walk implementation in R (it's slow...) ########################################

euclideanDistance <- function(p1, p2){
  sqrt( (p1[1] - p2[1])^2 + (p1[2]- p2[2])^2 )
}


adjCells <- function(r, xy){
  #r: resolution of raster
  #xy: xy location of center cell

  d <- c(rep(xy[,1]-r[1], 3), rep(xy[,1]+r[1],3), xy[,1], xy[,1],
         rep(c(xy[,2]+r[2], xy[,2], xy[,2]-r[2]), 2), xy[,2]+r[2], xy[,2]-r[2])

  d <- matrix(d, ncol=2)

  #ngh_cells <- cellFromXY(dem, d)

}

randomWalk <- function(dem, source_point, slp_thresh = 30, exp_div = 3, per_fct = 2, reps = 100){

  # get initial cell center
  xy = coordinates(source_point)
  cntr_cell <- cellFromXY(dem, xy = xy)
  i = 0
  n = 0
  sim_paths <- vector(mode = "list", length = reps)
  cell_pos <- 9999
  r = res(dem)
  v_dem <- getValues(dem)
  v_blank <- rep(0, ncol(dem)*nrow(dem))


  # Get distance to each cell
  #ngh_cells <- adjacent(dem, cntr_cell, directions = 8, pairs = FALSE, id = TRUE)
  d <- adjCells(r, xy)
  ngh_cells <- cellFromXY(dem, d)

  ngh_points <- xyFromCell(dem, ngh_cells)
  cntr_point <- xyFromCell(dem, cntr_cell)

  cell_dist <- apply(ngh_points, 1, euclideanDistance, p2 = cntr_point)

  # Start repeated simulations

  for(k in 1:reps){

    #print(k)
    path_cells <- list()
    xy_cell = coordinates(source_point)
    cntr_cell <- cellFromXY(dem, xy = xy_cell)

    i = 0
    n = 0
    prv_pos <- 9999


    while(n == 0){

      i = i + 1

      #ngh_cells <- adjacent(dem, cntr_cell, directions = 8, pairs = FALSE, id = TRUE)
      d <- adjCells(r, xy_cell)
      ngh_cells <- cellFromXY(dem, d)


      if(length(ngh_cells) < 8){
        break
      }

      elv_values <- v_dem[c(ngh_cells, cntr_cell)]
      elv_ngh <- elv_values[1:8]
      elv_cntr <- elv_values[9]

      lower_elv <- elv_ngh <= elv_cntr

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
      xy_cell <- xyFromCell(dem, cntr_cell)


    }

    sim_paths[[k]] <- unlist(path_cells)

  }

  # Merge paths

  v_sims <- v_blank

  for(k in 1:reps){
    v_path <- v_blank
    v_path[sim_paths[[k]]] <- 1

    v_sims <- v_path + v_sims
  }

  return(setValues(dem, v_sims))

}

start_time <- Sys.time()
r_sims <- randomWalk(crop_dem, source_point, reps = 1000)
print(Sys.time() - start_time )

library(microbenchmark)
microbenchmark(randomWalk(crop_dem, source_point, reps = 10))

profvis(randomWalk(crop_dem, source_point, reps = 100))

r_sims[r_sims == 0] <- NA
plot(r_sims)
plot(source_point, add = TRUE)
plot(runout_polygon, add = TRUE)

#plot(rasterCdf(r_sims))
#plot(source_point, add = TRUE)
#plot(runout_polygon, add = TRUE)


# Testing code speed ##################
library(microbenchmark)


dem = crop_dem
slp_thresh = 30
exp_div = 3
per_fct = 2
reps = 1


# get initial cell center
cntr_cell <- cellFromXY(dem, xy = coordinates(source_point))
i = 0
n = 0
sim_paths <- vector(mode = "list", length = reps)
cell_pos <- 9999
r = res(dem)
xy = coordinates(source_point)

# Get distance to each cell
ngh_cells <- adjacent(dem, cntr_cell, directions = 8, pairs = FALSE, id = TRUE)
ngh_points <- xyFromCell(dem, ngh_cells)
cntr_point <- xyFromCell(dem, cntr_cell)


# test strt
adj_cells <- function(r, xy){


  d <- c(rep(xy[,1]-r[1], 3), rep(xy[,1]+r[1],3), xy[,1], xy[,1],
         rep(c(xy[,2]+r[2], xy[,2], xy[,2]-r[2]), 2), xy[,2]+r[2], xy[,2]-r[2])


  d <- matrix(d, ncol=2)

  cellFromXY(dem, d)

}

adj_cells(r, xy)

microbenchmark(adjacent(dem, cntr_cell, directions = 8, pairs = FALSE, id = TRUE), adj_cells(r, xy))


#test end


cell_dist <- apply(ngh_points, 1, euclideanDistance, p2 = cntr_point)



# Start repeated simulations

for(k in 1:reps){

  #print(k)

  path_cells <- list()
  cntr_cell <- cellFromXY(dem, xy = coordinates(source_point))

  row_col <- rowColFromCell(dem, cntr_cell) # testing

  i = 0
  n = 0
  prv_pos <- 9999


  while(n == 0){

    i = i + 1


    ngh_cells <- adjacent(dem, cntr_cell, directions = 8, pairs = FALSE, id = TRUE)
    ngh_cells <- adj_cells(r, xy)

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



