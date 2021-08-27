# Notes
# - Much faster using matrix indicies instead of raster... work with matrix then convert back to raster...
# - Need to make pcmRun faster... start w supresswarnings call, find alnterantive to deal with NAN
# PCM can be made faster by making a quicker NaN handling than supress...


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


dem = crop_dem
slp_thresh = 30
exp_div = 3
per_fct = 2
reps = 1

# random walk implementation in R (it's slow...) ########################################

library(compiler)

pcmRun <- function(mu = 0.1, md = 40, v_p = 1, theta_p = 30, theta_i= 20, l = 12.5){
  #vp = v i-1
  #l depend on direction of movement from previous cell

  g = 9.80665

  alpha <-  g*(sin(theta_i*pi/180) - mu*cos(theta_i*pi/180))
  beta <- -2*l / (md)

  if(theta_p > theta_i){
    delta_theta = theta_p - theta_i
  } else {
    delta_theta = 0
  }

  if(alpha*(md)*(1-exp(beta))+(v_p)^2*exp(beta)*cos(delta_theta*pi/180) < 0){
    v_i = NaN
  } else {

    v_i <- sqrt(alpha*(md)*(1-exp(beta))+(v_p)^2*exp(beta)*cos(delta_theta*pi/180))

  }

  return(v_i)
}

#pcmRun(mu = 0.1, md = 40, v_p = 1, theta_p = 30, theta_i = 20, l = 12.5)

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

adjRowCol <- function(rowcol){

  d <- c(rowcol[1] - 1, rowcol[2] -1,
         rowcol[1], rowcol[2] -1,
         rowcol[1] + 1, rowcol[2] -1,

         rowcol[1] - 1, rowcol[2] + 1,
         rowcol[1], rowcol[2] + 1,
         rowcol[1] + 1, rowcol[2]  + 1,

         rowcol[1] - 1, rowcol[2],
         rowcol[1] + 1, rowcol[2])

  d <- matrix(d, ncol = 2, byrow = TRUE)

  return(d)
}



pcmRunComp <- cmpfun(pcmRun)
adjRowColComp <- cmpfun(adjRowCol)


pcmRW <- function(dem, source_point, mu = 0.1, md = 40, slp_thresh = 30, exp_div = 3, per_fct = 2, reps = 100){

  # get initial cell center
  xy = coordinates(source_point)
  cntr_cell <- cellFromXY(dem, xy = xy)
  i = 0
  n = 0
  sim_paths <- vector(mode = "list", length = reps)
  sim_vels <- vector(mode = "list", length = reps)
  cell_pos <- 9999
  r = res(dem)
  v_dem <- getValues(dem)
  v_blank <- rep(0, ncol(dem)*nrow(dem))

  rw <- rowFromY(dem, xy[2])
  cl <- colFromX(dem, xy[1])
  rowcol <- c(rw, cl)

  # Get distance to each cell
  #ngh_cells <- adjacent(dem, cntr_cell, directions = 8, pairs = FALSE, id = TRUE)
  d <- adjCells(r, xy)
  ngh_cells <- cellFromXY(dem, d)

  ngh_points <- xyFromCell(dem, ngh_cells)
  cntr_point <- xyFromCell(dem, cntr_cell)

  cell_dist <- apply(ngh_points, 1, euclideanDistance, p2 = cntr_point)

  m_dem <- as.matrix(dem)
  m_blank <- as.matrix(setValues(dem, 0))

  # Start repeated simulations

  for(k in 1:reps){

    #print(k)
    path_cells <- list()
    vel_cells <- list()
    cntr_cell <- rowcol

    i = 0
    n = 0
    prv_pos <- 9999
    v_p = 1
    theta_p = 1


    while(n == 0){

      i = i + 1

      #ngh_cells <- adjacent(dem, cntr_cell, directions = 8, pairs = FALSE, id = TRUE)
      #d <- adjCells(r, xy_cell)
      #ngh_cells <- cellFromXY(dem, d)

      d <- adjRowColComp(cntr_cell)


      if(length(d) < 16){
        break
      }

      #elv_values <- v_dem[c(ngh_cells, cntr_cell)]

      elv_values <- m_dem[d]
      elv_ngh <- elv_values[1:8]
      cell_id <- 1:8
      elv_cntr <- m_dem[cntr_cell[1], cntr_cell[2]]

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

      #cells <- ngh_cells[lower_elv]
      cells <- cell_id[lower_elv]

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
         # prv_pos <- which(nxt_cell == cell_id)
        } else {
          nxt_cell <- N
         # prv_pos <- which(nxt_cell == cell_id)
        }

      } else {
        # otherwise mfdf criterion
        N <- cells[gamma_i >= gamma_max^exp_div]
        trans_prob <- prob[gamma_i >= gamma_max^exp_div]

        if(length(N) > 1){
          nxt_cell <- sample(N, size = 1, prob = trans_prob)
         # prv_pos <- which(nxt_cell == cell_id)
        } else {
          nxt_cell <- N
         # prv_pos <- which(nxt_cell == cell_id)
        }

      }

      prv_pos <- nxt_cell
      # apply PCM model to determine velocity


      v_i <- pcmRunComp(mu = mu, md = md, v_p = v_p, theta_p = theta_p, theta_i = beta_ngh[nxt_cell == cells], l = cell_dist[nxt_cell])

      if(is.nan(v_i)){
        break
      }

      v_p <- v_i
      theta_p = beta_ngh[nxt_cell == cells]

      path_cells[[i]] <- d[nxt_cell,]
      vel_cells[[i]] <- v_i
      cntr_cell <- d[nxt_cell,]
    }

    sim_paths[[k]] <- matrix(unlist(path_cells), ncol = 2, byrow = TRUE)
    sim_vels[[k]] <- unlist(vel_cells)

  }

  # Merge paths

  m_sims <- m_blank
  m_svel <- m_blank

  for(k in 1:reps){
    m_path <- m_blank
    m_path[sim_paths[[k]]] <- 1
    m_sims <- m_path + m_sims

    m_vel <- m_blank
    m_vel[sim_paths[[k]]] <- sim_vels[[k]]
    m_svel <- m_vel + m_svel

  }

  m_svel <- m_svel / k

  r_sim <- raster::raster(m_sims, xmn=dem@extent@xmin,
                          xmx=dem@extent@xmax,
                          ymn=dem@extent@ymin,
                          ymx=dem@extent@ymax,
                          crs=crs(dem))

  r_svel <- raster::raster(m_svel, xmn=dem@extent@xmin,
                          xmx=dem@extent@xmax,
                          ymn=dem@extent@ymin,
                          ymx=dem@extent@ymax,
                          crs=crs(dem))

  return(list(parea = r_sim, avg_vel = r_svel))

}

start_time <- Sys.time()
r_sims <- pcmRW(crop_dem, source_point, mu = 0.1, md = 40, reps = 1000)
print(Sys.time() - start_time )

plot(r_sims$avg_vel)

r_sims$parea[r_sims$parea == 0] <- NA
plot(r_sims$parea)
plot(source_point, add = TRUE)
plot(runout_polygon, add = TRUE)

profvis(pcmRW(crop_dem, source_point, mu = 0.1, md = 40, reps = 1000))

r_sims[r_sims == 0] <- NA
plot(r_sims)
plot(source_point, add = TRUE)
plot(runout_polygon, add = TRUE)

#plot(rasterCdf(r_sims))
#plot(source_point, add = TRUE)
#plot(runout_polygon, add = TRUE)


