rasterThreshold <- function(x, prob_range = c(0.9, Inf)){
  if(prob_range[2] != Inf){

    m <- c(-Inf, prob_range[1], NA, prob_range, 1, prob_range[2], Inf, NA)
    rclmat <- matrix(m, ncol=3, byrow=TRUE)

  } else {

    m <- c(-Inf, prob_range[1], NA, prob_range, 1)
    rclmat <- matrix(m, ncol=3, byrow=TRUE)

  }

  rc <- reclassify(x, rclmat)
  rc

}
