
prep_dist <- function(xy, method = c("euclid", "k-NN"), k = 1) {
  
  method <- method[1]
  if(!method %in% c("euclid", "k-NN")) stop("Unavailable method")
  
  if(is.null(rownames(xy))) stop("xy must have names for columns") 
  if (!inherits(xy, "matrix") || !inherits(xy, "data.frame")) stop("Non convenient xy")
  if(ncol(xy) < 2) stop("xy must contain at list two values for plots coordinates")
  if(nrow(xy) < 2) stop("xy must contain at list two plots")
  
  if(method == "euclid") {
    
    dist_xy <- dist(xy)
    
  }
  
  else if(method == "k-NN") {
    
    knn <- kNNdist(dist(xy), k, all = T)
    dist_xy <- as.dist(knn)
    
  }
  
  return(dist_xy)
  
}