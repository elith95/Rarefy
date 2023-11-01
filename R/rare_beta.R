
rare_beta <- function(comm, dist_xy = NULL, method = c("whittaker", "jaccard", "bray", "cody"), random=99, spatial = FALSE) {
  
  method <- method[1]
  if(!method %in% c("whittaker", "jaccard", "bray", "cody")) stop("Unavailable method")
  if(method == 'cody' && !spatial) stop("cody's index can be used only with spatial = TRUE")
  
  if(is.null(colnames(comm))) stop("comm must have names for columns") 
  if (!inherits(comm, "matrix") && !inherits(comm, "data.frame")) stop("Non convenient comm")
  if (any(comm < 0)) stop("Negative value in comm")
  if(suppressWarnings(any(rowSums(comm)) == 0)) stop("Empty row")
  if(suppressWarnings(any(colSums(comm)) == 0)) {
    v <- apply(comm, 2, function(col) any(col != 0 ))
    comm <- comm[, v]
  }
  
  if(!spatial && !is.null(dist_xy)) warning("Spatial is FALSE, so classic rarefaction will be calculated")
  if(spatial && is.null(dist_xy)) stop("Object of class 'dist' expected for dist_xy")
  if(!is.null(dist_xy)) {
    if (!inherits(dist_xy, "dist")) stop("Object of class 'dist' expected for dist_xy") 
    dist_xy <- as.matrix(dist_xy)
    if (nrow(comm) != nrow(dist_xy)) stop("comm and dist_xy don't have the same number of plots")
    if(!is.null(rownames(dist_xy)) && !is.null(rownames(comm)) ) {
      if(any(!rownames(comm) %in% rownames(dist_xy))) stop("comm and dist_xy must have the same names for the plots")
    } else if(!is.null(rownames(dist_xy)) && is.null(rownames(comm))) {
      rownames(comm) <- rownames(dist_xy)
      warning("comm has no row names")
      warning("row names of dist_xy set as row names of comm")
    } 
    if (any(dist_xy < 0)) stop("Negative value in dist_xy") 
  }
  
  if(method == 'whittaker') {
    r_fin <- array(dim = c(ifelse(spatial, nrow(comm), random), nrow(comm)))
    nami <- rownames(comm)
    
    for(i in 1:ifelse(spatial, nrow(comm), random)) {
      if(spatial) v <- nami[order(dist_xy[, i])]
      else v <- sample(1:nrow(comm), nrow(comm))
      x <- comm[v,]
      sub <- specnumber(x)
      sub <- cummean(sub)
      x <- apply(x, 2, cumsum)
      g <- specnumber(x)
      r_fin[i,] <- (g - sub) / sub
    }
  }
  
  else if(method == 'jaccard' || method == 'bray') {
    r_fin <- array(dim = c(ifelse(spatial, nrow(comm), random), nrow(comm)-1))
    nami <- rownames(comm)
    
    for(i in 1:ifelse(spatial, nrow(comm), random)) {
      if(spatial) v <- nami[order(dist_xy[, i])]
      else v <- sample(1:nrow(comm), nrow(comm))
      x <- comm[v,]
      beta <- as.matrix(vegdist(x, method = method, binary = ifelse(method == 'jaccard', TRUE, FALSE)))
      r_fin[i,] <- unlist(lapply(2:nrow(comm), function(j) {
        m <- beta[1:j,1:j]
        return(mean(m))
      }))
    }
  }
  
  else if(method == 'cody') {
    r_fin <- array(dim = c(nrow(comm), nrow(comm)-1))
    nami <- rownames(comm)
    for(i in 1:nrow(comm)) {
      v <- nami[order(dist_xy[, i])]
      x <- comm[v,]
      p1 <- x[1,]
      x1 <- apply(x, 2, cumsum)
      g <- apply(x1[-1,], 1, function(x) length(x[x[p1 == 0] > 0]))
      l <- apply(x[-1,], 1, function(x) length(x[x == 0]))
      r_fin[i,] <- (g + l) / 2
    }
  }
  
  rare <- colMeans(r_fin)
  IC_up <- rare + (1.96 * (sd(r_fin) / sqrt(ifelse(spatial, nrow(comm), random))))
  IC_low <- rare - (1.96 * (sd(r_fin) / sqrt(ifelse(spatial, nrow(comm), random))))
  rare <- c(NA, rare)
  IC_up <- c(NA, IC_up)
  IC_low <- c(NA, IC_low)
  df <- data.frame(rare, IC_up, IC_low)
  colnames(df) <- c('Rarefaction', 'IC_up', 'IC_low')
  
  return(df)
  
}

