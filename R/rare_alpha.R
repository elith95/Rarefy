
rare_alpha <- function(comm, dist_xy = NULL, method = c("HCDT", "hill"), q = 0, random = 99, spatial = FALSE) {
  
  method <- method[1]
  if(!method %in% c("HCDT", "hill")) stop("Unavailable method")
  
  if (!inherits(comm, "matrix") && !inherits(comm, "data.frame")) stop("Non convenient comm")
  if (any(comm<0)) stop("Negative value in comm")
  if(suppressWarnings(any(rowSums(comm)) == 0)) stop("Empty row")
  if(suppressWarnings(any(colSums(comm)) == 0)) {
    v <- apply(comm, 2, function(col) any(col != 0))
    comm <- comm[,v]
  }
  
  if(!spatial && !is.null(dist_xy)) warning("Spatial is FALSE, so classic rarefaction will be calculated")
  if(spatial && is.null(dist_xy)) stop("Object of class 'dist' expected for dist_xy")
  if(!is.null(dist_xy)){
    if (!inherits(dist_xy, "dist")) stop("Object of class 'dist' expected for dist_xy") 
    dist_xy<-as.matrix(dist_xy)
    if (any(dist_xy < 0)) stop("Negative value in dist_xy")
    if (nrow(comm) != nrow(dist_xy)) stop("comm and dist_xy don't have the same number of plots")
    if(!is.null(rownames(dist_xy)) && !is.null(rownames(comm)) ) {
      if(any(!rownames(comm) %in% rownames(dist_xy))) stop("comm and dist_xy must have the same names for the plots")
    } else if(!is.null(rownames(dist_xy)) && is.null(rownames(comm))) {
      rownames(comm) <- rownames(dist_xy)
      warning("comm has no rownames")
      warning("rownames of dist_xy set as rownames of comm")
    } 
  }
  
  if(method == 'HCDT') {
    r_fin <- array(dim = c(ifelse(spatial, nrow(comm), random), nrow(comm)))
    nami <- rownames(comm)
    
    for(i in 1:ifelse(spatial, nrow(comm), random)) {
      if(spatial) v <- nami[order(dist_xy[,i])]
      else v <- sample(1:nrow(comm), nrow(comm))
      x <- comm[v,]
      x <- apply(x,2,cumsum)
      x <- sweep(x,1,rowSums(x),"/")
      if(q == 1) r_fin[i,] <- apply(x,1, function(x) -sum(x[x>0]*log(x[x>0])))
      else r_fin[i,] <- apply(x, 1, function(x) (1-sum(x[x>0]**q))/(q-1))
    }
  }
  
  else if(method == 'hill') {
    r_fin <- array(dim = c(ifelse(spatial, nrow(comm), random), nrow(comm)))
    nami <- rownames(comm)
    for(i in 1:ifelse(spatial, nrow(comm), random)) {
      if(spatial) v <- nami[order(dist_xy[,i])]
      else v <- sample(1:nrow(comm), nrow(comm))
      x <- comm[v,]
      x <- sweep(x, 1, rowSums(x),"/")
      x <- apply(x, 2, cummean) 
      if(q == 1) r_fin[i,] <- apply(x,1, function(x) exp(-sum(x[x>0]*log(x[x>0]))))
      else r_fin[i,] <- apply(x,1, function(x) sum((x[x>0]**q))**(1/(1-q)))
    }
  }
  
  rare <- colMeans(r_fin,na.rm = TRUE)
  IC_up <- rare + (1.96*(sd(r_fin,na.rm = TRUE)/sqrt(ifelse(spatial,nrow(comm),random))))
  IC_low <- rare - (1.96*(sd(r_fin,na.rm=TRUE)/sqrt(ifelse(spatial,nrow(comm),random))))
  df <- data.frame(rare, IC_up, IC_low)
  colnames(df) <- c('Rarefaction', 'IC_up', 'IC_low')
  
  return(df)
}


