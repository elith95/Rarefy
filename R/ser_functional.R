

ser_functional <- function(comm, dist_f = NULL, dist_xy, method=c('rao','chao'), tau = NA, q = 0, comparison = FALSE, resampling = 99) {
  
  #Selection of the method
  method <- method[1]
  if (!method %in% c("rao", "chao")) stop("Unavailable method")
  
  #Community matrix control
  if (!inherits(comm, "matrix") && !inherits(comm, "data.frame")) stop("Non convenient comm")
  if (any(comm < 0)) stop("Negative value in comm")
  if(suppressWarnings(any(rowSums(comm)) == 0)) stop("Empty row")
  if(suppressWarnings(any(colSums(comm)) == 0)) {
    v<-apply(comm, 2, function(col) any(col !=0 ))
    comm<-comm[,v]
  }
    
  #Distance matrix control
  if (!inherits(dist_xy, "dist")){
    if(!is.vector(dist_xy)){
      if(!is.matrix(dist_xy)){
        if(!is.data.frame(dist_xy))
          stop("Object dist_xy must be of class vector, matrix, data.frame or dist") }}}
  
  if(!inherits(dist_xy, "dist")) {
    if(is.vector(dist_xy)){
      if(length(dist_xy) != nrow(comm)) stop("Incorrect definition of object gradient")
      if(!is.null(names(dist_xy))){
        if(any(!rownames(comm) %in% names(dist_xy))) stop("Names in gradient must be the same as row names in community")
        dist_xy <- dist_xy[rownames(comm)] }}
    else { if (nrow(comm) != nrow(dist_xy)) stop("comm and dist_xy don't have the same number of plots") }
    dist_xy <- dist(scale(dist_xy)) 
    warning("dist_xy is used to calculate a distance matrix between sampling units using the method euclidean of the function dist()")
  }
  
  dist_xy <- as.matrix(dist_xy)
  if (nrow(comm) != nrow(dist_xy)) stop("comm and dist_xy don't have the same number of plots")
  if (any(dist_xy<0)) stop("Negative value in dist_xy")
  if(!is.null(rownames(dist_xy)) && !is.null(rownames(comm)) ) {
    if(any(!rownames(comm) %in% rownames(dist_xy))) stop("comm and dist_xy must have the same names for the plots")
  } else if(!is.null(rownames(dist_xy)) && is.null(rownames(comm))) {
    rownames(comm)<-rownames(dist_xy)
    warning("comm has no row names")
    warning("row names of dist_xy set as row names of comm")
  } 
  
  #Functional traits distance matrix control
  if(is.null(dist_f)) stop("dist_f must have a value")
  if(!is.null(dist_f)) {
    if (!inherits(dist_f, "dist")) stop("Object of class 'dist' expected for dist_f")
    if (!is.euclid(sqrt(dist_f))) warning("Squared Euclidean or Euclidean property expected for dist_f")
    dist_f <- as.matrix(dist_f)
    if (ncol(comm) != nrow(dist_f)) stop("comm and dist_f don't have the same number of species")
    if(!is.null(rownames(dist_f)) && !is.null(colnames(comm)) && any(!colnames(comm) %in% rownames(dist_f))) stop("Names in dist_f must be the same as col names in comm") 
    if (any(dist_f < 0)) stop("Negative value in dist_f") 
  }
  
  #Rao index
  if(method == 'rao') {
    r_fin <- array(dim = c(nrow(comm), nrow(comm)))
    com <- sweep(comm, 1, rowSums(comm), "/")
    nami <- rownames(comm)
    
    for(i in 1:nrow(comm)) {
      v <- nami[order(dist_xy[,i])]
      x <- com[v,]
      x <- apply(x,2,cummean)
      r_fin[i,] <- apply(x,1,function(newp) t(as.vector(newp)) %*% dist_f %*% as.vector(newp))
    }
    
    fin <- colMeans(r_fin, na.rm=TRUE)
    IC_up <- fin + (1.96*(sd(r_fin,na.rm=TRUE)/sqrt(ncol(comm))))
    IC_low <- fin - (1.96*(sd(r_fin,na.rm=TRUE)/sqrt(ncol(comm))))
    df <- data.frame(fin, IC_up, IC_low)
    colnames(df) <- c('Rarefaction','IC_up','IC_low')
    
    if(comparison) {
      nsr <- rare_Rao(comm, as.dist(dist_f), resampling=resampling)
      df <- data.frame(df, nsr)
    }
  }
  
  #Chao index
  if(method == 'chao') {
    dist_f[which(dist_f > tau,arr.ind = T)] <- tau
    com <- sweep(comm, 1, rowSums(comm), "/")
    nami <- rownames(comm)
    r_fin <- array(dim = c(nrow(comm), nrow(comm)))
    
    for(i in 1:nrow(comm)) {
      v <- nami[order(dist_xy[, i])]
      x <- com[v,]
      x <- apply(x,2,cummean)
      for(j in 1:nrow(comm)) {
        v1 <- as.matrix(x[j,])
        a <- as.vector((1 - dist_f/tau) %*% as.vector(v1))
        d <- x[j,]
        d <- d[a!=0]
        nplus <- sum(d)
        a <- a[a!=0]
        d <- d/a
        if(q==1){
          r_fin[i,j] <- exp(sum(-d*a/nplus*log(a/nplus)))
        }else{
          r_fin[i,j] <- (sum(d*(a/nplus)^q))^(1/(1-q))
        }
      } 
    }
    
    fin <- colMeans(r_fin, na.rm=TRUE)
    IC_up <- fin + (1.96*(sd(r_fin, na.rm = TRUE)/sqrt(nrow(comm))))
    IC_low <- fin - (1.96*(sd(r_fin, na.rm = TRUE)/sqrt(nrow(comm))))
    df <- data.frame(fin, IC_up, IC_low)
    colnames(df) <- c('Rarefaction', 'IC_up', 'IC_low')
    
    if(comparison) {
      r_fin <- do.call(rbind, lapply(1:resampling, function(i) {
        v <- sample(1:nrow(com), nrow(com))
        ch <- NA
        x <- com[v,]
        x <- apply(x, 2, cummean)
        
        for(j in 1:nrow(com)) {
          v1 <- as.matrix(x[j,])
          a <- as.vector((1 - dist_f/tau) %*% as.vector(v1))
          d <- x[j,]
          d <- d[a!=0]
          nplus <- sum(d)
          a <- a[a!=0]
          d <- d/a
          if(q == 1){
            ch[j] <- exp(sum(-d*a/nplus*log(a/nplus)))
          }else{
            ch[j] <- (sum(d*(a/nplus)^q))^(1 / (1-q))
          }
        } 
        
        return(ch)
      }))
      
      fin <- colMeans(r_fin, na.rm = TRUE)
      IC_up <- fin + (1.96*(sd(r_fin, na.rm = TRUE)/sqrt(nrow(comm))))
      IC_low <- fin - (1.96*(sd(r_fin, na.rm = TRUE)/sqrt(nrow(comm))))
      df1 <- data.frame(fin, IC_up, IC_low)
      colnames(df1) <- c('Classic Rarefaction', 'IC_up', 'IC_low')
      df <- data.frame(df, df1)
    }
  }
  
  return(df)
  
}

