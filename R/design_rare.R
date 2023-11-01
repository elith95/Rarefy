
design_rare <- function(comm, dist_xy = NULL, method = c("fun_div", "equation"), args, fun = NULL, formula = NULL, resampling = 99, spatial = FALSE, mean = F, beta = F) {
  
  method <- method[1]
  if(!method %in% c("fun_div", "equation")) stop("Unavailable method")
  
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
  
  if(!is.list(args)) stop("args must be a list")
  if(is.null(names(args))) stop("args must have names")
  
  if(method == "fun_div") {
    if(!inherits(fun, 'character')) stop("fun_div must be a character") 
    if(!exists(fun)) stop("the function doesn't exist") 
    
    f <- match.fun(fun)
    arg <- as.list(args(f)) 
    v <- names(arg)
    
    if(!all(names(args) %in% v)) stop("The arguments must be the ones specified by the function")
    
    ind <- match(NA, unlist(args))
    ind <- match(names(unlist(args)[ind]), names(args))
    
    nami <- rownames(comm)
    r_fin <- array(dim = c(ifelse(spatial, nrow(comm), resampling), nrow(comm)))
    
    if(beta) {
      
      for(i in 1:ifelse(spatial, nrow(comm), resampling)) {
        if(spatial) v <- nami[order(dist_xy[, i])]
        else v <- sample(1:nrow(comm), nrow(comm))
        x <- comm[v,]
        args[[ind]] <- x
        beta <- as.matrix(do.call(f, args))
        r_fin[i,] <- unlist(lapply(2:nrow(comm), function(j) {
          m <- beta[1:j,1:j]
          return(mean(m))
        }))
      }
    }
    
    else {
      
      for(i in 1:ifelse(spatial, nrow(comm), resampling)) {
        if(spatial) v <- nami[order(dist_xy[,i])]
        else v <- sample(1:nrow(comm), nrow(comm))
        x <- comm[v,]
        x <- apply(x, 2, ifelse(mean, cummean, cumsum))
        args[[ind]] <- x
        r_fin[i,] <- as.matrix(do.call(f,args))
      }
      
    }
  }
  
  else if(method == "equation") {
    
    if(!is.character(formula)) stop("formula must be a string")
    str <- strsplit(gsub("[^a-zA-Z]", " ", formula), " ")[[1]]
    if(!names(args) %in% str[nzchar(str)]) stop("The names in args must correspond to the formula terms")
    if(beta) stop("With the method equation in this package version you can't produce a distance matrix as output of the formula")
    
    nami <- rownames(comm)
    r_fin <- array(dim = c(ifelse(spatial, nrow(comm), resampling), nrow(comm)))
    
    for(i in 1:ifelse(spatial, nrow(comm), resampling)) {
      if(spatial) v <- nami[order(dist_xy[,i])]
      else v <- sample(1:nrow(comm), nrow(comm))
      x <- comm[v,]
      x <- apply(x, 2, ifelse(mean, cummean, cumsum))
      args[[ind]] <- x
      r_fin[i,] <- unlist(lapply(2:nrow(comm), function(j) {
        res <- eval(parse(text = formula), args)
        return(res)
      }))
    }
  }
  
  fin <- colMeans(r_fin, na.rm=TRUE)
  IC_up <- fin + (1.96*(sd(r_fin,na.rm=TRUE)/sqrt(ifelse(spatial,nrow(comm),resampling))))
  IC_low <- fin - (1.96*(sd(r_fin,na.rm=TRUE)/sqrt(ifelse(spatial,nrow(comm),resampling))))
  df <- data.frame(fin, IC_up, IC_low)
  colnames(df) <- c('Rarefaction', 'IC_up', 'IC_low')
  
  return(df)
  
}



