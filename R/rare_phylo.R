
rare_phylo <- function(comm, tree, method = c("faith", "barker", "Ia", "hill", "tsallis", "renyi"), exp = 0, resampling = 99) {
  
  method <- method[1]
  if(!method %in% c("faith", "barker", "Ia", "tsallis", "hill", "renyi")) stop("Unavailable method")
  
  if(is.null(colnames(comm))) stop("comm must have names for columns") 
  if (!inherits(comm, "matrix") && !inherits(comm, "data.frame")) stop("Non convenient comm")
  if (any(comm < 0)) stop("Negative value in comm")
  if(suppressWarnings(any(rowSums(comm)) == 0)) stop("Empty row")
  if(suppressWarnings(any(colSums(comm)) == 0)) {
    v <- apply(comm, 2, function(col) any(col != 0 ))
    comm <- comm[, v]
  }
  
  if(is.null(tree)) stop("tree must have a value")
  if(!is.null(tree)) {
    if(!inherits(tree, "phylo") && !inherits(tree, "phylo4")) stop("tree must be of class phylo or phylo4") 
    if(inherits(tree, "phylo4")) suppressWarnings(tree<-as(tree, "phylo"))
    if(any(!colnames(comm) %in% tree$tip.label)) stop("comm contains tip names that are not available in tree")
    if(is.null(tree$edge.length)) stop("edge lengths are required for tree") 
  }
  
  if(method == 'faith') {
    r_fin <- array(dim = c(resampling, nrow(comm)))
    
    for(i in 1:resampling) {
      v <- sample(1:nrow(comm), nrow(comm))
      x <- comm[v,]
      x <- apply(x, 2, cumsum)
      sub_tree <- apply(x, 1, function(x) drop.tip(tree, which(!(tree$tip.label %in% colnames(comm)[x > 0]))))
      r_fin[i,] <- unlist(lapply(sub_tree, function(x) sum(x$edge.length)))
    }
    
    rare <- colMeans(r_fin)
  }
  
  else if(method == 'barker') {
    com <- sweep(comm, 2, colSums(comm), "/")
    res <- NA
    
    get.leaves <- function(x, st){
      leaves.node <- tips(st, x[2])
    }
    
    r_fin <- do.call(rbind, lapply(1:resampling, function(i) {
      v <- sample(1:nrow(com), nrow(com))
      x <- com[v,]
      x <- apply(x, 2, cumsum)
      sub_tree <- apply(x, 1, function(x) drop.tip(tree, which(!(tree$tip.label %in% colnames(com)[x > 0]))))
      lv <- lapply(sub_tree, function(y) apply(y$edge, 1, function(z) get.leaves(z, y)))
      for(j in 1:nrow(x)) {
        v <- x[j,]
        rel <- unlist(lapply(lv[[j]], function(l) mean(v[l])))
        res[j] <- nrow(sub_tree[[j]]$edge) * ((sum(sub_tree[[j]]$edge.length * rel)) / sum(rel))
      }
      
      return(res)
    }))
    
    rare <- colMeans(r_fin)
  }
  
  else if(method == 'Ia') {
    r_fin <- array(dim = c(resampling, nrow(comm)))
    
    for(i in 1:resampling) {
      v <- sample(1:nrow(comm),nrow(comm))
      x <- comm[v,]
      x <- apply(x, 2, cumsum)
      r_fin[i,] <- pIa(tree, x, exponent = exp)[[1]]
    }
    
    rare <- colMeans(r_fin)
  }
  
  else if(method == 'hill' || method == 'tsallis' || method == 'renyi') {
    r_fin <- array(dim = c(resampling, nrow(comm)))
    
    for(i in 1:resampling) {
      v <- sample(1:nrow(comm), nrow(comm))
      x <- comm[v,]
      x <- apply(x, 2, cumsum)
      r_fin[i,] <- evodivparam(tree, x, method = method, q = exp)
    }
    
    rare <- colMeans(r_fin)
  }
  
  IC_plus <- rare + (1.96 * (sd(r_fin) / sqrt(resampling)))
  IC_neg <- rare - (1.96 * (sd(r_fin) / sqrt(resampling)))
  df <- data.frame(rare, IC_plus, IC_neg)
  colnames(df) <- c('Rarefaction', 'IC_plus', 'IC_neg')
  
  return(df)
  
}

