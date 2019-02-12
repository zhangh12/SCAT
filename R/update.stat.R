
update.stat <- function(meta.stat, stat, lambda, nsamples, snp.id){
  
  if(!is.null(meta.stat)){
    id <- which(meta.stat$SNP.ID %in% snp.id)
    meta.stat <- meta.stat[id, ]
  }
  
  nstudy <- length(stat)
  
  lam <- NULL
  nsam <- list()
  st <- list()
  k <- 0
  for(i in 1:length(stat)){
    if(!any(snp.id %in% stat[[i]][, 'SNP.ID'])){
      next
    }
    lam <- c(lam, lambda[i])
    k <- k + 1
    nsam[[k]] <- nsamples[[k]]
    id <- which(stat[[i]]$SNP.ID %in% snp.id)
    st[[i]] <- stat[[i]][id, ]
    st[[i]]$SNP.ID <- NULL
  }
  
  if(k == 0){
    stop()
  }
  
  list(meta.stat = meta.stat, 
       stat = st, 
       lambda = lam, 
       nsamples = nsam)
  
}
