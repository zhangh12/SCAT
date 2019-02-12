
recover.stat <- function(meta.stat, stat, nsamples, cond, test, ref.geno, ref.info, options){

  snps <- c(cond, test)
  ref.geno <- ref.geno[, snps, drop = FALSE]

  maf <- colMeans(ref.geno[, test, drop = FALSE], na.rm = TRUE)/2
  maf <- pmin(maf, 1-maf)
  maf <- formatC(maf, format = 'f', digits = 2)
  names(maf) <- NULL
  
  ref.cor <- cor(ref.geno, use = "pairwise.complete.obs", method = "pearson")
  ref.cor[is.na(ref.cor)] <- 0

  nsnps <- length(snps)
  score0 <- rep(0, nsnps)
  names(score0) <- snps

  rownames(meta.stat) <- meta.stat$SNP.ID

  for(rs in snps){
    score0[rs] <- meta.stat[rs, 'BETA']/meta.stat[rs, 'SE']^2
    if(ref.info[rs, 'RefAllele'] != meta.stat[rs, 'RefAllele']){
      ref.cor[rs, ] <- -ref.cor[rs, ]
      ref.cor[, rs] <- -ref.cor[, rs]
    }
  }

  foo <- function(mat, rid, nsam){
    rid <- intersect(rid, rownames(mat))
    if(length(rid) == 0){
      return(NULL)
    }
    tmp <- sapply(mat[rid, 'Direction'], function(u){ifelse(unlist(strsplit(u, ''), use.names = FALSE) == '?', 0, 1)})
    if(is.vector(tmp)){
      tmp <- matrix(tmp, nrow = 1)
    }
    colnames(tmp) <- rid
    
    es <- t(tmp) %*% (tmp * nsam)
    se <- mat[rid, 'SE']
    wt <- es * outer(diag(es), diag(es), '*')^(-.5) / outer(se, se, '*')
    rownames(wt) <- rid
    colnames(wt) <- rid
    wt
    
  }
  
  nstudy <- length(stat)
  wt <- matrix(0, length(snps), length(snps))
  rownames(wt) <- snps
  colnames(wt) <- snps
  max.total.N <- 0
  kk = 0
  for(i in 1:nstudy){
    kk <- kk+1
    max.total.N <- max.total.N + sum(nsamples[[i]])
    tmp <- foo(stat[[i]], snps, nsamples[[i]])
    if(is.null(tmp)){
      next
    }
    rs <- rownames(tmp)
    wt[rs, rs] <- wt[rs, rs] + tmp
  }

  V <- wt * ref.cor
  rs <- sort(names(score0))
  score0 <- score0[rs] / sqrt(max.total.N)
  V <- V[rs, rs, drop = FALSE] / max.total.N
  ref.cor <- ref.cor[rs, rs, drop = FALSE]

  test.id <- which(names(score0) == test)

  re <- try(inv <- solve(V[-test.id, -test.id, drop = FALSE]))

  if(any(ref.cor[test.id, -test.id]^2 > options$R2) || error.try(re)){
    p <- NA
  }else{
    U <- V[test.id, -test.id] %*% inv
    mu <- U %*% score0[-test.id]
    s <- V[test.id, test.id] - U %*% V[-test.id, test.id]
    p <- pchisq((score0[test.id] - mu)^2/s, df = 1, lower.tail = FALSE)
  }
  
  cond.pval <- meta.stat[cond, 'P']
  test.pval <- meta.stat[test, 'P']
  
  cond.dir <- meta.stat[cond, 'Direction']
  test.dir <- meta.stat[test, 'Direction']

  cond.rs <- ref.info[cond, 'SNP']
  test.rs <- ref.info[test, 'SNP']

  rho <- as.vector(ref.cor[test.id, -test.id])
  max.rho <- max(abs(rho))
  sgn.rho <- ifelse(rho[which.max(abs(rho))[1]] > 0, '+', '-')
  
  reformat <- function(x){
    formatC(x, format = 'e', digits = 1)
  }
  
  # use Idx (index) SNP to refer the SNP being conditioned on
  data.frame(Idx.SNP = paste(cond.rs, collapse = ','), Test.SNP = test.rs, 
             Idx.Pos = paste(cond, collapse = ','), Test.Pos = test, 
             Idx.Dir = paste(cond.dir, collapse = '/'), Test.Dir = test.dir, 
             Idx.P = paste(reformat(cond.pval), collapse = ','), Test.P = reformat(test.pval), 
             Max.R2 = formatC(max.rho^2, format = 'f', digits = 2), 
             Cor.Dir = sgn.rho, 
             Cond.P = p, 
             MAF = maf, 
             stringsAsFactors = FALSE)
  
}

