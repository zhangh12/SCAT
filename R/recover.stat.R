
recover.stat <- function(meta.stat, nsamples, cond, test, ref.geno, ref.info, options){

  snps <- c(cond, test)
  ref.geno <- ref.geno[, snps, drop = FALSE]

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

  tmp <- sapply(meta.stat[snps, 'Direction'], function(u){ifelse(unlist(strsplit(u, ''), use.names = FALSE) == '?', 0, 1)})
  if(is.vector(tmp)){
    tmp <- matrix(tmp, nrow = 1)
  }
  colnames(tmp) <- snps

  es1 <- t(tmp) %*% (tmp * nsamples)
  es2 <- matrix(sum(nsamples), nsnps, nsnps)

  n <- colSums(tmp * nsamples)
  es3 <- outer(n, n, Vectorize(function(a1, a2){min(a1, a2)}))
  se <- meta.stat[snps, 'SE']
  wt1 <- es1 * outer(diag(es1), diag(es1), '*')^(-.5) / outer(se, se, '*')
  wt2 <- es2 * outer(diag(es2), diag(es2), '*')^(-.5) / outer(se, se, '*')
  wt3 <- es3 * outer(diag(es3), diag(es3), '*')^(-.5) / outer(se, se, '*')

  rownames(wt1) <- snps
  colnames(wt1) <- snps
  rownames(wt2) <- snps
  colnames(wt2) <- snps
  rownames(wt3) <- snps
  colnames(wt3) <- snps

  V1 <- wt1 * ref.cor
  V2 <- wt2 * ref.cor
  V3 <- wt3 * ref.cor

  max.total.N <- sum(nsamples)
  rs <- sort(names(score0))
  score0 <- score0[rs] / sqrt(max.total.N)
  V1 <- V1[rs, rs, drop = FALSE] / max.total.N
  V2 <- V2[rs, rs, drop = FALSE] / max.total.N
  V3 <- V3[rs, rs, drop = FALSE] / max.total.N
  ref.cor <- ref.cor[rs, rs, drop = FALSE]

  test.id <- which(names(score0) == test)

  re1 <- try(inv1 <- solve(V1[-test.id, -test.id, drop = FALSE]))
  re2 <- try(inv2 <- solve(V2[-test.id, -test.id, drop = FALSE]))
  re3 <- try(inv3 <- solve(V3[-test.id, -test.id, drop = FALSE]))

  if(any(ref.cor[test.id, -test.id]^2 > options$R2) || error.try(re1)){
    p1 <- NA
  }else{
    U1 <- V1[test.id, -test.id] %*% inv1
    mu1 <- U1 %*% score0[-test.id]
    s1 <- V1[test.id, test.id] - U1 %*% V1[-test.id, test.id]
    p1 <- pchisq((score0[test.id] - mu1)^2/s1, df = 1, lower.tail = FALSE)
  }

  if(any(ref.cor[test.id, -test.id]^2 > options$R2) || error.try(re2)){
    p2 <- NA
  }else{
    U2 <- V2[test.id, -test.id] %*% inv2
    mu2 <- U2 %*% score0[-test.id]
    s2 <- V2[test.id, test.id] - U2 %*% V2[-test.id, test.id]
    p2 <- pchisq((score0[test.id] - mu2)^2/s2, df = 1, lower.tail = FALSE)
  }

  if(any(ref.cor[test.id, -test.id]^2 > options$R2) || error.try(re3)){
    p3 <- NA
  }else{
    U3 <- V3[test.id, -test.id] %*% inv3
    mu3 <- U3 %*% score0[-test.id]
    s3 <- V3[test.id, test.id] - U3 %*% V3[-test.id, test.id]
    p3 <- pchisq((score0[test.id] - mu3)^2/s3, df = 1, lower.tail = FALSE)
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
  # data.frame(cond.rs = paste(cond.rs, collapse = ','), test.rs = test.rs,
  #            cond = paste(cond, collapse = ','), test = test,
  #            p1 = p1, p2 = p2, p3 = p3,
  #            cond.dir = paste(cond.dir, collapse = '/'), test.dir = test.dir,
  #            rho = max.rho,
  #            stringsAsFactors = FALSE)
  
  # data.frame(Cond.SNP = paste(cond.rs, collapse = ','), Test.SNP = test.rs,
  #            Cond.Pos = paste(cond, collapse = ','), Test.Pos = test,
  #            Cond.Dir = paste(cond.dir, collapse = '/'), Test.Dir = test.dir,
  #            Max.R2 = max.rho^2, 
  #            Cond.P = p1, 
  #            stringsAsFactors = FALSE)
  
  
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
             Cond.P = p1, 
             stringsAsFactors = FALSE)

}

