
cond.test <- function(meta.stat, reference, nsamples, model, ref.info, options){

  res <- NULL
  for(i in 1:nrow(model)){
    msg <- paste0("Performing conditional test (", i, '/', nrow(model), '): ', date())
    message(msg)

    cond <- unlist(strsplit(model$cond[i], ','), use.names = FALSE)
    test <- unlist(strsplit(model$test[i], ','), use.names = FALSE)
    snps <- c(cond, test)
    ref <- ref.info[snps, ]
    ref.id <- unique(ref$Reference.ID)
    ref.geno <- NULL
    for(j in ref.id){
      if(is.null(ref.geno)){
        ref.geno <- read.bed(reference$bed[j], reference$bim[j], reference$fam[j], sel.snps = ref$SNP[ref$Reference.ID == j])
      }else{
        ref.geno <- cbind(ref.geno, read.bed(reference$bed[j], reference$bim[j], reference$fam[j], sel.snps = ref$SNP[ref$Reference.ID == j]))
      }
    }

    names(snps) <- ref$SNP
    colnames(ref.geno) <- snps[colnames(ref.geno)]
    maf <- apply(ref.geno, 2, function(u){m <- mean(u, na.rm = TRUE)/2; min(m, 1-m)})
    rm.snps <- colnames(ref.geno)[which(maf < options$maf)]
    if(length(rm.snps) > 0){
      cond <- setdiff(cond, rm.snps)
      test <- setdiff(test, rm.snps)
      if(length(cond) == 0 || length(test) == 0){
        next
      }
      model$cond[i] <- paste(cond, sep = '', collapse = ',')
      model$test[i] <- paste(test, sep = '', collapse = ',')
      snps <- c(cond, test)
      ref.geno <- ref.geno[, snps, drop = FALSE]
    }

    for(s in test){
      res <- rbind(res, recover.stat(meta.stat, nsamples, cond, s, ref.geno, ref.info, options))
    }
    #print(dim(res))
  }

  res

}
