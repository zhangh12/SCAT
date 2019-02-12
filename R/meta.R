
meta <- function(summary.files, model, lambda, nsamples){

  validate.summary.files(summary.files)

  lambda <- validate.lambda(summary.files, lambda)

  sf <- load.summary.files(summary.files, model, lambda, nsamples)

  ref.allele <- extract.reference.allele(sf$stat)

  conf.snps <- extract.conflictive.snps(sf$stat, ref.allele)

  rcs <- remove.conflictive.snps(sf$stat, ref.allele, conf.snps)

  m <- merge.stat(rcs$stat, rcs$ref.allele, conf.snps)
  
  list(meta.stat = m$meta.stat, stat = m$stat, lambda = sf$lambda, nsamples = sf$nsamples)
  

}


