
meta <- function(summary.files, nsamples, lambda = NULL, sel.snps = NULL){

  validate.summary.files(summary.files)

  lambda <- validate.lambda(summary.files, lambda)

  sf <- load.summary.files(summary.files, lambda, nsamples, sel.snps)

  ref.allele <- extract.reference.allele(sf$stat)

  conf.snps <- extract.conflictive.snps(sf$stat, ref.allele)

  rcs <- remove.conflictive.snps(sf$stat, ref.allele, conf.snps)

  meta.stat <- merge.stat(rcs$stat, rcs$ref.allele, conf.snps, sf$lambda)

  list(meta.stat = meta.stat, nsamples = sf$nsamples)

}


