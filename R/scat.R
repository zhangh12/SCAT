

scat <- function(summary.files, model, reference, lambda, nsamples, min.maf = 0.05, max.R2 = 0.90){

  options <- list(maf = min.maf, R2 = max.R2)
  
  # we don't check format at this stage, will add it later
  validate.input(summary.files, reference, lambda, nsamples)

  # SNPs may be totally missing in some summary files, 
  # so lambda and nsamples need to be updated (discard some elements)
  m <- meta(summary.files, model, lambda, nsamples)

  model <- parse.model(model, m$meta.stat$SNP.ID)

  ref.info <- load.reference.info(reference, model)

  model <- update.model(model, ref.info$SNP.ID)

  m <- update.stat(m$meta.stat, m$stat, m$lambda, m$nsamples, ref.info$SNP.ID)

  cond.test(m$meta.stat, m$stat, reference, m$nsamples, model, ref.info, options)

}


