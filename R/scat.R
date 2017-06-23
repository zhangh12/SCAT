

scat <- function(summary.files, model, reference, lambda, nsamples, min.maf = 0.05, max.R2 = 0.90){

  options <- list(maf = min.maf, R2 = max.R2)
  
  # we don't check format at this stage, will add it later
  validate.input(summary.files, reference, lambda, nsamples)

  m <- meta(summary.files, nsamples, lambda)
  meta.stat <- m$meta.stat
  nsamples <- m$nsamples

  model <- parse.model(model, meta.stat$SNP.ID)

  ref.info <- load.reference.info(reference, model)

  model <- update.model(model, ref.info$SNP.ID)

  meta.stat <- meta.stat[meta.stat$SNP.ID %in% ref.info$SNP.ID, ]

  cond.test(meta.stat, reference, nsamples, model, ref.info, options)

}


