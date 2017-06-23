

scat <- function(summary.files, model, reference, lambda, nsamples, options = NULL){

  # we don't check format at this stage
  validate.input(summary.files, reference, lambda, nsamples)

  m <- meta(summary.files, nsamples, lambda)
  meta.stat <- m$meta.stat
  nsamples <- m$nsamples

  model <- parse.model(model, meta.stat$SNP.ID)

  ref.info <- load.reference.info(reference, model)

  model <- update.model(model, ref.info$SNP.ID)

  meta.stat <- meta.stat[meta.stat$SNP.ID %in% ref.info$SNP.ID, ]

  #model <- convert.model(ref.info, model)

  print('0.1.2')
  cond.test(meta.stat, reference, nsamples, model, ref.info, options)

}


