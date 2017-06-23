
validate.input <- function(summary.files, reference, lambda, nsamples){

  if(!is.null(summary.files)){
    validate.summary.files(summary.files)
  }

  if(!is.null(reference)){
    validate.reference(reference)
  }

  if(!is.null(summary.files) && !is.null(lambda)){
    validate.lambda(summary.files, lambda)
  }

  if(!is.null(lambda) && !is.null(nsamples)){
    validate.sample.size(lambda, nsamples)
  }

}
