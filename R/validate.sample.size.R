
validate.sample.size <- function(lambda, nsamples){

  if(!is.list(nsamples)){
    msg <- 'nsamples should be a list'
    stop(msg)
  }

  if(length(nsamples) != length(lambda)){
    msg <- 'Length of nsamples and lambda should be equal'
    stop(msg)
  }

}
