# only use (test or condition on) SNPs available in PLINK files
update.model <- function(model, snps){

  for(i in 1:nrow(model)){
    cond <- unlist(strsplit(model$cond[i], ','), use.names = FALSE)
    cond <- intersect(cond, snps)
    model$cond[i] <- paste(cond, collapse = ',', sep = '')

    test <- unlist(strsplit(model$test[i], ','), use.names = FALSE)
    test <- intersect(test, snps)
    model$test[i] <- paste(test, collapse = ',', sep = '')
  }

  id <- apply(model, 1, function(r){all(r != '')})

  if(!any(id)){
    msg <- 'No model can be tested'
    stop(msg)
  }

  model <- model[id, ]

  model

}
