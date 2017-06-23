
convert.model <- function(ref.info, model){

  for(i in 1:nrow(model)){
    cond <- unlist(strsplit(model$cond[i], ','), use.names = FALSE)
    cond <- intersect(cond, ref.info$SNP.ID)
    model$cond[i] <- paste(ref.info[cond, 'SNP'], collapse = ',', sep = '')

    test <- unlist(strsplit(model$test[i], ','), use.names = FALSE)
    test <- intersect(test, ref.info$SNP.ID)
    model$test[i] <- paste(ref.info[test, 'SNP'], collapse = ',', sep = '')
  }

  id <- apply(model, 1, function(r){all(r != '')})

  if(!any(id)){
    msg <- 'No model can be tested'
    stop(msg)
  }

  model <- model[id, ]

  model

}

