
parse.model <- function(model, snp.id){

  msg <- paste("Parsing models:", date())
  message(msg)

  foo <- function(u){
    u <- unlist(strsplit(u, ':'))
    chr <- as.integer(u[1])
    pos <- as.integer(unlist(strsplit(u[2], '-')))
    if(length(pos) == 1){
      pos <- c(pos, pos)
    }
    pos <- c(min(pos), max(pos))
    paste(chr, ':', pos[1]:pos[2], sep = '')
  }

  cond <- NULL
  test <- NULL
  model$cond <- gsub(' ', '', model$cond)
  model$test <- gsub(' ', '', model$test)
  for(i in 1:nrow(model)){
    cd <- unlist(strsplit(model$cond[i], ','), use.names = FALSE)
    cd <- unlist(sapply(cd, foo), use.names = FALSE)
    cd <- intersect(cd, snp.id)
    if(length(cd) == 0){
      next
    }
    #model$cond[i] <- paste(cond, collapse = ',', sep = '')
    tt <- unlist(strsplit(model$test[i], ','), use.names = FALSE)
    tt <- as.vector(unlist(sapply(tt, foo), use.names = FALSE))
    tt <- setdiff(tt, cd)
    tt <- intersect(tt, snp.id)
    while(length(tt) > 0){
      ntt <- min(1e3, length(tt))
      test <- c(test, paste(tt[1:ntt], collapse = ',', sep = ''))
      cond <- c(cond, paste(cd, collapse = ',', sep=''))
      tt <- tt[-c(1:ntt)]
    }
    #model$test[i] <- paste(test, collapse = ',', sep = '')
  }

  model <- data.frame(cond, test, stringsAsFactors = FALSE)

  model

}
