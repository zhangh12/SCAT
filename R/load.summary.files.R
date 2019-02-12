
# inflation lambda is adjusted here
load.summary.files <- function(summary.files, model, lambda, nsamples){

  msg <- paste("Loading summary files:", date())
  message(msg)

  header <- c('SNP', 'Chr', 'Pos', 'RefAllele', 'EffectAllele', 'BETA') # columns that must be provided by users
  opt.header <- c('P', 'SE')

  complete.header <- c(header, opt.header, 'Direction')

  nfiles <- length(summary.files)
  stat <- list()

  snp.id <- NULL
  fid <- 0
  for(i in 1:nfiles){
    st <- read.table(summary.files[i], header = TRUE, as.is = TRUE, nrows = 1e4)
    header.map <- colnames(st)
    colnames(st) <- convert.header(colnames(st), complete.header)
    tmp <- (header %in% colnames(st))
    if(!all(tmp)){
      msg <- paste0("Columns below were not found in ", summary.files[i], ":\n", paste(header[!tmp], collapse = " "))
      stop(msg)
    }
    names(header.map) <- colnames(st)

    col.class <- sapply(st, class)
    col.id <- which(colnames(st) %in% complete.header)
    col.class[-col.id] <- "NULL"
    col.class[c('SNP', 'RefAllele', 'EffectAllele')] <- 'character'
    col.class[c('Chr', 'Pos')] <- 'integer'
    names(col.class) <- header.map[names(col.class)]
    #st <- read.table(summary.files[i], header = TRUE, as.is = TRUE, colClasses = col.class)
    st <- data.table::setDF(data.table::fread(summary.files[i], header = TRUE, showProgress = FALSE, verbose = FALSE, select = which(col.class != 'NULL')))
    colnames(st) <- convert.header(colnames(st), complete.header)

    if(!any(opt.header %in% colnames(st))){
      msg <- paste0("Neither SE nor P is not provided in ", summary.files[i])
      stop(msg)
    }


    if(!('P' %in% colnames(st))){
      st$P <- NA
    }

    if(!('SE' %in% colnames(st))){
      st$SE <- NA
    }

    if(!('Direction' %in% colnames(st))){
      msg <- paste0('Direction is absent in ', summary.files[i], '. Function scat() therefore assumed equal sample sizes for all SNPs in that file. Please verify if this assumption is reasonable in your data. Violation of this assumption can lead to false positive in conditional analysis. If this assumption is valid, feel free to ignore this warning. ')
      warning(msg, immediate. = TRUE)
      st$Direction <- ifelse(st$BETA == 0, '0', ifelse(st$BETA > 0, '+', '-'))
    }

    nc <- unique(nchar(st$Direction))
    if(length(nc) != 1){
      msg <- paste0('String lengths of Direction are unequal in ', summary.files[i])
      stop(msg)
    }

    st <- st[, which(toupper(colnames(st)) %in% toupper(complete.header))]

    dup <- duplicated(st$SNP)
    if(any(dup)){
      dup.snps <- unique(st$SNP[dup])
      msg <- paste("SNPs below are duplicated: ", paste(dup.snps, collapse = " "))
      stop(msg)
    }

    id.no.SE.P <- which(is.na(st$SE) & is.na(st$P))
    if(length(id.no.SE.P) > 0){
      msg <- paste("For SNPs below, neither SE nor P is not provided in", summary.files[i], ":\n", paste(st$SNP[id.no.SE.P], collapse = " "))
      stop(msg)
    }

    st$RefAllele <- toupper(st$RefAllele)
    st$EffectAllele <- toupper(st$EffectAllele)

    id.no.SE <- which(is.na(st$SE))
    id.no.P <- which(is.na(st$P))

    if(length(id.no.SE) > 0){
      z2 <- qchisq(st$P[id.no.SE], df = 1, lower.tail = FALSE)
      st$SE[id.no.SE] <- abs(st$BETA[id.no.SE]/sqrt(z2))
    }
    
    st$SE <- st$SE * sqrt(lambda[i])
    st$P <- pchisq((st$BETA/st$SE)^2, df = 1, lower.tail = FALSE)

    fid <- fid + 1
    rownames(st) <- st$SNP
    st <- st[complete.cases(st), ]

    st$SNP.ID <- paste0(st$Chr, ':', st$Pos)
    stat[[fid]] <- st
    snp.id <- unique(c(snp.id, st$SNP.ID))
    rm(st)
    gc()

  }

  if(length(stat) == 0){
    msg <- "No SNPs to be included in analysis"
    stop(msg)
  }

  # snp.id <- unique(snp.id)
  
  model <- parse.model(model, snp.id)

  snp.id <- unique(unlist(strsplit(c(model$cond, model$test), ',')))
  
  m <- update.stat(NULL, stat, lambda, nsamples, snp.id)
  rm(snp.id)
  gc()
  
  list(stat = m$stat, lambda = m$lambda, nsamples = m$nsamples)

}

