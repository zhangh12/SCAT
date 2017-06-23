

load.reference.info <- function(reference, model){

  msg <- paste("Loading reference information:", date())
  message(msg)

  cond <- unlist(sapply(model$cond, function(u){unlist(strsplit(u, ','), use.names = FALSE)}), use.names = FALSE)
  test <- unlist(sapply(model$test, function(u){unlist(strsplit(u, ','), use.names = FALSE)}), use.names = FALSE)
  snp.id <- unique(c(cond, test))

  if("matrix" %in% class(reference)){
    reference <- as.data.frame(reference)
  }

  ref.info <- NULL
  nref <- nrow(reference)
  col.class <- rep('NULL', 6)
  col.class[c(2, 5, 6)] <- 'character'
  col.class[c(1, 4)] <- 'integer'
  bim.files <- reference$bim
  for(i in 1:nref){
    tmp <- try(bim <- read.table(bim.files[i], header = FALSE, as.is = TRUE, colClasses = col.class), silent = TRUE)
    if(error.try(tmp)){
      msg <- paste0('Cannot load ', bim.files[i])
      stop(msg)
    }

    colnames(bim) <- c("Chr", "SNP", "Pos", "RefAllele", "EffectAllele")
    bim$SNP.ID <- paste(bim$Chr, ':', bim$Pos, sep = '')
    bim$Reference.ID <- i
    bim$RefAllele <- toupper(bim$RefAllele)
    bim$EffectAllele <- toupper(bim$EffectAllele)
    bim <- bim[bim$SNP.ID %in% snp.id, ]
    ref.info <- rbind(ref.info, bim)

    rm(bim)
    gc()
  }

  if(is.null(ref.info)){
    msg <- "No SNPs are found in bim files"
    stop(msg)
  }

  ref.info <- ref.info[order(ref.info$Reference.ID), ]
  rownames(ref.info) <- ref.info$SNP.ID

  ref.info

}

