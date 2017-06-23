

merge.stat <- function(stat, ref.allele, conf.snps, lambda){

  msg <- paste("Merging summary statistics:", date())
  message(msg)
  
  stat <- update.direction(stat, ref.allele)

  RefAllele <- ref.allele$RefAllel
  EffectAllele <- ref.allele$EffectAllele

  snps <- names(RefAllele)
  nsnp <- length(snps)

  BETA <- rep(0, nsnp)
  SE <- rep(0, nsnp)
  Chr <- rep(0, nsnp)
  Pos <- rep(0, nsnp)
  names(BETA) <- snps
  names(SE) <- snps
  names(Chr) <- snps
  names(Pos) <- snps

  nstudy <- length(stat)
  for(i in 1:nstudy){
    s <- stat[[i]][, 'SNP']
    stat[[i]]$sgn <- ifelse(stat[[i]][, 'RefAllele'] == RefAllele[s] & stat[[i]][, 'EffectAllele'] == EffectAllele[s], 1, -1)
    stat[[i]][, 'SE'] <- stat[[i]][, 'SE'] * sqrt(lambda[i])
    stat[[i]][, 'P'] <- pchisq(stat[[i]][, 'BETA']^2/stat[[i]][, 'SE']^2, df = 1, lower.tail = FALSE)
    BETA[s] <- BETA[s] + stat[[i]][, 'sgn'] * stat[[i]][, 'BETA']/stat[[i]][, 'SE']^2
    SE[s] <- SE[s] + 1/stat[[i]][, 'SE']^2
    Chr[s] <- stat[[i]][, 'Chr']
    Pos[s] <- stat[[i]][, 'Pos']
  }

  SE <- sqrt(1/SE)
  BETA <- BETA * SE^2
  P <- pchisq(BETA^2/SE^2, df = 1, lower.tail = FALSE)

  SNP <- as.character(names(BETA))
  meta.stat <- data.frame(SNP = SNP, Chr = Chr, Pos = Pos,
                          RefAllele = RefAllele[SNP], EffectAllele = EffectAllele[SNP],
                          BETA = BETA[SNP], SE = SE[SNP], P = P[SNP], stringsAsFactors = FALSE)
  rownames(meta.stat) <- NULL

  header <- c('SNP', 'Chr', 'Pos', 'RefAllele', 'EffectAllele', 'BETA', 'SE', 'P')
  meta.stat <- meta.stat[, header]

  Direction <- rep('', nsnp)
  names(Direction) <- meta.stat$SNP
  for(i in 1:nstudy){
    nc <- nchar(stat[[i]][1, 'Direction'])
    d <- rep(paste0(rep('?', nc), collapse = ''), nsnp)
    names(d) <- meta.stat$SNP
    d[stat[[i]][, 'SNP']] <- stat[[i]][, 'Direction']
    Direction <- paste(Direction, d, sep = '')
    names(Direction) <- meta.stat$SNP
  }
  meta.stat$Direction <- Direction

  if(!is.null(conf.snps)){
    meta.stat <- meta.stat[!(meta.stat$SNP %in% conf.snps), ]
  }

  meta.stat <- meta.stat[order(meta.stat$P), ]
  meta.stat$SNP.ID <- paste(meta.stat$Chr, ':', meta.stat$Pos, sep = '')
  rownames(meta.stat) <- meta.stat$SNP.ID

  # list(meta.stat = meta.stat, conf.snps = conf.snps)
  meta.stat

}

