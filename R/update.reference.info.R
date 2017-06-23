
update.reference.info <- function(ref.info, snps){

  ref.info <- ref.info[ref.info$SNP %in% snps, ]
  if(nrow(ref.info) == 0){
    msg <- 'No SNP left'
    stop(msg)
  }

  rownames(ref.info) <- ref.info$SNP

  ref.info

}
