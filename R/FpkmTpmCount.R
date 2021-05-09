#' TPM FPKM COUNT conversion
#'
#' TPM FPKM COUNT conversion
#'
#' @param counts raw count of RNAseq data.
#' @param fpkm FPKM.
#' @return generate a expression matrix (TPM or FPKM).
#' @export
#' @keywords dataConversion
#' @examples
#' tpm <- countToTpm(counts = counts)
#' fpkm <- countToFpkm(counts = counts)
#' tpm <- fpkmToTpm(counts = fpkm)


countToTpm <- function(counts, effLen=5)
{
    rate <- log(counts) - log(effLen)
    denom <- log(sum(exp(rate)))
    exp(rate - denom + log(1e6))
}
 
countToFpkm <- function(counts, effLen=5)
{
    N <- sum(counts)
    exp( log(counts) + log(1e9) - log(effLen) - log(N) )
}
 
fpkmToTpm <- function(fpkm)
{
    exp(log(fpkm) - log(sum(fpkm)) + log(1e6))
}

fff <- function(x, r=6, c=6)
{ 
  x[1:r,1:c]
}

ttt <- function(x, r=6, c=6)
{ 
  x[1:r,(ncol(x)-c):ncol(x)]
}

ftr <- function(x)
{ 
  a <- x[,2:ncol(x)]
  rownames(a) <- x[,1]
  return(a)
}