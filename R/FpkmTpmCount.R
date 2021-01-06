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

ff <- function(x)
{ 
  x[1:5,1:5]
}