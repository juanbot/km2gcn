#' Function that generates the correlation for distance, in line
#' with a signed network
#' @param a numerical vector with gene expression for gene a
#' @param b numerical vector with gene expression for gene a
#' @param signed is the network signed?
#' @return The signed correlation if signed is TRUE
#' @export
corDistance = function(a,b,signed=TRUE){
  if(signed)
    return(0.5 * (1 + WGCNA::corFast(a,b)))
  return(abs(WGCNA::corFast(a,b)))
}
