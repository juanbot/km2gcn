#' The function selects the best module for a gene on the basis of the highest
#' correlation of that gene's expression with the eigengenes
#' @param gene numerical vector with the gene expression
#' @param centroids Matrix with many rows as \code{gene} components and a column for each module eigengene
#' @param signed is the network signed?
#' @return The index of the eigengene within the matrix passed as argument
#' @export
getBestModuleCor <- function(gene,centroids,signed=TRUE){

  return(which.max(corDistance(centroids,gene,signed)))
}
