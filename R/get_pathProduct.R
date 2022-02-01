#' @title Compute the product of all correlations along a path
#' @description This function returns the product of all correlations along a path between two nodes in the correlations
#' graph.
#' @param corrMat Matrix Object. It represent the correlation matrix
#' @param nodes Vector. It contains all the nodes of the path that links \code{nodes[1]} and \code{nodes[length(nodes)-1]}
#' @return Scalar double which represents the product of all correlations along a path between two nodes in the
#' correlations graph.
#' @author Alessandro De Carlo \email{alessandro.decarlo01@@universitadipavia.it}


get_pathProduct <- function(corrMat,nodes){

    corr <- array(0L,dim = length(nodes)-1)
    for (i in 1:length(nodes)-1) {
        corr[i] <- corrMat[nodes[i],nodes[i+1]]
    }
    return(prod(corr))
}
