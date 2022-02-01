#' @title Validate input covariance matrix
#' @description function used for check that input matrix is a matrix object which represents a valid covariance
#' matrix. If these requisites are not respected, excecution will be stopped.
#' @param covMatrix object that represents input covariance matrix to validate.
#' @note Eigenvalues are non considered because the package can handle not PD matrices too
#' by approximating them to nearest PD.
#' @author Alessandro De Carlo \email{alessandro.decarlo01@@universitadipavia.it}
#'  @seealso \code{\link[Matrix]{nearPD}}
#' @seealso \code{\link{validate_array}}
#' @seealso \code{\link{validate_corrMatrix}}
#'


validate_covMatrix <- function(covMatrix){
    if(!is.matrix(covMatrix)){
        stop("covMatrix must be a matrix object.")
    }

    if(!isSymmetric(covMatrix)){
        stop("covMatrix must be symmetric.")
    }

}
