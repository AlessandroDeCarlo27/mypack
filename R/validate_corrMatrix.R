#' @title Validate input correlation matrix
#' @description function used for check that input matrix is a matrix object which represents a valid correlation
#' matrix. If these requisites are not respected, excecution will be stopped.
#' @param corrMatStart object that represents input correlation matrix to be validated.
#' @note Eigenvalues are non considered because the functions in which \emph{validate_corrMatrix} is called
#' can handle not PD matrices too by approximating them to nearest PD.
#' @author Alessandro De Carlo \email{alessandro.decarlo01@@universitadipavia.it}
#' @seealso \code{\link[Matrix]{nearPD}}
#' @seealso \code{\link{validate_array}}
#' @seealso \code{\link{validate_covMatrix}}


validate_corrMatrix <- function(corrMatStart){
    #corrMatStart must be a matrix
    if (!is.matrix(corrMatStart)) {
        stop("A matrix object must be passed as input.")
    }
    #corrMatStart must be symmetric
    if(!isSymmetric(corrMatStart)){
        stop("Input matrix must be symmetric.")
    }

    #trace of corrMatStart must be equal to #column of matrix
    if(sum(diag(corrMatStart))!=ncol(corrMatStart)){
        stop("Invalid correlation matrix. Elements along diagonal must be equal to 1")
    }

    if(any(corrMatStart[upper.tri(corrMatStart,diag = F)]>1,na.rm = T)){
        stop("Invalid correlation matrix. Correlations >1 were found")
    }

}
