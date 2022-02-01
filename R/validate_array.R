#' @title Validate vectors of statistical parameters (mean, sd, cv)
#' @description This function is used to validate vectors of statistical parameters (mean, sd, cv) which are typically
#' the input of the functions in this package. They must be array objects without NA, NaN or infinite values and
#' their length must be equal to the size of the associated covariance matrix/correlation matrix (square matrix).
#' If these requisites are not respected, excecution will be stopped.
#' @param name_arr string with name of the object that represent input array to validate
#' @param arr object that represent input array to validate
#' @param matRef matrix object that represent correlation matrix/covariance matrix
#' @param name_mat string with name of the object that represent correlation matrix/covariance matrix
#' @author Alessandro De Carlo \email{alessandro.decarlo01@@universitadipavia.it}
#' @seealso  \code{\link{validate_corrMatrix}}
#' @seealso  \code{\link{validate_covMatrix}}
#'


validate_array <- function(name_arr,arr,matRef,name_mat){
    if(!is.array(arr)){
        stop(paste(name_arr,"must be an array object",sep=" "))
    }

    if(length(arr)!=ncol(matRef)|length(arr)!=nrow(matRef)){
        stop(paste("Wrong dimensions: number of elements in", name_arr,"must match the dimensions of", name_mat,".",sep=" "))
    }

    if(any(is.na(arr))||any(is.nan(arr))||any(is.infinite(arr))){
        stop(paste(name_arr,"contains illegal values (NA,NaN Inf)."))
    }


}
