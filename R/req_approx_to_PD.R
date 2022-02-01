#' @title Required approximation to nearest PD for Coviarance/Correlation Matrix
#' @description Function that test if input matrix object is PD or not. If it is PD, then input matrix is returned.
#' Otherwise, if \emph{flag_force} is \emph{TRUE}, input matrix is approximated to nearest PD, new matrix and
#' normalized Frobenius and Infinity norms are returned. If \emph{flag_force} is \emph{FALSE}, then a warning message is
#' displayed and input matrix is returned as output.
#' @importFrom Matrix nearPD norm
#' @param matr matrix object representing a covariance/correlation matrix
#' @param name string with the name of matrix object representing a covariance/correlation matrix
#' @param flag_force boolean flag for managing input matrices not PD. if \emph{TRUE}, input matrix is approximated
#' to nearest PD if necessary; if \emph{FALSE} approximation to nearest PD is not performed, a warning message is
#' displayed and input matrix is returned as output.
#' @param is_corr boolean flag used to distinguish Correlation to Covariance matrices.
#' @return A list containing:\tabular{ll}{
#'    \code{matr} \tab Matrix object. It contains input correlation/covariance matrix, if the one passed in
#'                     input is PD or if \emph{flag_force} is set to \emph{FALSE}. If input matrix object is not
#'                     PD and  \emph{flag_force} is set to \emph{FALSE}, \code{matr} contains the nearest PD matrix
#'                     to input one.\cr
#'    \tab \cr
#'    \code{normF} \tab Scalar double. It represents the Frobenius Norm of the difference between input matrix and
#'                      its nearest PD, normalized for  the Frobenius Norm of input matrix. If input matrix is PD
#'                      or \emph{flag_force} is set to \emph{FALSE}, this element is not returned in output list.\cr
#'    \tab \cr
#'    \code{normInf} \tab Scalar double. It represents the Infinity norm of the difference between input matrix and
#'                      its nearest PD, normalized for  the Infinity norm of input matrix. If input matrix is PD
#'                      or \emph{flag_force} is set to \emph{FALSE}, this element is not returned in output list.\cr
#' }
#'@author Alessandro De Carlo \email{alessandro.decarlo01@@universitadipavia.it}
#' @seealso \code{\link[Matrix]{nearPD}}
#' @seealso \code{\link[Matrix]{norm}}
#'
#'



req_approx_to_PD<-function(matr,name,flag_force,is_corr){

    out_list <- list()

    if(min(eigen(matr)$values)<0){

        if(flag_force){
            warning(paste("Input",name,"is not PD. It will be approximated to nearest PD. Set force_toPD=FALSE to avoid",sep=" "))
            new_matPD<- Matrix::nearPD(matr,corr = is_corr)
            out_list[["matr"]] <- as.matrix(new_matPD$mat)
            out_list[["normF"]] <- new_matPD$normF/Matrix::norm(matr,"f")
            out_list[["normInf"]] <- Matrix::norm(out_list[["matr"]]-matr,"i")/
                Matrix::norm(matr,"i")
        }else{
            warning(paste("Input",name,"is not PD. Set force_toPD=TRUE to approximate it to the nearest PD",sep=" "))
            out_list[["matr"]] <- matr
        }

    }else{
        out_list[["matr"]] <- matr
    }

    return(out_list)
}
