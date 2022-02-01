#' @title Convert Normal parameters to Log-Normal parameters
#' @description This function is used to convert mean vector and covariance matrix of a multivariate Normal
#' distribution to mean vector and covariance matrix of the associated multivariate Log-Normal distribution.
#' @param mu Array object containing means of multivariate Normal distribution
#' @param covMatrix Matrix object containing covariance matrix of multivariate Normal distribution
#' @export
#' @return A list containing:\tabular{ll}{
#' \code{muLn} \tab Array object. It contains mean of multivariate Log-Normal distribution associated to
#'                    Normal one.\cr
#'    \tab \cr
#'    \code{sigmaLn} \tab Matrix object. It contains covariance matrix of multivariate Log-Normal distribution
#'    associated to Normal one.\cr
#' }
#' @author Alessandro De Carlo \email{alessandro.decarlo01@@universitadipavia.it}
#' @seealso \code{\link{logn_to_normal}}
#' @examples
#' #define correlations
#' corr<- diag(rep(1,4))
#' corr[1,4] <- 0.9
#' corr[4,1]<-corr[1,4]
#' corr[2,4] <- -0.3
#' corr[4,2] <- corr[2,4]
#' corr[3,2] <- -0.2
#' corr[2,3] <- corr[3,2]


#' #define sd of variables
#' sd2 <- array(c(rep(1,4)))
#' #obtain covariance matrix
#' covMatrix2 <- sd2%*%t(sd2)*corr
#' #define mean vector
#' mu2 <- array(rep(2.5,4))
#' normal_to_logn(mu2,covMatrix2)
#' #output:
#'
#' # $muLn
#' # [1] 20.08554 20.08554 20.08554 20.08554
#' #
#' # $sigmaLn
#' #           [,1]   [,2]      [,3]      [,4]
#' # [1,] 693.2044    0.00000   0.00000  588.8459
#' # [2,]   0.0000  693.20436 -73.12923 -104.5614
#' # [3,]   0.0000  -73.12923 693.20436    0.0000
#' # [4,] 588.8459 -104.56139   0.00000  693.2044




normal_to_logn <- function(mu,covMatrix){

    validate_covMatrix(covMatrix)
    validate_array("mu",mu,covMatrix,"covMatrix")

    muLn <- array(c(rep(0,length(mu))))
    names(muLn) <- names(mu)
    sigmaLn <- matrix(0L, ncol = length(mu), nrow = length(mu))
    colnames(sigmaLn) <- colnames(covMatrix)
    row.names(sigmaLn) <- row.names(covMatrix)

    for (i in 1:length(mu)) {
        for (j in 1:length(mu)) {
            sigmaLn[i,j] <- (exp(mu[i]+mu[j]+0.5*covMatrix[i,i]+0.5*covMatrix[j,j]))*(exp(covMatrix[i,j])-1)
        }
    }

    for (i in 1:length(mu)) {
        muLn[i] <- exp(mu[i]+0.5*covMatrix[i,i])
    }

    out_list <- list()
    out_list[["muLn"]] <-muLn
    out_list[["sigmaLn"]] <- sigmaLn
    return(out_list)
}
