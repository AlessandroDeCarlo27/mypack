#' @title Sampling random numbers from a multivariate Log-Normal distribution
#' @description This function allows to generate random samples from a Log-Normal distribution using the parameters
#' of Log-Normal distribution instead of the ones of underlying Normal.
#' @importFrom stats qnorm runif
#' @importFrom Matrix nearPD norm
#' @export
#' @param Nsamples Number of random samples to extract from multivariate Log-Normal distribution
#' @param force_toPD boolean flag. If input correlation/covariance matrix is not PD, it is approximated to its
#'                   nearest PD when this flag is set to \emph{TRUE} (default value). If it is set to \emph{FALSE},
#'                   approximation to nearest PD is not performed and a warning message is printed.
#' @param full_output boolean flag. If it is set to \emph{TRUE} (default value), in output list is returned also
#' the covariance matrix of Normal distribution underlying the Log-Normal one. If the covariance matrix of Normal
#' distribution is not PD and this flag is activated, output list will contain the initial version of Normal
#' covariance matrix, its nearest PD and the normalized Frobenius and Infinity norms of the difference.
#' @param ... One of the following five combination of named arguments:
#' \enumerate{
#' \item \code{mu=input_mean, covMatrix = input_covariance}
#' \item \code{mu=input_mean, sd=input_sd, corrMatrix=input_correlations}
#' \item \code{mu=input_mean, cv=input_cv, corrMatrix=input_correlations}
#' \item \code{sd=input_sd, cv_input_cv, corrMatrix=input_correlations}
#' \item \code{sd=input_sd, cv_input_cv, covMatrix = input_covariance}
#' }
#' where:
#' \itemize{
#' \item \code{mu}: Array object containing means of multivariate Log-Normal distribution
#' \item \code{sd}: Array object containing standard deviations of multivariate Log-Normal distribution
#' \item \code{cv}: Array object containing CVs of multivariate Log-Normal distribution. CV are \bold{not} expressed
#' with percentage.
#' \item \code{covMatrix}: Matrix object containing covariance matrix of multivariate Log-Normal distribution
#' \item \code{corrMatrix}: Matrix object containing correlation matrix of multivariate Log-Normal distribution
#' }
#'
#'
#' @details
#' It is necessary to specify as input of this function the correlation/covariance structure of Log-Normnal
#' variables. If this matrix is not PD, it will be approximated to its nearest PD using \code{nearPD} function
#' of \code{Matrix} package. Details about the approximation performed will be provided to user in output list
#' with the nearest PD of input matrix and Forbenius and Infinity normalized norms of their differences.
#'
#' This function evaluate two conditions before computing the parameters of Normal distribution:
#' \enumerate{
#' \item For each couple of variables \eqn{X1}, \eqn{X2} is tested if \eqn{\rho12} is within a certain range
#' whose bounds are computed according to the CVs of \eqn{X1} and \eqn{X2}. This is a necessary condition (becomes
#' sufficent too when correlation matrix is a 2x2) for obtaining a PD correlation/covariance matrix for the
#' Normal distribution underlying the Log-Normal one. If it is not satisfied, a warning message is displayed
#' but covariance matrix of Normal distribution is still computed even if it is not PD.
#' \item For each couple of variables \eqn{X1}, \eqn{X2} is tested if \eqn{\rho12*cv1*cv2 > -1}. If exists a couple
#' which does not respect this condition, it means that covariance matrix of underlying Normal distribution cannot
#' be computed. In fact, the covariance between \eqn{X1}, \eqn{X2} of Normal distribution, \eqn{\sigma^2}12,
#' is defined as \eqn{ln((\rho12*cv1*cv2 )+1)}. In this case, a warning message will be displayed and a recap of
#' validation performed is returned but the sampling will not be performed.
#' }
#'
#' If condition 2 is fulfilled, independently by contition 1, sampling is guaranteed. If Normal covariance matrix
#' is not PD it is approximated to its nearest PD using \code{nearPD} function of \code{Matrix} package. It is
#' necessary having Normal covariance matrix PD since this sampling strategy is based on Cholesky decomposition
#' of it. Sampling strategy can be summed up with the following steps:
#' \enumerate{
#' \item Generate a random uniform matrix \code{Nsamples}x\code{Nvariables}
#' \item Compute the inverse CDF of standard Normal distribution for these samples
#' \item Multiply the matrix obtained in previous step for Cholesky decomposition of Normal covariance matrix
#' \item Sum the mean values of each variables
#' \item Perform exponential transformation
#' }
#'
#' The approximation to the nearest PD of initial covariance/correlation matrix and/or Normal covariance will alter
#' correlation/covariance structure among variables. In order to appreciate this difference, this function returns both
#' Frobenius and Infinity norms. See the section below for further details.
#'
#' @return A list with variable length containing:\tabular{ll}{
#' \code{input_covMat_adj} \tab Matrix object with the nearest PD matrix to input covariance matrix. This object
#' is returned only if input covariance matrix is not PD and \code{flag_force=TRUE}\cr
#' \tab \cr
#' \code{input_corrMat_adj} \tab Matrix object with the nearest PD matrix to input correlation matrix. This object
#' is returned only if input covariance matrix is not PD and \code{flag_force=TRUE}\cr
#' \tab \cr
#' \code{normF_input_adj} \tab Scalar double. It represents the Frobenius Norm of the difference between input covariace/correlation matrix and
#' its nearest PD, normalized for  the Frobenius Norm of input matrix. If input matrix is PD
#' or \emph{flag_force} is set to \emph{FALSE}, this element is not returned in output list.\cr
#' \tab \cr
#' \code{normInf_input_adj} \tab Scalar double. It represents the Infinity norm of the difference between input covariace/correlation matrix and
#' its nearest PD, normalized for  the Infinity norm of input matrix. If input matrix is PD
#' or \emph{flag_force} is set to \emph{FALSE}, this element is not returned in output list.\cr
#' \code{validation_res} \tab Matrix object which contains a brief recap of the tested conditions on the input
#'    correlation matrix. It is composed by N(=number of Log-Normal variables couples)  and 8 columns:\cr
#'    \tab \itemize{
#'    \item{\code{var1}:} {numerical index of the first Log-Normal variable of the couple}
#'    \item{\code{var2}:} {numerical index of the second Log-Normal variable of the couple}
#'    \item{\code{lower}:} {lower bound for the range of Log-Normal correlation between \code{var1} and \code{var2}}
#'    \item{\code{upper}:} {upper bound for the range of Log-Normal correlation between \code{var1} and \code{var2}}
#'    \item{\code{tested_corr}:} {input correlation value between \code{var1} and \code{var2}}
#'    \item{\code{is_valid_bound}:} {numerical flag. If \emph{1} \code{tested_corr} is inside the range, \emph{0}
#'    otherwise}
#'    \item{\code{tested_for_logTransf}:} {value computed for testing condition 2 (\eqn{\rho12*cv1*cv2})}
#'    \item{\code{can_logTransf}:} {numerical flag. If \emph{1} covariance of normal distribution between
#'    \code{var1} and \code{var2} can be computed.}
#'    } \cr
#' \tab \cr
#'    \code{normal_cov} \tab Matrix object. It contains covariance matrix of multivariate Normal distribution
#'    associated to Log-Normal one when Normal covariance is not approximated to its nearest PD. This object is
#'    returned only if \code{full_output=TRUE}.\cr
#'    \tab \cr
#'    \code{samples} \tab Matrix object. It contains \code{Nsamples} rows while the number of columns is equal to
#'    the number of Log-Normal variables.\cr
#'    \tab \cr
#'    \code{normal_cov_not_adjusted} \tab Matrix object. It contains covariance matrix of multivariate Normal distribution
#'    associated to Log-Normal one when it is not PD. This object is returned only if \code{full_output=TRUE} and
#'    Normal covariance matrix is approximated to its nearest PD.\cr
#'    \tab \cr
#'    \code{normal_cov_adjusted} \tab Matrix object. It contains the nearest PD to covariance matrix of
#'     multivariate Normal distribution associated to Log-Normal. This object is returned only if \code{full_output=TRUE}
#'     and Normal covariance matrix is approximated to its nearest PD.\cr
#'    \tab \cr
#'    \code{logN_corr_initial} \tab Matrix object. It contains the initial correlation structure of multivariate
#'    Log-Normal distribution. This object is returned only if Normal covariance matrix is approximated to its
#'    nearest PD.\cr
#'    \tab \cr
#'    \code{logNcorr_adj} \tab Matrix object. It contains the correlation structure of multivariate
#'    Log-Normal distribution considering the approximation of Normal covariance matrix to its nearest PD.
#'    This object is returned only if Normal covariance matrix is approximated to its nearest PD.\cr
#'    \tab \cr
#'    \code{normF_logN_corrMat} \tab Scalar double. It represents the Frobenius norm of the difference between
#'    initial correlation matrix and the one obtained after the approximation of Normal covariance matrix,
#'    normalized for the norm of the first. It is returned only if Normal covariance matrix is approximated to
#'    its nearest PD. \cr
#'    \tab \cr
#'    \code{normInf_logN_corrMat} \tab Scalar double. It represents the Infinity norm of the difference between
#'    initial correlation matrix and the one obtained after the approximation of Normal covariance matrix,
#'    normalized for the norm of the first. It is returned only if Normal covariance matrix is approximated to
#'    its nearest PD.\cr
#'    \tab \cr
#'    \code{logNcov_initial} \tab Matrix object. It contains the initial covariance matrix of multivariate
#'    Log-Normal distribution. This object is returned only if Normal covariance matrix is approximated to its
#'    nearest PD.\cr
#'    \tab \cr
#'    \code{logNcov_adj} \tab Matrix object. It contains the covariance matrix of multivariate
#'    Log-Normal distribution considering the approximation of Normal covariance matrix to its nearest PD.
#'    This object is returned only if Normal covariance matrix is approximated to its nearest PD.\cr
#'    \code{normF_logN_covMat} \tab Scalar double. It represents the Frobenius norm of the difference between
#'    initial covariance matrix and the one obtained after the approximation of Normal covariance matrix,
#'    normalized for the norm of the first. It is returned only if Normal covariance matrix is approximated to
#'    its nearest PD. \cr
#'    \code{normF_logN_covMat} \tab Scalar double. It represents the Infinity norm of the difference between
#'    initial covariance matrix and the one obtained after the approximation of Normal covariance matrix,
#'    normalized for the norm of the first. It is returned only if Normal covariance matrix is approximated to
#'    its nearest PD. \cr
#' }
#' @author Alessandro De Carlo \email{alessandro.decarlo01@@universitadipavia.it}
#' @seealso \code{\link{req_approx_to_PD}}
#' @seealso \code{\link{get_logNcorr_bounds}}
#' @seealso \code{\link{validate_logN_corrMatrix}}
#' @seealso \code{\link{logn_to_normal}}
#' @seealso \code{\link[Matrix]{nearPD}}
#' @seealso \code{\link[Matrix]{norm}}
#'
#' @examples
#' #different ways to run this function
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
#' cv2 <- sd2/mu2
#' Nsamples <- 1000
#' #call 1
#' mvlogn(Nsamples,mu=mu2,covMatrix=covMatrix2)
#' #call 2
#' mvlogn(Nsamples,mu=mu2, sd=sd2, corrMatrix=corr)
#' #call 3
#' mvlogn(Nsamples, mu=mu2, cv=cv2, corrMatrix=corr)
#' # call 4
#' mvlogn(Nsamples, sd = sd2, cv=cv2, corrMatrix=corr)
#' #call 5
#' mvlogn(Nsamples, sd=sd2, cv=cv2, covMatrix=covMatrix2)
#'
#' #all 5 calls allows to sample from the same distribution


mvlogn <- function(Nsamples,force_toPD=TRUE,full_output=TRUE,...){

    #get extra parameters
    plist <- list(...)
    #generate list of outputs
    out_list <- list()

    #pre-processing
    if(!is.null(plist$mu)&!is.null(plist$covMatrix)){
        #case 1: mu and covMatrix
        mu<-plist$mu
        covMatrix_in <- plist$covMatrix


        validate_covMatrix(covMatrix_in)
        validate_array("mu",mu,covMatrix_in,"covMatrix")
        test_approx <- req_approx_to_PD(covMatrix_in,"covMatrix",force_toPD,FALSE)
        covMatrix <- test_approx$matr
        if(length(test_approx)>1){
            out_list[["input_covMat_adj"]] <-test_approx$matr
            out_list[["normF_input_adj"]] <- test_approx$normF
            out_list[["normInf_input_adj"]] <- test_approx$normInf
        }

        sd <- as.array(sqrt(diag(covMatrix)))
        corrMatrix <- covMatrix/(sd%*%t(sd))


    }else if(!is.null(plist$mu)&!is.null(plist$corrMatrix)&!is.null(plist$sd)){
        #case 2: mu, sd and corrMatrix
        mu<- plist$mu
        corrMatrix_in <- plist$corrMatrix
        sd <- plist$sd

        validate_array("mu",mu,corrMatrix_in,"corrMatrix")
        validate_array("sd",sd,corrMatrix_in,"corrMatrix")
        validate_corrMatrix(corrMatrix_in)
        test_approx <- req_approx_to_PD(corrMatrix_in,"corrMatrix",force_toPD,TRUE)
        corrMatrix <- test_approx$matr
        if(length(test_approx)>1){
            out_list[["input_corrMat_adj"]] <- test_approx$matr
            out_list[["normF_input_adj"]] <- test_approx$normF
            out_list[["normInf_input_adj"]] <- test_approx$normInf
        }
        covMatrix<-sd%*%t(sd)*corrMatrix

    }else if(!is.null(plist$cv) & !is.null(plist$mu)& !is.null(plist$corrMatrix)){
        #case3: mu, cv, corrMatrix
        mu<- plist$mu
        corrMatrix_in <- plist$corrMatrix
        cv <- plist$cv

        validate_array("mu",mu,corrMatrix_in,"corrMatrix")
        validate_array("cv",cv,corrMatrix_in,"corrMatrix")
        validate_corrMatrix(corrMatrix_in)
        test_approx <- req_approx_to_PD(corrMatrix_in,"corrMatrix",force_toPD,TRUE)
        corrMatrix <- test_approx$matr
        if(length(test_approx)>1){
            out_list[["input_corrMat_adj"]] <- test_approx$matr
            out_list[["normF_input_adj"]] <- test_approx$normF
            out_list[["normInf_input_adj"]] <- test_approx$normInf
        }

        sd<- mu*cv
        covMatrix<-sd%*%t(sd)*corrMatrix

    }else if(!is.null(plist$cv)& !is.null(plist$sd) & !is.null(plist$corrMatrix)){
        #case4: sd, cv, corrMat
        sd<- plist$sd
        corrMatrix_in <- plist$corrMatrix
        cv <- plist$cv

        validate_array("sd",sd,corrMatrix_in,"corrMatrix")
        validate_array("cv",cv,corrMatrix_in,"corrMatrix")
        validate_corrMatrix(corrMatrix_in)
        test_approx <- req_approx_to_PD(corrMatrix_in,"corrMatrix",force_toPD,TRUE)
        corrMatrix <- test_approx$matr
        if(length(test_approx)>1){
            out_list[["input_corrMat_adj"]] <- test_approx$matr
            out_list[["normF_input_adj"]] <- test_approx$normF
            out_list[["normInf_input_adj"]] <- test_approx$normInf
        }

        mu<-sd/cv
        covMatrix<-sd%*%t(sd)*corrMatrix

    }else if(!is.null(plist$cv)&!is.null(plist$sd)&!is.null(plist$covMatrix)){
        #case 5: sd, cv, covMat
        sd <- plist$sd
        covMatrix_in <- plist$covMatrix
        cv <- plist$cv

        validate_array("sd",sd,covMatrix_in,"covMatrix")
        validate_array("cv",cv,covMatrix_in,"covMatrix")
        validate_covMatrix(covMatrix_in)
        test_approx <- req_approx_to_PD(covMatrix_in,"covMatrix",force_toPD,FALSE)
        covMatrix <- test_approx$matr
        if(length(test_approx)>1){
            out_list[["input_covMat_adj"]] <-test_approx$matr
            out_list[["normF_input_adj"]] <- test_approx$normF
            out_list[["normInf_input_adj"]] <- test_approx$normInf
        }
        mu<-sd/cv
        corrMatrix <- covMatrix/(sd%*%t(sd))

    }else{
        stop("Wrong list of inputs. Function can be used with the following sets of inputs
         mu, covMatrix
         mu, corrMatrix
         mu, cv, corrMatrix
         sd, cv, corrMatrix
         sd, cv, corrMatrix")
    }

    normal_params <- logn_to_normal(mu,covMatrix)

    muN <- normal_params$muN
    sigmaN <- normal_params$sigmaN

    if(is.null(muN)){
        #it means that input covariance matrix of normal distribution underlying normal one could not be computed
        new_out <- append(out_list,normal_params)
        return(new_out)
    }

    out_list[["validation_res"]] <-normal_params$validation_res

    if (min(eigen(sigmaN)$values)<0){
        spd_transf <- Matrix::nearPD(sigmaN)
        newSigma <- as.matrix(spd_transf$mat)
        if(full_output){
            out_list[["normal_cov_not_adjusted"]] <- sigmaN
            out_list[["normal_cov_adjusted"]] <- newSigma
            out_list[["normF_normal_adj"]] <- spd_transf$normF/Matrix::norm(sigmaN,"f")
            out_list[["normInf_normal_adj"]] <- Matrix::norm(sigmaN-newSigma)/Matrix::norm(sigmaN,"i")
        }
        sigmaN <-newSigma
        warning("Covariance matrix of Normal distribution associated to initial covariance matrix (lognormal distribution) is not PD. It is necessary to approximate it to the nearest PD. Input correlation structure may change.")
        #compute logN params from normal params after the approximation
        lognParam<- normal_to_logn(muN,sigmaN)
        newsd <- sqrt(diag(lognParam$sigmaLn))
        logNcorr <- lognParam$sigmaLn/(newsd%*%t(newsd))
        out_list[["logN_corr_initial"]] <-corrMatrix
        out_list[["logNcorr_adj"]] <- logNcorr
        out_list[["normF_logN_corrMat"]] <- Matrix::norm(logNcorr-corrMatrix,"F")/
            Matrix::norm(corrMatrix,"F")
        out_list[["normInf_logN_corrMat"]] <- Matrix::norm(logNcorr-corrMatrix,"I")/
            Matrix::norm(corrMatrix,"I")
        out_list[["logNcov_initial"]] <- covMatrix
        out_list[["logNcov_adj"]] <- lognParam$sigmaLn

        out_list[["normF_logN_covMat"]] <- Matrix::norm(lognParam$sigmaLn-covMatrix,"F")/
            Matrix::norm(covMatrix,"F")
        out_list[["normInf_logN_covMat"]] <- Matrix::norm(lognParam$sigmaLn-covMatrix,"I")/
            Matrix::norm(covMatrix,"I")
    }else{
        if(full_output){
            out_list[["normal_cov"]] <- sigmaN
        }
    }

    sigmaN_dec <- chol(sigmaN)
    U <- matrix(runif(length(mu)*Nsamples),ncol = length(mu))
    Zu <- qnorm(U)
    Zu_corr <- Zu%*%sigmaN_dec

    ts <- t(matrix(rep(muN,Nsamples),ncol  = Nsamples))

    samples <- exp(Zu_corr)*exp(ts)


    out_list[["samples"]] <- samples
    return(out_list)
}
