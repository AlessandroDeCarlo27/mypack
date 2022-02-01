#' @title Validate correlation matrix of Log-Normal distribution
#' @description Given a set of Log-Normal variables with their means, standard deviations and correlation matrix,
#' this function evaluates if the correlation structure respects two conditions:
#' \enumerate{
#' \item For each couple of variables \eqn{X1}, \eqn{X2} is tested if \eqn{\rho12} is within a certain range
#' whose bounds are computed according to the CVs of \eqn{X1} and \eqn{X2}. This is a necessary condition (becomes
#' sufficent too when correlation matrix is a 2x2) for obtaining a PD correlation/covariance matrix for the
#' Normal distribution underlying the Log-Normal one.
#' \item For each couple of variables \eqn{X1}, \eqn{X2} is tested if \eqn{\rho12*cv1*cv2 > -1}. If exists a couple
#' which does not respect this condition, it means that covariance matrix of underlying Normal distribution cannot
#' be computed. In fact, the covariance between \eqn{X1}, \eqn{X2} of Normal distribution, \eqn{\sigma^2}12,
#' is defined as \eqn{ln((\rho12*cv1*cv2 )+1)}.
#' }
#' @export
#' @param mu Array object which contains mean values of variables with a Log-Normal distribution
#' @param sd Array object which contains sd values of variables with a Log-Normal distribution
#' @param corrMatrix Matrix object which contains correlations of variables with a Log-Normal distribution
#' @return A list containing:\tabular{ll}{
#'    \code{is_valid_logN_corrMat} \tab Boolean flag which indicates if the input correlation matrix satisfies
#'    condition 1 \cr
#'    \tab \cr
#'    \code{can_log_transf_covMat} \tab Boolean flag which indicates if the input correlation matrix satisfies
#'    condition 2 \cr
#'    \tab \cr
#'    \code{validation_res} \tab Matrix object which contains a brief recap of the tested conditions on the input
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
#'    }
#'    \cr
#' }
#' @note
#' Within this package, function used for sampling from multivariate Log-Normal distribution (\code{mvlogn})
#' can overcome situations in which covariance matrix of the Normal distribution underlying the Log-Normal one is
#' not PD (i.e. condition 1 not satisfied) by approximating it to nearest PD. This may alter the desired
#' correlation structure but simulation is still allowed. On the contrary, if condition 2 is not satisfied,
#' simulation cannot be performed because computing Normal covariance matrix is impossible since the argument of
#' logarithmic trasformation would be negative for some couples of variables.
#' @author Alessandro De Carlo \email{alessandro.decarlo01@@universitadipavia.it}
#' @seealso \code{\link[Matrix]{nearPD}}
#' @seealso \code{\link{get_logNcorr_bounds}}
#' @seealso \code{\link{logn_to_normal}}
#' @examples
#'
#' #CONDITION 1 AND 2 SATISFIED
#' #input correlation matrix
#' corr<- diag(rep(1,4))
#' corr[1,4] <- 0.9
#' corr[4,1]<-corr[1,4]
#' corr[2,4] <- -0.2
#' corr[4,2] <- corr[2,4]
#' corr[3,2] <- -0.1
#' corr[2,3] <- corr[3,2]
#' #input standard deviations
#' sd2 <- array(c(rep(1,4)))
#' #input means
#' mu2 <- array(rep(2.5,4))
#' validate_logN_corrMatrix(mu2,sd2,corr)
#'
#' #output
#'
#' #$is_valid_logN_corrMat
#' #[1] TRUE
#' #
#' #$can_log_transf_covMat
#' #[1] TRUE
#'
#' #$validation_res
#' #      var1 var2 lower      upper tested_corr   is_valid_bound tested_for_logTransf  can_logTransf
#' #[1,]    1    2 -0.862069     1         0.0              1                0.000             1
#' #[2,]    1    3 -0.862069     1         0.0              1                0.000             1
#' #[3,]    1    4 -0.862069     1         0.9              1                0.144             1
#' #[4,]    2    3 -0.862069     1        -0.1              1               -0.016             1
#' #[5,]    2    4 -0.862069     1        -0.2              1               -0.032             1
#' #[6,]    3    4 -0.862069     1         0.0              1                0.000             1


validate_logN_corrMatrix <- function(mu,sd,corrMatrix){

    #parsing
    validate_corrMatrix(corrMatrix)
    validate_array("mu",mu,corrMatrix,"corrMatrix")
    validate_array("sd",sd,corrMatrix,"corrMatrix")
    #getting bounds
    bounds <- get_logNcorr_bounds(mu,sd)
    #test bounds and if a covariance matrix of underlying Normal distribution
    #can be computed
    cv <- sd/mu
    tested <- matrix(0L, ncol = 1, nrow = nrow(bounds))
    tested_logt <- matrix(0L, ncol = 1, nrow = nrow(bounds))
    test_res <- matrix(F, ncol = 2, nrow = nrow(bounds))
    for (j in 1:nrow(bounds)) {
        tested[j,1] <- corrMatrix[bounds[j,1],bounds[j,2]]
        test_res[j,1] <- corrMatrix[bounds[j,1],bounds[j,2]]>bounds[j,"lower"] &&
            corrMatrix[bounds[j,1],bounds[j,2]]<bounds[j,"upper"]
        tested_logt[j,1] <- (corrMatrix[bounds[j,1],bounds[j,2]]*cv[bounds[j,1]]*cv[bounds[j,2]])
        test_res[j,2] <- tested_logt[j,1] > -1
    }
    final_recap <- cbind(bounds,tested,test_res[,1],tested_logt,test_res[,2])
    colnames(final_recap) <- c("var1","var2","lower","upper","tested_corr","is_valid_bound",
                               "tested_for_logTransf","can_logTransf")
    out_list <- list()
    out_list[["is_valid_logN_corrMat"]] <- all(test_res[,1])
    out_list[["can_log_transf_covMat"]] <- all(test_res[,2])
    out_list[["validation_res"]] <- final_recap
    return(out_list)
}
