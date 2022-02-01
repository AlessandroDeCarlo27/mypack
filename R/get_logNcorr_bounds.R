#' @title Compute ranges of Correlations for Log-Normal distributions
#' @description Correlations between Log-Normal distributed variables should be in specific ranges
#' depending on their CVs.These bounds represent a necessary condition for obtaining a PD covariance/correlation
#' matrix for the normal distribution underlying.
#' @importFrom utils combn
#' @export
#' @param mu Array object which contains mean values of variables with a Log-Normal distribution
#' @param sd Array object which contains sd values of variables with a Log-Normal distribution
#' @return A Matrix object with \emph{nrow} equal to the number of Log-Normal variables couples and 4 columns:
#' \itemize{
#' \item{\code{var1}:} {numerical index of the first Log-Normal variable of the couple}
#' \item{\code{var2}:} {numerical index of the second Log-Normal variable of the couple}
#' \item{\code{lower}:} {lower bound for the range of Log-Normal correlation between \code{var1} and \code{var2}}
#' \item{\code{upper}:} {upper bound for the range of Log-Normal correlation between \code{var1} and \code{var2}}
#' }
#' @author Alessandro De Carlo \email{alessandro.decarlo01@@universitadipavia.it}
#' @examples
#' sd2 <- array(c(1,2,0.5,3)) #array with standard deviations of variables
#' mu2 <- array(c(3,0.2,1,8)) #array with means of variables
#' get_logNcorr_bounds(mu2,sd2)
#'
#' #returns
#' #      var1 var2   lower     upper
#' # [1,]    1    2 -0.1506242 0.3025073
#' # [2,]    1    3 -0.8529277 0.9942674
#' # [3,]    1    4 -0.8885903 0.9996221
#' # [4,]    2    3 -0.1275056 0.3517665
#' # [5,]    2    4 -0.1443341 0.3146269
#' # [6,]    3    4 -0.8398526 0.9968250
#'
#'


get_logNcorr_bounds <- function(mu,sd){

    # parsing
    if(!is.array(mu)){
        stop("mu must be an array object")
    }
    if(!is.array(sd)){
        stop("sd must be an array object")
    }
    if(length(mu)!=length(sd)){
        stop("mu and sd must have the same length")
    }

    cv <- sd/mu
    combos <- t(combn(1:length(mu),2))
    bounds <- matrix(0L,nrow = nrow(combos),ncol = 2)
    for (i in 1:nrow(combos)) {
        cv1 <- cv[combos[i,1]]
        cv2 <- cv[combos[i,2]]
        L <- sqrt(log(cv1^2+1)*log(cv2^2+1))
        bounds[i,1] <- (1/cv1)*(1/cv2)*(exp(-L)-1) #lower
        bounds[i,2] <- (1/cv1)*(1/cv2)*(exp(L)-1) #upper bound
    }
    bounds_recap <- cbind(combos,bounds)
    colnames(bounds_recap) <- c("var1","var2","lower","upper")
    return(bounds_recap)

}
