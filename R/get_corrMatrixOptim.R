#' @title Get final correlation matrix with estimated indirect correlations
#' @description This function returns the final correlations matrix with the estimated values of indirect correlations.
#' @param x_opt vector with estimated indirect correlations matrix.
#' @param cbase Matrix object. It represents the input correlation matrix. Indirect correlations to be estimated must be
#' explicitated with NA in this matrix.
#' @param var_optim Matrix object. It  has a row number equal to the number of indirect correlations to be estimated and
#' and 4 columns:
#' \itemize{
#'    \item{\code{var1}:} {numerical index of the first variable of the indirect correlation couple}
#'    \item{\code{var2}:} {numerical index of the second variable of the indirect correlation couple}
#'    \item{\code{lower}:} {lower bound for the range of indirect correlation between \code{var1} and \code{var2}}
#'    \item{\code{upper}:} {upper bound for the range of indirect correlation between \code{var1} and \code{var2}}
#'    }
#' @return Matrix object which contains the original correlation matrix with indirect correlations estimated.
#' @author Alessandro De Carlo \email{alessandro.decarlo01@@universitadipavia.it}
#' @seealso \code{\link{optim.fun}}



#function that returns the correlation matrix optimized
get_corrMatrixOptim <- function(x_opt,cbase,var_optim){
    cbase[lower.tri(cbase,diag = F)] <- 0
    corr_eval <- cbase
    counter <- 1
    for (i in 1:nrow(var_optim)) {
        if(is.na(var_optim[i,3])||is.na(var_optim[i,4])){
            str2eval <- paste("corr_eval[",var_optim[i,1],",",var_optim[i,2],"]<- 0",sep="")
            eval(parse(text=str2eval))
        }else{
            str2eval <- paste("corr_eval[",var_optim[i,1],",",var_optim[i,2],"]<- ",x_opt[counter],sep="")
            eval(parse(text=str2eval))
            counter <- counter+1
        }
    }
    corr_eval <- corr_eval+t(corr_eval)-diag(rep(1,ncol(corr_eval)))
    return(corr_eval)

}
