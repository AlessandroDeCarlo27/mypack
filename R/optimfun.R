#' @title Optimization function
#' @description This function provide a general implementation of the functional that has to be minimized in order to
#' estimate indierect correlations. Given the starting correlation matrix with the indirect correlations to be estimated
#' identified by NA, it automatically generate a matrix which depends by a vector, x, of unknown variables which
#' represents indirect correlations. During the optimizazion task, in order to obtain a PD correlation matrix, the task
#' is to minmize the negative of the minimum eigenvalue in this matrix.
#' @param x vector of unknown variables (indirect correlations)
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
#' @return The negative of the minimum eigenvalue of correlation matrix. This represents the quantity to minimize in
#' the optimization problem for estimating indirect correlations.
#' @author Alessandro De Carlo \email{alessandro.decarlo01@@universitadipavia.it}



optim.fun <- function(x,cbase,var_optim){
    cbase[lower.tri(cbase,diag = F)] <- 0
    corr_eval <- cbase
    counter <- 1
    for (i in 1:nrow(var_optim)) {
        if(is.na(var_optim[i,3])||is.na(var_optim[i,4])){
            str2eval <- paste("corr_eval[",var_optim[i,1],",",var_optim[i,2],"]<- 0",sep="")
            eval(parse(text=str2eval))
        }else{
            str2eval <- paste("corr_eval[",var_optim[i,1],",",var_optim[i,2],"]<- x[",counter,"]",sep="")
            eval(parse(text=str2eval))
            counter <- counter+1
        }

    }
    corr_eval <- corr_eval+t(corr_eval)-diag(rep(1,ncol(corr_eval)))
    return(-min(eigen((corr_eval))$values))
}
