#' @title Estimation of Indirect Correlations
#' @description This function estimates indirect correlations given an input correlation matrix. Indirect effect must
#' be specified by placing \code{NA} values inside input correlation matrix.
#' @export
#' @importFrom pracma fmincon
#' @importFrom Matrix nearPD
#' @param corrMatStart Matrix object. It contains correlation the matrix with fixed correlations values and indirect
#' correlations to be estimated which are represented by \code{NA} values inside the matrix.
#' @param force_estimate Boolean flag. If the optimization fails (i.e. final correlation matrix is not PD), when this
#' flag is \emph{TRUE} the matrix obtained is approximated to its nearest PD. This may alter fixed initial correlations.
#' If this flag is set to \emph{FALSE}, this further optimization is skipped and a warning message is showed.
#' @details Indirect correlations matrix are estimated solving a constrained optimization problem. Starting
#' from the fixed correlations, a correlation graph is built. Then, for each couple of variables whose indirect
#' correlation is of interest (i.e. \code{NA} values), all the possible paths among them are considered (without
#' visiting a node more than once). For each path its cost by multiplying the correlations along it.
#' The maximum and the minimum costs provide a reasonable range for the indirect correlation effect.
#' If does not exist any path betwen two nodes, that indirect correlation will not be estimated and it will be
#' automatically set to 0. If there is a unique path, the bound is computed considering \eqn{cost +- 0.5*cost}.
#'
#' Given the bounds of indirect correlations, a constrained optimization problem is solved by minimizing the negative
#' of minimum eigenvalue of correlation matrix. The starting value for the estimate of indirect correlations is provided
#' by the middle-point of the bounds computed. If this first optimization fails (i.e. final correlation matrix is not PD), user can force a second optimization
#' step in which the previously obtained correlation matrix is approximated to its nearest PD. Note that this step may
#' alter initial fixed correlations.
#'
#' Note that for estimating a indirect correlation effect between two variables, at least one path that joins them in
#' the correlation graph must exists. If for all indirect correlations declared does not exist any path, this function
#' prints a warning message and plots the graph of correlations providing a more direct contribute for debugging.
#'
#' @return A list containing \tabular{ll}{
#'   \code{corrMatfinal} \tab Matrix object that contains final correlation matrix with indirect correlations estimated \cr
#'    \tab \cr
#'    \code{optim}\tab List of objects containing the outputs provided by the solver (\code{fmincon} of \code{pracma}
#'    package) used for the constrained optimization. It is returned with this name in output list
#'    If the optimization step can converge to a PD correlation or if the optimization step fails and \code{force_estimate}
#'    is set to \code{FALSE}.
#'    \cr
#'    \tab \cr
#'    \code{optim1} \tab The same of \code{optim}, it is returned only when constrained optimization cannot converge to
#'    a PD correlation matrix and \code{force_estimate} is set to \code{TRUE}. \cr
#'    \tab \cr
#'    \code{optim2} \tab List of objects containing the outputs provided by the function \code{nearPD} of \code{Matrix}
#'    package used for approximating to nearest PD the correlation matrix obtained by solving the constrained optimization
#'    problem. It is returned only when constrained optimization cannot converge to a PD correlation matrix and \code{force_estimate} is set to \code{TRUE}.
#'    }
#'
#' @author Alessandro De Carlo \email{alessandro.decarlo01@@universitadipavia.it}
#' @seealso \code{\link[Matrix]{nearPD}}
#' @seealso \code{\link[pracma]{fmincon}}
#' @seealso \code{\link{estimate_corr_bounds}}
#' @seealso \code{\link{optim.fun}}
#'
#' @examples
#' #define initial correlation structure
#' c_start <- diag(rep(1,10))
#' c_start[1,2] <- -0.6
#' c_start[1,3] <- -0.75
#' c_start[2,3] <-0.95
#' c_start[2,4] <- 0.75
#' c_start[2,6] <- -0.6
#' c_start[2,8] <- 0.75
#' c_start[3,4] <- 0.6
#' c_start[3,8] <-0.75
#' c_start[4,7] <- 0.6
#' c_start[4,8]<-0.75
#' c_start[5,7] <- -0.95
#' c_start <- c_start+t(c_start)-diag(rep(1,ncol(c_start)))
#' #set to NA idnirect correlations to be estimated
#' c_start[c_start==0]<-NA
#'
#' colnames(c_start)<- paste(rep("X",10),1:10,sep = "")
#' rownames(c_start) <- paste(rep("X",10),1:10,sep = "")
#' # let's plot initial correlation matrix
#' plot_graph_corr(c_start,"Graph of Initial Correlation Matrix")
#' r<-estimate_indirect_corr(c_start)
#' #see final output
#' plot_graph_corr(r$corrMatfinal,'Graph of Final Correlation Matrix')
#'
#'







estimate_indirect_corr <- function(corrMatStart,force_estimate=FALSE){

    #VALIDATE INPUT
    validate_corrMatrix(corrMatStart)

    #indirect correlations to estimate must be declared with NA inside the input matrix
    if(!any(is.na(corrMatStart))){
        stop("No indirect correlations to estimate are declared")
    }

    #bounds of correlation matrix
    bounds <- estimate_corr_bounds(corrMatStart)

    #test if all bounds are NA: this implies that all indirect correlations estimates will be set to 0
    if(all(is.na(bounds[,3]))&&all(is.na(bounds[,4]))){
        list_variables <- paste(paste(bounds[,1],bounds[,2],sep="--"),collapse = " ")
        warning(paste("Cannot estimate indirect correlations, paths between all indicated variables are missing:\n",list_variables,sep=""))
        plot_graph_corr(corrMatStart,'Independet Variables')
        matOpt <- get_corrMatrixOptim(NULL,corrMatStart,bounds)
        output_optim <- list()
        output_optim[["corrMatfinal"]]<-matOpt
        return(output_optim)
    }

    #set to 0 correlations to be estimated
    corrMatStart[is.na(corrMatStart)]<-0
    # get indices of couples for which an indirect effect exists according to
    #graph path analysis
    notNa_idx <- !(is.na(bounds[,3])&is.na(bounds[,4]))
    x0 <- bounds[notNa_idx,3]+(sign(bounds[notNa_idx,3])*
                                   ((abs(bounds[notNa_idx,4])-abs(bounds[notNa_idx,3]))/2))
    #solve constrained optimization problem
    r <- pracma::fmincon(x0,optim.fun,lb=bounds[notNa_idx,3],
                         ub=bounds[notNa_idx,4],cbase=corrMatStart,
                         var_optim=bounds)

    #compute correlation matrix with estimated indirect correlations
    corrMatfinal <- get_corrMatrixOptim(r$par,corrMatStart,bounds)
    output_optim <- list()
    output_optim[["optimizationBounds"]] <-bounds
    if(r$val<0){
        #optimization successfull
        output_optim[["optim"]]<-r
        output_optim[["corrMatfinal"]] <- corrMatfinal
        return(output_optim)
    }else{
        if(force_estimate){
            warning("Optimization step-1 failed. Matrix obtained in not semidefinite-positive.
                Matrix of step-1 will be approximated to nearest correlation matrix.
                This may change fixed correlations.")
            output_optim[["optim1"]] <-r
            r2 <- Matrix::nearPD(corrMatfinal,corr = TRUE)
            output_optim[["optim2"]] <-r2
            output_optim[["corrMatfinal"]] <- as.matrix(r2$mat)
        }else{
            warning("Optimization step-1 failed. Matrix obtained in not semidefinite-positive.
                To approximate Matrix of step-1 to nearest correlation matrix set force_estimate=TRUE")
            output_optim[["optim"]]<-r
            output_optim[["corrMatfinal"]] <- corrMatfinal
        }
    }

    return(output_optim)

}
