#' @title Estimate bounds of indirect correlations
#' @description For each indirect correlation to be estimated, this function returns the range of its possible values.
#' Indirect correlations must be identified in input correlation matrix with \code{NA} values. From input
#' correlation structure is derived the associated graph. Then, for each indirect correlation between the generic
#' couple of variables (\eqn{X1},\eqn{X2}),  all the possible paths in the graph that links \eqn{X1} and \eqn{X2}
#' visiting one node at most once are considered.
#'
#' For each path is computed its cost by multiplying the correlations along it, range of indirect correlation values is
#' taken considering the minumum and the maximum costs. If does not exist any path betwen two nodes, maximum and minimum values will be set to \code{NA}. If
#' there is a unique path, the bound is computed considering \eqn{cost +- 0.5*cost}.
#' @export
#' @importFrom igraph graph_from_adjacency_matrix all_simple_paths
#' @param corrMat Matrix Object which represents a correlation matrix. Indirect correlations must be set to
#' \code{NA} within the matrix. This matrix must contain at least two \code{NA} values (it is symmetric).
#' @return A matrix object with N(= number of indirect correlations) rows and 4 columns:
#' \itemize{
#'    \item{\code{var1}:} {numerical index of the first variable of the indirect correlation couple}
#'    \item{\code{var2}:} {numerical index of the second variable of the indirect correlation couple}
#'    \item{\code{lower}:} {lower bound for the range of indirect correlation between \code{var1} and \code{var2}}
#'    \item{\code{upper}:} {upper bound for the range of indirect correlation between \code{var1} and \code{var2}}
#'    }
#' @author Alessandro De Carlo \email{alessandro.decarlo01@@universitadipavia.it}
#' @seealso \code{\link[igraph]{graph_from_adjacency_matrix}}
#' @seealso \code{\link[igraph]{all_simple_paths}}
#' @seealso \code{\link{get_pathProduct}}
#' @seealso \code{\link{get_bounds}}
#' @examples
#' # create a correlation matrix, define some correlations
#' c_start <- diag(rep(1,6))
#' c_start[1,2] <- -0.6
#' c_start[1,3] <- -0.75
#' c_start[2,3] <-0.95
#' c_start[2,4] <- 0.75
#' c_start[2,6] <- -0.6
#' c_start <- c_start+t(c_start)-diag(rep(1,ncol(c_start)))
#' #set to NA indirect correlations
#' c_start[c_start==0]<-NA
#' #plot correlation graph
#' plot_graph_corr(c_start,"Graph of Correlation Matrix")
#' #get bounds of correlations
#' estimate_corr_bounds(c_start)
#'
#' #output: variable 5 is not directly correlated with any of the others so it
#' #is impossible to establish a path to it. So, indirect correlations bounds with variable 5 are set to NA
#'
#' #       var1 var2     lower   upper
#' #[1,]    1    4 -0.534375 -0.4500
#' #[2,]    1    5        NA      NA
#' #[3,]    1    6  0.360000  0.4275
#' #[4,]    2    5        NA      NA
#' #[5,]    3    4  0.337500  0.7125
#' #[6,]    3    5        NA      NA
#' #[7,]    3    6 -0.570000 -0.2700
#' #[8,]    4    5        NA      NA
#' #[9,]    4    6 -0.540000 -0.3600
#' #[10,]    5    6        NA      NA




estimate_corr_bounds <- function(corrMat){

    if(!is.matrix(corrMat)){
        stop("A matrix object must be passed as input.")
    }

    if(ncol(corrMat)!=nrow(corrMat)){
        stop("A square matrix object must be passed as input.")
    }

    if(!any(is.na(corrMat))){
        stop("No indirect correlations to estimate are declared")
    }

    #build the graph associated to the corresponding correlation matrix.
    #correlations to be estimated (indirect) are fixed to NA in input correlation Matrix
    #then they are set to 0.
    #nodes: variables
    #edges: correlations !=0 in input correlation matrix
    conn <- corrMat
    conn[is.na(conn)] <-0
    conn<-conn-diag(rep(1,ncol(corrMat))) #don't consider loops
    conn[conn!=0]<-1
    adjm <- conn
    #adjm is the adjacency matrix in which:
    #1: a correlation between a couple of variable is known
    #0: a correlation between a couple of variable is unknown
    g1 <- igraph::graph_from_adjacency_matrix( adjm, mode="undirected")
    varIdx <- 1:ncol(corrMat)
    #get number of NA elements (is the number of correlations to estimate)
    logical_na <- is.na(corrMat[upper.tri(corrMat)])
    Ndeterm <- length(logical_na[logical_na==T])
    bounds <- matrix(0L,ncol = 4,nrow = Ndeterm)
    colnames(bounds) <- c("var1","var2","lower","upper")

    #get the couple of parameters for which indirect correlations must be estimed
    logical_na_matrix <- corrMat
    logical_na_matrix[lower.tri(logical_na_matrix,diag = T)] <- 0
    logical_na_matrix <- is.na(logical_na_matrix)
    counter <- 1
    for (i in 1:nrow(logical_na_matrix)) {
        for (j in 1:ncol(logical_na_matrix)) {
            if(logical_na_matrix[i,j]&j>i){
                bounds[counter,1]<-i
                bounds[counter,2]<-j
                counter <- counter+1
            }
        }
    }

    #computing bounds of indirect correlations to estimate
    for (i in 1:nrow(bounds)) {
        paths <- igraph::all_simple_paths(g1,bounds[i,1],bounds[i,2],cutoff=-1)
        bounds[i,3:4] <- get_bounds(corrmat = corrMat,paths)

    }
    return(bounds)

}
