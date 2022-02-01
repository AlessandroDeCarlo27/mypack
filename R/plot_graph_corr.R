#' @title Plot of correlation graph
#' @description Function that generate a of the graph associated to correlation matrix
#' \if{html}{
#'
#'   \out{<div style="text-align: center">}\figure{man/plotGcorr.png}{options: style="width:750px;max-width:75\%;"}\out{</div>}
#'
#' }
#' \if{latex}{
#'
#'   \out{\begin{center}}\figure{man/plotGcorr.png}\out{\end{center}}
#'
#' }
#' @importFrom qgraph qgraph
#' @export
#' @param corrMat Matrix object which represents correlation matrix
#' @param title string object which contains the title of the plot. Default values is "Graph of correlation matrix".
#' @param title_size double with the size of the title. Default value is 1.2
#' @param edge_width double containing the width of the edges. Default value is 0.2
#' @param edge_label_size double which contains the size of the labels plotted on the edges. They correnspond with
#' correlations. Default value is 1.2
#' @param edge_label_position double in [0,1] interval which indicates the relative position of the label on the graph.
#' Default value is 0.6.
#' @param custom_labels vector with custom node names. If it is \emph{NULL}, then column names of correlation matrix
#' are used. If correlation matrix has not column names, indices of column are taken. Default value is NULL.
#' @param bg_color string with hexadecimal code of background color. Default value is gray "#F8F8F8"
#' @param shade_factor double for tuning edge shades. Default value is 0.2.
#' @note No more than 3 letters can be showed in each node.
#' @author Alessandro De Carlo \email{alessandro.decarlo01@@universitadipavia.it}
#' @examples
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
#' plot_graph_corr(c_start,"Graph of Correlation Matrix")
#'
#'






plot_graph_corr <- function(corrMat,
                            title="Graph of Correlation Matrix",
                            title_size=1.2,
                            edge_width = 0.2,
                            edge_label_size = 1.2,
                            edge_label_position=0.6,
                            custom_labels = NULL,
                            bg_color ="#F8F8F8",
                            shade_factor = 0.2){

    if(!is.matrix(corrMat)){
        stop("A matrix object must be passed as input.")
    }

    if(ncol(corrMat)!=nrow(corrMat)){
        stop("A square matrix object must be passed as input.")
    }

    wei_conn <- corrMat
    wei_conn[is.na(wei_conn)]<-0
    wei_conn <- wei_conn-diag(rep(1,ncol(wei_conn)))
    qgraph::qgraph(wei_conn,
                   edge.labels=TRUE,
                   title=title,
                   title.cex = title_size,
                   edge.width = edge_width,
                   edge.label.cex=edge_label_size,
                   edge.label.position=edge_label_position,
                   labels=custom_labels,
                   bg=bg_color,
                   colFactor = shade_factor)
}
