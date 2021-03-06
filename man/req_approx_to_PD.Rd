% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/req_approx_to_PD.R
\name{req_approx_to_PD}
\alias{req_approx_to_PD}
\title{Required approximation to nearest PD for Coviarance/Correlation Matrix}
\usage{
req_approx_to_PD(matr, name, flag_force, is_corr)
}
\arguments{
\item{matr}{matrix object representing a covariance/correlation matrix}

\item{name}{string with the name of matrix object representing a covariance/correlation matrix}

\item{flag_force}{boolean flag for managing input matrices not PD. if \emph{TRUE}, input matrix is approximated
to nearest PD if necessary; if \emph{FALSE} approximation to nearest PD is not performed, a warning message is
displayed and input matrix is returned as output.}

\item{is_corr}{boolean flag used to distinguish Correlation to Covariance matrices.}
}
\value{
A list containing:\tabular{ll}{
   \code{matr} \tab Matrix object. It contains input correlation/covariance matrix, if the one passed in
                    input is PD or if \emph{flag_force} is set to \emph{FALSE}. If input matrix object is not
                    PD and  \emph{flag_force} is set to \emph{FALSE}, \code{matr} contains the nearest PD matrix
                    to input one.\cr
   \tab \cr
   \code{normF} \tab Scalar double. It represents the Frobenius Norm of the difference between input matrix and
                     its nearest PD, normalized for  the Frobenius Norm of input matrix. If input matrix is PD
                     or \emph{flag_force} is set to \emph{FALSE}, this element is not returned in output list.\cr
   \tab \cr
   \code{normInf} \tab Scalar double. It represents the Infinity norm of the difference between input matrix and
                     its nearest PD, normalized for  the Infinity norm of input matrix. If input matrix is PD
                     or \emph{flag_force} is set to \emph{FALSE}, this element is not returned in output list.\cr
}
}
\description{
Function that test if input matrix object is PD or not. If it is PD, then input matrix is returned.
Otherwise, if \emph{flag_force} is \emph{TRUE}, input matrix is approximated to nearest PD, new matrix and
normalized Frobenius and Infinity norms are returned. If \emph{flag_force} is \emph{FALSE}, then a warning message is
displayed and input matrix is returned as output.
}
\seealso{
\code{\link[Matrix]{nearPD}}

\code{\link[Matrix]{norm}}
}
\author{
Alessandro De Carlo \email{alessandro.decarlo01@universitadipavia.it}
}
