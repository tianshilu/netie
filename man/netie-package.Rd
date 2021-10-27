\name{netie-package}
\alias{netie-package}
\docType{package}
\title{
\packageTitle{netie}
}
\description{
\packageDescription{netie}
}
\details{

The DESCRIPTION file:
\packageDESCRIPTION{netie}
\packageIndices{netie}
~~ An overview of how to use the package, including the most important functions ~~
 netie(input_data,sigma_square = 100000,alpha = 10,beta = 2,sigma_p_sqr = 0.1,sigma_a_sqr = NULL,max_iter =1000,multi_sample = T) 
Please refer to https://github.com/tianshilu/Netie for more details.
}
\author{
\packageAuthor{netie}

Maintainer: \packageMaintainer{netie}
}
\references{
https://github.com/tianshilu/Netie
}
\keyword{ package neoantigen}

\examples{
data(input_data)
netie(input_data,sigma_square=100000,alpha=10,beta=2,
sigma_p_sqr=0.1,max_iter=1000,multi_sample=TRUE)
}
