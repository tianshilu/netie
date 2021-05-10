\name{netie-package}
\alias{netie-package}
\alias{netie}
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
 netie(input_data,sigma_square = 100000 ,
                                   alpha = 10,beta = 2,sigma_p_sqr = 0.1,sigma_a_sqr = NULL,max_iter =100000,multi_sample = T) 
Please refer to https://github.com/tianshilu/Netie for more details.
}
\author{
\Tianshi Lu{netie}

Maintainer: \Tianshi Lu{netie}
}
\references{
~~ Literature or other references for background information ~~
https://github.com/tianshilu/Netie
}
~~ Optionally other standard keywords, one per line, from file KEYWORDS in the R documentation directory ~~
\keyword{ package }
\seealso{
~~ Optional links to other man pages, e.g. ~~
~~ \code{\link[<pkg>:<pkg>-package]{<pkg>}} ~~
}
\examples{
netie(input_data,sigma_square = 100000 ,
                                   alpha = 10,beta = 2,sigma_p_sqr = 0.1,sigma_a_sqr = NULL,max_iter =100000,multi_sample = T)
}
