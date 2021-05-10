\name{input_data}
\alias{input_data}
\docType{data}
\title{
%%   ~~ data name/kind ... ~~
input_data
}
\description{
%%  ~~ A concise (1-5 lines) description of the dataset. ~~
one kidney cancer patient with 2 samples
}
\usage{data("input_data")}
\format{
  A data frame with 297 observations on the following 7 variables.
  \describe{
    \item{\code{mutation_id}}{a character vector }
    \item{\code{sample_id}}{a character vector}
    \item{\code{cluster_id}}{a numeric vector}
    \item{\code{cellular_prevalence}}{a numeric vector}
    \item{\code{cellular_prevalence_std}}{a numeric vector}
    \item{\code{variant_allele_frequency}}{a numeric vector}
    \item{\code{neo_load}}{a numeric vector}
  }
}
\details{
%%  ~~ If necessary, more details than the __description__ above ~~
}
\source{
%%  ~~ reference to a publication or URL from which the data were obtained ~~
}
\references{
%%  ~~ possibly secondary sources and usages ~~
}
\examples{
data(input_data)
## maybe str(input_data) ; plot(input_data) ...
}
\keyword{datasets}
