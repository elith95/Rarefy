\name{prep_dist}
\alias{prep_dist}
\title{
Distance matrix preparation for Rarefy
}
\description{
The function prepares the distance matrix to use to calculate spatially explicit rarefaction curves.
}
\usage{
prep_dist(xy, method = c("euclid", "k-NN"), k = 1)
}

\arguments{
  \item{xy}{data.frame or matrix object containing the coordinates to calculate a gradient or the distance between sampling units.}
  \item{method}{the method used to compute the distances, one among "euclid" and "k-NN".}
  \item{k}{number of neighbors if "k-NN" is selected as method.}
}
\details{
.
}
\value{
Object of class dist.
}
\references{
.
}
\author{
Elisa Thouverai \email{elisa.th95@gmail.com}
}

\examples{
\dontrun{.}
}
