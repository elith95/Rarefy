\name{design_rare}
\alias{design_rare}
\title{
Customized Rarefaction Curve
}
\description{
.
}
\usage{
design_rare(comm, dist_xy = NULL, 
  method = c("fun_div", "equation"), args, 
  fun = NULL, formula = NULL, resampling = 99, 
  spatial = FALSE, mean = F, beta = F)
}

\arguments{
   \item{comm}{a community dataframe or matrix with N plots as rows, S species as columns. Both the presence/absence (1/0) or the abundances of species in plots are allowed as entries. Plot and species names should be provided as row names and column names.}
  \item{dist_xy}{an object of class \code{dist} containing the pairwise geographic distances among the plots or an object of class, an object of class \code{vector} containing the order of the sampling units along a gradient or an object of class \code{matrix} or \code{data.frame} where each column represents a gradient along which the sampling units are ordered. The names of the labels must be the same as the rows of \code{comm}. If \code{dist_xy} is not of class dist, then a distance matrix is calculated using the Euclidean distance coefficient.}
  \item{method}{Can be "fun_div" if a user-defined function will be used for the calculation of the rarefaction curve; "equation" if a formula written by the user will be used.}
  \item{args}{a list with the arguments for user-defined function or the terms of the equation. The value NA should be given in place of the community matrix in the list. The names of the elements must correspond to the names of the arguments of the function passed or to the names of the terms used in the equation.}
  \item{fun}{a string with the name of the user-defined function for the diversity index in the rarefaction. The function must calculate the value of the chosen diversity index per plot and return a numeric vector with the values calculated.}
  \item{formula}{string containing the formula of the equation to calculate the index. The terms of the equation must correspond to the names of the args list. Also functions of other packages can be included in the formula.}
   \item{resampling}{number of times plots (rows) are randomly resampled from comm to calculate the mean accumulation curve for the non-spatially-explicit rarefaction.}
   \item{spatial}{logical. If TRUE, calculates spatially explicit rarefaction.}
   \item{mean}{logical. If TRUE, the total values for species abundance or presence among the sample size is calculated as their mean instead of the sum.}
   \item{beta}{ligical. If TRUE, the result of the function provided with fun argument will be a distance matrix.}
}

\value{
An object of class \code{data.frame} with 3 columns is returned: 

  - Rarefaction: : mean of the values of the accumulation curves for all the sampling dimensions;
  
  - IC_up: upper confidence interval;
  
  - IC_low: lower confidence interval.
  
If comparison is TRUE, the \code{data.frame} object will have six columns, with the values of the accumulation curve and confidence intervals for both spatially explicit and non-spatially explicit rarefaction.
}


\author{
Elisa Thouverai \email{elisa.th95@gmail.com}
}

\examples{

\dontrun{

library(adiv)

#Data
data(duneFVG)
data(duneFVG.xy)

#Spatial distance between sampling units
dist_sp <- dist(duneFVG.xy$tot.xy)

#List with the arguments
a <- list(NA,'Shannon')
names(a) <- c('comm','method')

#Classic Rarefaction
rare_shannon <- design_rare(duneFVG$total, 
method = "fun_div", fun = 'speciesdiv', args = a, mean = TRUE)

#Spatially explicit rarefaction
rare_sp_shannon <- design_rare(duneFVG$total, dist_sp, method = "fun_div", 
fun = 'speciesdiv', args = a, mean = TRUE, spatial = TRUE)
}


}

