\name{rare_phylo}
\alias{rare_phylo}
\alias{ser_phylo}

\title{Phylogenetic Diversity Rarefaction Curves}

\description{
\code{rare_phylo} calculates classic rarefaction curves using different indexes of phylogenetic diversity.

\code{ser_phylo} calculates spatially-explicit rarefaction curves using different indexes of phylogenetic diversity.
}
\usage{
rare_phylo(comm, tree, 
method = c("faith","barker","Ia","hill","tsallis","renyi"), 
exp = 0, resampling = 99)

ser_phylo(comm, tree, dist_xy, 
method = c("faith","barker","Ia","hill","tsallis","renyi"), 
exp=0, comparison = FALSE, resampling = 99)
}


\arguments{
  \item{comm}{a community dataframe or matrix with N plots as rows, S species as columns. Both the presence/absence (1/0) or the abundances of species in plots are allowed as entries. Plot and species names should be provided as row names and column names.}
  \item{tree}{an object of class \code{phylo} or \code{phylo4} containing the phylogenetic tree of the species which are present in \code{comm}.}
  \item{dist_xy}{an object of class \code{dist} containing the pairwise geographic distances among the plots or an object of class, an object of class \code{vector} containing the order of the sampling units along a gradient or an object of class \code{matrix} or \code{data.frame} where each column represents a gradient along which the sampling units are ordered. The names of the labels must be the same as the rows of \code{comm}. If \code{dist_xy} is not of class dist, then a distance matrix is calculated using the Euclidean distance coefficient.}
  \item{method}{the diversity index for the calculation of the rarefaction curve, one among "faith", "barker", "Ia", "hill", "tsallis" or "renyi". See details.}
  \item{exp}{parameter that determines the sensitivity of the measure to the relative abundance of the species for "Ia", "hill", "tsallis" and "renyi" indexes.}
  \item{resampling}{number of times plots (rows) are randomly resampled from comm to calculate the mean accumulation curve for the non-spatially-explicit rarefaction.}
  \item{comparison}{if TRUE, both non-spatially explicit and spatially explicit phylogenetic rarefactions are calculated.}
}

\details{ 
The available methods are:

\code{faith}: Faith's phylogenetic diversity (PD) is defined as the sum of branch lengths in a phylogenetic tree for the assemblage of species (Faith 1992): \deqn{PD= \sum_{i \in B}L_{i}}{PD=sum_i L_i} where \eqn{L_i}{L_i} is the branch length of the branch \eqn{i} and \eqn{B} is the number of branches in the tree.

\code{barker}: Barker's index is the abundance weighted Faith's PD. The number of branches (B) is multiplied by the weighted mean branch length, with weights equal to the average abundance of species sharing that branch (Vellend et al. 2010): \deqn{PDw= B \times \frac{\sum_{i}^{B}L_{i} A_{i}}{\sum_{i}^{B}A_{i}} }{PDw=B*((sum_i L_i A_i)/(sum_i Ai))} where \eqn{L_i}{L_i} is the branch length of the branch \eqn{i}, and \eqn{A_i}{A_i} is the average abundance of the species sharing the branch \eqn{i}. \eqn{B} is the number of branches in the tree.

\code{Ia}: Ia index, by Pavoine et al. (2009), calculates PD partitioned between evolutionary periods and between plots defined in terms of spatial and time units. Tsallis or HCDT entropy (Harvda and Charvat 1967; Daroczy 1970; Tsallis 1988)  (that measures diversity by regrouping individuals into categories) is computed for each period of the phylogenetic tree, from the number of lineages that descend from the period and from the relative abundances summed within these lineages within the focal community. With \eqn{exp=0}, HCDT is the richness (number of species) minus one and Ia is Faith's PD minus the height of the phylogenetic tree; with \eqn{exp} tending to 1 HCDT is a generalization of the Shannon index while with \eqn{exp=2} HCDT is the Simpson index and Ia is Rao's QE applied to phylogenetic distances between species. To apply Ia, the phylogeny must be ultrametric: \deqn{Ha= \frac{(1-\sum_{i=1}^{n}p_{i}^{a})}{(a-1)}}{Ha= (1-sum_i p_i^a)/(a-1)} the equation for the HCDT entropy, where \eqn{p_i}{p_i} is the relative abundance of the species \eqn{i} and \eqn{a} is the scaling constant that weights the importance of rarity of the species. \deqn{I_{a}= \sum_{K=1}^{N}(t_{K}-t_{K-1})H_{a,K}}{Ia= sum_K ((t_K-t_(K-1))H_a,K)} where \eqn{H_{a,K}}{H_a,K} is Ha applied to the period \eqn{K} and \eqn{t_K-t_{K-1}}{t_K-t_(K-1)} is the length of the period \eqn{K}.

Hill index (\code{hill}) and the HCDT (\code{tsallis}) and Renyi (\code{renyi}) entropies are adapted for the calculation of the phylogenetic diversity replacing the species with the units of the branch length in the phylogenetic tree (Pavoine & Ricotta 2019): \deqn{Hill=\left [ \sum_{i \in B} L_{i}(p_{i})^{q}]\right ]^{\frac{1}{1-q}}}{Hill= [sum_i Li(pi)^q]^(1/(1-q))}\deqn{HCDT= \frac{1-\sum_{i \in B} L_{i}(p_{i})^q}{q-1}}{HCDT= (1-sum_i L_i(p_i)^q)/(q-1)}\deqn{Renyi= \frac{1}{1-q}log\left [ \sum_{i \in B} L_{i}(p_{i})^q\right ]}{Renyi= (1/(1-q))log[sum_i L_i(p_i)^q]} where \eqn{L_i}{L_i} is the branch length of the branch \eqn{i}, \eqn{p_i}{p_i} is the relative abundance of the species sharing the branch \eqn{i} and \eqn{q} is the scaling constant that weights the importance of rarity of the species. \eqn{B} is the number of branches in the tree.

}

\value{
An object of class \code{data.frame} with 3 columns is returned: 

  - Rarefaction: : mean of the values of the accumulation curves for all the sampling dimensions;
  
  - IC_up: upper confidence interval;
  
  - IC_low: lower confidence interval.
  
If comparison is TRUE, the \code{data.frame} object will have six columns, with the values of the accumulation curve and confidence intervals for both spatially explicit and non-spatially explicit rarefaction.
}

\references{
Chao, A., Chiu, C.-H., Hsieh, T.C., Davis, T., Nipperess, D.A., Faith, D.P. (2014) Rarefaction and extrapolation of phylogenetic diversity.  \emph{Methods in Ecology and Evolution}, \bold{6}, 380--388.

Daroczy, Z. (1970) Generalized information functions. \emph{Information and Control}, \bold{16}, 36--51.

Faith, D.P. (1992) Conservation evaluation and phylogenetic diversity. \emph{Biological Conservation}, \bold{61}, 1--10.

Havrda, M., Charvat F. (1967) Quantification method of classification processes: concept of structural alpha-entropy. \emph{Kybernetik}, \bold{3}, 30--35

Pavoine, S., Love, M., Bonsall, M.B. (2009) Hierarchical partitioning of evolutionary and ecological patterns in the organization of phylogenetically-structured species assemblages: application to rockfish (genus: Sebastes) in the Southern California Bight. \emph{Ecology Letters}, \bold{12}, 898--908.

Pavoine, S., Ricotta, C. (2019) A simple translation from indices of species diversity to indices of phylogenetic diversity. \emph{Ecological Indicators}, \bold{101}, 552--561.

Swenson, N.G. (2014) Functional and Phylogenetic Ecology in R. Springer UseR! Series, Springer, New York, New York, U.S.A.

Tsallis, C. (1988) Possible generalization of Boltzmann-Gibbs statistics. \emph{Journal of Statistical Physics}, \bold{52}, 480--487.

Vellend, M., Cornwell, W.K., Magnuson-Ford, K., Mooers, A. (2010) Measuring phylogenetic biodiversity. \emph{Magurran & McGill}, 194--207.
}
\author{
Elisa Thouverai \email{elisa.th95@gmail.com}

with contributions of Sandrine Pavoine.
}

\examples{
\dontrun{
#Time consuming
require(picante)
require(geiger)
data(phylocom)
phylo<-treedata(phylocom$phylo,phylocom$sample[1,],warnings = FALSE)$phy
## Non-spatially explicit rarefaction 
raref<-rare_phylo(phylocom$sample,phylo,resampling=999) ##Faith index
plot(raref [,1], ylab="Faith", xlab="Number of plots", type="l", ylim=range(raref, na.rm
=TRUE))
lines(raref[,2], lty=2)
lines(raref[,3], lty=2)
rareb<-rare_phylo(phylocom$sample,phylo,method='barker',resampling=999) ##Barker index
plot(rareb [,1], ylab="Barker", xlab="Number of plots", type="l", ylim=range(rareb, na.rm
=TRUE))
lines(rareb[,2], lty=2)
lines(rareb[,3], lty=2)
rareia<-rare_phylo(phylocom$sample,phylo,method='Ia',resampling=999,exp=2) ##Ia index
plot(rareia [,1], ylab="Ia", xlab="Number of plots", type="l", ylim=range(rareia, na.rm
=TRUE))
lines(rareia[,2], lty=2)
lines(rareia[,3], lty=2)
}}

\keyword{models}
