\name{maphylo_bootstrap}
\alias{maphylo_bootstrap}
\title{
Bootstrapping of Expression Data
}
\description{
This function takes an ExpressionSet, generates bootstrap replicates and
calculates the pairwise distances for each bootstrap replicate. Returns a
list of distance matrices.
}
\usage{
maphylo_bootstrap(M, group = c(), bootstrap = 100, bootstrap_groups = FALSE, r
= seq(0.5, 1.4, by = 0.1), dm = "pearson")
}
\arguments{
  \item{M}{
    An ExpressionSet.
}
  \item{group}{
   A factor, can be used to group samples. Expression values are then
   averaged over the groups. 
}
  \item{bootstrap}{
  The number of bootstrap replicates.
}
  \item{bootstrap_groups}{
  If TRUE, then samples in groups are also sampled with replacement. Will
  result in much longer calculation times.
}
  \item{r}{
       numeric vector which specifies the relative sample sizes of
       bootstrap replications. For original sample size n and
       bootstrap sample size n', this is defined as
       r=n'/n.
}
  \item{dm}{
    The distance metric. Valid are "euclidean", "spearman", "tau" and
    "pearson".
}
}
\value{
A list of distances. The used probe sets are also attached for every bootstrap
replicate.
}
\references{
Riester M, Stephan-Otto Attolini C, Downey RJ, Singer S, Michor F (2010) A
Differentiation-Based Phylogeny of Cancer Subtypes. PLoS Comput Biol 6(5):
e1000777. doi:10.1371/journal.pcbi.1000777
}
\author{
M. Riester
}

\examples{
library(golubEsets)
data(Golub_Merge)

# sampleNames are the labels of the leaves
sampleNames(Golub_Merge) = paste(sampleNames(Golub_Merge),
    Golub_Merge$ALL.AML) 

# bootstrap
dists = maphylo_bootstrap(Golub_Merge,r=c(1))

# use Neighbor-Joining to reconstruct trees
trees = maphylo_reconstruct(dists)

# plot the first tree
plot(trees[[1]])

}
