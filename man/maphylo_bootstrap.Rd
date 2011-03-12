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
= seq(0.5, 1.4, by = 0.1), dm = "pearson", snow = FALSE)
}
\arguments{
  \item{M}{
    An ExpressionSet or matrix.
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
  \item{snow}{
    If specified, use the \code{snow} package for bootstrapping on a
    cluster. 
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

\dontrun{
# Typically it is recommended to filter uninformative genes, for
# example with a low variance filter.
library(genefilter)
Golub_Merge = nsFilter(Golub_Merge)$eset
}

# bootstrap
dists = maphylo_bootstrap(Golub_Merge,r=c(1))

\dontrun{
library(snow)
# bootstrap on 10 cluster nodes
cl = makeMPIcluster(10)
# use a thread-safe random number generator
clusterSetupRNG(cl)

dists = maphylo_bootstrap(Golub_Merge,r=c(1),bootstrap=1000,snow=cl)
stopCluster(cl)
}

# use Neighbor-Joining to reconstruct trees
trees = maphylo_reconstruct(dists)

\dontrun{
# it is also possible to reconstruct trees in parallel, especially useful 
# for slower methods and large trees
trees = maphylo_reconstruct(dists,snow=cl,pm=fastme.bal)
}

# plot the first tree
plot(trees[[1]])

# Write the trees in NEXUS format. Now analyze the tree (in Dendroscope for
# example). 
write.nexus(trees, file="maphylo.nxs")

}
