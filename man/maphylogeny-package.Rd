\name{maphylogeny-package}
\alias{maphylogeny-package}
\alias{maphylogeny}
\docType{package}
\title{
maphylogeny - Gene Expression Phylogenies
}
\description{
Fancy hierarchical clustering of microarray data based on phylogenetic
methods. Supports the \code{snow} package for bootstrapping on a cluster,
}
\details{
\tabular{ll}{
Package: \tab maphylogeny\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2011-02-23\cr
License: \tab GPL\cr
LazyLoad: \tab yes\cr
}
}
\author{
M. Riester

Maintainer: Who to complain to <markus@jimmy.harvard.edu>
}
\references{
Riester M, Stephan-Otto Attolini C, Downey RJ, Singer S, Michor F (2010) A
Differentiation-Based Phylogeny of Cancer Subtypes. PLoS Comput Biol 6(5):
e1000777. doi:10.1371/journal.pcbi.1000777
}
\keyword{ package }
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
