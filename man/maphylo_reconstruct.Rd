\name{maphylo_reconstruct}
\alias{maphylo_reconstruct}
\title{
Reconstruct Phylogenies
}
\description{
Takes a list of distance matrices and reconstructs trees with the specified
method.
}
\usage{
maphylo_reconstruct(dists, snow = FALSE, pm = nj, ...)
}
\arguments{
  \item{dists}{
    The list of distance matrices, typically obtained with
    \code{\link{maphylo_bootstrap}}.
}
  \item{snow}{
    If specified, use the \code{snow} package for bootstrapping on a
    cluster. 
}    
  \item{pm}{
    The pylogenetic reconstruction method. First argument takes a distance
    matrix. See the \code{ape} package.
}
  \item{\dots}{
    Further arguments for the pm method.
}
}
\value{
A list of trees.
}
\details{
See \code{\link{maphylo_bootstrap}} for detailed examples.
}
\references{
Citation: Riester M, Stephan-Otto Attolini C, Downey RJ, Singer S, Michor F
(2010) A Differentiation-Based Phylogeny of Cancer Subtypes. PLoS Comput Biol
6(5): e1000777. doi:10.1371/journal.pcbi.1000777

Paradis E., Claude J. & Strimmer K. 2004. APE: analyses of
phylogenetics and evolution in R language. Bioinformatics 20:
289-290.
      
}
\author{
M. Riester
}
\examples{

## The function is currently defined as
function(dists,pm=nj, ...) {
    lapply(dists, function(dist) pm(dist$dist, ...))
  }

\dontrun{
    # use Neighbor Joining (NJ)
    trees = maphylo_reconstruct(dists)

    # FastMe
    trees = maphylo_reconstruct(dists, pm = fastme.bal)
    
    # MPI
    cl = makeMPIcluster(10)
    trees = maphylo_reconstruct(dists, snow=cl, pm = fastme.bal)
}
}
