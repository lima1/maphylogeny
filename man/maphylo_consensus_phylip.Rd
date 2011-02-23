\name{maphylo_consensus_phylip}
\alias{maphylo_consensus_phylip}
\title{
PHYLIP fconsense
}
\description{
Uses PHYLIP fconsense (EMBOSS) to construct the majority consensus tree. 
}
\usage{
maphylo_consensus_phylip(trees, outfile = "phy", outgroup = 1)
}
\arguments{
  \item{trees}{
    A list of trees.
}
  \item{outfile}{
    The outfile prefix.
}
  \item{outgroup}{
    The outgroup.
}
}
\references{
Riester M, Stephan-Otto Attolini C, Downey RJ, Singer S, Michor F (2010) A
Differentiation-Based Phylogeny of Cancer Subtypes. PLoS Comput Biol 6(5):
e1000777. doi:10.1371/journal.pcbi.1000777
}
\author{
M. Riester
}
