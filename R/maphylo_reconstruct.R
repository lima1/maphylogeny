maphylo_reconstruct <-
function(dists,snow=FALSE,pm=nj, ...) {

    pmw <- function(dist, ...) pm(dist$dist, ...)
    
    if ("cluster" %in% class(snow)) {
        clusterEvalQ(snow, library(ape))
        return(parLapply(snow, dists, pmw, ...))
    }    
    
    lapply(dists, pmw, ...)
}

