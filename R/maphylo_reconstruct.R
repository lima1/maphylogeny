maphylo_reconstruct <-
function(dists,pm=nj, ...) {
    lapply(dists, function(dist) pm(dist$dist, ...))
}

