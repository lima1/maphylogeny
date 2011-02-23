maphylo_consensus_phylip <-
function(trees, outfile="phy", outgroup=1) {
    outfileall = paste(outfile,".trees",sep="")
    sapply(1:length(trees), function(i) write.tree(trees[[i]],
    file=outfileall, append=TRUE))
    system(paste("fconsense -intreefile ", outfileall, " -outgrno ",
    outgroup, 
    " -outfile ", outfile, ".fconsensus -outtreefile ",
    outfile, ".ftree", sep = ""))
}

