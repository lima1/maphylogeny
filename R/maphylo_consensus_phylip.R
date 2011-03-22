maphylo_consensus_phylip <-
function(trees, outfile="phy", outgroup=1) {
    outfileall = paste(outfile,".trees",sep="")
    if (file.exists(outfileall)) file.remove(outfileall)
    sapply(trees, write.tree, file=outfileall, append=TRUE)
    system(paste("fconsense -intreefile ", outfileall, " -outgrno ",
    outgroup, 
    " -outfile ", outfile, ".fconsensus -outtreefile ",
    outfile, ".ftree", sep = ""))
}

