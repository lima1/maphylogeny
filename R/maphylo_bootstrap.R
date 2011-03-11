maphylo_bootstrap <-
function(M, group=c(), bootstrap = 100,
bootstrap_groups=FALSE,r=seq(.5,1.4,by=.1), dm="pearson", snow=FALSE)
{
    if (length(group) == 0) group = as.factor(1:ncol(M))
    if (!is.factor(group)) stop("group not a factor")

    Mss = exprs(M)
    
    combine_matrix <- function(Mx) {
        t = sapply(levels(group), function(x) sampleNames(Mx)[group %in%
        x])
        My = sapply(t, function(x) { 
            if (length(x) > 1) { 
                y=x 
                if (bootstrap_groups) y = sample(x,replace=TRUE)
                return(apply(exprs(Mx[,y]),1,mean)) 
             }
             else return(exprs(Mx[,x]))} )
        rownames(My) = featureNames(Mx)     
        My
    }        
    if (length(group) != length(levels(group))) Mss = combine_matrix(M)
    
    do_test <- function(ri) {
        mysample = Mss

        if (bootstrap_groups) 
            mysample = combine_matrix(M)

        if (bootstrap > 0)
            mysample = mysample[sample(nrow(mysample), ri, replace=TRUE),]

        d = switch(dm,
            euclidean = dist(t(mysample), method="euclidean"),
            spearman  = spearman.dist(t(mysample)),
            tau       = tau.dist(t(mysample)),
            pearson   = as.dist((1-cor(mysample))/2),
            stop("Unknown distance metric"))

        list(dist=d,genes=row.names(mysample));
    }

    reps = as.vector(sapply(r,function(x) rep(nrow(M)*x,max(bootstrap,1),)))
    if ("cluster" %in% class(snow)) 
        dists = parLapply(snow, reps, do_test)
    else    
        dists = lapply(reps, do_test )

    dists
}

