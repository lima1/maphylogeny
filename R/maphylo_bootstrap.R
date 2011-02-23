maphylo_bootstrap <-
function(M, group=c(), bootstrap = 100,
bootstrap_groups=FALSE,r=seq(.5,1.4,by=.1), dm="spearman")
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
        if (bootstrap_groups) 
            mysample = combine_matrix(M[sample(nrow(Mss), ri, replace=TRUE),])
        else 
            mysample = Mss[sample(nrow(Mss), ri, replace=TRUE),]

        if (dm == "euclidean") {
            d = dist(t(mysample), method="euclidean")
        } else if (dm == "spearman") {
            d = spearman.dist(t(mysample))
        } else if (dm == "tau") {
            d = tau.dist(t(mysample))
        }
        else        
        d = as.dist((1-cor(mysample))/2)

        list(dist=d,genes=row.names(mysample));
    }

    reps = as.vector(sapply(r,function(x) rep(nrow(M)*x,max(bootstrap,1),)))
    dists = lapply(reps, function(i) do_test(i) )
}

