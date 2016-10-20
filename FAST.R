#title: FAST
#help: Extended Fourier Amplitude Sensitivity Test
#type: sensitivity
#author: Gilles Pujol
#require: sensitivity
#options: n='100',M='4',q='qunif',q.arg='list(0,1)'
#options.help: n=samples isze: length of the discretization of the s-space,M=number of harmonics to sum in the Fourier series decomposition,q=quantile functions names corresponding to wanted factors distributions,q.arg=quantile functions parameters

FAST <- function(options) {
    library(sensitivity)
    
    # all parameters are initialy strings, so you have to put as global non-string values
    options$n <- as.integer(options$n)
    options$M <- as.integer(options$M)
    options$q <- as.character(options$q)
    options$q.arg <- as.list(eval(parse(text = options$q.arg)))
    
    fast = new.env()
    lapply(names(options),function(x) assign(x,options[[x]],fast))
    return(fast)
}

getInitialDesign <- function(algorithm, d) {
    set.seed(1)
    print(d)
    algorithm$m <- fast99(model=NULL,factors=d,n=algorithm$n,M=algorithm$M,q=algorithm$q,q.arg=algorithm$q.arg)
    return(algorithm$m$X)
}

getNextDesign <- function(algorithm, X, Y) {
    return()
}

displayResults <- function(algorithm, X, Y) {
    algorithm$Yi = Y[,1]
    eval(expr=parse(text="tell(m,Yi)"),envir=algorithm) #tell(algorithm$m,Y)
    
    algorithm$files = "plot.png"
    png(file=algorithm$files,bg="transparent",height=600,width=600)
    plot(algorithm$m)
    dev.off()
    
    sobol_mat=data.frame('First order'=algorithm$m$D1,'Total order'=algorithm$m$Dt)
    rownames(sobol_mat) <- names(X)
    
    html=paste0("<HTML name='sensitivity'>",paste0(collapse="\n",sub("<table",replacement = "<table border=1",capture.output(print(xtable::xtable(  sobol_mat),type="html")))),"<br/><img src='",algorithm$files,"'/></HTML>")
    return(html)
}
