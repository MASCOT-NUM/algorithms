#title: PCC
#help: Partial Correlation Coefficients
#type: sensitivity
#author: Bertrand Iooss
#require: sensitivity
#options: rank='false|true',nboot='100',conf='0.95',n='100'
#options.help: rank=analysis done on the ranks ?,nboot=number of bootstrap replicates,conf=confidence level of the bootstrap confidence intervals,n=sample size

PCC <- function(options) {
    library(sensitivity)
    
    options$n <- as.integer(options$n)
    options$rank <- as.logical(options$rank)
    options$nboot <- as.integer(options$nboot)
    options$conf <- as.numeric(options$conf)

    pcc = new.env()
    lapply(names(options), function(x) assign(x, options[[x]], pcc))
    return(pcc)
}

getInitialDesign <- function(algorithm, d) {
    set.seed(1)
    
    X = matrix(runif(d*algorithm$n),ncol=d,nrow=algorithm$n)
    
    colnames(X) = paste0("x",1:d)
    return(data.frame(X))
}

getNextDesign <- function(algorithm, X, Y) {
    return()
}

displayResults <- function(algorithm, X, Y) {
    algorithm$Y = Y[, 1]
    algorithm$X = X
    
    algorithm$m = eval(expr = parse(text = "pcc(X,y=Y,nboot=nboot,rank=rank,conf=conf)"), envir = algorithm)
    
    algorithm$files = "plot.png"
    png(file = algorithm$files, bg = "transparent", height = 600, width = 600)
    plot(algorithm$m)
    dev.off()
    
    p = (print(algorithm$m))
    
    html = paste0("<HTML name='sensitivity'>", 
    paste0(collapse = "\n",gsub("\\\"", replacement = "'",sub("<table", replacement = "<table border=1", capture.output(print(xtable::xtable(p),type = "html"))))) , 
    "<br/><img src='", algorithm$files, "'/></HTML>")

    return(html)
}## here you can optionaly add some dependecy code, external functions, 
## ... needed in getInitialDesign, getNextDesign or displayResults
