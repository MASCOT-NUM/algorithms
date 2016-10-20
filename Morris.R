#title: Morris
#help: Morris's Elementary Effects Screening Method
#type: sensitivity
#author: Gilles Pujol
#require: sensitivity
#options: r='10',levels='4'
#options.help: r=integer giving the number of repetitions of the design,levels=integer specifying the number of levels of the design

Morris <- function(options) {
    library(sensitivity)
    options$r <- as.integer(options$r)
    options$levels <- as.integer(options$levels)
    morris = new.env()
    lapply(names(options), function(x) assign(x, options[[x]], 
        morris))
    return(morris)
}

getInitialDesign <- function(algorithm, d) {
    set.seed(1)
    algorithm$m <- morris(model = NULL, factors = d, r = algorithm$r, 
        design = list(type = "oat", levels = algorithm$levels))
    return(algorithm$m$X)
}

getNextDesign <- function(algorithm, X, Y) {
    return()
}

displayResults <- function(algorithm, X, Y) {
    algorithm$Yi = Y[, 1]
    eval(expr = parse(text = "tell(m,Yi)"), envir = algorithm)
    
    algorithm$files = "plot.png"
    png(file = algorithm$files, bg = "transparent", height = 600, width = 600)
    plot(algorithm$m)
    dev.off()
    
    mu = colMeans(algorithm$m$ee)
    mu.star = colMeans(abs(algorithm$m$ee))
    sig = sqrt((colSums(algorithm$m$ee^2) - mu^2 * 4)/(algorithm$r - 1))
    effect = rep("Weak effect", length(mu))
    for (i in 1:length(mu)) if (mu.star[i] > 0.2 * max(mu.star)) 
        if (sig[i] > 0.2 * max(sig)) 
            effect[i] = "Non linear or interactive effect"
        else effect[i] = "Linear effect"
        
    html = paste0("<HTML name='sensitivity'>", paste0(collapse = "\n", 
        sub("<table", replacement = "<table border=1", capture.output(print(xtable::xtable(cbind(mu, 
            mu.star, sig, effect)), type = "html")))), "<br/><img src='", 
        algorithm$files, "'/></HTML>")
        
    return(html)
}
