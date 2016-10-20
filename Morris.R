#title: Morris
#help: Morris's Elementary Effects Screening Method
#author: Gilles Pujol
#type: Sensitivity analysis
#output: Sensitivities
#options: r=10,levels=4
#require: sensitivity

Morris <- function(options) {
    library(sensitivity)
    
    # all parameters are initialy strings, so you have to put as global non-string values
    options$r <- as.integer(options$r)
    options$levels <- as.integer(options$levels)
    
    morris = new.env()
    lapply(names(options),function(x) assign(x,options[[x]],morris))
    return(morris)
}

getInitialDesign <- function(morris,d) {
    set.seed(1)
    morris$m <- morris(model=NULL,factors=d,r=morris$r,design=list(type="oat",levels=morris$levels))
    return(morris$m$X)
}

getNextDesign <- function(morris,X,Y) {
    return()
}

displayResults <- function(morris,X,Y) {
    morris$Yi = Y[,1]
    eval(expr=parse(text="tell(m,Yi)"),envir=morris) #tell(morris$m,Y)
    
    morris$files = "plot.png"
    png(file=morris$files,bg="transparent",height=600,width=600)
    plot(morris$m)
    dev.off()
    
    mu=colMeans(morris$m$ee)
    mu.star=colMeans(abs(morris$m$ee))
    sig=sqrt((colSums(morris$m$ee^2)-mu^2*4)/(morris$r-1))
    effect=rep("Weak effect",length(mu))
    for (i in 1:length(mu)) if (mu.star[i]>0.2*max(mu.star)) if (sig[i]>0.2*max(sig)) effect[i]="Non linear or interactive effect" else effect[i]="Linear effect"
    
    html=paste0("<HTML name='sensitivity'>",paste0(collapse="\n",sub("<table",replacement = "<table border=1",capture.output(print(xtable::xtable(  cbind(mu,mu.star,sig,effect)),type="html")))),"<br/><img src='",morris$files,"'/></HTML>")
    return(html)
}

