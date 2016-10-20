#title: GradientDescent
#help: First-order local optimization algorithm<br/>http://en.wikipedia.org/wiki/Gradient_descent
#type: optimization
#author: yann.richet@irsn.fr
#require: 
#options: nmax='100',delta='0.1',epsilon='0.01',target='-Inf'
#options.help: nmax=Maximum number of iterations,delta=initial step (in Y) to descend,epsilon=step (in X) to caculate finite differences derivatives,target=Shortcut to terminate if min<target

GradientDescent <- function(options) {
    options$nmax <- as.integer(options$nmax)
    options$delta <- as.numeric(options$delta)
    options$epsilon <- as.numeric(options$epsilon)
    options$target <- as.numeric(options$target)

    gradientdescent = new.env()
    gradientdescent$i = 0
    lapply(names(options), function(x)
        assign(x, options[[x]], gradientdescent))
    return(gradientdescent)
}

getInitialDesign <- function(algorithm, d) {
    algorithm$i = 0
    return(askfinitedifferences(rep(0.5,d),algorithm$epsilon));
}

getNextDesign <- function(algorithm, X, Y) {
    if (algorithm$i>algorithm$nmax) return();
    if (min(Y[,1])<algorithm$target) return();

    d = ncol(X)
    n = nrow(X)

    prevXn = as.matrix(X[(n-d):n,])
    prevYn = as.matrix(Y[(n-d):n,1])

    if (algorithm$i > 1)
        if (Y[n-d,1] > Y[n-d-d,1]) {
            algorithm$delta <- algorithm$delta / 2
            prevXn = as.matrix(X[(n-d-d-1):(n-d-1),])
            prevYn = as.matrix(Y[(n-d-d-1):(n-d-1),1])
        }

    grad = gradient(prevXn,prevYn)
    grad = grad / sqrt(sum(grad * grad))
    xnext = t(prevXn[1,] - (grad * algorithm$delta))
    for (t in 1:d) {
        if (xnext[t] > 1.0) {
            xnext[t] = 1.0;
        }
        if (xnext[t] < 0.0) {
            xnext[t] = 0.0;
        }
    }

    algorithm$i <- algorithm$i+1

    return(askfinitedifferences(xnext,algorithm$epsilon))
}

displayResults <- function(algorithm, X, Y) {
    m = min(Y)
    m.ix=which(Y==m)
    x = as.matrix(X)[m.ix[1],]

    resolution <- 600
    d = dim(X)[2]

    if(d>1) {
        algorithm$files <- paste("pairs_",algorithm$i-1,".png",sep="")
        png(file=algorithm$files,bg="transparent",height=resolution,width = resolution)
        red = (as.matrix(Y)-min(Y))/(max(Y)-min(Y))
        pairs(X,col=rgb(r=red,g=0,b=1-red),Y=Y[[1]],d=d,panel=panel.vec)
        dev.off()
    } else {
        algorithm$files <- paste("plot_",algorithm$i-1,".png",sep="")
        png(file=algorithm$files,bg="transparent",height=resolution,width = resolution)
        red = (as.matrix(Y)-min(Y))/(max(Y)-min(Y))
        plot(x=X[,1],y=Y[,1],xlab=names(X),ylab=names(Y),col=rgb(r=red,g=0,b=1-red))
        dev.off()
    }

    html=paste(sep="<br/>",paste("<HTML name='minimum'>minimum is ",m),paste(sep="","found at ",paste(collapse="= ",capture.output(x)),"<br/><img src='",algorithm$files,"' width='",resolution,"' height='",resolution,"'/></HTML>"))
    plotmin=paste("<Plot1D name='min'>",m,"</Plot1D>")

    if (d == 1) {
        plotx=paste("<Plot1D name='argmin'>",paste(x),"</Plot1D>")
    } else if (d == 2) {
        plotx=paste("<Plot2D name='argmin'>[",paste(collapse=",",x),"]</Plot2D>")
    } else {
        plotx=paste("<PlotnD name='argmin'>[",paste(collapse=",",x),"]</PlotnD>")
    }

    return(paste(html,plotmin,plotx))
}

askfinitedifferences <- function(x,epsilon) {
    xd <- as.array(x);
    for (i in 1:length(x)) {
        xdi <- as.array(x);
        if (xdi[i] + epsilon > 1.0) {
            xdi[i] <- xdi[i] - epsilon;
        } else {
            xdi[i] <- xdi[i] + epsilon;
        }
        xd <- rbind(xd,xdi,deparse.level = 0)
    }
    xd
}

displayResultsTmp <- function(g,X,Y) {
    displayResults(g,X,Y)
}

gradient <- function(xd,yd) {
    d = ncol(xd)
    grad = rep(0,d)
    for (i in 1:d) {
        grad[i] = (yd[i+1] - yd[1]) / (xd[i+1,i] - xd[1,i])
    }
    grad
}

panel.vec <- function(x, y , col, Y, d, ...) {
    #points(x,y,col=col)
    for (i in 1:(length(x)/(d+1))) {
        n0 = 1+(i-1)*(d+1)
        x0 = x[n0]
        y0 = y[n0]
        for (j in 1:d) {
            if (x[n0+j] != x0) {
                dx = (Y[n0]-Y[n0+j])/(max(Y)-min(Y))
                #break;
            }
        }
        for (j in 1:d) {
            if (y[n0+j] != y0) {
                dy = (Y[n0]-Y[n0+j])/(max(Y)-min(Y))
                #break;
            }
        }
        points(x=x0,y=y0,col=col[n0],pch=20)
        lines(x=c(x0,x0+dx),y=c(y0,y0+dy),col=col[n0])
        if (exists("x0p")) {
            lines(x=c(x0p,x0),y=c(y0p,y0),col=col[n0],lty=3)
        }
        x0p=x0
        y0p=y0
    }

}
