#title: REMBO
#help: Random Embedding Bayesian Optimization (REMBO)
#type: optimization
#author: Mickael Binois
#require: lhs,DiceKriging,DiceView,pso,MASS
#options: initBatchSize='4',batchSize='4',iterations='10',bounds='true|false',trend='y~1',covtype='matern3_2|matern5_2|gauss|powexp|exp',liar='max',search_min='false|true',d='3',low_size='sqrt(d)',kernel='Warped',A_norm='true'
#options.help: initBatchSize=Initial batch size,batchSize=iterations batch size,iterations=number of iterations,bounds=add input variables bounding values (2^d combinations),trend=Universal) kriging trend,covtype=Kriging covariance kernel,liar=liar value for in-batch loop (when batchsize>1),search_min=maximize or minimize,d=Low dimension (integer),low_size=?,kernel=?,A_norm=?

REMBO <- function(options) {
    library(lhs)
    library(DiceKriging)
    library(DiceView)
    library(pso)
    library(MASS)

    # Rembo options: d (low dimension for Rembo),
    #                D (high dimension/ Initial problem dimension)
    #                low_size=[-sqrt(d), sqrt(d)]^d (default box size for the low dimensional search)
    #                kernel = "Warped" (kernel for Rembo: "lowDim", "highDim", "Warped" (default))
    #                A_norm = true (should the row of the embedding matrix be normalized)


    # all parameters are initialy strings, so you have to put as global non-string values
    options$initBatchSize <- as.integer(options$initBatchSize)
    options$batchSize <- as.integer(options$batchSize)
    options$iterations <- as.integer(options$iterations)
    options$bounds <- as.logical(options$bounds)
    options$trend <- as.formula(options$trend)
    options$search_min <- as.logical(options$search_min)
    options$d <- as.integer(options$d)
    options$kernel <- as.character(options$kernel)
    options$A_norm <- as.logical(options$A_norm)

    rembo = new.env()
    rembo$i = 0
    lapply(names(options),function(x) assign(x,options[[x]],rembo))

    rembo$low_size <- eval(parse(text=options$low_size),rembo)

    print(rembo$low_size)

    return(rembo)
}

getInitialDesign <- function(algorithm, d) {
    set.seed(1)

    algorithm$D <- d
    A <- matrix(rnorm(algorithm$D*algorithm$d), algorithm$D, algorithm$d)
    if(algorithm$A_norm)
        A <- A/sqrt(rowSums(A^2))

    algorithm$A <- A

    if(algorithm$kernel == "Warped"){
        algorithm$pA <- A %*% ginv(t(A) %*% A) %*% t(A)
        algorithm$WObs <- NULL # for storing warped dimensional points
    }

    algorithm$lowObs <- numeric(0)  # for storing low dimensional points



    if(algorithm$initBatchSize < 100){
        lhs <- optimumLHS(n = max(algorithm$D+1,algorithm$initBatchSize), k = algorithm$d)
    }else{
        lhs <- maximinLHS(n = max(algorithm$D+1,algorithm$initBatchSize), k = algorithm$d)
    }
    if (algorithm$bounds) {
        e=c(0,1)
        id=1
        while(id<algorithm$d){
            e=rbind(cbind(e,0),cbind(e,1))
            id=id+1
        }
        Xinit=rbind(as.matrix(e),as.matrix(lhs))
    } else {
        Xinit=as.matrix(lhs)
    }

    ## Specific algorithm part
    #1) First resize to [-boxsize,boxsize]^d
    Yinit <- 2*Xinit - 1# design in the low dimensional space
    Xinit <- t(apply(Yinit, 1, mapping_to_X, algorithm$A)) # design is the high dimensional one

    #2) Check that no replicates are present (coming from the convex projection) and eventually replace them
    if(any(duplicated(Xinit))){
        size <- nrow(Xinit) # required size
        Ytemp <- Yinit[!duplicated(Xinit),]
        while(nrow(Ytemp) < size){
            tmp <- size - nrow(Ytemp)
            Yinit <- augmentLHS(Yinit, tmp)
            Xinit <- rbind(Xinit, t(apply(Yinit[-c(1:nrow(Xinit)),,drop = FALSE], 1, mapping_to_X, algorithm$A)))
            Ytemp <- Yinit[!duplicated(Xinit),]
        }
        Yinit <- Ytemp
        Xinit <- Xinit[!duplicated(Xinit),]
    }

    # assign("lowObs", Yinit, envir = algorithm)
    algorithm$lowObs <- Yinit
    if(algorithm$kernel=="Warped"){
        algorithm$WObs <- t(apply(Yinit, 1, Psi, A = algorithm$A, pA = algorithm$pA))
        # assign("WObs", WObs, envir = algorithm)
    }
    return(Xinit)
}

getNextDesign <- function(algorithm, X, Y) {
    if (algorithm$i > algorithm$iterations) return();

    d = algorithm$d
    if (dim(Y)[2] == 2) {
        noise.var <- as.array(Y[,2])^2
    } else {
        noise.var <- NULL
    }

    if (algorithm$search_min) {y=Y[,1]} else {y=-Y[,1]}

    # Different models depending on the kernel
    if(algorithm$kernel == "lowDim"){
        algorithm$kmi <- km(control = list(trace=FALSE), formula = algorithm$trend, optim.method='BFGS',
                        covtype = algorithm$covtype, noise.var = noise.var, design = algorithm$lowObs, response=y,
                        iso = T)
    }else{
        if(algorithm$kernel == "highDim"){
            algorithm$kmi <- km(control = list(trace=FALSE), formula = algorithm$trend, optim.method='BFGS',
                            covtype = algorithm$covtype, noise.var = noise.var, design = X, response=y,
                            iso = T)
        }
        if(algorithm$kernel == "Warped"){
            algorithm$kmi <- km(control = list(trace=FALSE), formula = algorithm$trend, optim.method='BFGS',
                            covtype = algorithm$covtype, noise.var = noise.var, design = algorithm$WObs, response=y,
                            iso = T)
        }
    }


    EGOi <- max_qEI_REMBO(model = algorithm$kmi, npoints = algorithm$batchSize, rembo = algorithm, L = algorithm$liar,
                          lower = -rep(algorithm$low_size, d),
                          upper = rep(algorithm$low_size, d), control=list(trace=FALSE))
    if (is.null(EGOi)) return()

    algorithm$lowObs <- rbind(algorithm$lowObs, EGOi$par)
    if(algorithm$kernel == "Warped"){
        algorithm$WObs <- rbind(algorithm$WObs, t(apply(EGOi$par, 1, Psi, A = algorithm$A, pA = algorithm$pA)))
    }

    Xnext <- t(apply(EGOi$par, 1, mapping_to_X, A = algorithm$A))

    algorithm$i <- algorithm$i + 1
    return(as.matrix(Xnext))
}

displayResults <- function(algorithm, X, Y) {
    algorithm$files <- paste("sectionview_",algorithm$i-1,".png",sep="")
    resolution <- 600

    if (dim(Y)[2] == 2) {
        noise.var <<- as.array(Y[,2])^2
        yname=paste0("N(",names(Y)[1],",",names(Y)[2])
    } else {
        noise.var <<- NULL
        yname=names(Y)
    }

    if (algorithm$search_min) {
        m = min(Y[,1])
        x = as.matrix(X)[which(Y[,1]==m),]
        if(algorithm$kernel == "lowDim"){
            xr = as.matrix(algorithm$lowObs)[which(Y[,1]==m),]
        }
        if(algorithm$kernel == "highDim"){
            xr = as.matrix(X)[which(Y[,1]==m),]
        }
        if(algorithm$kernel == "Warped"){
            xr = as.matrix(algorithm$WObs)[which(Y[,1]==m),]
        }

        html=paste(sep="<br/>",paste("<HTML>minimum is ",m),paste(sep="","found at ",paste(collapse="<br/>",paste(sep="= ",names(x),x)),"<br/><img src='",algorithm$files,"' width='",resolution,"' height='",resolution,"'/></HTML>"))
    } else {
        m = max(Y[,1])
        x = as.matrix(X)[which(Y[,1]==m),]
        if(algorithm$kernel == "lowDim"){
            xr = as.matrix(algorithm$lowObs)[which(Y[,1]==m),]
        }
        if(algorithm$kernel == "highDim"){
            xr = as.matrix(X)[which(Y[,1]==m),]
        }
        if(algorithm$kernel == "Warped"){
            xr = as.matrix(algorithm$WObs)[which(Y[,1]==m),]
        }
        html=paste(sep="<br/>",paste("<HTML>maximum is ",m),paste(sep="","found at ",paste(collapse="<br/>",paste(sep="= ",names(x),x)),"<br/><img src='",algorithm$files,"' width='",resolution,"' height='",resolution,"'/></HTML>"))
    }

    png(file=algorithm$files,bg="transparent",height=resolution,width = resolution)

    try(sectionview.km(algorithm$kmi,center=xr,Xname=names(algorithm$kmi@X),yname=yname))
    dev.off()

    return(html)
}

distXmin <- function(x,Xmin) {
    return(min(sqrt(rowSums((Xmin-matrix(x,nrow=nrow(Xmin),ncol=ncol(Xmin),byrow=TRUE))^2))))
}

EI <- function (x, model, plugin=NULL) {
    if (is.null(plugin)){ if (model@noise.flag) plugin <- min(model@y-2*sqrt(model@noise.var)) else plugin <- min(model@y) }
    m <- plugin

    ########################################################################################
    # Convert x in proper format(s)
    if (!is.matrix(x)) x <- matrix(x,ncol= model@d)
    d <- ncol(x)
    if (d != model@d){ stop("x does not have the right number of columns (",d," instead of ",model@d,")") }
    newdata <- x
    colnames(newdata) = colnames(model@X)

    ########################################################################################
    #cat("predict...")
    predx <- predict.km(object=model, newdata=newdata, type="UK", checkNames = FALSE)
    #cat(" done.\n")
    kriging.mean <- predx$mean
    kriging.sd   <- predx$sd

    xcr <- (m - kriging.mean)/kriging.sd

    xcr.prob <- pnorm(xcr)
    xcr.dens <- dnorm(xcr)
    res <- (m - kriging.mean) * xcr.prob + kriging.sd * xcr.dens

    too.close = which(kriging.sd/sqrt(model@covariance@sd2) < 1e-06)
    res[too.close] <- max(0,m - kriging.mean)

    return(res)
}

mapping_to_X <- function(y, A){
    Xmap <- A %*% y
    Xmap = pmin(Xmap, 1)
    Xmap = pmax(Xmap, -1)
    Xmap <- (Xmap + 1)/2
    return(Xmap)
}

max_EI_REMBO <- function(model, rembo, lower, upper, control=NULL) {

    d <- ncol(rembo$lowObs)

    if (is.null(control$print.level)) control$print.level <- 1
    if (is.null(control$max.parinit.iter)) control$max.parinit.iter <- 10^d
    if(d<=6) N <- 10*2^d else N <- 100*d
    if (is.null(control$pop.size))  control$pop.size <- N
    if (is.null(control$solution.tolerance))  control$solution.tolerance <- 1e-15

    pars=NULL
    for (i in 1:d) pars=cbind(pars,matrix(runif(N,lower[i],upper[i]),ncol=1))

    #t=Sys.time()
    if(rembo$kernel == "highDim"){
        par_W <- t(apply(pars, 1, mapping_to_X, rembo$A))
    }
    if(rembo$kernel == "Warped"){
        par_W <- t(apply(pars, 1, Psi, A = rembo$A, pA = rembo$pA))
    }
    if(rembo$kernel == "lowDim")
        par_W <- pars
    ei <- EI(par_W, model)

    #print(capture.output(Sys.time()-t))
    # print(cbind(pars,ei))

    good_start = which(ei==max(ei,na.rm=T))
    par0=matrix(pars[good_start[sample(1:length(good_start),1)],],nrow=1)

    o <- psoptim(par=par0, fn=function(x){
        if(rembo$kernel == "highDim"){
            x <- t(mapping_to_X(x, rembo$A))
        }
        if(rembo$kernel == "Warped"){
            x <- Psi(x, A = rembo$A, pA = rembo$pA)
        }
        EI(x,model)
    },lower = lower, upper = upper,
    control = list(fnscale = -1, trace = control$print.level, maxit = 10*d))

    o$par <- t(as.matrix(o$par))
    colnames(o$par) <- colnames(rembo$lowObs)
    o$value <- as.matrix(o$value)
    colnames(o$value) <- "EI"

    return(list(par=o$par, value=o$value, counts=o$counts,par.all=o$par.all))
}

max_qEI_REMBO <- function(model, npoints, L, rembo,  lower, upper,  control=NULL, ...) {
    n1 <- nrow(model@X)

    newlowDimPoints <- NULL

    for (s in 1:npoints) {
        oEGO <- max_EI_REMBO(model=model, rembo = rembo, lower=lower, upper=upper, control, ...)

        newlowDimPoints <- rbind(newlowDimPoints, oEGO$par)
        newX <- t(mapping_to_X(as.vector(oEGO$par), rembo$A)) # replace by a apply?

        if(rembo$kernel == "lowDim"){
            if(distXmin(oEGO$par,model@X)<=prod(upper-lower)*1E-10) {
                warning("Proposed a point already in design !");
                npoints=s-1;
                break;
            }
            newPointKrig <- oEGO$par
        }

        if(rembo$kernel == "highDim"){
            if (distXmin(newX,model@X)<=prod(upper-lower)*1E-10) {
                warning("Proposed a point already in design !");
                npoints=s-1;
                break;
            }
            newPointKrig <- newX
        }

        if(rembo$kernel == "Warped"){
            newW <- Psi(as.vector(oEGO$par), A = rembo$A, rembo$pA)
            if (distXmin(newW,model@X)<=prod(upper-lower)*1E-10) {
                warning("Proposed a point already in design !");
                npoints=s-1;
                break;
            }
            newPointKrig <- newW
        }

        model@X <- rbind(model@X, newPointKrig)

        if (L=="min")
            l = min(model@y)
        else if (L=="max")
            l = max(model@y)
        else if (L=="upper95")
            l = predict.km(object = model,newdata = newPointKrig,type="UK",light.return = TRUE)$upper95
        else if (L=="lower95")
            l = predict.km(object = model,newdata = newPointKrig,type="UK",light.return = TRUE)$lower95
        else l = L

        model@y <- rbind(model@y, l, deparse.level=0)

        model@F <- trendMatrix.update(model, Xnew=data.frame(newPointKrig))
        if (model@noise.flag) {
            model@noise.var = c(model@noise.var, 0) # here is the fix!
        }
        newmodel = NULL
        try(newmodel <- computeAuxVariables(model))
        if (is.null(newmodel)) {warning("Unable to update model !");npoints=s-1;break;}
        model = newmodel
    }

    if (npoints==0) return()
    return(list(par = newlowDimPoints, value = model@y[(n1+1):(n1+npoints),, drop=FALSE]))
}

Psi <- function(y, A, pA){

    Xmap <- t(A %*% y)
    Xmap <- pmin(Xmap, 1)
    Xmap <- pmax(Xmap, -1) # convex projection

    Xo <- t(pA %*% t(Xmap)) # orthogonal projection

    if(max(abs(Xo)) > 1){
        pivot <- Xo/max(abs(Xo))
        tmp <- dist(rbind(Xo, Xmap))[1] # distance between Xo and Xmap
        tmp2 <- sqrt(sum(pivot^2)) # norm of pivot
        Xo <- Xo * (tmp2 + tmp)/tmp2 # distortion
    }
    return(Xo)
}

warping <- function(y, A, type){
    if(type == 'none')
        return(y)
    if(type == 'convex')
        return(mapping_to_X(y, A))
    if(type== 'warped')
        return(Psi(y, A))
}
