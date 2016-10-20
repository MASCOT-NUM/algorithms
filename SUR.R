#title: SUR
#help: Global inversion using uncertainty reduction criterion (SUR) algorithm
#type: inversion
#output: Inversion set
#author: clement.chevalier@unine.ch 
#require: lhs,MASS,rgenoud,DiceKriging,KrigInv,DiceView,plotrix,randtoolbox,pbivnorm
#options: initBatchSize='4',batchSize='4',iterations='10',bounds='true|false',trend='y~1',covtype='matern3_2|matern5_2|gauss|exp|powexp',Tlim='0.0'
#options.help: initBatchSize=initial LHS size,batchSize=iterations sample size,iterations=number of iterations,bounds=Add bounding values to initial sample ?,trend=(Universal) kriging trend,covtype=Kriging covariance type,Tlim=Targeted value to inverse

SUR <- function(options) {
    
    library(lhs)
    library(MASS)
    library(rgenoud)
    library(DiceKriging)
    library(DiceView)
    library(plotrix)
    library(randtoolbox)
    library(pbivnorm)
    library(KrigInv)
    
    #all parameters are initialy strings, so you have to put as global non-string values
    options$initBatchSize <- as.integer(options$initBatchSize)
    options$batchSize <- as.integer(options$batchSize)
    options$iterations <- as.integer(options$iterations)
    options$bounds <- as.logical(options$bounds)
    
    options$trend <- as.formula(options$trend)
    
    options$Tlim <- as.numeric(options$Tlim)
    
    
    sur = new.env()
    sur$i <- 0
    
    lapply(names(options),function(x) assign(x,options[[x]],sur))
    return(sur)
}

getInitialDesign <- function(algorithm, d) {
    set.seed(1)
    if(algorithm$initBatchSize < 100){
        lhs <- optimumLHS(n=algorithm$initBatchSize,k=d)
    }else{
        lhs <- maximinLHS(n=algorithm$initBatchSize,k=d)
    }
    if (algorithm$bounds) {
        e=c(0,1)
        id=1
        while(id<d){
            e=rbind(cbind(e,0),cbind(e,1))
            id=id+1
        }
        Xinit=rbind(as.matrix(lhs),as.matrix(e))
    } else {
        Xinit=as.matrix(lhs)
    }
    return(Xinit)
}

getNextDesign <- function(algorithm, X, Y) {
    
    d = dim(X)[2]
    
    if (dim(Y)[2] == 2) {
        algorithm$noise.var <- as.array(Y[,2])^2
    } else {
        algorithm$noise.var <- NULL
    }
    
    if (algorithm$i==1){#initialize VectModel, first time we build a km model
        algorithm$model <- km(formula=algorithm$trend, design = X, response = Y[,1],noise.var = algorithm$noise.var,
                        covtype=algorithm$covtype,control=list(pop.size=50*d,wait.generations=4,
                                                         BFGSburnin=2,max.generations=10*d))
        algorithm$VectModel <- c(algorithm$model)
    }else{   
        #update model
        d = dim(X)[2]
        X.new <- X[(nrow(X)-algorithm$batchSize+1):nrow(X),];
        X.new<-data.frame(X.new)
        y.new <- Y[(nrow(X)-algorithm$batchSize+1):nrow(X),1]
        
        kmcontrol <- NULL
        aSize<-length(algorithm$VectModel)
        #model <- update_km(model=VectModel[[aSize]],NewX=X.new,NewY=y.new,CovReEstimate=TRUE,new.noise.var=new.noise.var,kmcontrol=kmcontrol)
        algorithm$model <- km(formula=algorithm$trend, design = X, response = Y[,1],noise.var = algorithm$noise.var,covtype=algorithm$covtype,control=list(pop.size=50*d,wait.generations=4,BFGSburnin=2,max.generations=10*d))
        algorithm$VectModel <- c(algorithm$VectModel,algorithm$model)
    }
    
    if (algorithm$i > algorithm$iterations) return();
    
    lower<-rep(0,d)
    upper<-rep(1,d)
    integcontrol<-list(n.points=300*d,n.points.among=3000*d,distrib="sur",init.distrib="sobol")
    optimcontrol<-list(method="genoud",pop.size=25*d*algorithm$batchSize,max.generations=10*d,wait.generations=2,unif.seed=1,int.seed=1,optim.option=2)	
    integration.param <- integration_design(integcontrol,d=d,lower=lower,upper=upper,model=algorithm$model,T=algorithm$Tlim)
    
    new.noise.var = ifelse(is.null(algorithm$noise.var),0,mean(algorithm$noise.var))
    
    oEGI <- max_sur_parallel(lower=lower,upper=upper,optimcontrol=optimcontrol,batchsize=algorithm$batchSize,integration.param=integration.param,T=algorithm$Tlim,model=algorithm$model,new.noise.var=new.noise.var)
    
    Xnext <- oEGI$par
    algorithm$uncertainty <- oEGI$value
    
    algorithm$i <- algorithm$i+1
    return(as.matrix(Xnext))
}

displayResults <- function(algorithm, X, Y) {
    
    if(!exists("VectModel",envir=algorithm)) return("initialising")
    
    model<-algorithm$VectModel[[length(algorithm$VectModel)]]
    d <- ncol(X)
    
    if (dim(Y)[2] == 2) {
        noise.var <- as.array(Y[,2])^2
    } else {
        noise.var <- NULL
    }
    
    if(all(model@covariance@range.val<1e-6)){
        #the kriging failed !
        html <- paste(sep="","<HTML>in iteration number ",algorithm$i-2,".<br/> Threshold level is: ",algorithm$Tlim,".<br/> Unable to estimate covariance parameters. The algorithm will evaluate random points...<br/> Use a less smooth covariance kernel, like 'exp', or wait for the next iteration</HTML>")
        return(html)
    }
    
    model2 <- km(formula=algorithm$trend, design = X, response = Y[,1],
                 noise.var = algorithm$noise.var,covtype=algorithm$covtype,
                 control=list(pop.size=50*d,wait.generations=4,BFGSburnin=2,max.generations=10*d))
    
    algorithm$files <- paste("result",algorithm$i-1,".png",sep="")
    height <- 500;width <- 500
    
    png(file=algorithm$files,height=height,width = width)
    TotalUncertainity <- 0
    
    if (d==2) {
        main <- paste("pn(x) after",algorithm$i-2,"iterations of SUR")
        xlab <- names(X)[1]			
        ylab <- names(X)[2]
        
        #calculates scaling
        Axscale<-(model2@X[2,1]-model2@X[1,1])/(model@X[2,1]-model@X[1,1])
        Ayscale<-(model2@X[2,2]-model2@X[1,2])/(model@X[2,2]-model@X[1,2])
        Bxscale <- model2@X[1,1]-Axscale*model@X[1,1]
        Byscale <- model2@X[1,2]-Ayscale*model@X[1,2]
        
        xscale <- c(Bxscale,Bxscale+Axscale)
        yscale <- c(Byscale,Byscale+Ayscale)
        lower <- c(xscale[1],yscale[1])
        upper <- c(xscale[2],yscale[2])
        
        res <- print_uncertainty_2d(model=model2,T=algorithm$Tlim,lower=lower,upper=upper,
                                    new.points=(algorithm$i-2)*algorithm$batchSize,main=main,xlab=xlab,ylab=ylab,
                                    col.points.end="blue",pch.points.end=8,
                                    cex.points=1.5,krigmeanplot=TRUE)
        
        TotalUncertainity <- algorithm$uncertainty
    }
    if(d==1){
        main <- paste("pn(x) after",algorithm$i-2,"iterations of SUR")
        xscale <- c(0,1)
        xlab <- names(X)[1]						
        
        Axscale<-(model2@X[2,1]-model2@X[1,1])/(model@X[2,1]-model@X[1,1])
        Bxscale <- model2@X[1,1]-Axscale*model@X[1,1]
        xscale <- c(Bxscale,Bxscale+Axscale)
        
        res <- print_uncertainty_1d(model=model2,T=algorithm$Tlim,lower=xscale[1],upper=xscale[2],
                                    new.points=(algorithm$i-2)*algorithm$batchSize,xlab=xlab,
                                    main=main,cex.points=1.5,pch.points.end=8,
                                    col.points.end="blue")
        
        TotalUncertainity <- algorithm$uncertainty
    }
    
    if(d > 2){
        #calculates scaling
        lower <- NULL
        upper <- NULL
        for(j in 1:d){
            Ajscale <- (model2@X[2,j]-model2@X[1,j])/(model@X[2,j]-model@X[1,j])
            Bjscale <- model2@X[1,j]-Ajscale*model@X[1,j]
            lower <- c(lower,Bjscale)
            upper <- c(upper,Bjscale + Ajscale)
        }
        res <- print_uncertainty_nd(model=model2,T=algorithm$Tlim,type="pn",lower=lower,upper=upper,
                                    nintegpoints=100,main="maximum probability of excursion",
                                    option="max")
        
        TotalUncertainity <- algorithm$uncertainty
    }
    
    dev.off()
    
    html <- paste(sep="","<HTML>in iteration number ",i-2,".<br/> Threshold level is: ",algorithm$Tlim,".<br/> Remaining volume uncertainity is: ",round(100*TotalUncertainity,2),"%<br/><img src='",algorithm$files,"' width='",width,"' height='",height,"'/></HTML>")
    return(html)
}

displayResultsTmp <- function(sur,X,Y) {
    displayResults(sur,X,Y)
}
