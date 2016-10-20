#title: Brent
#help: Brent method for 1D root finding
#type: inversion
#author: miguel.munoz-zuniga@ifpen.fr
#require: 
#options: Target='0.0',EPS='3.e-8',tol='1.e-8',iMax='100'
#options.help: Target=Targeted value to inverse,EPS=?,tol=tolerance for root estimation,iMax=Maximum iterations

Brent <- function(options) {
    options$EPS <- as.numeric(options$EPS)
    options$tol <- as.numeric(options$tol)
    options$Target <- as.numeric(options$Target)
    options$iMax <- as.integer(options$iMax)
    
    brent = new.env()
    brent$i = 0
    lapply(names(options), function(x)
        assign(x, options[[x]], brent))
    return(brent)
}

getInitialDesign <- function(algorithm, d) {
    if (d!=1) stop("Cannot find root of >1D function")
    
    algorithm$i <- 0
    algorithm$exit <- -1    # Reason end of algo
    
    x = c(0, 1)
    return(matrix(x))
}

getNextDesign <- function(algorithm, X, Y) {
    # to remain consistent with other iterations, as we expect to have 3 last values of X and Y
    if (algorithm$i <= 0) {
        X = rbind(X,X[2])
        Y = rbind(Y,Y[2])
    }
    
    if (algorithm$i >= algorithm$iMax) {
        algorithm$exit <- 2
        return(NULL)
    }
    algorithm$i <- algorithm$i + 1
    
    Y = Y - algorithm$Target

    a <- as.numeric(X[length(X) - 2, 1])
    b <- as.numeric(X[length(X) - 1, 1])
    c <- as.numeric(X[length(X), 1])
    fa <- as.numeric(Y[length(Y) - 2, 1])
    fb <- as.numeric(Y[length(Y) - 1, 1])
    fc <- as.numeric(Y[length(Y), 1])
    if (algorithm$i == 1 &
        fa * fb > 0) {
        # root must be bracketed for algorithm
        algorithm$exit <- 1
        return(NULL)
    }
    
    if (fb * fc > 0) {
        #Rename a, b, c and adjust bounding interval d
        c <- a
        fc <- fa
        d <- b - a
        e <- d
    } else {
        d = c-b
        e = d
    }

    if (abs(fc) < abs(fb)) {
        # b stand for the best approx of the root which will lie between b and c
        a = b
        b = c
        c = a
        fa = fb
        fb = fc
        fc = fa
    }
    
    tol1 = 2. * algorithm$EPS * abs(b) + 0.5 * algorithm$tol # Convergence check tolerance.
    xm = .5 * (c - b)
    if (abs(xm) <= tol1 | fb == 0) {
        # stop if fb = 0 return root b or tolerance reached
        Xnext = NULL
        algorithm$exit <- 0
        return(Xnext)
    }
    if ((abs(e) >= tol1) & (abs(fa) > abs(fb))) {
        s = fb / fa
        if (a == c) {
            #Attempt linear interpolation
            #print("Alinear")
            p = 2. * xm * s
            q = 1. - s
        } else {
            #Attempt inverse quadratic interpolation.
            #print("Aquadratic")
            q = fa / fc
            r = fb / fc
            p = s * (2. * xm * q * (q - r) - (b - a) * (r - 1.))
            q = (q - 1.) * (r - 1.) * (s - 1.)
        }
        
        if (p > 0) {
            q = -q # Check whether in bounds.
        }
        p = abs(p)
        if (2. * p < min(3. * xm * q - abs(tol1 * q), abs(e * q))) {
            #print("confirmInterpol")
            e <- d #Accept interpolation.
            d <- p / q
        } else {
            #print("bisection1")
            d <- xm #Interpolation failed, use bisection.
            e <- d
        }
    } else {
        # Bounds decreasing too slowly, use bisection.
        #print("bisection2")
        d = xm
        e <- d
    }
    a = b #Move last best guess to a.
    fa = fb
    if (abs(d) > tol1) {
        #then Evaluate new trial root.
        b = b + d
    } else {
        b = b + sign(xm) * tol1
    }
    Xnext = c(a, b, c)
    return(matrix(Xnext, ncol = 1))
}

displayResults <- function(algorithm, X, Y) {
    if (algorithm$exit == 1)
        exit.txt = "root not bracketed"
    else if (algorithm$exit == 2)
        exit.txt = "maximum iteration reached"
    else if (algorithm$exit == 0)
        exit.txt = "algorithm converged"
    else
        exit.txt = paste("error code", algorithm$exit)
    
    algorithm$files <- paste("result", algorithm$i, ".png", sep = "")
    height <- 500
    width <- 500
    
    png(file = algorithm$files,
        height = height,
        width = width)
    plot(X,
         Y,
         pch = 20,
         col = "grey70")
    #plot(as.matrix(X[3*i-1,1]),as.matrix(Y[3*i-1,1]),pch=20,col="grey70")
    abline(h = algorithm$Target,
           lty = 2,
           col = "grey70")
    dev.off()
    
    html <-
        paste0(
            " <HTML name='Root'>In iteration number ",
            algorithm$i,": ",exit.txt,
            ".<br/>",
            "the root approximation is ",
            X[3 * algorithm$i - 1, 1],
            ".<br/>",
            "corresponding to the value ",
            Y[3 * algorithm$i - 1, 1],
            "<br/><img src='",
            algorithm$files,
            "' width='",
            width,
            "' height='",
            height,
            "'/>",
            "</HTML>"
        )
    return(html)
}

