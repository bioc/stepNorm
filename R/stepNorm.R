
###########################################################
###   fitMethods.R
### A list  of fitting methods
###########################################################
require(MASS)
loessfit <- function(x, y, span=0.4, subset=TRUE,
                     degree=1, family="symmetric",
                     control=loess.control(trace.hat="approximate", iteration=5, surface="direct"),...)
  loess(y~x, span=span, degree=degree, family=family, control=control, subset=subset, ...)

medfit <- function(x, y, subset=TRUE)   rlm(y ~ 1, subset=subset)


rlmfit <- function(x, y, subset=TRUE)  {
  x <- x[subset]
  y <- y[subset]
  rlm(y ~ x)
}

loess2Dfit <- function(x1, x2, y, span=0.2, subset=TRUE,
                     degree=1, family="symmetric",
                     control=loess.control(trace.hat="approximate", iteration=5, surface="direct"),...)
  loess(y~x1 * x2, span=span, degree=degree, family=family, control=control, subset=subset, ...)

aov2Dfit <- function(x1, x2, y, subset=TRUE)     
  lm(y ~ as.factor(x1) + as.factor(x2), subset=subset)

rlm2Dfit <- function(x1, x2, y, subset=TRUE)   rlm(y ~ x1 * x2, subset=subset)

### modified from the tRNA library
### ref: Bioinformatics 2003 July 22;19(11):1325-32
spatialMedfit <- function(x1, x2, y, subset=TRUE, width = 11, height = width)
  {
    xyloc <- cbind(x1, x2)
    m <- matrix(NA, nrow = max(xyloc[,1]), ncol = max(xyloc[,2]))
    m[xyloc][subset] <- as.vector(y)[subset]

    m <- MedianSmooth(m, width, height)

    fitted <- m[xyloc]
    resid <- y - as.vector(m[xyloc])

    res <- list(fitted=fitted, resid=resid, fitted.w=m, width=width, height=height)
    class(res) <- "spatialMed"
    return(res)
  }


MedianSmooth <- function(m, width, height=width)
  {
        W <- c(height, width)
        w <- floor(W/2)
        if(!all(W==1+2*w)) stop("only odd dimensions supported")

        M <- matrix(nrow=nrow(m) + 2*w[1], ncol=ncol(m) + 2*w[2])
        i <- c(row(m) + w[1])
        j <- c(col(m) + w[2])
        ij <- cbind(i,j)
        M[ij] <- m

        coords.r <- (-w[1]):(w[1])
        coords.c <- (-w[2]):(w[2])
        win.med <- function(ij)
          {
            data <- M[ij[1] + coords.r, ij[2] + coords.c]
            median(data, na.rm = TRUE)
          }
        x <- apply(ij, 1, win.med)
        dim(x) <- dim(m)
        return(x)
 }


####################################################################
###	Miscs.R
### Miscellaneous functions

require(marray)
###### if no no.plates  given, then assume no empty spots on the slide
maCompPlate2 <- function(no.plates = NULL,  n = 384) 
{
  function(x)
    {
      is.int <- function(z) trunc(z)==z
      totalPlate <- maNspots(x)/n
      if (is.null(no.plates))
        {
          if (is.int(totalPlate))
            {
              tmp <- n/(maNgr(x) * maNgc(x))
              return(factor(rep(rep(1:totalPlate, rep(tmp, totalPlate)), (maNgr(x) * 
                                                             maNgc(x))))[maSub(x)])
            }
          else
            {
              totalPlate.c <- ceiling(totalPlate)
              totalPlate.f <- floor(totalPlate)
              tmp1 <- n/(maNgr(x) * maNgc(x))
              tmp2 <- ((totalPlate-totalPlate.f) * n ) / (maNgr(x) * maNgc(x))
              return(factor(rep(rep(1:totalPlate.c, c(rep(tmp1, totalPlate.f),tmp2)), (maNgr(x) * 
                                  maNgc(x))))[maSub(x)])
            }
        }
      else
        {
          tmp1 <- n / (maNgr(x) * maNgc(x))
          tmp2 <- maNsr(x) * maNsc(x) - tmp1 * no.plates
          return(factor(rep(c(rep(1:no.plates, rep(tmp1, no.plates)),
                rep(NA, tmp2)), maNgr(x) * maNgc(x))[maSub(x)]))
        } 
    }
}

square<- function (x) x^2

##########################################################
###	WithinNorm.R	
### Within-slide Normalization
##########################################################


noFit <- function(y.fun="maM", x.fun="maA")
  {
    function(marray)
      {
        y <- eval(call(y.fun, marray))
        x <- eval(call(x.fun, marray))
        fit <- list(varfun = c(x=x.fun, y=y.fun),
                   x = as.matrix(x),
                   y = as.matrix(y),
                   residuals = as.matrix(y), 
                   fitted = matrix(NA), 
                   fun = "noFit",
                   enp = 0,
                   df.residual= length(y))
        Res <- new("marrayFit", unclass(fit))
        return(Res)
      }
  }


fitWithin <- function(x.fun="maA", y.fun="maM", z.fun=TRUE, subset=TRUE, fun="medfit", ...)
{
  function(marray)
    {
      y <- eval(call(y.fun, marray))
      x <- eval(call(x.fun, marray))
      ifelse (is.character(z.fun), z <- as.factor(eval(call(z.fun, marray))), z<- as.factor(TRUE))
      subset <- maNum2Logic(length(y), subset)
      good <- !(is.infinite(y) | is.na(y) | is.infinite(x) | is.na(x))
      yfit <- rep(NA, length(y))
      resid <- rep(NA, length(y))
      N <- length(y)
      enp <- 0
      info <- z
      Models <- list()

      dots <- list(...)
      
     
      for (i in levels(info))
        {
          which <- (z == i) & !(is.na(z)) & good
          m <- do.call(fun, c(list(x=x, y=y, subset=which&subset), dots))
          yfit[which] <- predict(m, data.frame(x=x[which]))
          resid[which] <- y[which]-yfit[which]
          enp <- enp + calcEnp(m)
          }
      
      fit <- list(varfun = c(x=x.fun, y=y.fun, z=z.fun),
                  x = as.matrix(x),
                  y = as.matrix(y),
                  residuals = as.matrix(resid), 
                  fitted = as.matrix(yfit), 
                  enp = enp,
                  df.residual= N - enp,
                  fun = fun)
      Res <- new("marrayFit", unclass(fit))
      return(Res)
    }
}

fit2DWithin <- function(x1.fun="maSpotRow", x2.fun="maSpotCol",
                        y.fun="maM", subset=TRUE, fun="aov2Dfit", ...)
{
  function(marray)
    {
      y <- eval(call(y.fun, marray))
      x1 <- (maGridRow(marray) -1) * maNsr(marray) + eval(call(x1.fun, marray))
      x2 <- (maGridCol(marray) -1) * maNsc(marray) + eval(call(x2.fun, marray))
      good <- !(is.infinite(y) | is.na(y) | is.infinite(x1) | is.na(x1) |
                is.infinite(x2) | is.na(x2))
      subset <- maNum2Logic(length(y), subset)
      yfit <- rep(NA, length(y))
      resid <- rep(NA, length(y))
      N <- length(y)
      dots <- list(...)
                
      m <- do.call(fun, c(list(x1=x1, x2=x2, y=y, subset=good&subset),dots))
      if(fun=="spatialMedfit") yfit[good] <- m$fitted[good]
      else yfit[good] <- predict(m, data.frame(x1=x1[good],x2=x2[good]))
      resid[good] <- y[good]-yfit[good]
      enp <- calcEnp(m)
     
      fit <- list(
                 varfun = c(x1=x1.fun,x2=x2.fun,y=y.fun,z="TRUE"),
                 x = cbind(x1,x2),
                 y = as.matrix(y),
                 residuals = as.matrix(resid), 
                 fitted = as.matrix(yfit), 
                 enp = enp,
                 df.residual= N - enp,
                 fun = fun)
                
      Res <- new("marrayFit", unclass(fit))

      return(Res)
    }
}

##########################################################
###	stepNormClass.R	
### Class definition for stepNorm
##########################################################

###########################################################################
## Set Class
require(methods)
.initMarrayFit <- function(where) {
setClass("marrayFit", representation("list"), where=where)
 }

#etClass("marrayFit", representation("list"))

##########################################################
###	Criteria.R	
### Functions to calculate AIC and BIC values
##########################################################

#### calculates AIC values
calcAIC <- function (fit, subset=TRUE, scale = 0, enp, loss.fun=square) 
{

  if (!missing(enp)) enp <- enp
  else enp <- fit$enp
  residuals <- if(is.null(tryCatch(residuals(fit), error = function(e) NULL))) fit[subset]
               else residuals(fit)[subset]
  n <- length(residuals)
  k <- 2*enp
  RSS <- sum.na(loss.fun(residuals))
  dev <- if (scale > 0) RSS/scale - n
         else n * log(RSS/n)
  return(c(Dev=dev, enp=enp, penalty=2, Criterion=dev+k) )
}

#### calculates BIC values
calcBIC <- function (fit, subset=TRUE, scale = 0, enp, loss.fun=square) 
{

  if (!missing(enp)) enp <- enp
  else enp <- fit$enp
  residuals <- if(is.null(tryCatch(residuals(fit), error = function(e) NULL))) fit[subset]
               else residuals(fit)[subset]
  n <- length(residuals)
  k <- enp*log(n)
  RSS <- sum.na(loss.fun(residuals))
  dev <- if (scale > 0) RSS/scale - n
         else n * log(RSS/n)
  return(c(Dev=dev, enp=enp, penalty=log(n), Criterion=dev+k) )
}

#### extracts degrees of freedom
calcEnp <- function(X)
  {
    if((data.class(X) == "rlm") |(data.class(X)=="lm")) return(summary(X)$df[1])
    else if(data.class(X) == "loess") return(X$enp)
    else if (data.class(X) == "spatialMed") return(length(X$resid)/(X$width * X$height))
    else stop("input must be of class rlm or loess...")
  }


##########################################################
###	stepNorm.R	
### main function for stepNorm
##########################################################

## if using the default step procedures and BIC as the criterion,
## the user can just pass in the data, which should be of class
## "marrayRaw" or "marrayNorm"
## for example: dd <- stepWithinNorm(swirl)

makeStepList <- function(A=c("median","rlm","loess"), PT=c("median","rlm","loess"),
                      PL=c("median","rlm","loess"),Spatial2D=c("rlm2D","loess2D","aov2D","spatialMedian")) {
  Bias1 <- c("median","rlm","loess")
  Bias2 <- c("rlm2D","loess2D","aov2D","spatialMedian")
  if(!is.null(A)) {
    As <- pmatch(A, Bias1)
    As <- Bias1[As]
    nA <- length(As)
  }
  else nA <- 0
  if(!is.null(PT)) {
    PTs <- pmatch(PT, Bias1)
    PTs <- Bias1[PTs]
    nPT <- length(PTs)
  }
  else nPT <- 0
  if(!is.null(PL)) {
    PLs <- pmatch(PL, Bias1)
    PLs <- Bias1[PLs]
    nPL<- length(PL)
  }
  else nPL <- 0
  if(!is.null(Spatial2D)) {
    SPs <- pmatch(Spatial2D, Bias2)
    SPs <- Bias2[SPs]
    nSP <- length(SPs)
  }
  else nSP <- 0
  
  NormSteps <- list()
  if (nA>0){
    A.list <- list()
    for (i in 1:nA) {
      if (As[i]=="median") A.list <- c(A.list, list(median=fitWithin(fun="medfit")))
      if (As[i]=="rlm") A.list <- c(A.list, list(rlm=fitWithin(fun="rlmfit")))
      if (As[i]=="loess") A.list <- c(A.list, list(loess=fitWithin(fun="loessfit")))
    }
   NormSteps <- c(NormSteps, list(WholeChipA=A.list)) 
  }

 
  if (nPT>0){
    PT.list <- list()
    for (i in 1:nPT) {
      if (PTs[i]=="median") PT.list <- c(PT.list, list(median=fitWithin(z.fun="maPrintTip",fun="medfit")))
      if (PTs[i]=="rlm") PT.list <- c(PT.list, list(rlm=fitWithin(z.fun="maPrintTip",fun="rlmfit")))
      if (PTs[i]=="loess") PT.list <- c(PT.list, list(loess=fitWithin(z.fun="maPrintTip",fun="loessfit")))
    }
    NormSteps <- c(NormSteps, list(PrintTipA=PT.list))
  }

   if (nPL>0){
    PL.list <- list()
    for (i in 1:nPL) {
      if (PLs[i]=="median") PL.list <- c(PL.list, list(median=fitWithin(z.fun="maCompPlate",fun="medfit")))
      if (PLs[i]=="rlm") PL.list <- c(PL.list, list(rlm=fitWithin(z.fun="maCompPlate",fun="rlmfit")))
      if (PLs[i]=="loess") PL.list <- c(PL.list, list(loess=fitWithin(z.fun="maCompPlate",fun="loessfit")))
    }
    NormSteps <- c(NormSteps, list(PlateA=PL.list))
  }

   if (nSP>0){
    SP.list <- list()
    for (i in 1:nSP) {
      if (SPs[i]=="rlm2D") SP.list <- c(SP.list, list(rlm2D=fit2DWithin(fun="rlm2Dfit")))
      if (SPs[i]=="loess2D") SP.list <- c(SP.list, list(loess2D=fit2DWithin(fun="loess2Dfit")))
      if (SPs[i]=="aov2D") SP.list <- c(SP.list, list(aov2D=fit2DWithin(fun="aov2Dfit")))
      if (SPs[i]=="spatialMedian") SP.list <- c(SP.list, list(spatialMedian=fit2DWithin(fun="spatialMedfit")))
    }
    NormSteps <- c(NormSteps, list(Spatial2D=SP.list))
  }

  if (length(NormSteps)==0) stop("no normalization method selected")
  return(NormSteps)
}
  

stepWithinNorm <- function(marraySet, subset=TRUE,
			   wf.loc, 
			   criterion = c("BIC", "AIC"),
                           loss.fun = square){
  nSlide <- maNsamples(marraySet)
  nGene <- maNspots(marraySet)
  subset <- maNum2Logic(nGene, subset)
  nSub <- length(c(1:nGene)[subset])
  mnormSet <- as(marraySet, 'marrayNorm')
  slot(mnormSet, "maNormCall") <- match.call()
  criterion <- match.arg(criterion)
  crtFun <- switch(criterion, AIC = calcAIC,
                              BIC = calcBIC)
  if (missing(wf.loc)) {
    wf.loc <- list(
            wholeChipA = list(med = fitWithin(fun="medfit", subset=subset),
                              rlm = fitWithin(fun="rlmfit", subset=subset),
                              loess = fitWithin(fun="loessfit", subset=subset)),
            printTipA = list(med = fitWithin(z.fun="maPrintTip", fun="medfit", subset=subset),
                             rlm = fitWithin(z.fun="maPrintTip", fun="rlmfit", subset=subset),
                             loess = fitWithin(z.fun="maPrintTip", fun="loessfit", subset=subset)),
	    plateA = list(med = fitWithin(z.fun="maCompPlate", fun="medfit", subset=subset),
                          rlm = fitWithin(z.fun="maCompPlate", fun="rlmfit", subset=subset),
                          loess = fitWithin(z.fun="maCompPlate", fun="loessfit", subset=subset)),
            wholeChipSpatial = list(rlm2D = fit2DWithin(fun="rlm2Dfit", subset=subset),
                                    loess2D = fit2DWithin(fun="loess2Dfit", span=0.2, subset=subset),
                                    aov2D = fit2DWithin(fun="aov2Dfit", subset=subset),
                                    spatialMed = fit2DWithin(fun="spatialMedfit",width=11, subset=subset)))
  }
  
  nStep <- length(wf.loc)
        
  # give biases and norm methods names if necessary
  if(is.null(names(wf.loc))) names(wf.loc) <- paste("Bias", LETTERS[c(1:nStep)], sep=" ")
  sapply(wf.loc, function(z) {
    if (is.null(names(z))) names(z) <- paste("Model", c(1:length(z)), sep=" ")})

  # within norm
  Res <- list()
  for (i in 1:nSlide) {
    cat('Normalizing slide ', i, '...\n\n')
    marray <- marraySet[,i]
    mnorm <- as(marray, 'marrayNorm')
    count <- 0

    # null model -- no fitting
    fit0 <- noFit()(marray[subset])
    bCrt <-  nSub*log(sum.na(loss.fun(fit0$resid[subset]))/nSub)
    cat(criterion, "of null model: ", bCrt, "\n\n")

    # initialize values
    cCrt <- bCrt
    enp <- 0

    # vectors that contain results
    normStep <- c()
    From <- " "
    To <- " "
    Dev <- bCrt
    Enp <- "0.00"
    Penalty <- " "
    CRT <- bCrt

    while(count < nStep) {
      fits <- list()
      mCrtM <- c()
      count <- count + 1
      step <- wf.loc[[count]]

      msg <- paste(count, names(wf.loc)[count], sep=" -- ")
      cat("step ", msg, ":\n")
                  
      for (funName in step) {
        fit <- funName(mnorm)
        fits <- c(fits, list(fit))
	mCrtM <- rbind(mCrtM, crtFun(fit, subset=subset, loss.fun=loss.fun)[c("Dev", "enp", "penalty")])
      }
      mCrt <- mCrtM[,"Dev"] + (enp + mCrtM[,"enp"]) * mCrtM[,"penalty"]
      metNames <- names(step)
      cat(criterion, "of methods (", metNames, ") are: ", mCrt, "\n") 
      chosen <- which.min(mCrt)
                       
      if (mCrt[chosen] < cCrt)	{
        normStep <- c(normStep,paste(names(wf.loc)[count], "-", metNames[chosen], sep=""))
        maM(mnorm) <- as.matrix(fits[[chosen]]$resid)
        cCrt <- mCrt[chosen]
        enp <- enp + mCrtM[chosen, "enp"]
	From <- c(From, " ")
        To <- c(To, metNames[chosen])
        Dev <- round(c(Dev, mCrtM[chosen, "Dev"]),2)
        Enp <- c(Enp, paste("+", round(mCrtM[chosen,"enp"],2), sep=""))
        Penalty <- c(Penalty, round(mCrtM[chosen,"penalty"],2))
        CRT <- round(c(CRT, cCrt),2)
        cat("chosen : ", metNames[chosen], "\n\n")
      }
      else {
        From <- c(From, " ")
        To <- c(To, " ")
        Dev <- c(Dev, Dev[length(Dev)])
        Enp <- c(Enp, "+0.00")
        Penalty <- c(Penalty, " ")
        CRT <- c(CRT, CRT[length(CRT)])
        cat("this normalization step is not necessary\n\n")
      }
    }
                        
    maM(mnormSet)[,i] <- maM(mnorm)
    aovTable <- data.frame(From, To, Dev, Enp, Penalty, CRT)
    dimnames(aovTable) <- list(c("null", names(wf.loc)), c("From", "To", "Deviance", "Enp", "Penalty", "Criterion"))
    Res <- c(Res, list(aovTable))
    cat("Slide", i, "normalization steps: ", paste(normStep, collapse="-> "), "\n\n")
  }
  return(list(normdata=mnormSet, res=Res))
}

withinNorm <- function(marraySet, y="maM", subset=TRUE,
                       norm=c("none", "median", "rlm", "loess",
                         "medianPrintTip", "rlmPrintTip", "loessPrintTip",
                         "medianPlate", "rlmPlate", "loessPlate",
                         "aov2D", "rlm2D", "loess2D", "spatialMedian"), ...)  {
    norm.method <- match.arg(norm)
    cat(norm.method)
    func <- switch(norm.method,
            none=noFit(y.fun=y),
                   
            median = fitWithin(y.fun=y, z.fun=TRUE, subset=subset, fun="medfit"),
            rlm = fitWithin(y.fun=y, z.fun=TRUE, subset=subset, fun="rlmfit"),
            loess = fitWithin(y.fun=y, z.fun=TRUE, subset=subset, fun="loessfit", ...),
                   
            medianPrintTip = fitWithin(y.fun=y, z.fun="maPrintTip", subset=subset, fun="medfit"),
            rlmPrintTip = fitWithin(y.fun=y, z.fun="maPrintTip", subset=subset, fun="rlmfit"),
            loessPrintTip = fitWithin(y.fun=y, z.fun="maPrintTip", subset=subset, fun="loessfit",...),
                   
            medianPlate = fitWithin(y.fun=y, z.fun="maCompPlate", subset=subset, fun="medfit"),
            rlmPlate = fitWithin(y.fun=y, z.fun="maCompPlate", fun="rlmfit"),
            loessPlate = fitWithin(y.fun=y, z.fun="maCompPlate", fun="loessfit", ...),

            aov2D = fit2DWithin(y.fun=y, subset=subset, fun="aov2Dfit"),
            rlm2D = fit2DWithin(y.fun=y, subset=subset, fun="rlm2Dfit"),  
            loess2D = fit2DWithin(y.fun=y, subset=subset, fun="loess2Dfit",...),
            spatialMedian = fit2DWithin(y.fun=y, subset=subset, fun="spatialMedfit", ...))
                   
    mnormSet <- as(marraySet, 'marrayNorm')
    slot(mnormSet, "maNormCall") <- match.call()
    nSlide <- maNsamples(marraySet)
    for (i in 1:nSlide){
      marray <- marraySet[,i]
      fit <- func(marray)
      maM(mnormSet)[,i] <- fit$resid
    }
    return(mnormSet)
  }

seqWithinNorm <- function(marraySet, y="maM", subset=TRUE, loss.fun=square,
                          A=c("loess","rlm","median","none"), PT=c("median","rlm","loess","none"),
                          PL=c("median","rlm","loess","none"),
                          Spatial2D=c("none","aov2D","rlm2D","loess2D","spatialMedian"),
                          criterion=c("BIC","AIC")) {
  nSlide <- maNsamples(marraySet)
  nGene <- maNspots(marraySet)
  mnormSet <- as(marraySet, 'marrayNorm')
  slot(mnormSet, "maNormCall") <- match.call()
  Res <- list()
  criterion <- match.arg(criterion)
  crtFun <- switch(criterion, AIC = calcAIC,
                              BIC = calcBIC)
  A <- match.arg(A, c("loess","rlm","median","none"))
  PT <- match.arg(PT, c("median","rlm","loess","none"))
  PL <- match.arg(PL, c("median","rlm","loess","none"))
  Spatial2D <- match.arg(Spatial2D, c("none","aov2D","rlm2D","loess2D","spatialMedian"))
  A.fun <- switch(A,
                    loess = fitWithin(y.fun=y, z.fun=TRUE, subset=subset, fun="loessfit"),
                    rlm = fitWithin(y.fun=y, z.fun=TRUE, subset=subset, fun="rlmfit"),
                    median = fitWithin(y.fun=y, z.fun=TRUE, subset=subset, fun="medfit"),
                    none = noFit(y.fun=y))
  PT.fun <- switch(PT,
                     median = fitWithin(y.fun=y, z.fun="maPrintTip", subset=subset, fun="medfit"),
                     rlm = fitWithin(y.fun=y, z.fun="maPrintTip", subset=subset, fun="rlmfit"),
                     loess = fitWithin(y.fun=y, z.fun="maPrintTip", subset=subset, fun="loessfit"),
                     none =  noFit(y.fun=y))
  PL.fun <- switch(PL,
                     median = fitWithin(y.fun=y, z.fun="maCompPlate", subset=subset, fun="medfit"),
                     rlm = fitWithin(y.fun=y, z.fun="maCompPlate", subset=subset, fun="rlmfit"),
                     loess = fitWithin(y.fun=y, z.fun="maCompPlate", subset=subset, fun="loessfit"),
                     none = noFit(y.fun=y))
  Spatial2D.fun <- switch(Spatial2D,
                   none = noFit(y.fun=y),
                   aov2D = fit2DWithin(y.fun=y, subset=subset, fun="aov2Dfit"),
                   rlm2D = fit2DWithin(y.fun=y, subset=subset, fun="rlm2Dfit"),  
                   loess2D = fit2DWithin(y.fun=y, subset=subset, fun="loess2Dfit"),
                   spatialMedian = fit2DWithin(y.fun=y, subset=subset, fun="spatialMedfit"))

  wf.loc <- list(A=A.fun,PT=PT.fun,PL=PL.fun,Spatial2D=Spatial2D.fun)

  calcDev <- function(X, subset=TRUE, loss.fun=square){
    return(length(X[subset])*log(sum.na(loss.fun(X[subset]))/length(X[subset])))
  }
  
  for (i in 1:nSlide) {
    cat('\n\nNormalizing slide ', i, '...\n\n')
    marray <- marraySet[,i]
    mnorm <- as(marray, 'marrayNorm')
    count <- 0
    Dev <- calcDev(maM(mnorm), subset=subset)
    enp <- 0
    BIC <- Dev
    cat("Deviance of null model:", Dev, "\n")
    cat("enp", "\t", "Deviance", "\t", "BIC", "\n")
    for (j in 1:length(wf.loc)) {
      func <- wf.loc[[j]]
      fit <- func(mnorm)
      maM(mnorm) <- as.matrix(fit$residuals)
      enp <- c(enp, fit$enp)
      Dev <- c(Dev, calcDev(maM(mnorm), subset=subset))
      BIC <- c(BIC, crtFun(maM(mnorm), subset=subset, enp=sum(enp)))
      cat(round(enp[length(enp)],2),"\t", round(Dev[length(Dev)],2), "\t", round(BIC[length(BIC)],2), "\n")
    }
    maM(mnormSet)[,i] <- maM(mnorm)
    tmp <- list(c("null", names(wf.loc)), enp, Dev, BIC)
    names(tmp) <- c("Bias", "enp", "deviance", criterion)
    Res <- c(Res, list(tmp))
  }
      
  return(list(normdata=mnormSet, res=Res))
}


##########################################################################
.First.lib <- function(libname, pkgname, where) {
    require(methods)
    #require(Biobase) || stop("can't load without Biobase")
    require(marray)
    require(MASS)
    if(missing(where)) {
        where <- match(paste("package:", pkgname, sep=""), search())
        if(is.na(where)) {
            warning(paste("Not a package name: ",pkgname))
            return()
        }
    where <- pos.to.env(where)
    }
    .initMarrayFit(where)
    cacheMetaData(as.environment(where))
  }





