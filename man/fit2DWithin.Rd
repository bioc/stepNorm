\name{fit2DWithin}
\alias{fit2DWithin}

\title{Bivariate location normalization function for cDNA microarray data}
\description{
  This function performs 2D location normalization on cDNA
  micoroarray. It operates on class
  \code{\link[marray:marrayRaw-class]{marrayRaw}} or class
  \code{\link[marray:marrayNorm-class]{marrayNorm}}. It allows the user
  to choose from a set of four basic normalization procedures.
}

\usage{
fit2DWithin(x1.fun = "maSpotRow", x2.fun = "maSpotCol", y.fun = "maM",
subset=TRUE, fun = aov2Dfit, ...)
}

\arguments{
  \item{x1.fun}{Name of accessor method for spot row coordinates, usually \code{maSpotRow}.}
  \item{x2.fun}{Name of accessor method for spot column coordinates, usually \code{maSpotCol}.}
  \item{y.fun}{Name of accessor method for spot statistics, usually the
    log-ratio \code{maM}.}
  \item{subset}{A "logical" or "numeric" vector indicating the subset of
    points used to compute the  normalization values.}
  \item{fun}{Character string specifying the normalization procedures:
    \describe{
      \item{rlm2Dfit}{for robust linear regression using the
	\code{\link[MASS]{rlm}} function}
      \item{loess2Dfit}{for robust local regression using the
	\code{\link[stats]{loess}} function}
      \item{aov2Dfit}{for linear regression using the \code{\link{lm}}
	function}
      \item{spatialMedfit}{for spatial median normalization}
  }}
  \item{\dots}{Misc arguments for \code{fun}}
}
\details{
  The spot statistic named in \code{y} is regressed on spot row and
  column coordinates, using the function specified by the argument
  \code{fun}. Typically, \code{rlm2Dfit} and \code{loess2Dfit}, which
  treat row and column coordinates as numeric vectors, require a lot fewer parameters than
  \code{aov2Dfit} which specifies these two variables as
  categorical. \code{spatialMedfit} could yet fit the most complicated
  model, depending on size of the smoothing window specified; details
  see Wison et al (2003).
}

\value{
  The function \code{fit2DWithin} returns a function (\eqn{F}) with
  bindings for \code{x1.fun}, \code{x2.fun}, \code{y.fun}, \code{subset}
  and \code{fun}. When the function \eqn{F} is evaluated with an object
  of class \code{\link[marray:marrayNorm-class]{marrayNorm}} or
  \code{\link[marray:marrayNorm-class]{marrayRaw}}, it carries out
  normalization and returns an object of class \code{\link{marrayFit}}
  that contains the normalization information as a list with the following
  components: 
    \item{varfun}{: A character vector of names of predictor variables.}
    \item{x}{: A numeric matrix of predictor variables.}
    \item{y}{: A numeric matrix of responses.}
    \item{residuals}{: A numeric matrix of normalized values (typically
      log ratios (\eqn{M})).}
    \item{fitted}{: A numeric matrix of the fitted values.}
    \item{enp}{: The equivalent number of parameters; see
      \code{\link[stats]{loess}}.} 
    \item{df.residual}{: The residual degrees of freedom.}
    \item{fun}{: A character string indicating the name of the function
      used for normalization.}
  Note that the \code{residuals} component stores the normalized ratios.
}

\references{
  Y. H. Yang, S. Dudoit, P. Luu, and T. P. Speed (2001). Normalization
  for cDNA microarray data. In M. L. Bittner, Y. Chen, A. N. Dorsel, and
  E. R. Dougherty (eds), \emph{Microarrays: Optical Technologies and
    Informatics}, Vol. 4266 of \emph{Proceedings of SPIE}.
  
  D. L. Wilson, M. J. Buckley, C. A. Helliwell and I. W. Wilson
  (2003). New normalization methods for cDNA microarray
  data. \emph{Bioinformatics}, Vol. 19, pp. 1325-1332.
}

\author{
  Yuanyuan Xiao, \email{yxiao@itsa.ucsf.edu}, \cr
  Jean Yee Hwa Yang, \email{jean@biostat.ucsf.edu}
}

\seealso{\code{\link{fitWithin}}}
\examples{
## use the swirl data as example
data(swirl)

## 2D rlm normalization
rlm2D <- fit2DWithin(fun="rlm2Dfit")
swirl1.rlm <- rlm2D(swirl[,1])
norm.M <- swirl1.rlm$residuals ## matrix of normalized ratios

## 2D loess normalization, default span=0.2
loess2D <- fit2DWithin(fun="loess2Dfit")
swirl1.loess <- loess2D(swirl[,1])
## 2D loess normalization, span=0.4
\dontrun{
loess2D.1 <- fit2DWithin(fun="loess2Dfit", span=0.4)
swirl1.loess.1 <- loess2D.1(swirl[,1])}

## 2D aov normalization
aov2D <- fit2DWithin(fun="aov2Dfit")
swirl1.aov <- aov2D(swirl[,1])

## 2D spatial median normalization, default window width=3
spatialMed2D <- fit2DWithin(fun="spatialMedfit")
swirl1.spatialMed <- spatialMed2D(swirl[,1])
## 2D loess normalization, window width=9
\dontrun{
spatialMed2D.1 <- fit2DWithin(fun="spatialMedfit", width=9)
swirl1.spatialMed.1 <- spatialMed2D.1(swirl[,1])}
}

\keyword{models}
