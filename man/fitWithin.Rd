\name{fitWithin}
\alias{fitWithin}

\title{Simple location normalization function for cDNA microarray data}
\description{
  This function performs location normalization on cDNA micoroarray. It
  operates on class \code{\link[marray:marrayRaw-class]{marrayRaw}} or
  class \code{\link[marray:marrayNorm-class]{marrayNorm}}. It allows the
  user to choose from a set of three basic normalization procedures.
}

\usage{
fitWithin(x.fun = "maA", y.fun = "maM", z.fun = TRUE, subset=TRUE, fun = "medfit", ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x.fun}{Name of accessor method for spot intensity, usually \code{maA}.}
  \item{y.fun}{Name of accessor method for spot statistics, usually the log-ratio \code{maM}.}
  \item{z.fun}{Name of accessor method for spot statistic used to
    stratify the data, usually a layout parameter,
    e.g. \code{\link[marray:marrayLayout-class]{maPrintTip}} or
    \code{\link[marray]{maCompPlate}}. If \code{z} is not a
    character, e.g. NULL,  the data are not stratified.}
  \item{subset}{A "logical" or "numeric" vector indicating the subset of
          points used to compute the  normalization values.}
  \item{fun}{Character string specifying the normalization procedure:
    \describe{
      \item{medfit}{for global median location normalization}
      \item{rlmfit}{for global intensity or A-dependent location normalization using 
the \code{\link[MASS]{rlm}} function}
      \item{loessfit}{for global intensity or A-dependent location normalization using 
	the \code{\link[stats]{loess}} function}
      }}
  \item{\dots}{Miscs arguments to be passed in \code{fun}}
}

\details{
  Normalization is typically performed on the expression ratios of cDNA
  microarray data, using the function specified by argument
  \code{fun}. Currently, this function is to be chosen from:
  \code{medfit} (median), \code{rlmfit} (rlm) and
  \code{loessfit}(loess). When \code{z.fun} is provided as a character
  string, for example, \code{maPrintTip}, the normalization procedure is
  operated within each print-tip of the slide.
}

\value{
  The function \code{fitWithin} returns a function(\eqn{F}) with
  bindings for \code{x.fun}, \code{y.fun}, \code{z.fun}, \code{subset}
  and \code{fun}. When the function \eqn{F} is evaluated with an object
  of class \code{\link[marray:marrayNorm-class]{marrayNorm}} or
  \code{\link[marray:marrayRaw-class]{marrayRaw}}, it carries out
  normalization and returns an object of class \code{\link{marrayFit}}
  that contains the normalization information as a list with the following 
  list components:
    \item{varfun}{: A character vector of names of predictor variables.}
    \item{x}{: A numeric matrix of predictor variables.}
    \item{y}{: A numeric matrix of repsonses.}
    \item{residuals}{: A numeric matrix of normalized values (typically
      log ratios (\eqn{M})).}
    \item{fitted}{: A numeric matrix of the fitted values.}
    \item{enp}{: The equivalent number of parameters; see \code{\link[stats]{loess}}.}
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
}

\author{
  Yuanyuan Xiao, \email{yxiao@itsa.ucsf.edu}, \cr
  Jean Yee Hwa Yang, \email{jean@biostat.ucsf.edu}
}

\seealso{\code{\link{fit2DWithin}}}
\examples{
## using the swirl data as example
data(swirl)

## median normalization
med <- fitWithin(fun="medfit")
swirl1.med <- med(swirl[,1])
norm.M <- swirl1.med$residuals ## matrix of normalized ratios

## rlm normalization
rlmF <- fitWithin(fun="rlmfit")
swirl1.rlm <- rlmF(swirl[,1])

## loess normalization, default span=0.4
loessF <- fitWithin(fun="loessfit")
swirl1.loess <- loessF(swirl[,1])
## loess normalization, span=0.2
loessF.1 <- fitWithin(fun="loessfit", span=0.2)
swirl1.loess.1 <- loessF.1(swirl[,1])


## within-printtip loess normalization
loessP <- fitWithin(z.fun="maPrintTip", fun="loessfit")
swirl1.loessP <- loessP(swirl[,1])
}

\keyword{models}

