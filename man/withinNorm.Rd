\name{withinNorm}
\alias{withinNorm}

\title{Within-slide normalization function for cDNA spotted microarrays}
\description{
  This function is a wrapper function around \code{fitWtihin} and
  \code{fit2DWithin}. It allows the user to choose from a set of
  thirteen basic location normalization procedures. The function
  operates on an object of class \code{\link[marrayClasses]{marrayRaw}}
  or \code{\link[marrayClasses]{marrayNorm}} and returns an object of
  class \code{\link[marrayClasses]{Norm}}.}
}

\usage{
withinNorm(marraySet,  y = "maM", subset = TRUE, norm = c("none", 
    "median", "rlm", "loess", "medianPrintTip", "rlmPrintTip", 
    "loessPrintTip", "medianPlate", "rlmPlate", "loessPlate", 
    "aov2D", "rlm2D", "loess2D", "spatialMedian"), ...)

}
\arguments{
  \item{marraySet}{Object of class \code{\link[marrayClasses]{marrayRaw}} or class
    \code{\link[marrayClasses]{marrayNorm}}, containing intensity data for the batch of
    arrays to be normalized.}
  \item{y}{Name of accessor method for spot statistics, usually the
    log-ratio \code{maM}.}
  \item{subset}{A "logical" or "numeric" vector indicating the subset of
    points used to compute the  normalization values.}
  \item{norm}{A character string specifying the normalization
    procedures:
    \describe{
      \item{none: }{no normalization}
      \item{median: }{global median location normalization}
      \item{rlm: }{global intensity or A-dependent robust linear
	normalization using the \code{\link[MASS]{rlm}} function}
      \item{loess: }{global intensity or A-dependent robust nonlinear
	normalization using the \code{\link[MASS]{loess}} function}
      \item{medianPrintTip: }{within-print-tip-group median
	normalization}
      \item{rlmPrintTip: }{within-print-tip-group intensity or
	A-dependent robust linear normalization using the
	\code{\link[MASS]{rlm}} function}
      \item{loessPrintTip: }{within-print-tip-group intensity or
	A-dependent robust nonlinear normalization using the
	\code{\link[MASS]{loess}} function}
      \item{medianPlate: }{within-well-plate-group median
	normalization}
      \item{rlmPlate: }{within-well-plate-group intensity or
	A-dependent robust linear normalization using the
	\code{\link[MASS]{rlm}} function}
      \item{loessPlate: }{within-well-plate-group intensity or
	A-dependent robust nonlinear normalization using the
	\code{\link[MASS]{loess}} function}
      \item{aov2D: }{spatial bivariate location normalization using
	ANOVA}
      \item{rlm2D: }{spatial bivariate location normalization using the
	\code{\link[MASS]{rlm}} function}
      \item{loess2D: }{spatial bivariate location normalization using the
	\code{\link[MASS]{loess}} function}
      \item{spatialMedian: }{spatial location normalization using a
	spatial median approach (see Wilson et al. (2003) in
	reference)}
  }}
  \item{\dots}{Misc arguements for the specified \code{norm} function}
}

\details{
  The function \code{withinNorm} dispatches to the function
  \code{\link{fitWithin}} or \code{\link{fit2DWithin}} with specified
  arguments according to the choice of \code{norm}. For instance,
  when \code{norm="loess"} for global intenstiy dependent robust
  nonlinear normalization, \code{withinNorm} calls
  \code{fitWithin(fun="loess")} with the default span parameter set
  at 0.4. If a different span is preferred, it should be input by
  \code{span=0.2} through the argument \code{\dots} in the
  \code{withinNorm} function (see example below). For more details see
  \code{\link{fitWithin}}, \code{\link{fit2DWithin}} and individual
  fitting functions such as \code{\link{loessfit}}.
}
     

\value{
  An object of class \code{\link[marrayClasses]{marrayNorm}}, containing
  the normalized intensity data.
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

\seealso{\code{\link{seqWithinNorm}}, \code{\link{stepWithinNorm}},
  \code{\link{fitWithin}}, \code{\link{fit2DWithin}},
  \code{\link{loessfit}}, \code{\link{rlmfit}}.}

\examples{
# Examples use swirl dataset, for description type ? swirl
data(swirl)

# Apply loess normalization for the first slide, span=0.4
\dontrun{
res.swirl1 <- withinNorm(swirl[,1], norm="loess")}

# Apply loess normalization for the first slide, span=0.2
\dontrun{
res.swirl1 <- withinNorm(swirl[,1], norm="loess", span=0.2)}
}      
     
\keyword{models}

