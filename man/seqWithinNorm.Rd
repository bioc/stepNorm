\name{seqWithinNorm}
\alias{seqWithinNorm}

\title{Sequential within-slide normalization function}
\description{
  This function conducts cDNA microarray normalization in a sequential
  fashion. In a two-color cDNA array setting, within-slide normalization calibrates signals
  from the two channels to remove non-biological variation introduced by
  various processing steps.
}

\usage{
  seqWithinNorm(marraySet, y = "maM", subset = TRUE, loss.fun = square,
  A = c("loess", "rlm", "median", "none"),
  PT = c("median", "rlm", "loess", "none"),
  PL = c("median", "rlm", "loess", "none"),
  Spatial2D = c("none", "aov2D", "rlm2D", "loess2D", "spatialMedian"),
  criterion = c("BIC", "AIC"))
}
\arguments{
  \item{marraySet}{Object of class \code{\link[marray:marrayRaw-class]{marrayRaw}}
    or class \code{\link[marray:marrayNorm-class]{marrayNorm}}, containing
    intensity data for the batch of arrays to be normalized.}
  \item{y}{Name of accessor method for spot statistics, usually the
    log-ratio \code{maM}.}
  \item{subset}{A "logical" or "numeric" vector indicating the subset of
    points used to compute the  normalization values.}
  \item{loss.fun}{The loss function used in calculating deviance, the
    default uses squared sum of residuals; for absolute sum of
    residuals, use \code{\link{abs}}}
  \item{A}{A character string specifying the normalization model
    for the adjustment of intensity or \eqn{A} bias:
    \describe{
      \item{loess: }{global intensity or \eqn{A}-dependent robust nonlinear
	normalization using the \code{\link[stats]{loess}} function}
      \item{rlm: }{global intensity or \eqn{A}-dependent robust linear
	normalization using the \code{\link[MASS]{rlm}} function}
      \item{median: }{global median location normalization}
      \item{none: }{no normalization for the \eqn{A} bias}
    }
    If not specified, \code{loess} normalization will be applied.
  }
  \item{PT}{A character string specifying the normalization model
    for the adjustment of print-tip or \eqn{PT} bias:
    \describe{
      \item{median: }{within-print-tip-group median normalization}
      \item{rlm: }{within-print-tip-group robust linear normalization
	using the \code{\link[MASS]{rlm}} function}
      \item{loess: }{within-print-tip-group robust nonlinear
	normalization using the \code{\link[stats]{loess}} function}
      \item{none: }{no normalization for the \eqn{PT} bias}
    }
    If not specified, \code{median} normalization within print-tip will be applied.
  }
  \item{PL}{A character string specifying the normalization model
    for the adjustment of well-plate or \eqn{PL} bias:
    \describe{
      \item{median: }{within-well-plate median normalization}
      \item{rlm: }{within-well-plate robust linear normalization
	using the \code{\link[MASS]{rlm}} function}
      \item{loess: }{within-well-plate robust nonlinear
	normalization using the \code{\link[stats]{loess}} function}
      \item{none: }{no normalization for the \eqn{PL} bias}
    }
    If not specified, \code{median} normalization within well-plate will be applied.
  }
   \item{Spatial2D}{A character string specifying the normalization model
     for the adjustment of spatial 2D bias:
     \describe{
       \item{none: }{no normalization for the spatial 2D bias}
       \item{aov2D: }{spatial bivariate location normalization using
	 ANOVA}
       \item{rlm2D: }{spatial bivariate location normalization using the
	 \code{\link[MASS]{rlm}} function}
       \item{loess2D: }{spatial bivariate location normalization using
	 the \code{\link[stats]{loess}} function}
       \item{spatialMedian: }{spatial location normalization using a
	spatial median approach (see Wilson et al. (2003) in
	reference)}
    }
    If not specified, no normalization will be carried out in this step.
  }
  \item{criterion}{Character string specifying the criterion:
    \describe{
      \item{AIC:}{the AIC criterion is used; see \code{\link{calcAIC}}.}
      \item{BIC:}{the BIC criterion is used; see \code{\link{calcBIC}}.}
    }
    If no specification, \code{BIC} is used. Note that here we don't use
    the criterion to choose normalization model in each step. Criterion
    is calculated solely for informaion purpose.}
}

\details{
  Typical systematic non-biological variations of a two-color cDNA
  microarray include the dependence of ratio measurements (M) on
  intensity (\eqn{A}), print-tip IDs (\eqn{PT}), plate IDs (\eqn{PL}) and spatial
  heterogeneity of the slide (Spatial 2D). The sequential normalization
  procedure in \code{seqWithinNorm} normalizes a slide in a sequential
  fashion: \eqn{A} -> \eqn{PT} -> \eqn{PL} -> Spatial2D. In each step
  one kind of variation is targeted for correction, and the user chooses
  the normalization method as desired. We calculate the AIC/BIC
  criterion along the normalization steps, but they are not used for
  selection of models. 
}
     

\value{
  An object of class "list":
  \item{normdata}{an object of class
    \code{\link[marray:marrayNorm-class]{marrayNorm}}, containing the
    normalized intensity data.}
  \item{res}{a list of the sequential normalization result for each
    slide within the marray dataset. Each list component is also a list
    containing the name of the biases, deviance, equivalent number of
    parameters, AIC/BIC value for a certain slide.}
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

\seealso{\code{\link{stepWithinNorm}}, \code{\link{withinNorm}},
  \code{\link{fitWithin}}, \code{\link{fit2DWithin}},
  \code{\link{calcAIC}}, \code{\link{calcBIC}}.}

\examples{
# Examples use swirl dataset, for description type ? swirl
data(swirl)

# Apply sequential normalization for the first slide
# default: loess(A)->median(PT)->median(PL)-> none (Spatial2D)
\dontrun{
res.swirl1 <- seqWithinNorm(swirl[,1])

# normalized data
norm.swirl <- res.swirl1[[1]]

# sequential normalization information
step.swirl <- res.swirl1[[2]]

}
# median(A)->median(PT)->median(PL)->none(Spatial2D)
res.swirl <- seqWithinNorm(swirl[,1], A="median",PT="median",PL="median",Spatial2D="none")

}      
     
\keyword{models}
