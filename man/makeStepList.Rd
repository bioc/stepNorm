\name{makeStepList}
\alias{makeStepList}

\title{Construction of a stepwise normalization list}
\description{
  This function provides a user friendly way to construct a list for input
  to the function \code{\link{stepWithinNorm}}. The list indicates
  intended biases for correction and models for stepwise normalization.
}

\usage{
  makeStepList(A = c("median", "rlm", "loess"), PT = c("median", "rlm",
"loess"), PL = c("median", "rlm", "loess"), Spatial2D = c("rlm2D",
"loess2D", "aov2D", "spatialMedian"))
}

\arguments{
  \item{A}{A character string specifying the normalization models
    for the adjustment of intensity or \eqn{A} bias:
    \describe{
      \item{median: }{global median location normalization}
      \item{rlm: }{global intensity or \eqn{A}-dependent robust linear
	normalization using the \code{\link[MASS]{rlm}} function}
      \item{loess: }{global intensity or \eqn{A}-dependent robust nonlinear
	normalization using the \code{\link[MASS]{loess}} function}
    }
    The user can specify any of these three choices and the selected
    model will be compared based the goodness fit and model parsimony;
    If the correction of the $A$ bias is not desired, the user can set
    \code{A = NULL}.
  }
  \item{PT}{A character string specifying the normalization models
    for the adjustment of print-tip or \eqn{PT} bias:
    \describe{
      \item{median: }{within-print-tip-group median normalization}
      \item{rlm: }{within-print-tip-group robust linear normalization
	using the \code{\link[MASS]{rlm}} function}
      \item{loess: }{within-print-tip-group robust nonlinear
	normalization using the \code{\link[MASS]{loess}} function}
      \item{none: }{no normalization for the \eqn{PT} bias}
    }
    If the correction of the $PT$ bias is not desired, the user can set
    \code{PT = NULL}.
  }
  \item{PL}{A character string specifying the normalization models
    for the adjustment of well-plate or \eqn{PL} bias:
    \describe{
      \item{median: }{within-well-plate median normalization}
      \item{rlm: }{within-well-plate robust linear normalization
	using the \code{\link[MASS]{rlm}} function}
      \item{loess: }{within-well-plate robust nonlinear
	normalization using the \code{\link[MASS]{loess}} function}
      \item{none: }{no normalization for the \eqn{PL} bias}
    }
    If the correction of the $PL$ bias is not desired, the user can set
    \code{PL = NULL}.
   }
   \item{Spatial2D}{A character string specifying the normalization models
     for the adjustment of spatial 2D bias:
     \describe{
       \item{none: }{no normalization for the spatial 2D bias}
       \item{aov2D: }{spatial bivariate location normalization using
	 ANOVA}
       \item{rlm2D: }{spatial bivariate location normalization using the
	 \code{\link[MASS]{rlm}} function}
       \item{loess2D: }{spatial bivariate location normalization using
	 the \code{\link[MASS]{loess}} function}
       \item{spatialMedian: }{spatial location normalization using a
	spatial median approach (see Wilson et al. (2003) in
	reference)}
    }
    If the correction of the $PL$ bias is not desired, the user can set
    \code{Spatial2D = NULL}.
  }
}

\details{
  This function provides a user friendly way to specify the parameter
  \code{wf.loc} for the main stepwise normalization function
  \code{stepWithinNorm}; see examples for details.
}
     

\value{
  An object of class "list" for input to the \code{stepWithinNorm} function.
}

\author{
  Yuanyuan Xiao, \email{yxiao@itsa.ucsf.edu}, \cr
  Jean Yee Hwa Yang, \email{jean@biostat.ucsf.edu}
}

\seealso{\code{\link{stepWithinNorm}}.}

\examples{
# Examples use swirl dataset, for description type ? swirl
data(swirl)

# To use the default parameters, which adjusts A, PT, PL and Spatial2D
# biases sequentially and compares all models available, simple type

wf.loc <- makeStepList()

# To apply loess for the A bias, and to omit the Spatial2D step
wf.loc <- makeStepList(A=("loess"), Spatial2D=NULL)

# To compare only rlm and loess in the A bias step, and other biases as default
wf.loc <- makeStepList(A=c("rlm","loess"))

# input to the stepWithinNorm function
\dontrun{
step.swirl1 <- stepWithinNorm(swirl[,1],wf.loc=wf.loc)}
}      
     
\keyword{models}
