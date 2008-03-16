\name{marrayFit-class}
\docType{class}

\alias{marrayFit-class}
\alias{marrayFit}

\title{
  Class "marrayFit", storing parameters and results of
  post-normalization cDNA microarray data }

\description{
  A simple list-based class for the storage of parameters
  and results of normalization of cDNA microarray data.}

\section{Creating Objects from the Class}{
Objects can be created by calls of the form \code{new('marrayFit', fit)} where
\code{fit} is a list. Objects of \code{marrayFit} in the \code{StepNorm}
package are typically created by functions \code{\link{fitWithin}} and
\code{\link{fit2DWithin}}.
}

\section{List Components}{
This class contains no slots, but objects should contain the following
list components:
\describe{
  \item{varfun}{: A character vector of names of predictor variables.}
  \item{x}{: A numeric matrix of predictor variables.}
  \item{y}{: A numeric matrix of repsonses.}
  \item{residuals}{: A numeric matrix of normalized values (typically
    log ratios (\eqn{M}))}.
  \item{fitted}{: A numeric matrix of the fitted values.}
  \item{enp}{: The equivalent number of parameters; see \code{\link{loess}}.}
  \item{df.residual}{: The residual degrees of freedom.}
  \item{fun}{: A character string indicating the name of the function
    used for normalization.}
}
}


\section{Methods}{
  This class inherits directly from class \code{list} so any operation
  appropriate for lists will work on objects of this class.
}
    
\author{
   Yuanyuan Xiao, \email{yxiao@itsa.ucsf.edu}, \cr
   Jean Yee Hwa Yang, \email{jean@biostat.ucsf.edu}
}

\seealso{\code{\link{fitWithin}}, \code{\link{fit2DWithin}}.}
\examples{

## load in swirl data
data(swirl)

## median normalization for the first slide of the swirl data
medWithin <- fitWithin(fun="medfit")
## medFit is an object of class marrayFit
medFit <- medWithin(swirl[,1])
## normalized ratios is stored in:
norm.M <- medFit$residuals

}

\keyword{classes}



