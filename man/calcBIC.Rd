\name{calcBIC}
\alias{calcBIC}

\title{Extract BIC from a Fitted Model}
\description{
  Computes the Bayesian Information Criterion for a fitted parametric model.
}
\usage{
calcBIC(fit, subset=TRUE, scale = 0, enp, loss.fun = square)
}

\arguments{
  \item{fit}{fitted model; see details below}
  \item{subset}{A "logical" or "numeric" vector indicating the subset of
          points used to compute the fitted model.}
  \item{scale}{optional numeric specifying the scale parameter of the
    model; see \code{scale} in \code{\link{step}}.}
  \item{enp}{equivalent number of parameters in the fitted model. If
    missing, the \code{enp} component from fit will be used.}
  \item{loss.fun}{the loss function used to calculate deviance; the
    default uses the squared deviation from the fitted  values;
    one could also use absolute deviations (\code{\link{abs}}).}
}

\details{
  The argument \code{fit} can be an object of class
  \code{\link{marrayFit}}, in which case the \code{residuals} component
  from the \code{\link{marrayFit}} object will be extracted to calculate
  the deviance; the user can also pass in a numeric vector, in which
  case it will be interpreted as the residuals and the user needs to
  specify the argument \code{enp}.

  The criterion used is
  \deqn{BIC = -2*log{L} + k * enp,}
  where L is the likelihood and \code{enp} the equivalent number of
  parameters of \code{fit}. For linear models (as in \code{\link{marrayFit}}),
  \eqn{-2log{L}} is computed from the deviance.

  \code{k = log(n)} corresponds to the BIC and is the penalty for
  the number of parameters.
}

\value{
  A numeric vector of length 4, giving
  \item{Dev}{the deviance of the \code{fit}.}
  \item{enp}{the equivalent number of parameters of the
    \code{fit}.}
  \item{penalty}{the penalty for number of parameters.}
  \item{Criterion}{the Akaike Information Criterion for \code{fit}.}
}

\author{
  Yuanyuan Xiao, \email{yxiao@itsa.ucsf.edu}, \cr
  Jean Yee Hwa Yang, \email{jean@biostat.ucsf.edu}
}

\seealso{\code{\link{AIC}}, \code{\link{deviance}}, \code{\link{calcAIC}}.}

\examples{
## load in swirl data
data(swirl)

## fit a model
fit <- fitWithin(fun="medfit")
## res is an object of class marrayFit
res <- fit(swirl[,1])

## calculate BIC
calcBIC(res)
## or could pass in the residual vector, but then argument "enp" needs to be specified
calcBIC(res$residual, enp=1) 
}

\keyword{manip}
