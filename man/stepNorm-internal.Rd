\name{stepNorm-internal}
\alias{calcEnp}
\alias{noFit}
\alias{medfit}
\alias{loessfit}
\alias{rlmfit}
\alias{loess2Dfit}
\alias{aov2Dfit}
\alias{rlm2Dfit}
\alias{spatialMedfit}
\alias{MedianSmooth}
\alias{square}

\title{Internal stepNorm functions}
\description{
 Internal stepNorm functions.
}
\usage{
calcEnp(X)
noFit(y.fun="maM", x.fun="maA")
medfit(x, y, subset=TRUE)
loessfit(x, y, span=0.4, subset=TRUE, degree=1, family="symmetric", 
		control=loess.control(trace.hat="approximate", iteration=5, surface="direct"),...)
rlmfit(x, y, subset=TRUE)
loess2Dfit(x1, x2, y, span=0.2, subset=TRUE, degree=1, family="symmetric",
            control=loess.control(trace.hat="approximate", iteration=5, surface="direct"),...)
aov2Dfit(x1, x2, y, subset=TRUE)
rlm2Dfit(x1, x2, y, subset=TRUE)
spatialMedfit(x1, x2, y, subset=TRUE, width = 11,  height = width)
MedianSmooth(m, width, height=width)
square{x}
}

\details{
  These are not to be called directly by the user.
}

\keyword{internal}

