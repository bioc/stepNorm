\name{maCompPlate2}
\alias{maCompPlate2}

\title{Generate plate IDs}
\description{
  This function is a modification of the \code{\link[marray]{maCompPlate}}
  function in the \code{marray} library. It generates plate IDs
  from the dimensions of the grid and spot matrices. Unlike the
  \code{\link[marray]{maCompPlate}} function, the number of spots is not
  necessarily a multiple of the number of wells on a plate, therefore
  this function allows empty spots on the slide. 
}

\usage{
maCompPlate2(no.plates = NULL, n = 384)
}

\arguments{
  \item{no.plates}{object of class "numeric", number of plates used
    specified by the user. If a number is not specified, then it is
    assumed that there are no empty spots on the slide.}
  \item{n}{object of class "numeric", number of wells in each plate, usually 384 or 96.}
}

\details{
  This function can be used to handle three cases: 1) the number of
  spots is a multiple of the number of wells on a plate (usually 96 or
  384); 2) the number of spots is not a multiple of the number of wells
  on a plate, and several of spots on the slide are therefore left
  empty. In this case, the user needs to specify the number of total
  plates used; plate IDs of empty spots will be NAs; 3)the number of
  spots is not a multiple of the number of wells on a plate, but all
  spots on the slide are spotted, therefore there is one plate not fully
  used. In this case, the user does not need to specify the number of
  total plates (as this will not be an integer), the function assumes no
  empty spots on the slide automatically. See Examples below.  
}
\value{
   The function \code{maCompPlate2} returns a function with bindings for
   \code{no.plates} and \code{n}, which when receiving a object of 
   \code{\link[marray:marrayRaw-class]{marrayRaw}},
   \code{\link[marray:marrayNorm-class]{marrayNorm}}
   or \code{\link[marray:marrayLayout-class]{marrayLayout}} class,
   it returns a vector of plate IDs (\code{\link{factor}}).
}

\author{Yuanyuan Xiao}

\seealso{\code{\link[marray]{maCompPlate}}, \code{\link[marray:marrayLayout-class]{marrayLayout}}.}
\examples{
####### case 1: no empty spots on the slide, full plates used
L<-new("marrayLayout", maNgr=4, maNgc=4, maNsr=22, maNsc=24)
### "compPlate" is a function
compPlate <- maCompPlate2(n=384)
plate <- compPlate(L)
table(plate)
### can also use:
plate<-maCompPlate(L,384)
table(plate)

####### case 2: with empty spots on the slide, full plates used
L<-new("marrayLayout", maNgr=4, maNgc=4, maNsr=22, maNsc=26)
### "compPlate" is a function
compPlate <- maCompPlate2(no.plates=22,n=384)
plate <- compPlate(L)
table(plate)
### empty spots are NAs
unique(plate)

###### case 3: no empty spots on the slide, one plate not full
L<-new("marrayLayout", maNgr=4, maNgc=4, maNsr=22, maNsc=26)
### argument no.plates not specified, the function assumes no empty spots
compPlate <- maCompPlate2(n=384)
plate <- compPlate(L)
### 23 full plates (384), the 24th not full (304)
table(plate)
### no NAs, no empty spots
unique(plate)
}


\keyword{manip}


