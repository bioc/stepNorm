\name{stepWithinNorm}
\alias{stepWithinNorm}

\title{Stepwise within-slide normalization function}
\description{
  This function conducts cDNA microarray normalization in a stepwise
  fashion. In a two-color cDNA array setting, within-slide normalization calibrates signals
  from the two channels to remove non-biological variation introduced by
  various processing steps.
}

\usage{
stepWithinNorm(marraySet, subset=TRUE, wf.loc, criterion = c("BIC", "AIC"), loss.fun = square)
}
\arguments{
  \item{marraySet}{Object of class \code{\link[marrayClasses]{marrayRaw}} or class
    \code{\link[marrayClasses]{marrayNorm}}, containing intensity data for the batch of
    arrays to be normalized.}
  \item{subset}{A "logical" or "numeric" vector indicating the subset of
          points used to compute the  normalization values.}
  \item{wf.loc}{Object of class \code{\link{list}}, each component is a step for the removal of a
    particular systematic varation. Typically each step is also a list
    of several candidate models of different complexity, the best model
    will be chosen by the \code{criterion} specified. For a user
    friendly way of constructing such a list, consult the function
    \code{\link{makeStepList}}.If missing, the
    default procedure will be used, which we consider appropriate for
    most slides. See details for how to specify the list and how it is used.}
  \item{criterion}{Character string specifying the criterion used for the selection of the 
	best normalization procedure in each step. This argument can be
	specified using the  first letter of each method; if no
	specification is made, the default is \code{BIC}:
	\describe{
		\item{AIC:}{the AIC criterion is used}
		\item{BIC:}{the BIC criterion is used.}
	} }
  \item{loss.fun}{loss function; default set at using residual sum of squares.}
}
\details{
  Typical systematic non-biological variations of a two-color cDNA
  microarray include the dependence of ratio measurements (M) on
  intensity (A), print-tip IDs (PT), plate IDs (PL) and spatial
  heterogeneity of the slide (SP). The stepwise normalization procedure
  normalizes a slide in a stepwise fashion. In each step one kind of
  variation is targeted for correction. Within each step, various
  candidate models are assessed for their adequacy with respect to the
  observed data. The assessment is made based on a common model
  selection criterion, AIC (see \code{\link{calcAIC}}) or BIC (see
  \code{\link{calcBIC}}), and the best model is then chosen for the specified
  step.

  The argument \code{wf.loc} is a list of steps. Each step is also a
  list of models. The user uses the function \code{\link{fitWithin}} or
  \code{\link{fit2DWithin}} to specify a model. Below is a table of how
  to do so:
  \tabular{lll}{
    systematic variation \tab model \tab function \cr
    intenstiy (A) \tab median \tab fitWithin(fun="medfit") \cr
    A\tab robust linear \tab fitWithin(fun="rlmfit") \cr
    A\tab robust nonlinear \tab fitWithin(fun="loessfit") \cr
    print-tip (PT) \tab median \tab fitWithin(z.fun="maPrintTip", fun="medfit") \cr
    PT\tab robust linear \tab fitWithin(z.fun="maPrintTip", fun="rlmfit") \cr
    PT\tab robust nonlinear \tab fitWithin(z.fun="maPrintTip",fun="loessfit") \cr
    plate (PL) \tab median \tab fitWithin(z.fun="maCompPlate", fun="medfit") \cr
    PL\tab robust linear \tab fitWithin(z.fun="maComplate", fun="rlmfit") \cr
    PL\tab robust nonlinear \tab fitWithin(z.fun="maCompPlate", fun="loessfit") \cr
    spatial (SP) \tab robust linear \tab fit2DWithin(fun="rlm2Dfit") \cr
    SP\tab robust nonlinear(span=0.2) \tab fit2DWithin(fun="loess2Dfit", span=0.2) \cr
    SP\tab anova \tab fit2DWithin(fun="aov2Dfit") \cr
    SP\tab spatial median (11X11) \tab fit2DWithin(fun="spatialMedfit", width=11)}

    If the \code{wf.loc} is not specified by the user, the default
    procedure conducts normalization in four steps: A -> PT -> PL -> SP
    and models are as described in the table above. The user can choose
    not to follow such a procedure by passing in a different list,
    however we advocate normalizing the intensity (A) variation first as
    it is usually the source of most variation in most slides. The list
    can be easier specified using the function \code{\link{makeStepList}} by
    inputing models as character strings, see \code{\link{makeStepList}}
    for details.
}
     

\value{
  An object of class "list":
  \item{normdata}{an object of class \code{\link[marrayClasses]{marrayNorm}}, containing
  the normalized intensity data.}
  \item{res}{a dataframe of the stepwise normalization result,
    containing the name of the model chosen for each step, deviance,
    equivalent number of parameters, AIC/BIC value.}
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

\seealso{\code{\link{seqWithinNorm}}, \code{\link{withinNorm}},
  \code{\link{fitWithin}}, \code{\link{fit2DWithin}},
  \code{\link{calcAIC}}, \code{\link{calcBIC}}.}

\examples{
# Examples use swirl dataset, for description type ? swirl
data(swirl)

# Apply stepwise normalization for the first slide
res.swirl1 <- stepWithinNorm(swirl[,1])

# normalized data
norm.swirl <- res.swirl1[[1]]

# stepwise procedure
step.swirl <- res.swirl1[[2]]

# using a stepwise procedure different than the default
# corrects intensity (A) and print-tip (PT), this can be
# carried out in two ways:
# 1)
steps <- list(
	    wholeChipA = list(med = fitWithin(fun="medfit"),
                              rlm = fitWithin(fun="rlmfit"),
                              loess = fitWithin(fun="loessfit")),
            printTipA = list(med = fitWithin(z.fun="maPrintTip", fun="medfit"),
                             rlm = fitWithin(z.fun="maPrintTip", fun="rlmfit"),
                             loess = fitWithin(z.fun="maPrintTip",fun="loessfit")))
			     
#2)			     
steps <- makeStepList(PL=NULL, Spatial2D=NULL)
res.swirl <- stepWithinNorm(swirl[,1], wf.loc=steps)


# using AIC criterion for the first slide
res.swirl <- stepWithinNorm(swirl[,1], criterion="A")
}      
     
\keyword{models}
