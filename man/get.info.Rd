\name{get.info}
\alias{get.info}
\alias{classifier.info}
\alias{print.classifier.info}
\alias{.EMC92_params}
\alias{.HM19_params}
\alias{.IFM15_params}
\alias{.MILLENNIUM100_params}
\alias{.MRCIX6_params}
\alias{.UAMS17_params}
\alias{.UAMS70_params}
\alias{.UAMS80_params}
\title{
Obtain classifier information.
}
\description{
Retrieve the underlying data used for a specific classifier.
}
\usage{
get.info ( classifier )
}

\arguments{
    \item{classifier}{
    A string referring to a classifier in \code{get.classifiers()}.    
    }
}

\value{
    A call to this function gives an object of class \code{classifier.info}. It contains information which is used to determine classifier scores.
    \item{weights}{probes and corresponding weights}
    \item{intercept}{the intercept used}
    \item{decisionBoundaries}{decision boundaries used for classification}
    \item{data}{A Biobase Expression set containing the training data}
    \item{means}{A vector of means per probe set as observed in the training data}
    \item{sds}{A vector of standard deviations per probe set as observed in the training data}
    \item{do.run}{The function called after preprocessing for the actual score calculation}
    \item{citation}{A reference to the original description of the classifier}

}
\seealso{
    \code{\link{applyClassifier}}, \code{\link{Classifier}},  \code{\link{get.classifiers}}, \code{\link[Biobase]{ExpressionSet}}
}
\examples{
get.classifiers()        #obtain all classifier names

emc92info<-get.info("EMC92")
emc92info$weights #obtain all probe sets and corresponding weights.

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
