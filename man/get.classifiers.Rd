\name{get.classifiers}
\alias{get.classifiers}

\title{
Obtain classifier names.
}
\description{
The \code{\link{applyClassifier}} function has a number of build-in classifiers which can be used. This function provides a list of them.
}
\usage{
get.classifiers ( normalization )
}

\arguments{
    \item{normalization}{
    an optional string refering to the normalization methods. One of: \code{mas5},\code{gcrma},\code{frma},\code{rma}. If 'NULL' or empty, all classifiers are obtained.
}
}

\value{
    A character vector of classifier names
}
\seealso{
\code{\link{get.info}}, \code{\link{applyClassifier}}, \code{\link{Classifier}}, \code{\link[Biobase]{ExpressionSet}}
}
\examples{
get.classifiers()        #obtain all classifiers
get.classifiers('mas5')    #obtain all classifiers which use mas5 normalized input
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
