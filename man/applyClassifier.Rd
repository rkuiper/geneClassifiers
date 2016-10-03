\name{applyClassifier}
\alias{applyClassifier}
\alias{applyClassifiers}
\alias{applyClassifier.mas5}
\alias{applyClassifier.gcrma}
\alias{applyClassifier.frma}
\alias{applyClassifier.rma}
\alias{geneClassifiers}


\title{
Apply a classifier
}
\description{
Note: Results of this package are allowed to be used for research only!
Apply a classifier to a gene expression set.
}
\usage{
applyClassifier.mas5 (which, eset, dataislog2 = FALSE,used.MAS5TargetValue=NULL,...)
applyClassifier.gcrma(which, eset, ...)
applyClassifier.frma (which, eset, ...)
applyClassifier.rma  (which, eset, ...)
}
\arguments{
    \item{which}{
a string refering to the classifier of interest. Valid classifiers can be shown with e.g. \code{get.classifiers('mas5')} for mas5 classifiers.
}
    \item{eset}{
    a Bioconductor Biobase \code{\link[Biobase]{ExpressionSet}}. This set should be normalized with the appropriate method (i.e. mas5,gcrma,frma or rma) depending on the classifier to use.
}
\item{dataislog2}{ 
    If the expression set you use for input is already log2 transformed you should set this argument to \code{TRUE}.
}
\item{used.MAS5TargetValue}{
    In case of not applying batch normalization (see details), the MAS5 target value used for the expression set must be given.
}
    \item{\dots}{
    A list of additional arguments. See details for a description.
}
}
\details{
    A list of possible classifiers is obtained by \code{get.classifiers()}. By default the input dataset (i.e. eset) is assumed to contain many (n>=25) samples with corresponding probe-sets needed for classification. If one of these conditions is not met, a classifier outcome might be seriously affected. By default an error is given. Although strongly discouraged, it is possible to circumvent the security checks. If not all required probe-sets are included in the input set, you can explicitly pass the parameter \code{allow.reweighted = TRUE} to the \code{applyClassifier} function in order to determine the classifier outcome using less probe-sets (e.g. possible if the missing probe-sets are known to have minimal contribution).See \code{vignette("MissingCovariates")} for more information. If the input expression set has a small number of samples, the default batch correction becomes ineffective. If you are aware of the possible negative effects you can force to not use batch correction by passing the parameter \code{do.batchcorrection=FALSE}. 
}
\value{
    A list with object(s) of class Classifier.
}
\seealso{
\code{\link{Classifier}}, \code{\link{get.classifiers}}, \code{\link{get.info}}, \code{\link[Biobase]{ExpressionSet}},\code{vignette("MissingCovariates")}
}
\examples{
data(eset.mas5) 
##Apply the EMC92 classifier
classification = applyClassifier.mas5("EMC92",eset.mas5,used.MAS5TargetValue=500)

continuous_score = as.numeric( classification )
classification_score = as.factor( classification)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
