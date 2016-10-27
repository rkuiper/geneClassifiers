setOldClass("package_version")

#' @name ClassifierResults
#' @title An S4 class to store classifier results.
#'
#'
#' @description  This class stores classifier results as obtained after running
#' the \code{\link{runClassifier}} function.
#'
#' @slot classifierParameters An object of class \code{\link{ClassifierParameters}}
#' in which the applied classifier parameters are stored.
#' @slot score A numeric vector of resulting classifier scores
#' @slot batchCorrection A character vector indicating wheter batch correction was applied
#' @slot weightingType A character string indicating wheter the weighting type was complete
#' (i.e. no missing data), reweighted (i.e. missing data was handled based on correction
#' using the covariance structure in the classifiers training data), or reduced (i.e. missing
#' data but not reweighting the original probeset weighting)
#' @slot .geneClassifierVersion An object of class \code{\link{package_version}}

setClass("ClassifierResults",
    representation=representation(
        classifierParameters="ClassifierParameters",
        score="numeric",
        batchCorrection="logical",
        weightingType='character',
        .geneClassifierVersion="package_version"
    )
)

#' @importFrom utils packageVersion packageName
setMethod("initialize", signature(.Object = "ClassifierResults"),
    function(.Object, ...){
        x <- list(...)
        if (is.null(x[["weightingType"]])) stop('weightingType must be given')
        if (is.null(x[["batchCorrection"]])) stop('batchCorrection must be given')
        if (is.null(x[["score"]])) stop('score must be given')
        if (!(x[["weightingType"]]%in%c("complete","reweighted"))) stop('weightingType should be one of: complete, reweighted')
        .Object@weightingType = x[["weightingType"]]
        .Object@score = x[["score"]]
        .Object@classifierParameters = x[["classifierParameters"]]
        .Object@batchCorrection = x[["batchCorrection"]]
        .Object@.geneClassifierVersion = packageVersion(packageName())
        return(.Object)
    }
)

#' @rdname getScores-methods
#' @aliases getScores,ClassifierResults-method
#' @export
setMethod("getScores",
    signature  = signature("ClassifierResults"),
    definition = function(object) {
        return(object@score)
    }
)

#' @rdname getClassifier-methods
#' @aliases getClassifier,ClassifierResults-method
#' @export
setMethod("getClassifier",
    signature  = signature(value = "ClassifierResults"),
    definition = function(value){
        return(value@classifierParameters)
    }
)

#' @rdname getName-methods
#' @aliases getName,ClassifierResults-method
#' @export
setMethod("getName",
    signature  = signature("ClassifierResults"),
    definition = function(object) {
        return(getName(getClassifier(object)))
    }
)

#' @rdname getClassifications-methods
#' @aliases getClassifications,ClassifierResults-method
#' @export
#' @importFrom utils as.roman
setMethod("getClassifications",
    signature  = signature("ClassifierResults"),
    definition = function(object) {
        boundaries<-getDecisionBoundaries( getClassifier(object) )
        classifications<-sapply(boundaries,FUN=function(x,y){y>x},y=getScores(object))
        classifications<-factor(
            rowSums(classifications),
            levels=seq(from=0, to=length(boundaries )),
            labels=paste("Risk",as.roman(1+seq(from=0,length(boundaries))),sep='-'),
            ordered=TRUE)
        return(classifications)
    }
)

setMethod("show",
    signature  = "ClassifierResults",
    definition = function(object){
        classifications<-getClassifications(object)
        riskGroup<-table(classifications)
        cat("Note: Research use only\n")
        cat("Classifier:" , getName(object),"\n")
        decisionBoundaries<-getDecisionBoundaries( getClassifier(object) )
        for (boundary in decisionBoundaries){
            cat("\t",paste(">",boundary ,":",names(riskGroup)[-1]),"\n")
        }
        print(riskGroup)
        cat("\n\tBatch corrected   :",c("no","yes")[object@batchCorrection+1],"\n")
        cat("\tweighting type      :",object@weightingType,"\n")
        cat("----------------\n")
    }
)
