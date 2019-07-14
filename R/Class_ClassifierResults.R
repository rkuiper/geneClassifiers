setOldClass("package_version")


############################################################## 
#' Obtain the weighting type used to obtain a classifier result.
#'
#' \code{getWeightingTypes} returns weigthing type 
#' 
#'
#' @return either 'complete' or 'reweighted'
#'
#' @family classifier results
#'
#' @export
#' @docType methods
#' @rdname getWeightingType-methods
#'
#' @examples
#' getWeightingTypes()

getWeightingTypes <- function() {
    return(c("complete", "reweighted"))
}

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

setClass("ClassifierResults", representation = representation(classifierParameters = "ClassifierParameters", 
    score = "numeric", batchCorrection = "logical", weightingType = "character", 
    .geneClassifierVersion = "package_version"), prototype = list(.geneClassifierVersion = packageVersionInternal()))

#' @importFrom methods validObject
setValidity("ClassifierResults", function(object) {
    errTxt <- vector()
    if (length(getWeightingType(object)) != 1) {
        errTxt <- c(errTxt, "More than one weightingTypes received")
    } else if (!getWeightingType(object) %in% getWeightingTypes()) {
        errTxt <- c(errTxt, paste0("Weigthingtypes must be one of: ", paste(getWeightingTypes(), 
            collapse = ", ")))
    }
    validObject(getClassifier(object), complete = TRUE)
    
    if (length(errTxt) > 0) {
        return(errTxt)
    }
    return(TRUE)
})


setMethod("ClassifierResults", signature = signature(weightingType = "character", 
    batchCorrection = "logical", score = "numeric", classifierParameters = "ClassifierParameters"), 
    definition = function(weightingType, batchCorrection, score, classifierParameters) {
        new("ClassifierResults", classifierParameters = classifierParameters, score = score, 
            batchCorrection = batchCorrection, weightingType = weightingType)
    })

#' @rdname getScores-methods
#' @aliases getScores,ClassifierResults-method
#' @export
setMethod("getScores", signature = signature("ClassifierResults"), definition = function(object) {
    return(object@score)
})

#' @rdname getClassifier-methods
#' @aliases getClassifier,ClassifierResults-method
#' @export
setMethod("getClassifier", signature = signature(value = "ClassifierResults"), definition = function(value) {
    return(value@classifierParameters)
})


#' @rdname getBatchCorrection-methods
#' @aliases getBatchCorrection,ClassifierResults-method
#' @export
setMethod("getBatchCorrection", signature = signature(object = "ClassifierResults"), 
    definition = function(object) {
        return(object@batchCorrection)
    })


#' @rdname getWeightingType-methods
#' @aliases getWeightingType,ClassifierResults-method
#' @export
setMethod("getWeightingType", signature = signature(object = "ClassifierResults"), 
    definition = function(object) {
        return(object@weightingType)
    })

#' @rdname getName-methods
#' @aliases getName,ClassifierResults-method
#' @export
setMethod("getName", signature = signature(object = "ClassifierResults"), definition = function(object) {
    return(getName(getClassifier(object)))
})

#' @rdname getClassifications-methods
#' @aliases getClassifications,ClassifierResults-method
#' @export
#' @importFrom utils as.roman
setMethod("getClassifications", signature = signature(object = "ClassifierResults"), 
    definition = function(object) {
        boundaries <- getDecisionBoundaries(getClassifier(object))
        y <- getScores(object)
        classifications <- vapply(X = boundaries, FUN = function(x, y) {
            y > x
        }, FUN.VALUE = rep(FALSE, length(y)), y = y)
        classifications <- factor(rowSums(classifications), levels = seq(from = 0, 
            to = length(boundaries)), labels = paste("Risk", as.roman(1 + seq(from = 0, 
            length(boundaries))), sep = "-"), ordered = TRUE)
        return(classifications)
    })

setMethod("show", signature = signature(object = "ClassifierResults"), definition = function(object) {
    classifications <- getClassifications(object)
    riskGroup <- table(classifications)
    cat("Note: Research use only\n")
    cat("Classifier:", getName(object), "\n")
    decisionBoundaries <- getDecisionBoundaries(getClassifier(object))
    for (boundary in decisionBoundaries) {
        cat("\t", paste(">", boundary, ":", names(riskGroup)[-1]), "\n")
    }
    print(riskGroup)
    cat("\n\tBatch corrected   :", c("no", "yes")[getBatchCorrection(object) + 1], 
        "\n")
    cat("\tweighting type      :", getWeightingType(object), "\n")
    cat("----------------\n")
})
