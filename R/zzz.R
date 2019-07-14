#' @aliases showClassifierList
#' @title Show classifier names and descriptions.
#'
#' @description \code{showClassifierList} gives a data.frame of all implemented classifiers.
#'
#' @param normalizations an optional text argument of one or more normalization methods in order to filter the classifiers to be shown.
#'
#' @return A data.frame with columns: 'name','normalizationMethod' and 'description'
#'
#' @details The names of the classifiers shown can be used as input for the \code{\link{runClassifier}} function and the \code{\link{getClassifier}} function.
#'
#' @family workflow functions
#'
#' @docType methods
#' @rdname showClassifierList-methods
#'
#' @examples
#' showClassifierList()
#' data(exampleMAS5)
#' myData <- setNormalizationMethod(exampleMAS5, 'MAS5.0',targetValue=500)
#' results <- runClassifier('UAMS70', myData)
#' getScores( results )
#' getClassifications( results )
#'
#' @export
showClassifierList <- function(normalizations) {
    registeredClassifiers <- lapply(.classifierHookList, function(x) {
        x()
    })
    registeredClassifiers <- t(simplify2array(registeredClassifiers))
    if (!missing(normalizations)) {
        normalizations <- match.arg(toupper(normalizations), getNormalizationMethods(), 
            several.ok = TRUE)
        registeredClassifiers <- registeredClassifiers[registeredClassifiers[, "normalizationMethod"] %in% 
            normalizations, , drop = FALSE]
    }
    return(registeredClassifiers)
}

#' @rdname getClassifier-methods
#' @aliases getClassifier,character-method
setMethod("getClassifier", signature = signature(value = "character"), definition = function(value) {
    classifierName <- match.arg(value, showClassifierList()[, "name"])
    idx <- which(showClassifierList()[, "name"] == classifierName)
    classifierParameters <- .classifierHookList[[idx]](infoOnly = FALSE)
    classifierParameters
})

#' @aliases setNormalizationMethod
#' @title Prepare data.
#'
#' @description \code{setNormalizationMethod} is to be called prior to running a classifier.
#'
#' @param expressionSet An object of class \code{\link{ExpressionSet}} containing the gene
#' expression data.
#' @param method A character string indicating the normalization that was applied to the data.
#' Possible values are given by \code{getNormalizationMethods()}.
#' @param ... see details.
#'
#' @return An object of class \code{\link{FixedExpressionData}}
#'
#' @details The \code{\link{FixedExpressionData}} class forms together with the
#' \code{\link{ClassifierParameters}} class the basis for input to the  \code{\link{runClassifier}}
#' function. The data inside the \code{\link{FixedExpressionData}}-class has to be stored as it is
#' right after normalization.  This function may require some additional arguments:
#' \itemize{
#' \item \code{isLog2Transformed = TRUE} Use this argument if the data already underwent a
#' log2transformation, as is common e.g. in case of MAS5.0 normalization.
#' \item \code{targetValue = value} This is a MAS5.0 specific argument. It is the sample intensity
#' mean when the lowest and highest 2\% of intensities are discarded. If only part of the original
#' expression set is given to this function, then this argument is required.
#' }
#'
#' @family workflow functions
#'
#' @docType methods
#' @rdname setNormalizationMethod-methods
#'
#' @examples
#' data(exampleMAS5)
#' myData <- setNormalizationMethod(exampleMAS5, 'MAS5.0',targetValue=500)
#' results <- runClassifier('EMC92', myData)
#' getScores( results )
#' getClassifications( results )
#'
#' @importFrom methods new
#' @export
setNormalizationMethod <- function(expressionSet, method, ...) {
    if (!inherits(expressionSet, "ExpressionSet")) {
        stop("Expect argument 'expressionSet' to be an object of type ExpressionSet")
    }
    method <- match.arg(toupper(method), getNormalizationMethods())
    dotList <- list(...)
    FixedExpressionData(normalizationMethod = method, expressionMatrix = exprs(expressionSet), 
        ...)
    
}

