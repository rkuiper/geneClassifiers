##############################################################
#' Obtain classifier score.
#'
#' \code{getScores} returns the resulting scores from a classifier
#' run
#'
#' @param object An object of class \code{\link{ClassifierResults}}
#'
#' @return A numeric vector with scores per sample
#'
#' @family classifier results
#'
#' @export
#' @docType methods
#' @rdname getScores-methods
#'
#' @examples
#' data(exampleMAS5)
#' myData <- setNormalizationMethod(exampleMAS5, "MAS5.0", targetValue=500)
#' results <- runClassifier('EMC92', myData)
#' getScores( results )
#' getClassifications( results )
setGeneric("getScores",function(object){standardGeneric ("getScores")})

##############################################################
#' Obtain classifier classifications.
#'
#' \code{getClassifications} returns the resulting classifications.
#'
#' @param object An object of class \code{\link{ClassifierResults}}
#'
#' @return A vector of orderd factors with classifications per sample
#'
#' @family classifier results
#'
#' @export
#' @docType methods
#' @rdname getClassifications-methods
#'
#' @examples
#' data(exampleMAS5)
#' myData <- setNormalizationMethod(exampleMAS5, "MAS5.0",targetValue=500)
#' results <- runClassifier('EMC92', myData)
#' getScores( results )
#' getClassifications( results )
setGeneric("getClassifications",function(object){standardGeneric ("getClassifications")})

setGeneric("runProcess",function(object,expressionMat){standardGeneric ("runProcess")})

setGeneric("getValues",function(object){standardGeneric ("getValues")})

##############################################################
#' Obtain object names.
#'
#' \code{getName} returns the name associated with the requested object.
#'
#' @param object The object to get the name of.
#'
#' @return The return value is a character string
#'
#' @seealso \code{\link{ClassifierParameters}}
#' @seealso \code{\link{ClassifierResults}}
#'
#' @export
#' @docType methods
#' @rdname getName-methods
#'
#' @examples
#' aClassifier <- getClassifier("EMC92")
#' getName( aClassifier )
setGeneric("getName",function(object){standardGeneric ("getName")})

setGeneric("explicitlyChangeExprs<-",function(object,value){standardGeneric ("explicitlyChangeExprs<-")})
setGeneric("getExpressionEnvironment",function(object,value){standardGeneric ("getExpressionEnvironment")})

##############################################################
#' Obtain normalization method
#'
#' \code{getNormalizationMethod} returns the normalization method
#' associated with the object
#'
#' @param object An object of class \code{\link{FixedExpressionData}} or
#' \code{\link{ClassifierParameters}}
#'
#' @return A character string indicating the normalization method.
#'
#' @seealso \code{\link{getNormalizationMethods}}
#' @family classifier information functions
#' @family fixed data information extraction functions
#'
#' @export
#' @docType methods
#' @rdname getNormalizationMethod-methods
#'
#' @examples
#' data(exampleMAS5)
#' myData <- setNormalizationMethod(exampleMAS5, "MAS5.0", targetValue=500)
#' aClassifier <- getClassifier("EMC92")
#' getNormalizationMethod( myData )
#' getNormalizationMethod( aClassifier )
setGeneric("getNormalizationMethod",function(object){standardGeneric ("getNormalizationMethod")})

##############################################################
#' Obtain the targetValue
#'
#' \code{getTargetValue} returns the current applied targetValue
#' in the MAS5.0 gene expression data.
#'
#' @param object An object of class \code{\link{FixedExpressionData}}
#'
#' @return A numeric value
#'
#' @family fixed data information extraction functions
#'
#' @export
#' @docType methods
#' @rdname getTargetValue-methods
#'
#' @examples
#' data(exampleMAS5)
#' myData <- setNormalizationMethod(exampleMAS5, "MAS5.0", targetValue=500)
#' getTargetValue( myData )

setGeneric("getTargetValue",function(object){standardGeneric ("getTargetValue")})
setGeneric("setTargetValue<-",function(object,value){standardGeneric ("setTargetValue<-")})

setGeneric("rawExprs",function(object){standardGeneric ("rawExprs")})

setGeneric("addTransformationProcess",function(object,name,values,...){standardGeneric ("addTransformationProcess")})
setGeneric("getTransformationProcesses",function(object){standardGeneric ("getTransformationProcesses")})
setGeneric("removeTransformationProcesses",function(object,n){standardGeneric ("removeTransformationProcesses")})

setGeneric("getNormalizationParameters",function(object){standardGeneric ("getNormalizationParameters")})

##############################################################
#' Obtain a classifier definition.
#'
#' \code{getClassifier} returns a requested classifier definition.
#'
#' @param value Either a text value indicating a classifier name
#' (see \code{\link{showClassifierList}}), or an object of
#' class \code{\link{ClassifierResults}} as returned by the
#' \code{\link{runClassifier}} function.
#'
#' @return The return value is a classifier definition which
#' is encoded in an object of class \code{\link{ClassifierParameters}}.
#' This can be used as input argument for the \code{\link{runClassifier}}
#' function.
#'
#' @family classifier information functions
#' @seealso \code{\link{ClassifierParameters}} and \code{\link{runClassifier}}
#'
#' @export
#' @docType methods
#' @rdname getClassifier-methods
#'
#' @examples
#' getClassifier("EMC92")
setGeneric("getClassifier",function(value){standardGeneric ("getClassifier")})

##############################################################
#' Obtain probe-set names.
#'
#' \code{getProbeNames} returns the probe names associated with
#' the requested classifier.
#'
#' @param object An object of class \code{\link{ClassifierParameters}}
#' as returned by \code{\link{getClassifier}}
#'
#' @return The return value is a character vector of probe-set names.
#'
#' @family classifier information functions
#' @seealso \code{\link{ClassifierParameters}}
#'
#' @export
#' @docType methods
#' @rdname getProbeNames-methods
#'
#' @examples
#' aClassifier <- getClassifier("EMC92")
#' getProbeNames( aClassifier )
setGeneric("getProbeNames",function(object){standardGeneric ("getProbeNames")})

##############################################################
#' Obtain classifier weights.
#'
#' \code{getWeights} returns the probe weights associated with
#' the classifier.
#'
#' @param object An object of class \code{\link{ClassifierParameters}}
#' as returned by \code{\link{getClassifier}}
#'
#' @return A numeric vector.
#'
#' @family classifier information functions
#'
#' @export
#' @docType methods
#' @rdname getWeights-methods
#'
#' @examples
#' aClassifier <- getClassifier("EMC92")
#' getWeights(aClassifier)
setGeneric("getWeights",function(object){standardGeneric ("getWeights")})

##############################################################
#' Obtain classifier training data.
#'
#' \code{getTrainingData} returns the training data that was used
#' for building the classifier.
#'
#' @param object An object of class \code{\link{ClassifierParameters}}
#' as returned by \code{\link{getClassifier}}
#'
#' @return An object of class \code{\link{ExpressionSet}}
#'
#' @family classifier information functions
#'
#' @export
#' @docType methods
#' @rdname getTrainingData-methods
#'
#' @examples
#' aClassifier <- getClassifier("EMC92")
#' getTrainingData(aClassifier)
setGeneric("getTrainingData",function(object){standardGeneric ("getTrainingData")})

##############################################################
#' Obtain citations to the classifier
#'
#' \code{getCitations} Obtain citations to the classifier
#'
#' @param object An object of class \code{\link{ClassifierParameters}}
#' as returned by \code{\link{getClassifier}}
#'
#' @return A character vector
#'
#' @family classifier information functions
#'
#' @export
#' @docType methods
#' @rdname getCitations-methods
#'
#' @examples
#' aClassifier <- getClassifier("EMC92")
#' getCitations(aClassifier)
setGeneric("getCitations",function(object){standardGeneric ("getCitations")})


##############################################################
#' Obtain the decision boundaries defined for the classifier.
#'
#' \code{getDecisionBoundaries} returns the a numeric vector
#'  of boundary values that separate the risk groups.
#'
#' @param object An object of class \code{\link{ClassifierParameters}}
#' as returned by \code{\link{getClassifier}}
#'
#' @return A numeric vector
#'
#' @family classifier information functions
#'
#' @export
#' @docType methods
#' @rdname getDecisionBoundaries-methods
#'
#' @examples
#' aClassifier <- getClassifier("EMC92")
#' getDecisionBoundaries(aClassifier)
setGeneric("getDecisionBoundaries",function(object ){standardGeneric ("getDecisionBoundaries")})


##############################################################
#' Obtain classifiers' intercept.
#'
#' \code{getIntercept} returns the numeric value of the
#' classifier's intercept.
#'
#' @param object An object of class \code{\link{ClassifierParameters}}
#' as returned by \code{\link{getClassifier}}
#'
#' @return A numeric value
#'
#' @family classifier information functions
#'
#' @export
#' @docType methods
#' @rdname getIntercept-methods
#'
#' @examples
#' aClassifier <- getClassifier("EMC92")
#' getIntercept(aClassifier)
setGeneric("getIntercept",function(object){standardGeneric ("getIntercept")})

##############################################################
#' Obtain classifiers' description.
#'
#' \code{getDescription} returns the descriptive text associated
#' with the classifier.
#'
#' @param object An object of class \code{\link{ClassifierParameters}}
#' as returned by \code{\link{getClassifier}}
#'
#' @return A character string describing the classifier
#'
#' @family classifier information functions
#'
#' @export
#' @docType methods
#' @rdname getDescription-methods
#'
#' @examples
#' aClassifier <- getClassifier("EMC92")
#' getDescription(aClassifier)
setGeneric("getDescription",function(object ){standardGeneric ("getDescription")})

##############################################################
#' Perform classification.
#'
#' \code{runClassifier} performs classification by applying a
#' classifier to gene expression data.
#'
#' @param classifierParameters Either a text value indicating a
#' classifier name (see \code{\link{showClassifierList}}), or an
#' object of class \link{ClassifierParameters} as returned by the
#' \code{\link{getClassifier}} function.
#'
#' @param fixedExpressionData The data to be classified in the
#' form of a \code{\link{FixedExpressionData}} object as returned
#' by the \code{\link{setNormalizationMethod}} function.
#'
#' @param ... see details
#'
#' @return The classification results as an object of class
#' \code{\link{ClassifierResults}}.
#'
#' @details A list of possible classifiers is obtained by
#' \code{\link{showClassifierList}}. The data to be classified
#' is first to be processed by the \link{setNormalizationMethod}
#' function. By default the data is assumed to contain many (n>=25)
#' samples with corresponding probe-sets needed for classification.
#' If one of these conditions is not met, a classifier outcome might
#' be seriously affected. By default an error is given. Although
#' strongly discouraged, it is possible to circumvent the security
#' checks. If not all required probe-sets are included in the input
#' set, you can explicitly pass the parameter \code{allow.reweighted = TRUE}
#' to the \code{runClassifier} function in order to determine the
#' classifier outcome using less probe-sets (e.g. possible if the
#' missing probe-sets are known to have minimal contribution).See
#' \code{vignette("MissingCovariates")} for more information. If
#' the input data has a small number of samples, the default batch
#' correction becomes ineffective. If you are aware of the possible
#' negative effects you can force to not use batch correction by
#' passing the parameter \code{do.batchcorrection=FALSE}.
#'
#' @family workflow functions
#'
#' @export
#' @docType methods
#' @rdname runClassifier-methods
#'
#' @examples
#' data(exampleMAS5)
#' myData<-setNormalizationMethod(exampleMAS5,"MAS5.0",targetValue=500)
#' runClassifier("EMC92",myData)
setGeneric("runClassifier",function(classifierParameters, fixedExpressionData,... ){standardGeneric ("runClassifier")})

setGeneric("getMeans",function(object ){standardGeneric ("getMeans")})
setGeneric("getSds",function(object ){standardGeneric ("getSds")})

##############################################################
#' Obtain classifiers' event chain.
#'
#' \code{getEventChain} returns the event chain encoded in the
#' in the classifier. The eventchain indicates what preprocessing
#' steps are performed by the \code{\link{runClassifier}} function
#' prior to classification.
#'
#' @param object An object of class \code{\link{ClassifierParameters}}
#' as returned by \code{\link{getClassifier}}
#'
#' @return Returns the event chain encoded in the
#' in the classifier encoded as a named list.
#'
#'
#' @family classifier information functions
#' @seealso showClassifierList getClassifier runClassifier
#'
#' @export
#' @docType methods
#' @rdname getEventChain-methods
#'
#' @examples
#' aClassifier <- getClassifier("EMC92")
#' getEventChain(aClassifier)
setGeneric("getEventChain",function(object ){standardGeneric ("getEventChain")})

setGeneric("reWeightClassifier",function(object,probeNames, intercept,weights){standardGeneric ("reWeightClassifier")})

.classifierHookList<-list()

