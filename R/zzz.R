

.onAttach<-function(libname, pkgname) {
    packageStartupMessage('See vignette("geneClassifiers") for help')
    checkClassifierHooks()
}
.onLoad<-function(libname, pkgname) {}

#' @aliases showClassifierList
#' @title Show classifier names and descriptions.
#'
#' @description \code{showClassifierList} gives a data.frame of all implemented classifiers.
#'
#' @param normalizations an optional text argument of one or more normalization methods in order to filter the classifiers to be shown.
#'
#' @return A data.frame with columns: "name","normalizationMethod" and "description"
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
#' myData <- setNormalizationMethod(exampleMAS5, "MAS5.0",targetValue=500)
#' results <- runClassifier('UAMS70', myData)
#' getScores( results )
#' getClassifications( results )
#'
#' @export
showClassifierList<-function(normalizations){
    registeredClassifiers<-lapply(.classifierHookList,function(x){x()})
    registeredClassifiers<-t(simplify2array(registeredClassifiers))
    if (!missing(normalizations)){
        normalizations<-match.arg(toupper(normalizations),getNormalizationMethods(),several.ok=TRUE)
        registeredClassifiers<-registeredClassifiers[registeredClassifiers[,"normalizationMethod"]%in%normalizations,,drop=FALSE]
    }
    return(registeredClassifiers)
}

#' @rdname getClassifier-methods
#' @aliases getClassifier,character-method
setMethod("getClassifier",
    signature  = signature(value = "character"),
    definition = function(value){
    classifierName<-match.arg(value,showClassifierList()[,"name"])
    idx<-which(showClassifierList()[,"name"]==classifierName)
    classifierParameters<- .classifierHookList[[idx]](infoOnly=FALSE)
    classifierParameters
    }
)


doCheckEventChains<-function(eventChain,normalizationMethod){
    normalizationMethod<-match.arg(normalizationMethod,getNormalizationMethods())
    allowedKeyWords<-c("targetValue","truncate","to.log","allow.reweighted","to.meancentering","to.unitvariance","to.referencemeanvar")


    if (!all(names(eventChain)%in%allowedKeyWords)){
        diff<-setdiff(names(eventChain),allowedKeyWords)
        stop("Eventchain structure error: Unrecognized keyword(s) in eventChain: ",paste(diff,collapse=","))
    }

    #$ 'targetValue' (Must be given for MAS5.0 classifiers only. Only used for mas5 classifiers. Indicating the required targetValue)
    if (normalizationMethod=="MAS5.0") {
        if (!"targetValue"%in%names(eventChain)){stop("Eventchain structure error: keyword 'targetValue' must be given for a mas5 classifier.")}
        if (!is.numeric(eventChain[["targetValue"]])){stop("Eventchain structure error: keyword 'targetValue' must have a numeric value.")}
        if (eventChain[["targetValue"]]<1){stop("Eventchain structure error: keyvalue 'targetValue' must be a positive real number.")}
    } else if ("targetValue"%in%names(eventChain)){stop("Eventchain structure error: keyword 'targetValue' can be given for a mas5 classifiers only.")
    }
    #$ 'truncate' (Must be given. Trunctates values below a given value (applies to untransformed/unlogged data.)


    if (!"truncate"%in%names(eventChain)){stop("Eventchain structure error: keyword 'truncate' must be given.")}
    if (!is.numeric(eventChain[["truncate"]])){stop("Eventchain structure error: keyword 'truncate' must have a numeric value.")}

    #$ 'to.log' (Although probably only logical to apply in MAS5.0 normalized data, always applicable)
        if (!is.null(eventChain[["to.log"]])){
        if (!is.numeric(eventChain[["to.log"]])){stop("Eventchain structure error: keyword 'to.log' must be a numerical value representing the base number for the log.")}
    }


    #$ 'allow.reweighted' Must be given (should reweigting be allowed in case of missing probe-set(s)?)
    if ("allow.reweighted"%in%names(eventChain)){
        if (!is.logical(eventChain[["allow.reweighted"]])){stop("Eventchain structure error: keyword 'allow.reweighted' must be a logical TRUE or FALSE.")}
    } else {
        stop("Eventchain structure error: keyword 'allow.reweighted' must be given.")
    }

    #$ 'to.meancentering' (per probe-set based centering to mean=0)
    if ("to.meancentering"%in%names(eventChain)){
        if (!is.logical(eventChain[["to.meancentering"]])){stop("Eventchain structure error: keyword 'to.meancentering' must be a logical TRUE or FALSE.")}
    }

    #$ 'to.unitvariance' (per probe-set based centering to mean=0)
    if ("to.unitvariance"%in%names(eventChain)){
        if (!is.logical(eventChain[["to.unitvariance"]])){stop("Eventchain structure error: keyword 'to.unitvariance' must be a logical TRUE or FALSE.")}
    }

    #$ 'to.referencemeanvar' (set probe-set means and variances to equal values as observed in training set)
    if ("to.referencemeanvar"%in%names(eventChain)){
        if (!is.logical(eventChain[["to.referencemeanvar"]])){stop("Eventchain structure error: keyword 'to.referencemeanvar' must be a logical TRUE or FALSE.")}
    }
}

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
#' myData <- setNormalizationMethod(exampleMAS5, "MAS5.0",targetValue=500)
#' results <- runClassifier('EMC92', myData)
#' getScores( results )
#' getClassifications( results )
#'
#' @importFrom methods new
#' @export
setNormalizationMethod<-function(expressionSet, method, ...){
    if (!inherits(expressionSet,"ExpressionSet")) {
        stop("Expect argument 'expressionSet' to be an object of type ExpressionSet")
    }
    method <- match.arg( toupper(method), getNormalizationMethods() )
    dotList <- list(...)
    return(new("FixedExpressionData",
        normalizationMethod=method,
        expressionMatrix=exprs(expressionSet),
        ...)
    )
}

checkClassifierHooks<-function(){
    funcDefCorrect<-unlist(lapply(.classifierHookList,function(x){
        fml<-formals(x)
        return(all(names(fml)==c("infoOnly","...")) & all(unlist(fml)==c(TRUE,"")))
    }))
    if (any(funcDefCorrect==FALSE)){
        stop("Found ",sum(funcDefCorrect),
        " classifierHooks that did not have the correct function defintion")
    }
    registeredClassifiers<-lapply(.classifierHookList,function(x){x()})
    argumentReturnOrder<-unlist(lapply(registeredClassifiers,function(x){
        all(names(x)==c("name","normalizationMethod","description"))
    }))
    namesCorrect<-unlist(lapply(registeredClassifiers,function(x){
        is.character(x["name"]) & nchar(x["name"])>2
    }))
    normalizationMethodCorrect<-unlist(lapply(registeredClassifiers,function(x){
        x["normalizationMethod"]%in%getNormalizationMethods()
    }))
    namesDuplicated<-duplicated(unlist(lapply(registeredClassifiers,function(x){
        x["name"]
    })))
    if (any(namesDuplicated==TRUE)){
        notOK<-registeredClassifiers[which(namesDuplicated==TRUE),drop=FALSE]
        txt<-paste(unlist(lapply(notOK,function(x){x["name"]})),collapse=",")
        stop("The following classifierHooks were found to have names that occured more than once: ",txt,"\n")
    }
    if (any(argumentReturnOrder==FALSE)){
        notOK<-registeredClassifiers[which(argumentReturnOrder==FALSE)]
        txt<-paste(unlist(lapply(notOK,function(x){x["name"]})),collapse=",")
        stop("The following classifierHooks did not have a correct argument return order: ",txt,"\n")
    }
    if (any(namesCorrect==FALSE)){
        notOK<-registeredClassifiers[which(namesCorrect==FALSE)]
        txt<-paste(unlist(lapply(notOK,function(x){x["name"]})),collapse=",")
        stop("Some classifierHooks did not have a correct classifier name (i.e. >2 characters)\n")
    }
    if (any(normalizationMethodCorrect==FALSE)){
        notOK<-registeredClassifiers[which(normalizationMethodCorrect==FALSE)]
        txt<-paste(unlist(lapply(notOK,function(x){x["name"]})),collapse=",")
        stop(
            "The following classifierHooks did not provided a known normalization method (i.e. one of ",
            paste(getNormalizationMethods(),collapse=","),"): ",txt,"\n"
        )
    }
    registeredClassifiers<-lapply(.classifierHookList,function(x){x(infoOnly=FALSE)})
}
