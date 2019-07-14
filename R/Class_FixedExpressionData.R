#' @aliases getNormalizationMethods
#' @title Show currenly implemented normalization methods.
#'
#' @description \code{getNormalizationMethods} returns a character vector of 
#' currenlty available normalization methods.
#'
#'
#' @details The given normlization methods can be used in the
#'
#' @family workflow functions
#'
#' @docType methods
#' @rdname getNormalizationMethod-methods
#'
#' @examples
#' data(exampleMAS5)
#'
#' showClassifierList()
#' getNormalizationMethods()
#'
#' myData <- setNormalizationMethod(exampleMAS5, 'MAS5.0',targetValue=500)
#' results <- runClassifier('UAMS70', myData)
#'
#' getScores( results )
#' getClassifications( results )
#'
#' @export
getNormalizationMethods <- function() {
    ## By convention uppercase
    return(c("MAS5.0", "GCRMA"))
}

#'Example MAS5.0 ExpressionSet
#'
#'An \code{\link{ExpressionSet}}. The data contains a sample of gene expression 
#'data from patients included in the HOVON65/GMMG-HD4 trial on multiple myeloma. 
#'The data was MAS5.0 normalized to a target value of 500.
#'
"exampleMAS5"

setOldClass("package_version")

#' @name FixedExpressionData
#' @title An S4 class to store classifier parameters.
#'
#' @description  This class stores gene expression data together with 
#' information on the normalization method and additional normalization related 
#' parameters. In order to ensure the data is not manipulated in unforeseen 
#' ways, manipulation is strictly controled through adding transformations 
#' which are predefined in the \code{TransformationProcess}-class. Upon reading
#' the data by the \code{exprs} function, the transformations areperformed in
#' the order the were added.
#'
#' @slot normalizationMethod A character string indicating the normalization
#' method that was applied to the data. Possible values are give by 
#' \code{\link{getNormalizationMethods}}.
#' @slot expressionEnvironment A locked environment in which the expression
#' matrix is stored.
#' @slot normalizationParameters A list with normalization specific values.
#' @slot transformationProcess A locked environment to which the transformation
#' processes are added.
#' @slot .geneClassifierVersion An object of class \code{\link{package_version}}
setClass("FixedExpressionData", representation = representation(normalizationMethod = "character", 
    expressionEnvironment = "environment", normalizationParameters = "list", transformationProcess = "environment", 
    .geneClassifierVersion = "package_version"), prototype = list(expressionEnvironment = as.environment(list(expressionMatrix = matrix(numeric(0), 
    ncol = 0, nrow = 0))), transformationProcess = as.environment(list(processes = list())), 
    .geneClassifierVersion = packageVersionInternal()))
setValidity("FixedExpressionData", function(object) {
    
    errTxt <- vector()
    if (!any(getNormalizationMethod(object) %in% getNormalizationMethods()) | length(getNormalizationMethod(object)) != 
        1) {
        errTxt <- c(errTxt, "Unknown normalization method")
    }
    if (!is.list(getTransformationProcesses(object))) {
        errTxt <- c(errTxt, "Unknown transformationProcess type")
    }
    if ("MAS5.0" %in% getNormalizationMethod(object)) {
        if (is.null(getNormalizationParameters(object)[["targetValue"]])) {
            errTxt <- c(errTxt, "targetValue must be given for MAS5.0")
        }
    }
    
    if (!is.matrix(rawExprs(object))) {
        errTxt <- c(errTxt, "'expressionMatrix' must be a numeric matrix.")
    }
    if (mode(rawExprs(object)) != "numeric") {
        errTxt <- c(errTxt, "'expressionMatrix' must be numeric.")
    }
    if (any(is.na(rawExprs(object)))) {
        errTxt <- c(errTxt, "NA values found in expression matrix.")
    }
    if (nrow(object) > 0) {
        if (is.null(rownames(object))) {
            errTxt <- c(errTxt, "'expressionMatrix' must have row names.")
        } else if (sum(duplicated(rownames(object))) > 0) {
            errTxt <- c(errTxt, "duplicated row names found.")
        }
    }
    if (ncol(object) > 0) {
        if (is.null(colnames(object))) {
            errTxt <- c(errTxt, "'expressionMatrix' must have column names.")
        } else if (sum(duplicated(colnames(object))) > 0) {
            errTxt <- c(errTxt, "duplicated column names found.")
        }
    }
    
    if (length(errTxt) > 0) {
        return(errTxt)
    }
    return(TRUE)
})

#' @importFrom utils packageVersion packageName
#' @importFrom stats quantile
setMethod("FixedExpressionData", signature = signature(normalizationMethod = "character", 
    expressionMatrix = "matrix"), definition = function(normalizationMethod, expressionMatrix, 
    ...) {
    isLog2Transformed <- FALSE  ##Assumption to be checked below
    targetValue <- NULL
    normalizationParameters <- list()
    
    ## process argument list
    dotList <- list(...)
    arg.isLog2Transformed <- NULL
    arg.targetValue <- NULL
    if (!is.null(dotList[["isLog2Transformed"]])) {
        if (!is.logical(dotList[["isLog2Transformed"]])) {
            stop("'isLog2Transformed' argument must be either TRUE of FALSE")
        }
        arg.isLog2Transformed <- dotList[["isLog2Transformed"]]
    }
    if (!is.null(dotList[["targetValue"]])) {
        if (!is.numeric(dotList[["targetValue"]])) {
            stop("'targetValue' argument must be numeric")
        }
        if (dotList[["targetValue"]] <= 0) {
            stop("'targetValue' argument must >0")
        }
        arg.targetValue <- dotList[["targetValue"]]
    }
    
    if (mode(expressionMatrix) != "numeric") {
        stop("'expressionMatrix' argument must be a numeric matrix.")
    }
    if (any(is.na(expressionMatrix))) {
        stop("NA values found in expression matrix.")
    }
    normalizationMethod <- match.arg(normalizationMethod, getNormalizationMethods())
    .getTargetValue <- function(expressionMatrix) {
        ## expressionMatrix must be not log2 transformed!
        e <- t(expressionMatrix)
        q <- apply(e, 1, quantile, c(0.02, 0.98), type = 2)
        naMat <- (e >= q[1, ] & e <= q[2, ])
        naMat[!naMat] <- NA
        targetValue <- unique(signif(rowMeans(naMat * e, na.rm = TRUE), 6))
        ## differs between samples -> return NULL
        if (length(targetValue) != 1) {
            targetValue <- NULL
        }
        return(targetValue)
    }
    .isAlreadylog2Transformed.check <- function(expressionMatrix) {
        storedSeed <- .Random.seed
        set.seed(1)
        aSample <- sample(expressionMatrix, 10000, replace = TRUE)
        .Random.seed <- storedSeed
        isLog2Transformed <- as.logical(quantile(aSample, 0.75, na.rm = TRUE) < 30 | 
            any(expressionMatrix < 0))
        return(isLog2Transformed)
    }
    
    ## Normalization method specific code:
    if (normalizationMethod == "MAS5.0") {
        isLog2Estimate <- .isAlreadylog2Transformed.check(expressionMatrix)
        ## if arg.isLog2Transformed given, warn if discordant with data
        if (!is.null(arg.isLog2Transformed)) {
            isLog2Transformed <- arg.isLog2Transformed
            if (arg.isLog2Transformed != isLog2Estimate) {
                warning("Reconsider argument isLog2Transformed = ", arg.isLog2Transformed)
            }
        } else {
            ## Assume estimate is correct and warn if estimate is TRUE
            isLog2Transformed <- isLog2Estimate
            if (isLog2Estimate == TRUE) {
                warning("The data seems to be log2 transformed. If this is incorrect, please provide the logical argument 'isLog2Transformed=FALSE'.")
            }
        }
        if (isLog2Transformed) 
            expressionMatrix[] <- 2^expressionMatrix
        
        calculatedTargetValue <- .getTargetValue(expressionMatrix)
        if (!is.null(arg.targetValue)) {
            targetValue <- arg.targetValue
            if (!is.null(calculatedTargetValue)) {
                if (arg.targetValue != calculatedTargetValue) {
                  warning("Using the given targetValue, but it does not correspond to the observed targetValue in the data of ", 
                    calculatedTargetValue, ".")
                }
            }
        } else if (is.null(arg.targetValue) & is.null(calculatedTargetValue)) {
            stop("Cannot determine the targetValue that was used for MAS5.0. Please provide the right 'targetValue=' as an argument.")
        } else if (is.null(arg.targetValue) & !is.null(calculatedTargetValue)) {
            targetValue <- calculatedTargetValue
        }
        normalizationParameters[["targetValue"]] <- targetValue
    }
    if (normalizationMethod == "GCRMA") {
    }
    object <- new("FixedExpressionData", normalizationMethod = normalizationMethod, 
        normalizationParameters = normalizationParameters)
    explicitlyChangeExprs(object) <- expressionMatrix
    return(object)
    
})

setMethod("show", signature = "FixedExpressionData", definition = function(object) {
    cat("____________________\n")
    cat("Fixed expression set\n\n")
    cat("Normalization method: ", getNormalizationMethod(object), "\n")
    cat("Number of samples   : ", dim(object)[2], "\n")
    cat("Number of features  : ", dim(object)[1], "\n")
    tp <- getTransformationProcesses(object)
    if (length(tp) > 0) {
        tpNames <- lapply(getTransformationProcesses(object), getName)
        tpNames <- paste(unlist(tpNames), collapse = "; ")
        cat("Applied transformation processes: ", tpNames, "\n")
    }
})

#' @importFrom methods new
#' @importFrom Biobase copyEnv
setMethod("addTransformationProcess", signature = signature(object = "FixedExpressionData", 
    name = "character"), definition = function(object, name, values, ...) {
    process <- TransformationProcess(name = name, values = values)
    aData <- copyEnv(object@transformationProcess)
    curLen <- length(aData[["processes"]])
    aData[["processes"]][[curLen + 1]] <- process
    lockEnvironment(aData, bindings = TRUE)
    object@transformationProcess <- aData
    return(object)
})

setMethod("getTransformationProcesses", signature = signature(object = "FixedExpressionData"), 
    definition = function(object) {
        return(object@transformationProcess[["processes"]])
    })

#' @importFrom Biobase copyEnv
setMethod("removeTransformationProcesses", signature = signature(object = "FixedExpressionData"), 
    definition = function(object, n = 1) {
        aData <- copyEnv(object@transformationProcess)
        curLen <- length(aData[["processes"]])
        n <- min(n, curLen)
        aData[["processes"]] <- aData[["processes"]][seq_len(max(0, curLen - n))]
        lockEnvironment(aData, bindings = TRUE)
        object@transformationProcess <- aData
        return(object)
    })

#' @rdname getNormalizationMethod-methods
#' @aliases getNormalizationMethod,FixedExpressionData-method
#' @export
setMethod("getNormalizationMethod", signature = signature("FixedExpressionData"), 
    definition = function(object) {
        object@normalizationMethod
    })

setMethod("getNormalizationParameters", signature = signature("FixedExpressionData"), 
    definition = function(object) {
        object@normalizationParameters
    })

#' @rdname getTargetValue-methods
#' @aliases getTargetValue,FixedExpressionData-method
#' @export
setMethod("getTargetValue", signature = signature("FixedExpressionData"), definition = function(object) {
    return(getNormalizationParameters(object)[["targetValue"]])
})

#' @importFrom methods validObject
setReplaceMethod("setTargetValue", signature = signature(object = "FixedExpressionData", 
    value = "numeric"), definition = function(object, value) {
    ## Replaces the rawexpression data with the new targetvalue
    if (!getNormalizationMethod(object) %in% c("MAS5.0")) {
        stop("Function 'setTargetValue' does only apply to MAS5.0 normalized data")
    }
    currentTargetValue <- getTargetValue(object)
    e <- rawExprs(object)
    explicitlyChangeExprs(object) <- e * value/currentTargetValue
    object@normalizationParameters[["targetValue"]] <- value
    validObject(object)
    return(object)
})

setMethod("getExpressionEnvironment", signature = signature("FixedExpressionData"), 
    definition = function(object) {
        return(object@expressionEnvironment)
    })

#' @importMethodsFrom Biobase exprs
setMethod("exprs", signature = signature("FixedExpressionData"), definition = function(object) {
    e <- rawExprs(object)
    processes <- getTransformationProcesses(object)
    for (process in processes) {
        e <- runProcess(process, e)
    }
    return(e)
})

setMethod("rawExprs", signature = signature("FixedExpressionData"), definition = function(object) {
    e <- getExpressionEnvironment(object)[["expressionMatrix"]]
    return(e)
})

############################################################## 
#' Dimensions of an Object
#'
#' Retrieve the dimension of an object.
#'
#' @param x an R object, for example a matrix, array or data frame.
#'
#' @return Retrieves the 'dim attribute of the object.  It is 'NULL' or a
#'    vector of mode 'integer'.
#'
#' @family fixed data information extraction functions
#'
#' @export
#' @docType methods
#' @rdname dim-methods
#' @aliases dim,FixedExpressionData-method
#' @export
#'
#' @examples
#' data(exampleMAS5)
#' myData <- setNormalizationMethod(exampleMAS5, 'MAS5.0', targetValue=500)
#' dim(myData)
#' dim(myData[1:10,1:3])
setMethod("dim", signature = signature("FixedExpressionData"), definition = function(x) {
    return(dim(rawExprs(x)))
})

#' @importMethodsFrom BiocGenerics rownames
setMethod("rownames", signature = signature(x = "FixedExpressionData"), definition = function(x) {
    return(rownames(rawExprs(x)))
})

#' @importMethodsFrom BiocGenerics colnames
setMethod("colnames", signature = signature(x = "FixedExpressionData"), definition = function(x) {
    return(colnames(rawExprs(x)))
})

############################################################## 
#' Extract
#'
#' Extract Parts of an Object
#'
#' @param x An object of class \code{\link{FixedExpressionData}}
#' @param i the rows index
#' @param j the column index
#' @param drop unused
#' @param ... unused
#' @return An object of class \code{\link{FixedExpressionData}}
#'
#' @family fixed data information extraction functions
#'
#' @export
#' @docType methods
#' @rdname bracket-methods
#'
#' @examples
#' data(exampleMAS5)
#' myData <- setNormalizationMethod(exampleMAS5, 'MAS5.0', targetValue=500)
#' dim(myData)
#' dim(myData[1:10,1:3])
#' dim(myData[[1:10,1:3]])
#' @export
setMethod("[", signature = signature(x = "FixedExpressionData", i = "ANY", j = "ANY"), 
    definition = function(x, i, j, ...) {
        return(x[[i, j, drop = FALSE]])
    })


#' @rdname bracket-methods
#' @aliases bracket,FixedExpressionData-method
#' @export
setMethod("[", signature = signature(x = "FixedExpressionData", i = "ANY", j = "missing"), 
    definition = function(x, i, j, ...) {
        return(x[[i, seq_len(ncol(x)), drop = FALSE]])
    })



#' @rdname bracket-methods
#' @aliases bracket,FixedExpressionData-method
#' @export
setMethod("[", signature = signature(x = "FixedExpressionData", i = "missing", j = "ANY"), 
    definition = function(x, i, j, ...) {
        return(x[[seq_len(nrow(x)), j, drop = FALSE]])
    })


#' @rdname bracket-methods
#' @aliases bracket,FixedExpressionData-method
#' @export
setMethod("[[", signature = signature(x = "FixedExpressionData", i = "ANY", j = "missing"), 
    definition = function(x, i, j, ...) {
        return(x[[i, seq_len(ncol(x)), drop = FALSE]])
    })

#' @rdname bracket-methods
#' @aliases bracket,FixedExpressionData-method
#' @export
setMethod("[[", signature = signature(x = "FixedExpressionData", i = "missing", j = "ANY"), 
    definition = function(x, i, j, ...) {
        return(x[[seq_len(nrow(x)), j, drop = FALSE]])
    })

#' @rdname bracket-methods
#' @aliases bracket,FixedExpressionData-method
#' @importFrom methods validObject
#' @export
setMethod("[[", signature = signature(x = "FixedExpressionData", i = "ANY", j = "ANY"), 
    definition = function(x, i, j, ...) {
        
        oldProcesses <- getTransformationProcesses(x)
        
        newObject <- new("FixedExpressionData", normalizationMethod = getNormalizationMethod(x), 
            normalizationParameters = getNormalizationParameters(x))
        explicitlyChangeExprs(newObject) <- rawExprs(x)[i, j, drop = FALSE]
        for (process in oldProcesses) {
            newObject <- addTransformationProcess(newObject, name = getName(process), 
                values = getValues(process))
        }
        validObject(newObject)
        return(newObject)
    })

#' @importFrom methods validObject
#' @importFrom Biobase copyEnv
setReplaceMethod("explicitlyChangeExprs", signature = signature(object = "FixedExpressionData", 
    value = "matrix"), definition = function(object, value) {
    
    # if (!all(dim(value)==dim(object))) stop('New and old expression matrices do not
    # have similar dimensions.')
    aData <- copyEnv(getExpressionEnvironment(object))
    aData[["expressionMatrix"]] <- value
    lockEnvironment(aData, bindings = TRUE)
    object@expressionEnvironment <- aData
    validObject(object)
    return(object)
})

