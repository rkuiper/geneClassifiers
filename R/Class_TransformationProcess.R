setOldClass("package_version")

setClass("TransformationProcess", representation(name = "character", values = "environment", 
    hook = "function", .geneClassifierVersion = "package_version"), prototype = list(.geneClassifierVersion = packageVersionInternal()))

#' @importFrom utils packageVersion packageName
setMethod("TransformationProcess", signature(name = "character", values = "ANY"), 
    definition = function(name, values, ...) {
        if (!name %in% c("LOG", "SHIFTPERPROBE", "LOWERTRUNCATE", "MULTIPLYPERPROBE")) {
            stop("Trying to initialze an unknown process:", name, ".")
        }
        
        procValues <- new.env()
        procValues[["processes"]] <- list()
        
        if (name == "LOG") {
            if (length(values) == 0) {
                stop("LOG TransformationProcess values argument must be given: the base number for the logarithm.")
            }
            if (length(values) != 1) {
                stop("LOG TransformationProcess only takes one value: the base number for the logarithm.")
            }
            if (!is.numeric(values)) {
                stop("LOG TransformationProcess requires a numeric base value")
            }
            if (values <= 0) {
                stop("LOG TransformationProcess requires a positive numeric base value")
            }
            procValues[["processes"]] <- values
            hook <- function(object, expressionMat) {
                expressionMat[] <- log(expressionMat, getValues(object))
                return(expressionMat)
            }
        }
        if (name == "SHIFTPERPROBE") {
            if (!is.numeric(values) | !is.vector(values)) {
                stop("SHIFTPERPROBE TransformationProcess requires a numeric value or numeric vector named by probe names when applied.")
            }
            if (length(values) > 1 & is.null(names(values))) {
                stop("SHIFTPERPROBE TransformationProcess requires a numeric value or numeric vector named by probe names when applied.")
            }
            procValues[["processes"]] <- values
            hook <- function(object, expressionMat) {
                values <- getValues(object)
                if (length(values) == 1 & all(is.null(names(values)))) {
                  values <- rep(values, nrow(expressionMat))
                  names(values) <- rownames(expressionMat)
                }
                intersectingNames <- intersect(names(values), rownames(expressionMat))
                expressionMat[intersectingNames, ] <- expressionMat[intersectingNames, 
                  ] + values[intersectingNames]
                return(expressionMat)
            }
        }
        if (name == "LOWERTRUNCATE") {
            if (length(values) != 1) {
                stop("LOWERTRUNCATE TransformationProcess only takes one value: the value below which to truncate.")
            }
            if (!is.numeric(values)) {
                stop("LOWERTRUNCATE TransformationProcess requires a numeric value")
            }
            procValues[["processes"]] <- values
            hook <- function(object, expressionMat) {
                expressionMat[expressionMat <= getValues(object)] <- getValues(object)
                return(expressionMat)
            }
        }
        if (name == "MULTIPLYPERPROBE") {
            if (!is.numeric(values) | !is.vector(values)) {
                stop("MULTIPLYPERPROBE TransformationProcess requires a numeric value or numeric vector named by probe names when applied.")
            }
            if (length(values) > 1 & is.null(names(values))) {
                stop("MULTIPLYPERPROBE TransformationProcess requires a numeric value or numeric vector named by probe names when applied.")
            }
            procValues[["processes"]] <- values
            hook <- function(object, expressionMat) {
                values <- getValues(object)
                if (length(values) == 1 & all(is.null(names(values)))) {
                  values <- rep(values, nrow(expressionMat))
                  names(values) <- rownames(expressionMat)
                }
                intersectingNames <- intersect(names(values), rownames(expressionMat))
                expressionMat[intersectingNames, ] <- expressionMat[intersectingNames, 
                  ] * values[intersectingNames]
                return(expressionMat)
            }
        }
        lockEnvironment(procValues, bindings = TRUE)
        new("TransformationProcess", name = name, values = procValues, hook = hook)
    })

setMethod("runProcess", signature = signature(object = "TransformationProcess", expressionMat = "matrix"), 
    definition = function(object, expressionMat) {
        return(object@hook(object, expressionMat))
    })

setMethod("getValues", signature = signature(object = "TransformationProcess"), definition = function(object) {
    return(object@values[["processes"]])
})

#' @rdname getName-methods
#' @aliases getName,TransformationProcess-method
setMethod("getName", signature = signature(object = "TransformationProcess"), definition = function(object) {
    return(object@name)
})
