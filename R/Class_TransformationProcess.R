setOldClass("package_version")

setClass("TransformationProcess", representation(
    name="character",
    values= "environment",
    hook="function",
    .geneClassifierVersion="package_version"
))

#' @importFrom utils packageVersion packageName
setMethod("initialize", signature(.Object = "TransformationProcess"),
    function(.Object,...){
        dotList <- list(...)
        if (!dotList[["name"]]%in%c("LOG","SHIFTPERPROBE","LOWERTRUNCATE","MULTIPLYPERPROBE")){stop("Trying to initialze an unknown process:", dotList[["name"]],".")}
        .Object@.geneClassifierVersion=packageVersion(packageName())
        .Object@name=dotList[["name"]]
        .Object@values<-new.env()
        .Object@values[["processes"]]<-list()
        if (.Object@name=="LOG"){
            if (length(dotList[["values"]])==0){stop("LOG TransformationProcess values argument must be given: the base number for the logarithm.")}
            if (length(dotList[["values"]])!=1){stop("LOG TransformationProcess only takes one value: the base number for the logarithm.")}
            if (!is.numeric(dotList[["values"]])){stop("LOG TransformationProcess requires a numeric base value")}
            if (dotList[["values"]]<=0){stop("LOG TransformationProcess requires a positive numeric base value")}
            .Object@values[["processes"]]<-dotList[["values"]]
            .Object@hook<-function(expressionMat){
                expressionMat[]<-log(expressionMat,getValues(.Object))
                return(expressionMat)
            }
        }
        if (.Object@name=="SHIFTPERPROBE"){
            if (!is.numeric(dotList[["values"]]) | !is.vector(dotList[["values"]])){stop("SHIFTPERPROBE TransformationProcess requires a numeric value or numeric vector named by probe names when applied.")}
            if (length(dotList[["values"]])>1 & is.null(names(dotList[["values"]]))){ stop("SHIFTPERPROBE TransformationProcess requires a numeric value or numeric vector named by probe names when applied.")}
            .Object@values[["processes"]]<-dotList[["values"]]
            .Object@hook<-function(expressionMat){
                values<-getValues(.Object)
                if (length(values)==1 & all(is.null(names(values)))){
                    values<-rep(values,nrow(expressionMat))
                    names(values)<-rownames(expressionMat)
                }
                intersectingNames<-intersect(names(values),rownames(expressionMat))
                expressionMat[intersectingNames,]<-expressionMat[intersectingNames,]+values[intersectingNames]
                return(expressionMat)
            }
        }
        if (.Object@name=="LOWERTRUNCATE"){
            if (length(dotList[["values"]])!=1){stop("LOWERTRUNCATE TransformationProcess only takes one value: the value below which to truncate.")}
            if (!is.numeric(dotList[["values"]])){stop("LOWERTRUNCATE TransformationProcess requires a numeric value")}
            .Object@values[["processes"]]<-dotList[["values"]]
            .Object@hook<-function(expressionMat){
                expressionMat[expressionMat<=getValues(.Object)]<-getValues(.Object)
                return(expressionMat)
            }
        }
        if (.Object@name=="MULTIPLYPERPROBE"){
            if (!is.numeric(dotList[["values"]]) | !is.vector(dotList[["values"]])){ stop("MULTIPLYPERPROBE TransformationProcess requires a numeric value or numeric vector named by probe names when applied.")}
            if (length(dotList[["values"]])>1 & is.null(names(dotList[["values"]]))){ stop("MULTIPLYPERPROBE TransformationProcess requires a numeric value or numeric vector named by probe names when applied.")}
            .Object@values[["processes"]]<-dotList[["values"]]
            .Object@hook<-function(expressionMat){
                values<-getValues(.Object)
                if (length(values)==1 & all(is.null(names(values)))){
                    values<-rep(values,nrow(expressionMat))
                    names(values)<-rownames(expressionMat)
                }
                intersectingNames<-intersect(names(values),rownames(expressionMat))
                expressionMat[intersectingNames,]<-expressionMat[intersectingNames,]*values[intersectingNames]
                return(expressionMat)
            }
        }
        lockEnvironment(.Object@values,bindings=TRUE)
        return(.Object)
    }
)

setMethod("runProcess",
    signature  = signature(object="TransformationProcess",expressionMat="matrix"),
    definition = function(object,expressionMat) {
        return(object@hook(expressionMat))
    }
)

setMethod("getValues",
    signature  = signature(object="TransformationProcess"),
    definition = function(object) {
        return( object@values[["processes"]])
    }
)

#' @rdname getName-methods
#' @aliases getName,TransformationProcess-method
setMethod("getName",
    signature  = signature(object="TransformationProcess"),
    definition = function(object) {
        return( object@name)
    }
)
