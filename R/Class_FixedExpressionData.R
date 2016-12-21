#' @aliases getNormalizationMethods
#' @title Show currenly implemented normalization methods.
#'
#' @description \code{getNormalizationMethods} returns a character vector of currenlty
#' available normalization methods.
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
#' myData <- setNormalizationMethod(exampleMAS5, "MAS5.0",targetValue=500)
#' results <- runClassifier('UAMS70', myData)
#'
#' getScores( results )
#' getClassifications( results )
#'
#' @export
getNormalizationMethods<-function(){
    ##By convention uppercase
    return(c("MAS5.0","GCRMA"))
}

#' Example MAS5.0 ExpressionSet
#'
#'An \code{\link{ExpressionSet}}. The data contains a sample of gene expression data from
#'patients included in the HOVON65/GMMG-HD4 trial on multiple myeloma. The data was MAS5.0
#'normalized to a target value of 500.
#'
"exampleMAS5"

setOldClass("package_version")

#' @name FixedExpressionData
#' @title An S4 class to store classifier parameters.
#'
#' @description  This class stores gene expression data together with information on
#' the normalization method and additional normalization related parameters. In order
#' to ensure the data is not manipulated in unforeseen ways, manipulation is strictly
#' controled through adding transformations which are predefined in the \code{TransformationProcess}-class.
#' Upon reading the data by the \code{exprs} function, the transformations are
#' performed in the order the were added.
#'
#' @slot normalizationMethod A character string indicating the normalization method that
#' was applied to the data. Possible values are give by \code{\link{getNormalizationMethods}}.
#' @slot expressionEnvironment A locked environment in which the expression matrix is stored.
#' @slot normalizationParameters A list with normalization specific values.
#' @slot transformationProcess A locked environment to which the transformation processes are added.
#' @slot .geneClassifierVersion An object of class \code{\link{package_version}}
setClass("FixedExpressionData",
    representation = representation(
        normalizationMethod="character",
        expressionEnvironment= "environment",
        normalizationParameters="list",
        transformationProcess="environment",
        .geneClassifierVersion="package_version"
    ),prototype=list(
    .geneClassifierVersion = packageVersionInternal()
    ),validity=function(object) {
        errTxt<-vector()
        if (!getNormalizationMethod(object)%in%getNormalizationMethods()){errTxt<-c(errTxt,"Unknown normalization method")}
        if (!is.list(getTransformationProcesses(object))){errTxt<-c(errTxt,"Unknown transformationProcess type")}

        if (getNormalizationMethod(object)=="MAS5.0") {
            if  (is.null(getNormalizationParameters(object)[["targetValue"]])) {errTxt<-c(errTxt,"targetValue must be given for MAS5.0") }
        }
        if (length(errTxt)>0) {return(errTxt)}
        return(TRUE)
    }
)

#' @importFrom utils packageVersion packageName
#' @importFrom stats quantile
setMethod("FixedExpressionData",
    signature  = signature(normalizationMethod="character",expressionMatrix="matrix"),
    definition = function(normalizationMethod,expressionMatrix,...){
        dotList<-list(...)
        normalizationMethod<-match.arg(normalizationMethod,getNormalizationMethods())
        .getTargetValue<-function(expressionMatrix){
            ##expressionMatrix must be not log2 transformed!
            e<-t(expressionMatrix)
            q<-apply(e,1,quantile,c(0.02,0.98),type=2)
            naMat<-(e>=q[1,] & e<=q[2,])
            naMat[!naMat]<-NA
            targetValue<-unique(signif(rowMeans(naMat*e,na.rm=TRUE),6))
            if(length(targetValue)!=1){targetValue<-NULL}
            return(targetValue)
        }
        .isAlreadylog2Transformed.check<-function(expressionMatrix){
            storedSeed<-.Random.seed
            set.seed(1)
            aSample<-sample(expressionMatrix,1e5,replace=TRUE)
            .Random.seed<-storedSeed
            isLog2Transformed<-as.logical(quantile(aSample,0.75,na.rm=TRUE)<30 | any(expressionMatrix<0))
            return(isLog2Transformed)
        }
        if (!is.matrix(expressionMatrix)){
            stop("'expressionMatrix' must be a numeric matrix.")
        }
        if (mode(expressionMatrix)!="numeric"){
            stop("'expressionMatrix' must be numeric.")
        }
        if (is.null(rownames(expressionMatrix)) | is.null(colnames(expressionMatrix))){
            stop("'expressionMatrix' must have row and column names.")
        }
        if (sum(duplicated(rownames(expressionMatrix)))>0 ) {
            stop("duplicated row names found.")
        }
        if (sum(duplicated(colnames(expressionMatrix)))>0) {
            stop("duplicated column names found.")
        }

        expressionEnvironment<-new.env()
        expressionEnvironment[["expressionMatrix"]]<-expressionMatrix
        normalizationParameters<-list()
        transformationProcess<-new.env()
        transformationProcess[["processes"]]<-list()

        ##Normalization method specific code:
        if (normalizationMethod == "MAS5.0"){
            calculatedEstimate<-.isAlreadylog2Transformed.check(expressionEnvironment[["expressionMatrix"]])
            if (!is.null(dotList[["isLog2Transformed"]])) {
                if(dotList[["isLog2Transformed"]]!=calculatedEstimate) {
                    warning("Reconsider argument isLog2Transformed = ",dotList[["isLog2Transformed"]])
                }
            } else {
                dotList[["isLog2Transformed"]]<-calculatedEstimate
                if (calculatedEstimate==TRUE) {
                     warning("The data seems to be log2 transformed. If this is incorrect, please provide the logical argument 'isLog2Transformed=FALSE'.")
                }
            }
            if (dotList[["isLog2Transformed"]]){
                 expressionEnvironment[["expressionMatrix"]]<-2^expressionEnvironment[["expressionMatrix"]]
            }

            calculatedTargetValue<-.getTargetValue(expressionEnvironment[["expressionMatrix"]])
            if (!is.null(dotList[["targetValue"]])){
                if (!is.numeric(dotList[["targetValue"]])){
                    stop("targetValue must be a positive numeric value")
                }
                if (dotList[["targetValue"]]<=0){
                    stop("targetValue must be a positive numeric value")
                }
            }
            if (!is.null(dotList[["targetValue"]]) & !is.numeric(dotList[["targetValue"]])){
                stop("targetValue must be a positive numeric value")
            }
            if (!is.null(dotList[["targetValue"]])&!is.null(calculatedTargetValue)){
                if(dotList[["targetValue"]]!=calculatedTargetValue) {
                    warning("Using the given targetValue, but it does not correspond to the observed targetValue in the data of ",calculatedTargetValue,".")
                }
            } else if (is.null(dotList[["targetValue"]]) & is.null(calculatedTargetValue)){
                stop("Cannot determine the targetValue that was used for MAS5.0. Please provide the right 'targetValue=' as an argument.")
            } else if (is.null(dotList[["targetValue"]]) & !is.null(calculatedTargetValue)) {
                dotList[["targetValue"]]<-calculatedTargetValue
            }
            normalizationParameters[["targetValue"]]<- dotList[["targetValue"]]
        }
        if (normalizationMethod == "GCRMA"){}
        lockEnvironment(expressionEnvironment,bindings=TRUE)
        lockEnvironment(transformationProcess,bindings=TRUE)
        new("FixedExpressionData",normalizationMethod=normalizationMethod, 
            expressionEnvironment= expressionEnvironment,
            normalizationParameters=normalizationParameters,
            transformationProcess=transformationProcess
        )
    }
)

setMethod("show",
    signature  = "FixedExpressionData",
    definition = function(object){
        cat("____________________\n")
        cat("Fixed expression set\n\n")
        cat("Normalization method: ",getNormalizationMethod(object),"\n")
        cat("Number of samples   : ",dim(object)[2],"\n")
        cat("Number of features  : ",dim(object)[1],"\n")
        tp<-getTransformationProcesses(object)
        if (length(tp)>0){
            tpNames<-lapply(getTransformationProcesses(object),getName)
            tpNames<-paste(unlist(tpNames),collapse="; ")
            cat("Applied transformation processes: ",tpNames,"\n")
        }
    }
)

#' @importFrom methods new
#' @importFrom Biobase copyEnv
setMethod("addTransformationProcess",
    signature  = signature(object="FixedExpressionData",name="character"),
    definition = function(object,name,values,...){
        process<-TransformationProcess(name = name, values = values )
        aData <- copyEnv( object@transformationProcess )
        curLen <- length(aData[["processes"]])
        aData[["processes"]][[curLen+1]] <- process
        lockEnvironment(aData,bindings=TRUE)
        object@transformationProcess<-aData
        return( object )
    }
)

setMethod("getTransformationProcesses",
    signature  = signature(object="FixedExpressionData"),
    definition = function(object){
        return(object@transformationProcess[["processes"]])
    }
)

#' @importFrom Biobase copyEnv
setMethod("removeTransformationProcesses",
    signature  = signature(object="FixedExpressionData"),
    definition = function(object,n=1){
        aData <- copyEnv(object@transformationProcess)
        curLen<-length(aData[["processes"]])
        n<-min(n,curLen)
        aData[["processes"]]<-aData[["processes"]][seq_len(max(0,curLen-n))]
        lockEnvironment(aData,bindings=TRUE)
        object@transformationProcess<-aData
        return(object)
    }
)

#' @rdname getNormalizationMethod-methods
#' @aliases getNormalizationMethod,FixedExpressionData-method
#' @export
setMethod("getNormalizationMethod",
    signature  = signature("FixedExpressionData"),
    definition = function(object) {
        object@normalizationMethod
    }
)

setMethod("getNormalizationParameters",
    signature  = signature("FixedExpressionData"),
    definition = function(object) {
        object@normalizationParameters
    }
)

#' @rdname getTargetValue-methods
#' @aliases getTargetValue,FixedExpressionData-method
#' @export
setMethod("getTargetValue",
    signature  = signature("FixedExpressionData"),
    definition = function(object) {
        return( getNormalizationParameters(object)[["targetValue"]] )
    }
)

setReplaceMethod("setTargetValue",
    signature  = signature(object="FixedExpressionData",value="numeric"),
    definition = function(object,value) {
        ##Replaces the rawexpression data  with the new targetvalue
        if (!getNormalizationMethod(object)%in%c("MAS5.0")){
            stop("Function 'setTargetValue' does only apply to MAS5.0 normalized data")
        }
        currentTargetValue<-getTargetValue(object)
        e<-rawExprs(object)
        explicitlyChangeExprs(object)<-e*value/currentTargetValue
        object@normalizationParameters[["targetValue"]]<-value
        return(object)
    }
)

setMethod("getExpressionEnvironment",
    signature  = signature("FixedExpressionData"),
    definition = function(object) {
        return( object@expressionEnvironment )
    }
)

#' @importMethodsFrom Biobase exprs
setMethod("exprs",
    signature  = signature("FixedExpressionData"),
    definition = function(object) {
        e<-rawExprs(object)
        processes<-getTransformationProcesses(object)
        for (process in processes){
            e<-runProcess(process,e)
        }
        return(e)
    }
)

setMethod("rawExprs",
    signature  = signature("FixedExpressionData"),
    definition = function(object) {
        e<-getExpressionEnvironment(object)[["expressionMatrix"]]
        return( e )
    }
)

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
#' myData <- setNormalizationMethod(exampleMAS5, "MAS5.0", targetValue=500)
#' dim(myData)
#' dim(myData[1:10,1:3])
setMethod("dim",
    signature  = signature("FixedExpressionData"),
    definition = function(x) {
        return( dim(rawExprs(x)) )
    }
)

#' @importMethodsFrom BiocGenerics rownames
setMethod("rownames",
    signature  = signature(x="FixedExpressionData"),
    definition = function(x) {
        return( rownames(rawExprs(x)) )
    }
)

#' @importMethodsFrom BiocGenerics colnames
setMethod("colnames",
    signature  = signature(x="FixedExpressionData"),
    definition = function(x) {
        return( colnames(rawExprs(x)) )
    }
)

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
#' myData <- setNormalizationMethod(exampleMAS5, "MAS5.0", targetValue=500)
#' dim(myData)
#' dim(myData[1:10,1:3])
#' dim(myData[[1:10,1:3]])
#' @export
setMethod("[",
    signature  = signature(x = "FixedExpressionData"),
    definition = function(x, i, j, ...) {
        return( x[[i,j,...]] )
    }
)

#' @rdname bracket-methods
#' @aliases bracket,FixedExpressionData-method
#' @export
setMethod("[[",
    signature  = signature(x="FixedExpressionData"),
    definition = function(x, i, j, ...) {
        dotList<-list(normalizationMethod=getNormalizationMethod(x),
        expressionMatrix=rawExprs(x)[i,j,...])
        dotList<-c(dotList,getNormalizationParameters(x))
        newObject<-do.call(FixedExpressionData,dotList)
        oldProcesses<-getTransformationProcesses(x)
        for (process in oldProcesses){
            newObject<-addTransformationProcess(newObject,name=getName(process),values=getValues(process))
        }
        return(newObject)
    }
)

#' @importFrom Biobase copyEnv
setReplaceMethod("explicitlyChangeExprs",
    signature  = signature(object="FixedExpressionData",value="matrix"),
    definition = function(object, value) {
        aData <- copyEnv(getExpressionEnvironment(object))
        aData[["expressionMatrix"]]<-value
        lockEnvironment(aData,bindings=TRUE)
        object@expressionEnvironment<-aData
        return(object)
    }
)

