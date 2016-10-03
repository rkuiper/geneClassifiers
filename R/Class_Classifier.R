setOldClass("package_version")

setGeneric("as.factor", function (x) {base::as.factor(x)})

setClass("Classifier",
representation(params="list",score="numeric",batchCorrection="logical", ins="character",weightingType='character',version="package_version"))

setMethod("initialize", signature(.Object = "Classifier"),
    function(.Object, ...){ 
        x <- list(...)
        if (is.null(x$weightingType))         stop('weightingType must be given');
        if (is.null(x$batchCorrection)) stop('batchCorrection must be given')
        if (is.null(x$ins))         stop('ins must be given');    
        if (is.null(x$score))         stop('score must be given');    
        if (is.null(x$params))         stop('params must be given');    
        if (!all(x$ins%in%names(x$params$weights))) stop('not all ins in the classifier');    
        if (!(x$weightingType%in%c("complete","reweighted"))) stop('weightingType should be one of: complete, reweighted');    
        .Object@weightingType = x$weightingType;
        .Object@ins = x$ins;
        .Object@score = x$score;
        .Object@params = x$params;
        .Object@batchCorrection = x$batchCorrection;
        .Object@version = packageVersion("geneClassifiers");
        return(.Object)
    }
)


setMethod("show", "Classifier",
    function(object){
        classifications<-sapply(object@params$decisionBoundaries,FUN=function(x,y){y>x},y=object@score)
        classifications<-matrix(classifications,ncol=1)
        classifications<-factor(rowSums(classifications),levels=c(0:length(object@params$decisionBoundaries)),labels=paste("Risk",as.roman(c(0:length(object@params$decisionBoundaries))+1),sep='-'),ordered=TRUE)

        riskGroup<-table(classifications)
        cat("Note: Research use only\n")
        cat("Classifier:" , object@params$name,"\n")
        for (i in paste(">",object@params$decisionBoundaries ,":",names(riskGroup)[-1])){cat("\t",i,"\n");}
        print(riskGroup)

        cat("\n\tBatch corrected   :",c("no","yes")[object@batchCorrection+1],"\n")
        cat("\tweighting type      :",object@weightingType,"\n")
        cat("----------------\n")
    }
)
setMethod("as.numeric",  signature("Classifier"), function(x) { x@score })
setMethod("as.factor",  signature("Classifier"), function(x) {    
    classifications<-sapply(x@params$decisionBoundaries,FUN=function(x,y){y>x},y=x@score)
    classifications<-factor(rowSums(classifications),levels=c(0:length(x@params$decisionBoundaries)))
    classifications
})

