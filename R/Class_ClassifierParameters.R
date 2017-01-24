setOldClass("package_version")

#' @name ClassifierParameters
#' @title An S4 class to store classifier parameters.
#'
#'
#' @description  This class stores classifier related information
#' This is information on probe-sets used and their weightings,
#' means, standard deviations and covariance structure as observed in
#' the classifiers training data, and the description of the procedure
#' on how to preprocess new data prior to application of the classifier.
#'
#' @slot name A character string indicating the name of the classifier
#' @slot description A short description of the classifier
#' @slot citations A character vector of citations to literature
#' @slot normalizationMethod A character string indicating the normalization
#' method to apply
#' @slot eventChain A list of preprocessing steps

#' @slot probeNames A character vector
#' @slot intercept A numeric value
#' @slot weights A numeric vector
#' @slot decisionBoundaries A numeric vector with values that separate
#' the risk-groups

#' @slot doRun A function which is called for the actual classification

#' @slot means A numeric vector of probe-set means as observed in the
#' trainingsset (if available)
#' @slot sds A numeric vector of probe-set standard deviations as observed
#' in the trainingsset (if available)

#' @slot .geneClassifierVersion An object of class \code{\link{package_version}}
setClass("ClassifierParameters",
    representation = representation(
        citations = "vector",
        name = "character",
        description="character",
        normalizationMethod="character",
        probeNames = "character",
        weights = "numeric",
        intercept="numeric",
        means = "numeric",
        sds = "numeric",
        decisionBoundaries="numeric", #One or more classification thresholds
        doRun="function",
        eventChain="list",
        hasTrainingData="logical",
        .geneClassifierVersion="package_version"
    ),prototype=list(
        .geneClassifierVersion = packageVersionInternal()
    ),validity=function(object) {
        errTxt<-vector()
        if (length(getName(object))!=1){errTxt<-c(errTxt,"Classifier name must be given as a single character string")} 
        if (length(getNormalizationMethod(object))!=1){errTxt<-c(errTxt,"Normalization method must be given as a single character string")} 
        else if (!getNormalizationMethod(object)%in%getNormalizationMethods()){errTxt<-c(errTxt,"Unknown normalization method")}
        if (length(getDescription(object))!=1){errTxt<-c(errTxt,"Classifier description must be given as a single character string")} 
        if (length(getProbeNames(object))==0) {errTxt<-c(errTxt,"Classifier probe names must be given.")}
        if (length(getWeights(object))!=length(getProbeNames(object))) {errTxt<-c(errTxt,"Classifier weights and probe names are of different length.")}
        if (length(hasTrainingData(object))!=1) {errTxt<-c(errTxt,"Classifier hasTrainingData indication must be either TRUE of FALSE.")}
        ##Do not check eventchain and dorun here. They were already part of the unit tests
        if (length(errTxt)>0) {return(c(getName(object),errTxt))}
        return(TRUE)
    }
)
#' @importMethodsFrom Biobase preproc experimentData
#' @importFrom utils packageVersion packageName data
setMethod("ClassifierParameters",
    definition = function(name, doRun, hasTrainingData, normalizationMethod, weights, description="",intercept=numeric(0),means=numeric(0),sds=numeric(0),decisionBoundaries=numeric(0),eventChain = list(),citations=character(0)){
        
        new("ClassifierParameters",
            name=name,
            normalizationMethod=normalizationMethod,
            description=description,
            citations=citations,
            probeNames=names(weights),
            intercept=intercept,
            weights=weights,
            means=means,
            sds=sds,
            decisionBoundaries=decisionBoundaries,
            doRun=doRun,
            eventChain=eventChain,
            hasTrainingData=hasTrainingData
        )
    }
)

setMethod("show",
    signature  = signature(object = "ClassifierParameters"),
    definition = function(object){
    cat("______________\n")
        cat("Classifier: ", getName(object),"\n")
        cat("Description: ", getDescription(object),"\n")
        cat("Based on n = ", length(getProbeNames(object))," probe sets\n")
    boundaries<-getDecisionBoundaries
    if (length(boundaries)==0){
        cat("Number of risk groups: none (only continuous risk score)")
    } else if (length(boundaries)>0){
        cat("Number of risk groups: n = ", length(boundaries)+1,"\n")
    }
        cat("To be used on ", getNormalizationMethod(object)," normalized data\n")
        allCitations<-paste(unlist(getCitations(object)),collapse="; ")
        regExPattern<-paste("(?<=.{",options("width"),"})",sep="")
        cat( strsplit(allCitations, regExPattern, perl = TRUE)[[1]],sep="\n")
    }
)


setMethod("hasTrainingData",
    signature  = signature(object = "ClassifierParameters"),
    definition = function(object){
        return(object@hasTrainingData)
    }
)


#' @rdname getEventChain-methods
#' @aliases getEventChain,ClassifierParameters-method
#' @export
setMethod("getEventChain",
    signature  = signature(object = "ClassifierParameters"),
    definition = function(object){
        return(object@eventChain)
    }
)

#' @rdname getDescription-methods
#' @aliases getDescription,ClassifierParameters-method
#' @export
setMethod("getDescription",
    signature  = signature(object = "ClassifierParameters"),
    definition = function(object){
        return(object@description)
    }
)


#' @rdname getNormalizationMethod-methods
#' @aliases getNormalizationMethod,ClassifierParameters-method
#' @export
setMethod("getNormalizationMethod",
    signature  = signature(object = "ClassifierParameters"),
    definition = function(object){
        return(object@normalizationMethod)
    }
)

#' @rdname getDecisionBoundaries-methods
#' @aliases getDecisionBoundaries,ClassifierParameters-method
#' @export
setMethod("getDecisionBoundaries",
    signature  = signature(object = "ClassifierParameters"),
    definition = function(object){
        return(object@decisionBoundaries)
    }
)

#' @rdname getIntercept-methods
#' @aliases getIntercept,ClassifierParameters-method
#' @export
setMethod("getIntercept",
    signature  = signature(object = "ClassifierParameters"),
    definition = function(object){
        return(object@intercept)
    }
)

#' @rdname getTrainingData-methods
#' @aliases getTrainingData,ClassifierParameters-method
#' @export
setMethod("getTrainingData",
    signature  = signature(object = "ClassifierParameters"),
    definition = function(object){
        if (!hasTrainingData(object)) return(NULL)
        orgClassifierName<-gsub("_reweighted","",getName(object))
        expectedDataName<-paste(orgClassifierName,"_training",sep="")
        trainingsData<- local({
            try(eval(parse(text=expectedDataName)),silent=TRUE)
        })
        if (!inherits(trainingsData,"ExpressionSet")) stop("Could not obtain training data (",expectedDataName,") for classifier", getName(object))
        return(trainingsData[getProbeNames(object),])
    }
)

#' @rdname getProbeNames-methods
#' @aliases getProbeNames,ClassifierParameters-method
#' @export
setMethod("getProbeNames",
    signature  = signature(object = "ClassifierParameters"),
    definition = function(object){
        return(object@probeNames)
    }
)

#' @rdname getWeights-methods
#' @aliases getWeights,ClassifierParameters-method
#' @export
setMethod("getWeights",
    signature  = signature(object = "ClassifierParameters"),
    definition = function(object){
        pn<-getProbeNames(object)
        w<-object@weights
        if (length(pn)==length(w))  names(w)<-pn
        return(w)
    }
)


#' @importFrom Biobase copyEnv
setMethod("reWeightClassifier",
    signature  = signature(object = "ClassifierParameters",fixedExpressionData ="FixedExpressionData"),
    definition = function(object,fixedExpressionData){
        ##Function reWeightClassifier is assumed to be called by the runClassifier function which performs the required checks
         .reweight<-function(expressionSetTrain,betas,probesPresent,probesAbsent){
            expressionSetTrain<-t(exprs(expressionSetTrain))
            mat<-cbind(1,expressionSetTrain[,c(probesPresent,probesAbsent),drop=FALSE]) #complete data
            mat0<-cbind(1,expressionSetTrain[,probesPresent,drop=FALSE]) #known data
            mat1<-expressionSetTrain[,probesAbsent,drop=FALSE] #unknown data
            Omega<-solve(t(mat0)%*%mat0)%*%t(mat0)%*%mat1
            fi<-c(betas[1],betas[probesPresent])+Omega%*%betas[probesAbsent]
            ##Assume normality to adjust new variance and bias toward the old values ensure compatability with classification thresholds
            sdNew<-sd(mat0%*%fi)
            sdOld<-sd(mat%*%c(betas[1],betas[c(probesPresent,probesAbsent)]))
            fi[-1]<-fi[-1]*sdOld/sdNew
            orgMu<-sum(colMeans(mat)*c(betas[1],betas[c(probesPresent,probesAbsent)]))-betas[1]
            newMu<-sum(colMeans(mat0)*fi)-fi[1]
            fi[1]<-betas[1] + orgMu-newMu ##reset intercept to original value
            names(fi)<-c("Intercept",probesPresent)
            betas<-fi
            return(betas)
        }

        probesPresent <-getProbeNames(object)[ which(getProbeNames(object)%in%rownames(fixedExpressionData))]
        probesAbsent  <-getProbeNames(object)[ which(!getProbeNames(object)%in%rownames(fixedExpressionData))]
        
         
           trainingData<-getTrainingData(object)
        
        idx<-match(probesPresent,getProbeNames(object))
        betas<-.reweight(
            expressionSetTrain = trainingData,
            betas = c(getIntercept(object),getWeights(object)),
            probesPresent,probesAbsent
        )
        newIntercept<-betas[1]
        newWeights<-betas[-1,drop=FALSE]
        newObject<- ClassifierParameters(
            name=paste(getName(object),"_reweighted",sep=""),
            description=getDescription(object),
            normalizationMethod=getNormalizationMethod(object),
            weights=newWeights,
            intercept=newIntercept,
            means = getMeans(object)[idx],
            sds = getSds(object)[idx],
            decisionBoundaries=getDecisionBoundaries(object),
            doRun=getDoRun(object),
            citations = getCitations(object),
            eventChain=getEventChain(object),
            hasTrainingData= hasTrainingData(object)
        )
        return(newObject)
    }
)


#' @rdname getCitations-methods
#' @aliases getCitations,ClassifierParameters-method
#' @export
setMethod("getCitations",
    signature  = signature(object = "ClassifierParameters"),
    definition = function(object){
        return(object@citations)
    }
)

#' @rdname getName-methods
#' @aliases getName,ClassifierParameters-method
#' @export
setMethod("getName",
    signature  = signature(object = "ClassifierParameters"),
    definition = function(object){
        return(object@name)
    }
)

#' @rdname getMeans-methods
#' @aliases getMeans,ClassifierParameters-method
#' @export
setMethod("getMeans",
    signature  = signature(object = "ClassifierParameters"),
    definition = function(object){
        pn<-getProbeNames(object)
        m<-object@means
        if (length(pn)==length(m))  names(m)<-pn
        return(m)
    }
)

#' @rdname getSds-methods
#' @aliases getSds,ClassifierParameters-method
#' @export
setMethod("getSds",
    signature  = signature(object = "ClassifierParameters"),
    definition = function(object){
        pn<-getProbeNames(object)
        s<-object@sds
        if (length(pn)==length(s))  names(s)<-pn
    return(s)
    }
)


setMethod("getDoRun",
    signature  = signature(object = "ClassifierParameters"),
    definition = function(object){
    return(object@doRun)
    }
)
#' @rdname runClassifier-methods
#' @aliases runClassifier,character,FixedExpressionData-method
#' @export
setMethod("runClassifier",
    signature  = signature(classifierParameters = "character", fixedExpressionData="FixedExpressionData"),
    definition = function(classifierParameters, fixedExpressionData,...){
        dotList<-list(...)
        classifierParameters<-getClassifier(classifierParameters)
        runClassifier(classifierParameters,fixedExpressionData,...)
    }
)

#' @rdname runClassifier-methods
#' @aliases runClassifier,ClassifierResults,FixedExpressionData-method
#' @export
#' @importFrom methods new
#' @importFrom stats sd
#' @importMethodsFrom Biobase exprs
setMethod("runClassifier",
    signature  = signature(classifierParameters = "ClassifierParameters",fixedExpressionData ="FixedExpressionData"),
    definition = function(classifierParameters, fixedExpressionData,...){
        dotList<-list(...)
        nTP<-length(getTransformationProcesses(fixedExpressionData))
        if (nTP>0){warning("The provided data already contains transformationProcesses. Consider removing them by using the function 'removeTransformationProcesses'.")}
        normalizationMethodClassifier<-getNormalizationMethod(classifierParameters)
        normalizationMethodData<-getNormalizationMethod(fixedExpressionData)
        if (normalizationMethodClassifier!=normalizationMethodData){
            stop("Classifier: ", getName(classifierParameters),
            " can only be applied to ", normalizationMethodClassifier,
            " normalized data. The given data is ",normalizationMethodData,
            " nomalized."
            )
        }
        eventChain = getEventChain(classifierParameters)
        allow.reweighted = FALSE
        do.batchcorrection = TRUE
        if (!is.null(dotList[["allow.reweighted"]])) allow.reweighted<-dotList[["allow.reweighted"]]
        if (!is.null(dotList[["do.batchcorrection"]])) do.batchcorrection<-dotList[["do.batchcorrection"]]
        if (do.batchcorrection & ncol(fixedExpressionData)<3)  {
            stop("Only a small number of samples available (n=",
            ncol(fixedExpressionData),
            "). Please consider setting 'do.batchcorrection = FALSE' but be aware this could severely induce biases!"
            )
        }
        if (do.batchcorrection & ncol(fixedExpressionData)<20) {
            warning("Small number of available samples (n=",
            ncol(fixedExpressionData),
            "). Include more samples or consider setting 'do.batchcorrection = FALSE' but be aware this could severely induces biases!"
            )
        }
        weightingType = "complete"
        batchCorrection = "yes"
        probesAbsent <- getProbeNames(classifierParameters)[ which(!getProbeNames(classifierParameters)%in%rownames(fixedExpressionData))]
        if (length(probesAbsent)>0) {
            allowed.by.classifier = FALSE
            if (!is.null(eventChain[["allow.reweighted"]])) {
                if (eventChain[["allow.reweighted"]]){ ##allowed by classifier?
                    allowed.by.classifier = TRUE
                    if (allow.reweighted) {##Explicitly allowed by user?
                        classifierParameters<-reWeightClassifier(  classifierParameters,fixedExpressionData)
                        weightingType="reweighted"
                    } else {
                        stop("Found ",length(probesAbsent),
                            " of " , length(getProbeNames(classifierParameters)),
                            " covariates to be missing for classifier ", getName(classifierParameters),
                            ". Consider using 'allow.reweighted = TRUE'"
                            )
                    }
                }
            } 
            if (! allowed.by.classifier) {
                stop("Probe sets missing (n = ",length(probesAbsent),
                    " of " , length(getProbeNames(classifierParameters)),") which are needed to apply classifier ",
                    getName(classifierParameters),
                    ".  Reweighting is not allowed for this classifier."
                )
            }
        }    
        fixedExpressionData<-fixedExpressionData[ getProbeNames(classifierParameters), ] ###Filter probes
        for (eventIdx in seq_len(length(eventChain))){
            keyword = names(eventChain[eventIdx])
            keyvalue = as.numeric(eventChain[eventIdx])
            if (keyword == "targetValue"){
                setTargetValue(fixedExpressionData)<-keyvalue
            } else if (keyword == "truncate"){
                fixedExpressionData<-addTransformationProcess(fixedExpressionData,"LOWERTRUNCATE",keyvalue)
            }
            else if (keyword == "to.log"){
                fixedExpressionData<-addTransformationProcess(fixedExpressionData,"LOG",keyvalue)
            } else if (keyword == "to.meancentering" & keyvalue>0){
                if (do.batchcorrection) {
                    meanValues<-apply(exprs(fixedExpressionData),1,mean,na.rm=TRUE)
                } else {
                    meanValues<-getMeans(classifierParameters)
                    if (length(meanValues)==0) {stop("Classifier",getName(classifierParameters), "has no reference means given and thus only batchwise approach (i.e. do.batchcorrection = TRUE) is possible.")}
                }
                fixedExpressionData<-addTransformationProcess(fixedExpressionData,"SHIFTPERPROBE",values=-1*meanValues)
            } else if (keyword == "to.unitvariance" & keyvalue>0 ){
              
                if (do.batchcorrection) {
                    sdValues<-apply(exprs(fixedExpressionData),1,sd)
                } else {
                    sdValues<-getSds(classifierParameters)
                    if (length(sdValues)==0) {stop("Classifier",getName(classifierParameters), "has no reference sd's given and thus only batchwise approach (i.e. do.batchcorrection = TRUE) is possible.")}
                }
                fixedExpressionData<-addTransformationProcess(fixedExpressionData,"MULTIPLYPERPROBE",values=1/sdValues)
            }
            else if (keyword == "to.referencemeanvar" & keyvalue>0 ) {
                if (do.batchcorrection) {
                    meanValues<-apply(exprs(fixedExpressionData),1,mean,na.rm=TRUE)
                    sdValues<-apply(exprs(fixedExpressionData),1,sd,na.rm=TRUE)
                    refMeanValues<-getMeans(classifierParameters)[rownames(fixedExpressionData)]
                    refSdValues<-getSds(classifierParameters)[rownames(fixedExpressionData)]
                    #Note that the required refmeans and refsds were already checked for in the unit tests!
                    fixedExpressionData<-addTransformationProcess(fixedExpressionData,"SHIFTPERPROBE",values=-meanValues)
                    fixedExpressionData<-addTransformationProcess(fixedExpressionData,"MULTIPLYPERPROBE",values=refSdValues/sdValues)
                    fixedExpressionData<-addTransformationProcess(fixedExpressionData,"SHIFTPERPROBE",values=refMeanValues)
                } else {
                    #Do nothing
                }
            }
            else if (keyword == "allow.reweighted") {}
            else {
                stop("Encountered unknown keyword (",keyword,") in the eventChain of ", getName(classifierParameters))
            }
        }
        score<-getDoRun(classifierParameters)(classifierParameters,fixedExpressionData)
        classifierResults<-ClassifierResults(classifierParameters=classifierParameters,score=score,batchCorrection=do.batchcorrection,weightingType=weightingType)
        return(classifierResults)
    }
)


