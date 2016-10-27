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

#' @slot data  A locked environment in which - if available - the expression
#' matrix is stored.
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
    data = "environment",
     probeNames = "character",
    weights = "numeric",
    intercept="numeric",
    means = "numeric",
    sds = "numeric",
    decisionBoundaries="numeric", #One or more classification thresholds
    doRun="function",
    eventChain="list",
       .geneClassifierVersion="package_version"
    )
)
#' @importMethodsFrom Biobase preproc experimentData
#' @importFrom utils packageVersion packageName data
setMethod("initialize",
    signature  = signature(.Object = "ClassifierParameters"),
    definition = function(.Object, ...){
        dotList <- list(...)
    .Object@.geneClassifierVersion=packageVersion(packageName())

    loadTrainingData<-function(classifierName){
        trainingsData<-NULL
        expectedDataName<-paste(classifierName,"_ExpressionSet_training",sep="")
        theData<-data(package = packageName())
        if (expectedDataName%in%theData[["results"]][,"Item"]){
            trainingsData<- local({
                txt<-paste("data(",expectedDataName,",envir=environment())",sep="")
                eval(parse(text=txt))
                eval(parse(text=paste("return(",expectedDataName,")",sep="")))
               })
            if (!inherits(trainingsData,"ExpressionSet")){
                stop("The trainings data available for classifier ",
                classifierName," is not of type ExpressionSet")
            }
        }


        return(trainingsData)
    }
    if (!is.character(dotList[["name"]]) | length(dotList[["name"]])!=1 ){
        stop("Classifier name must be given as a single character string")
    }
    if (!is.character(dotList[["normalizationMethod"]]) | length(dotList[["normalizationMethod"]])!=1 ){
        stop("normalizationMethod must be given as a single character string for classifier ",dotList[["name"]],".")
    }
    dotList[["normalizationMethod"]] <-match.arg(dotList[["normalizationMethod"]],getNormalizationMethods())

    if (is.null(dotList[["description"]])){
        dotList[["description"]]<-character(0)
    }
    if (!is.character(dotList[["description"]]) | length(dotList[["name"]])!=1 ){
        stop("Classifier description must be given as a single character string for classifier ",dotList[["name"]],".")
    }
    ##Try to load the training data
    trainingData<-loadTrainingData(dotList[["name"]])
    if (is.null(dotList[["probeNames"]])){
        stop("Please provide probeNames for classifier ",dotList[["name"]],".")
    }
    if (!is.null(trainingData)){
        intersection<-intersect(rownames(trainingData),dotList[["probeNames"]])
        if (    length(intersection)!=nrow(trainingData) |
            length(intersection)!=length(dotList[["probeNames"]])
        ) { stop("The trainings dataset was found to contain more (or less) probes than given by the probeNames argument for classifier ",dotList[["name"]],".")
        }
    }
    if (is.null(dotList[["intercept"]])){
        dotList[["intercept"]]<-numeric(0)
    }
    if (is.null(dotList[["weights"]])){
        stop("Please provide weights for classifier ",dotList[["name"]],".")
    }

    if (length(dotList[["weights"]])!=length(dotList[["probeNames"]])){
        stop("Found a different number of weights and probeNames for classifier ",dotList[["name"]],".")
    }



    if (is.null(dotList[["means"]]) ){
        dotList[["means"]]<-as.numeric(rep(NA,length(dotList[["probeNames"]])))
    }
    if (length(dotList[["means"]])!=length(dotList[["probeNames"]])){
        stop("Found a different number of means and probeNames for classifier ",dotList[["name"]]," .")
    }

    if (is.null(dotList[["sds"]]) ){
        dotList[["sds"]]<-as.numeric(rep(NA,length(dotList[["probeNames"]])))
    }
    if (length(dotList[["sds"]])!=length(dotList[["probeNames"]])){
        stop("Found a different number of sds values and probeNames for classifier ",dotList[["name"]],".")
    }

    if (is.null(dotList[["decisionBoundaries"]]) ){
        dotList[["decisionBoundaries"]]<-numeric(0)
    }

    if (!is.null(trainingData)){
        preprocData<-preproc(experimentData(trainingData))
        if (is.null(preprocData[["normalisationMethod"]])) {
            stop("Please provide a 'normalisationMethod' in the preproc annotation of the trainingset" )
        }
        if (!preprocData[["normalisationMethod"]]%in%getNormalizationMethods()){
            stop("The 'normalisationMethod' in the preproc annotation of the trainingset is not one of:",getNormalizationMethods() )
        }
        if (preprocData[["normalisationMethod"]]!=dotList[["normalizationMethod"]]){
            stop("The 'normalisationMethod' in the preproc annotation of the trainingset (",preprocData[["normalisationMethod"]],") does not correspond to the normalisationMethod (",dotList[["normalizationMethod"]],") given for for classifier ",dotList[["name"]]," .")
        }
        if (preprocData[["normalisationMethod"]]=="MAS5.0"){
            if (is.null(preprocData[["targetValue"]])) {
                stop("Please provide a 'targetValue' in the preproc annotation of the trainingset" )
            }
        }

    }
    if (is.null(dotList[["eventChain"]])){
        stop("No eventChain given for for classifier ",dotList[["name"]]," .")
    }
    doCheckEventChains(dotList[["eventChain"]],dotList[["normalizationMethod"]])

    .Object@name=dotList[["name"]]
    .Object@normalizationMethod=dotList[["normalizationMethod"]]
    .Object@description=dotList[["description"]]
    .Object@citations=dotList[["citations"]]
        .Object@data=new.env()
    .Object@data[["trainingData"]]<-trainingData
        .Object@probeNames=dotList[["probeNames"]]
        .Object@intercept=dotList[["intercept"]]
        .Object@weights=dotList[["weights"]]
        .Object@means=dotList[["means"]]
        .Object@sds=dotList[["sds"]]
        .Object@decisionBoundaries=dotList[["decisionBoundaries"]]
        .Object@doRun=dotList[["doRun"]]
    .Object@eventChain=dotList[["eventChain"]]

        return(.Object)
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
        return(object@data[["trainingData"]])
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
        w<-object@weights
    names(w)<-getProbeNames(object)
        return(w)
    }
)

#' @importFrom Biobase copyEnv
setMethod("reWeightClassifier",
    signature  = signature(object = "ClassifierParameters",probeNames="character", intercept="numeric",weights="numeric"),
    definition = function(object,probeNames, intercept,weights){
    if (length(probeNames)!=length(weights)){ stop("Length probeNames and weights do not correspond during reweighting")}
    curProbeNames<-getProbeNames(object)
    if (!all(probeNames%in%curProbeNames)){stop("Some given probe names are not part of the classifier")}
    idx<-match(probeNames,curProbeNames)

    object@weights<-weights
    object@intercept<-intercept
    object@probeNames<-probeNames
    object@means<-object@means[idx]
    object@sds<-object@sds[idx]
    if (!is.null(getTrainingData)){
        aData<-copyEnv(object@data)
        aData[["trainingData"]]<-aData[["trainingData"]][probeNames,]
        lockEnvironment(aData,bindings=TRUE)
        object@data<-aData
    }
    object@name=paste(object@name,"_reweighted",sep="")

        return(object)
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

setMethod("getMeans",
    signature  = signature(object = "ClassifierParameters"),
    definition = function(object){
    m<-object@means
    names(m)<-getProbeNames(object)
        return(m)
    }
)

setMethod("getSds",
    signature  = signature(object = "ClassifierParameters"),
    definition = function(object){
        s<-object@sds
    names(s)<-getProbeNames(object)
    return(s)
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

        .reweight<-function(expressionSetTrain,betas,probesPrensent,probesAbsent){
            expressionSetTrain<-t(exprs(expressionSetTrain))
            mat<-cbind(1,expressionSetTrain[,c(probesPrensent,probesAbsent),drop=FALSE]) #complete data
            mat0<-cbind(1,expressionSetTrain[,probesPrensent,drop=FALSE]) #known data
            mat1<-expressionSetTrain[,probesAbsent,drop=FALSE] #unknown data
            Omega<-solve(t(mat0)%*%mat0)%*%t(mat0)%*%mat1
            fi<-c(betas[1],betas[probesPrensent])+Omega%*%betas[probesAbsent]
            ##Assume normality to adjust new variance and bias toward the old values ensure compatability with classification thresholds
            sdNew<-sd(mat0%*%fi)
            sdOld<-sd(mat%*%c(betas[1],betas[c(probesPrensent,probesAbsent)]))
            fi[-1]<-fi[-1]*sdOld/sdNew
            orgMu<-sum(colMeans(mat)*c(betas[1],betas[c(probesPrensent,probesAbsent)]))-betas[1]
            newMu<-sum(colMeans(mat0)*fi)-fi[1]
            fi[1]<-betas[1] + orgMu-newMu ##reset intercept to original value
            names(fi)<-c("Intercept",probesPrensent)
            betas<-fi
            return(betas)
        }

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

        ##Check if all probe-sets needed are in fixedExpressionData...
        probesPrensent<-getProbeNames(classifierParameters)[ which(getProbeNames(classifierParameters)%in%rownames(fixedExpressionData))]
        probesAbsent  <-getProbeNames(classifierParameters)[ which(!getProbeNames(classifierParameters)%in%rownames(fixedExpressionData))]

        if (length(probesAbsent)>0) {
            if (length(probesPrensent)==0){
                stop("No classifier probesets found in provided data for ",getName(classifierParameters))
            }
            if (allow.reweighted & weightingType == "complete" & "allow.reweighted"%in%names(eventChain)) {
                if (as.numeric(eventChain["allow.reweighted"])>0){
                trainingData<-getTrainingData(classifierParameters)
                    if (!is.null(trainingData) ) { ##Training set is given
                        betas<-.reweight(
                            expressionSetTrain = trainingData,
                            betas = c(getIntercept(classifierParameters),getWeights(classifierParameters)),
                            probesPrensent,probesAbsent
                        )
                        classifierParameters<-reWeightClassifier(  classifierParameters,probeNames=probesPrensent,intercept=betas[1],weights=betas[-1,drop=FALSE])

                        weightingType = "reweighted"
                        warning("Using a reweighted classifier for ",
                        getName(classifierParameters),
                        " because ", length(probesAbsent),
                        " of " , length(probesPrensent)+length(probesAbsent),
                        " covariates are missing"
                        )
                    } else {
                        warning("Cannot create a reweighted classifier ",
                        getName(classifierParameters),
                        " because the trainingset is unknown. ",length(probesAbsent),
                        " of " , length(probesPrensent)+length(probesAbsent),
                        " covariates are missing."
                        )
                    }
                }
            } else if (weightingType == "complete") {
                if ("allow.reweighted"%in%names(eventChain)) {
                    if (as.numeric(eventChain["allow.reweighted"])>0){
                        stop("Found ",length(probesAbsent),
                        " of " , length(probesPrensent)+length(probesAbsent),
                        " covariates to be missing for classifier ", getName(classifierParameters),
                        ". Consider using 'allow.reweighted = TRUE'"
                        )
                    }
                }
                stop("Found ", length(probesAbsent),
                " of " , length(probesPrensent)+length(probesAbsent),
                " covariates to be missing for classifier ", getName(classifierParameters),
                ". Reweighting is not allowed for this classifier"
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
            meanValues<-getMeans(classifierParameters)
            if (do.batchcorrection) {
                meanValues<-apply(exprs(fixedExpressionData),1,mean,na.rm=TRUE)
            }
            fixedExpressionData<-addTransformationProcess(fixedExpressionData,"SHIFTPERPROBE",values=-1*meanValues)
        } else if (keyword == "to.unitvariance" & keyvalue>0 ){
            sdValues<-getSds(classifierParameters)
            if (do.batchcorrection) {
                sdValues<-apply(exprs(fixedExpressionData),1,sd)
            }
            fixedExpressionData<-addTransformationProcess(fixedExpressionData,"MULTIPLYPERPROBE",values=1/sdValues)
        }
        else if (keyword == "to.referencemeanvar" & keyvalue>0 ) {
            if (do.batchcorrection) {
                meanValues<-apply(exprs(fixedExpressionData),1,mean,na.rm=TRUE)
                sdValues<-apply(exprs(fixedExpressionData),1,sd,na.rm=TRUE)
                refMeanValues<-getMeans(classifierParameters)[rownames(fixedExpressionData)]
                refSdValues<-getSds(classifierParameters)[rownames(fixedExpressionData)]
                fixedExpressionData<-addTransformationProcess(fixedExpressionData,"SHIFTPERPROBE",values=-meanValues)
                fixedExpressionData<-addTransformationProcess(fixedExpressionData,"MULTIPLYPERPROBE",values=refSdValues/sdValues)
                fixedExpressionData<-addTransformationProcess(fixedExpressionData,"SHIFTPERPROBE",values=refMeanValues)
            }
        }
        else if (keyword == "allow.reweighted") {}
        else {
            stop("Encountered unknown keyword (",keyword,") in the eventChain of ", getName(classifierParameters))
        }
        }
        score<-classifierParameters@doRun(classifierParameters,fixedExpressionData)
        classifierResults<-new("ClassifierResults",classifierParameters=classifierParameters,score=score,batchCorrection=do.batchcorrection,weightingType=weightingType)
        return(classifierResults)

    }
)


