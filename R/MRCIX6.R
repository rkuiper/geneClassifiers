
.classifierHookList[[length(.classifierHookList)+1]]<-
function(infoOnly=TRUE,...){
    name="MRCIX6"
    normalizationMethod="MAS5.0"
    description="A risk classifier for multiple myeloma"
    hasTrainingData = TRUE
    #training set MAS5.0 log2 means
    means<-c(
        8.02738041819488,10.5981054650442,
        8.36951781518735,12.4309125118577,
        12.2597337824519,14.8329028643051
    )
    #training set MAS5.0 log2 sds
    sds<-c(
        1.47608610664399,0.4624421847404,
        1.52422557089425,0.747362639863772,
        0.62481474641846,0.79127962225654
    )
    weights<- rep(1,6)
    names(weights)<-c(
        "203755_at","216326_s_at",
        "203213_at","218034_at",
        "200608_s_at","217731_s_at"
    )
    intercept<-0
    decisionBoundaries<-0.999
    citations<-list("Dickens NJ, Walker BA, Leone PE, et al. Homozygous deletion mapping in myeloma samples identifies genes and an expression signature relevant to pathogenesis and outcome. Clin Cancer Res. 2010;16(6):1856-1864.")
    doRun=function(thisClassifier,fixedExpressionData){
        data<-t(exprs(fixedExpressionData)[getProbeNames(thisClassifier),])
        score<-0
        if (all(c("203755_at","216326_s_at")%in%colnames(data))){
            score<-score+1*((data[,"203755_at"]-data[,"216326_s_at"]) > 0)
        }
        if (all(c("203213_at","218034_at")%in%colnames(data))){
            score<-score+1*((data[,"203213_at"]-data[,"218034_at"]) > 0)
        }
        if (all(c("200608_s_at","217731_s_at")%in%colnames(data))){
            score<-score+1*((data[,"200608_s_at"]-data[,"217731_s_at"]) > 0)
        }
        names(score)<-rownames(data)
        return(score)
    }
    eventChain=list(
        "targetValue" = 500,
        "truncate" = -Inf,
        "to.log"=2,
        "to.referencemeanvar"=TRUE,
        "allow.reweighted"=FALSE
    )
    if(infoOnly) {
        return( c( "name"=name ,"normalizationMethod"=  normalizationMethod ,"description"=description) )
    }
    ClassifierParameters(
        name=name,
        description=description,
        normalizationMethod=normalizationMethod,
        weights=weights,
        intercept=intercept,
        means = means,
        sds = sds,
        decisionBoundaries=decisionBoundaries,
        doRun=doRun,
        citations = citations,
        eventChain=eventChain,
        hasTrainingData=hasTrainingData
    )
}
