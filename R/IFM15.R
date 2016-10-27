.classifierHookList[[length(.classifierHookList)+1]]<-
function(infoOnly=TRUE,...){
    name="IFM15"
    normalizationMethod="MAS5.0"
    description="A risk classifier for multiple myeloma"
    probeNames<-c(
        "209683_at"  ,"200779_at"  ,"203657_s_at" ,"201425_at" ,
        "217752_s_at","200783_s_at","202486_at"   ,"202951_at" ,
        "208644_at"  ,"202470_s_at","212098_at"   ,"228737_at" ,
        "204072_s_at","228677_s_at","1565162_s_at","231736_x_at"
    )
    #training set means are unknown
    #means<-
    #training set sds are unknown
    #sds<-
    weights<-c(
        -0.29483393,-0.2491429,-0.17822457,-0.13137903,
        0.38697337, 0.31490195, 0.30371178, 0.29530369,
        0.27578783, 0.26987655, 0.25043791, 0.21956366,
        0.21255699, 0.19243758, 0.08886402, 0.08886402
    )
    intercept<-0
    decisionBoundaries<-0.956
    citations<-list("Decaux O, Lode L, Magrangeas F, et al. Prediction of survival in multiple myeloma based on gene expression profiles reveals cell cycle and chromosomal instability signatures in high-risk patients and hyperdiploid signatures in low-risk patients: a study of the Intergroupe Francophone du Myelome. J Clin Oncol. 2008;26(29):4798-4805")
    doRun=function(thisClassifier,fixedExpressionData){
        e<-exprs(fixedExpressionData)[getProbeNames(thisClassifier),]
        (t(getWeights(thisClassifier))%*%e+getIntercept(thisClassifier))[1,]
    }
    eventChain=list(
        "targetValue" = 500,
        "truncate" = -Inf,
        "to.log"=2,
        "to.meancentering"=TRUE,
        "to.unitvariance"=TRUE,
        "allow.reweighted"=FALSE
         )
    if(infoOnly) {
        return( c( "name"=name ,"normalizationMethod"=  normalizationMethod ,"description"=description) )
    }
    new("ClassifierParameters",
        name=name,
        description=description,
        normalizationMethod=normalizationMethod,
        weights=weights,
        intercept=intercept,
        probeNames=probeNames,
#        means = means,
#        sds = sds,
        decisionBoundaries=decisionBoundaries,
        doRun=doRun,
        citations = citations,
        eventChain=eventChain
    )
}
