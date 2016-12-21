.classifierHookList[[length(.classifierHookList)+1]]<-
function(infoOnly=TRUE,...){
    name="UAMS17"
    normalizationMethod="MAS5.0"
    description="A risk classifier for multiple myeloma"
    hasTrainingData = TRUE

    weights<-c(
        0.283,-0.296,-0.208,+0.314,-0.287,
        +0.251,+0.193,+0.269,+0.375,+0.158,
        +0.316,+0.232,-0.251,-0.23,-0.402,
        +0.191,+0.148
    )
    names(weights)<-c(
        "200638_s_at","1557277_a_at","200850_s_at",
        "201897_s_at","202729_s_at","203432_at",
        "204016_at","205235_s_at","206364_at",
        "206513_at","211576_s_at","213607_x_at",
        "213628_at","218924_s_at","219918_s_at",
        "220789_s_at","242488_at"
    )
    intercept<-0
    #means are based on the mas5 log2 data
    means<-c(
        12.7066030077155,5.97008377969383,11.2490695936341,
        10.097036288508,10.2405530788169,10.1236916286819,
        6.94359447318795,8.63798709271333,8.03664149238594,
        11.2162558849659,7.47121611222391,8.91432869887769,
        10.9743537880763,10.6465992818962,8.61438301078837,
        8.05996950086986,7.09081971226595
    )
    #sds are based on the mas5 log2 data
    sds<-c(
        0.536455858686746,1.31675763096059,0.488914867823306,
        0.941803268650604,0.960269151350048,0.662456898945665,
        1.21502520635941,0.811635628837525,1.00113287806509,
        0.91766256010483,1.58684159508562,0.99951858573628,
        0.47602844762267,0.788790892324437,1.76361960133134,
        1.31800408948714,1.35791489054738
    )
    decisionBoundaries<-1.5
    citations<-list("John D. Shaughnessy, Fenghuang Zhan, Bart E. Burington, et. al.; Blood Mar 2007, 109 (6) 2276-2284; DOI: 10.1182/blood-2006-07-038430")
    doRun=function(thisClassifier,fixedExpressionData){
        e<-exprs(fixedExpressionData)[getProbeNames(thisClassifier),]
        (t(getWeights(thisClassifier))%*%e+getIntercept(thisClassifier))[1,]
    }
    eventChain=list(
        "targetValue" = 500,
        "truncate" = -Inf,
        "allow.reweighted"=TRUE,
        "to.log"=2,
        "to.meancentering"=TRUE,
        "to.unitvariance"=TRUE
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
