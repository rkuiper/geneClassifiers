.classifierHookList[[length(.classifierHookList)+1]]<-
function(infoOnly=TRUE,...){
    name="HM19"
    normalizationMethod="GCRMA"
    description="A risk classifier for multiple myeloma"
    hasTrainingData = TRUE
    #means are based on the mas5 log2 data
    means<-c(
        5.80843379667409,9.23160864112744,6.12584428860212,
        3.09461289710876,2.41294266158072,6.95247713002351,
        4.29846110127437,2.68195699105932,3.26351320694089,
        9.11614329187661,2.64178820366459,5.52467647461559,
        2.76447936071926,4.51820607118098,4.02797658562296,
        10.3711591498012,7.87487883163514,12.320981823562,
        8.66750196731196
    )
    #sds are based on the mas5 log2 data
    sds<-c(
        2.1538620919844,0.83775519834133,1.16525794621433,
        1.53141604660456,0.537115117544859,1.89816898681006,
        1.76958645739875,0.876960963486335,0.846985109747,
        1.41175412673134,1.25337802593257,1.80702259705485,
        1.06392790859626,1.40660194441016,1.08807945657841,
        0.908833675641071,1.53827225084067,0.813580000704891,
        2.03505832355478
    )
    weights<-c(rep(1,15), rep(-1,4))
    names(weights)<-c(
    '203755_at'  , '218460_at' , '226936_at' , '219855_at'  , '233660_at' ,
    '203358_s_at', '203764_at' , '218726_at' , '221520_s_at','234672_s_at',
    '214464_at'  , '226980_at' , '225687_at' , '229553_at'  ,'219978_s_at',
    '225272_at'  , '235353_at' , '204031_s_at','220945_x_at'
    )
    intercept<-0
    decisionBoundaries<-c(28.4,54.6)
    citations<-list("Reme T, Hose D, Theillet C, Klein B. Modeling risk stratification in human cancer. Bioinformatics. 2013;29(9):1149-1157")
    doRun=function(thisClassifier,fixedExpressionData){
        e<-exprs(fixedExpressionData)[getProbeNames(thisClassifier),]
        (t(getWeights(thisClassifier))%*%e+getIntercept(thisClassifier))[1,]
    }
    eventChain=list(
            "truncate" = -Inf,
            "allow.reweighted"=TRUE,
            "to.referencemeanvar"=TRUE
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
