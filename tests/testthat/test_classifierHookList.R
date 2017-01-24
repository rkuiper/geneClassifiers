
##All elements in .classifierHookList are of type ClassifierParameters?
test_that("returns ClassifierParameters Object",{
    isClassifierParametersObject<-lapply(.classifierHookList,function(x){
        inherits(x(FALSE),"ClassifierParameters")
    })

    allOK<-all(isClassifierParametersObject)
    expect_true(allOK)
})


#** Test info retrieval
test_that("classifierHookList Info",{
    structureOK<-all(unlist(lapply(geneClassifiers:::.classifierHookList,function(x){all(names(x())==c("name","normalizationMethod","description"))})))
    expect_equal(structureOK , TRUE)
})

test_that("functionDefs for hook",{
    #Test: classifierHooks have the correct function defintion
     funcDefCorrect<-unlist(lapply(.classifierHookList,function(x){
        fml<-formals(x)
        return(all(names(fml)==c("infoOnly","...")) & all(unlist(fml)==c(TRUE,"")))
    }))
    allOK<-all(funcDefCorrect)
    expect_true(allOK)
})


test_that("classifierHookList argumentReturnOrder",{
    registeredClassifiers<-lapply(.classifierHookList,function(x){x()})
    argumentReturnOrder<-unlist(lapply(registeredClassifiers,function(x){
        all(names(x)==c("name","normalizationMethod","description"))
    }))
    allOK<-all(argumentReturnOrder)
    expect_true(allOK)
})

test_that("classifierHookList namesCorrect",{
    registeredClassifiers<-lapply(.classifierHookList,function(x){x()})
    namesCorrect<-unlist(lapply(registeredClassifiers,function(x){
        is.character(x["name"]) & nchar(x["name"])>2
    }))
    allOK<-all(namesCorrect)
    expect_true(allOK)
})

test_that("classifierHookList normalizationMethodCorrect",{
    registeredClassifiers<-lapply(.classifierHookList,function(x){x()})
    normalizationMethodCorrect<-unlist(lapply(registeredClassifiers,function(x){
        x["normalizationMethod"]%in%getNormalizationMethods()
    }))
    allOK<-all(normalizationMethodCorrect)
    expect_true(allOK)
})

#** Test classifier retrieval





allowedKeyWords<-c("targetValue","truncate","to.log","allow.reweighted","to.meancentering","to.unitvariance","to.referencemeanvar")
##Loop through all classifiers
lapply(.classifierHookList,function(x){
    classifier<-x(FALSE)
    nm<-getNormalizationMethod(classifier)
    
    ec<-getEventChain(classifier)
    probeNames<-getProbeNames(classifier)
    
    #$ 'only allowed keywords in eventchain?
    testdescription<-paste("Eventchain keywords for",getName(classifier))
    test_that(testdescription,{
        correctEchainKeywords<-all(names(ec)%in%allowedKeyWords)
        expect_true(correctEchainKeywords)
    }) 


    #$ 'trainingData available?
    if (hasTrainingData(classifier)) {
        testdescription<-paste("is trainingData available for",getName(classifier))
        test_that(testdescription,{
            trainingData<-try(getTrainingData(classifier),silent=TRUE)
            expect_true(!inherits(trainingData,"try_error"))
        })
    }

    
    #$ 'probeNames given?
    testdescription<-paste("Probenames given for",getName(classifier))
    test_that(testdescription,expect_true(!is.null(probeNames) & all(!is.na(probeNames))))

    #$ 'weights given?
    testdescription<-paste("Weights given for",getName(classifier))
    test_that(testdescription,expect_true(!is.null(getWeights(classifier)) & all(!is.na(getWeights(classifier)))))

    #$ 'targetValue' (Must be given for MAS5.0 classifiers only. Only used for mas5 classifiers. Indicating the required targetValue)
    if (nm == "MAS5.0") {
        testdescription<-paste("Eventchain structure error: keyword 'targetValue' must be given for a mas5 classifier (",getName(classifier),")")
        test_that(testdescription, expect_true("targetValue"%in%names(ec)))

        testdescription<-paste("Eventchain structure error: keyword 'targetValue' must have a numeric value. (",getName(classifier),")")
        test_that(testdescription, expect_true(is.numeric(ec[["targetValue"]])))

        testdescription<-paste("Eventchain structure error: keyvalue 'targetValue' must be a positive real number. (",getName(classifier),")")
        test_that(testdescription, expect_true(ec[["targetValue"]]>=1))

    } else { ##Not MAS5.0
        testdescription<-paste("Targetvalue given for non MAS5.0. (",getName(classifier),")")
        test_that(testdescription, expect_true(!"targetValue"%in%names(ec))) 
    }

  
    #$ 'truncate' (Must be given. Truncates values below a given value (applies to untransformed/unlogged data.)
    testdescription<-paste("Eventchain structure error: keyword 'truncate' must be given as numeric value. (",getName(classifier),")")
    test_that(testdescription, expect_true("truncate"%in%names(ec) & is.numeric(ec[["truncate"]])))

    #$ 'to.log' (Although probably only logical to apply in MAS5.0 normalized data, always applicable)
    if (!is.null(ec[["to.log"]])){
        testdescription<-paste("Eventchain structure error: keyword 'to.log' must be a numerical value representing the base number for the log. (",getName(classifier),")")
        test_that(testdescription, expect_true(is.numeric(ec[["to.log"]])))
    }

    #$ 'allow.reweighted' Must be given (should reweigting be allowed in case of missing probe-set(s)?)
    testdescription<-paste("Eventchain structure error: keyword 'allow.reweighted' must be a logical TRUE or FALSE. (",getName(classifier),")")
    test_that(testdescription, expect_true("allow.reweighted"%in%names(ec) & is.logical(ec[["allow.reweighted"]])))


    #$ 'to.meancentering' (per probe-set based centering to mean=0)
    if ("to.meancentering"%in%names(ec)){ 
        testdescription<-paste("Eventchain structure error: keyword 'to.meancentering' must be a logical TRUE or FALSE. (",getName(classifier),")")
        test_that(testdescription, expect_true(is.logical(ec[["to.meancentering"]])))
        ##Note that under normal circumstances (i.e. batch correction) refrence means are not used thus not required!
    }

    #$ 'to.unitvariance' (per probe-set based scaling to sd=1)
    if ("to.unitvariance"%in%names(ec)){
        testdescription<-paste("Eventchain structure error: keyword 'to.unitvariance' must be a logical TRUE or FALSE. (",getName(classifier),")")
        test_that(testdescription, expect_true(is.logical(ec[["to.unitvariance"]])))
        ##Note that under normal circumstances (i.e. batch correction) refrence sds are not used thus not required!
    }

    #$ 'to.referencemeanvar' (set probe-set means and variances to equal values as observed in training set)
    if ("to.referencemeanvar"%in%names(ec)){
         testdescription<-paste("Eventchain structure error: keyword 'to.referencemeanvar' must be a logical TRUE or FALSE. (",getName(classifier),")")
        test_that(testdescription, expect_true(is.logical(ec[["to.referencemeanvar"]])))
        
        if (length(probeNames)>0) {
            if (ec[["to.referencemeanvar"]]==TRUE) {
                testdescription<-paste("Eventchain structure error: means and sds must be given when keyword 'to.referencemeanvar' == TRUE. (",getName(classifier),")")
                test_that(testdescription, expect_true(length(getMeans(classifier))>0 & length(getMeans(classifier))==length(probeNames) & length(getSds(classifier))==length(probeNames)))
            }
        }
    }


    #$ 'allow.reweighted' = TRUE (only sensical if trainingset is given)
    if ("allow.reweighted"%in%names(ec)){
        testdescription<-paste("Eventchain structure error: keyword 'allow.reweighted' must be a logical TRUE or FALSE. (",getName(classifier),")")
        test_that(testdescription, expect_true(is.logical(ec[["allow.reweighted"]])))     
        
        if (ec[["allow.reweighted"]]==TRUE) {
                testdescription<-paste("Eventchain error: 'allow.reweighted'=TRUE only meaningful if trainingData given. (",getName(classifier),")")
                test_that(testdescription, expect_true(hasTrainingData(classifier)))
        }
    }


    testdescription<-paste("Check normalization method for",getName(classifier))
    test_that(testdescription,{
        correctNormalizationMethod<-nm%in%getNormalizationMethods()
        expect_true(correctNormalizationMethod)
    })     
    
    
    if (hasTrainingData(classifier)){
        trainingData<-getTrainingData(classifier)
        intersection<-intersect(rownames(trainingData),probeNames)
        testdescription <- paste("trainings dataset contain more (or less) probes than probeNames for classifier ",getName(classifier),".")
        test_that(testdescription,expect_true(length(intersection)==nrow(trainingData) & length(intersection)==length(probeNames)))

        testdescription<-paste("Check normalization method training data for",getName(classifier))
        preprocData<-preproc(experimentData(trainingData))
        trainDatNormMethod <- preprocData[["normalisationMethod"]]

        test_that(testdescription,{
            isOK<-FALSE;
            if (!is.null(trainDatNormMethod)) {
                isOK<-trainDatNormMethod%in%getNormalizationMethods() & trainDatNormMethod==nm
                if (nm=="MAS5.0"){ 
                    isOK<- isOK & !is.null(preprocData[["targetValue"]]) ##Provide a target value in case of MAS5.0
            
                }
            }
            expect_true(isOK)
        })
         
    }
})
 


test_that("classifierHookList namesDuplicated",{
    registeredClassifiers<-lapply(.classifierHookList,function(x){x()})
    namesDuplicated<-duplicated(unlist(lapply(registeredClassifiers,function(x){
        x["name"]
    })))
    allOK<-!any(namesDuplicated)
    expect_true(allOK)
})





