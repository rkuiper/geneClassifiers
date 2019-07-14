

test_that("Checking ClassifierResults Class",{
    expect_true(is(new("ClassifierResults"),"ClassifierResults")) ##Empty unusable object    

    expect_true( is(
        new("ClassifierResults",
        score=seq_len(10)/5,
        classifierParameters=getClassifier("EMC92"),
        batchCorrection=TRUE,
        weightingType='complete'
        ),"ClassifierResults")) 
    expect_error(
        new("ClassifierResults",
        score=seq_len(10)/5,
        classifierParameters=getClassifier("EMC92"),
        batchCorrection=TRUE,
        weightingType='WrOnG ValUe'
        ),"invalid class") #

    # auto recursive testing classes
    tmpParams<-getClassifier("EMC92")
    tmpParams@name<-c("A","B")
    expect_error(
        new("ClassifierResults",
            score=seq_len(10)/5,
            classifierParameters=tmpParams,
            batchCorrection=TRUE,
            weightingType='complete'
        ),
        "invalid class"
    )
   

    
})
 



