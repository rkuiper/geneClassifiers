test_that("getNormalizationMethods() return valid methods",{
    expect_equal(getNormalizationMethods(),c("MAS5.0","GCRMA"));
}
)


test_that("showClassifierList returns valid matrix (with and without argument)",{
    ok<-TRUE;
    cl<-showClassifierList();
    ok<-ok & is.matrix(cl)
    if (ok) {
        ok<- ok & length(intersect(c("name","normalizationMethod","description"),colnames(cl)))==3
    }
    if (ok) {
      cl<-showClassifierList("MAS5.0");
      ok<- ok & length(intersect(c("name","normalizationMethod","description"),colnames(cl)))==3
    }
    expect_true(ok)
})




test_that("getClassifier function with signture ('character') must return object of class 'ClassifierParameters'",{
    myClassifier<-getClassifier("EMC92")
    expect_true(is(myClassifier,"ClassifierParameters"))
})

test_that("setNormalizationMethod function should return an object of 'FixedExpressionData'",{
    eset<-Biobase::ExpressionSet(matrix(seq(100)+1e6,ncol=10)) 
    eset.fixed <- setNormalizationMethod(eset, "MAS5.0",targetValue=500)
    expect_true(is(eset.fixed,"FixedExpressionData"))
})

test_that("setNormalizationMethod must receive an ExpressionSet",{
    mat<-matrix(seq(100)+1e6,ncol=10)
    expect_error(is(setNormalizationMethod(mat, "MAS5.0"),"to be an object of type ExpressionSet"))
})




