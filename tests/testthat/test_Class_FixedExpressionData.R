test_that("getNormalizationMethods() return valid methods",{
    expect_equal(length(intersect(c("MAS5.0","GCRMA"),getNormalizationMethods())),2);
}
)

test_that("Checking FixedExpressionData Class",{
    expect_true(is(new("FixedExpressionData"),"FixedExpressionData")) ##Empty unusable object    
})


test_that("Creating FixedExpressionData with constructor function signature(normalizationMethod='character',expressionMatrix='matrix')",{
    mat<-matrix(rep("A",100),ncol=10)
    expect_error(FixedExpressionData("MAS5",mat),"must be a numeric matrix")
    mat<-matrix(rep(1e16,100),ncol=10)
    expect_error(FixedExpressionData("MAS5",mat),"must have row names")
    expect_error(FixedExpressionData("MAS5",mat),"must have column names")
    mat<-matrix(rep(1e16,100),ncol=10);rownames(mat)<-rep("A",nrow(mat))
    expect_error(FixedExpressionData("MAS5",mat),"must have column names")
    mat<-matrix(rep(1e16,100),ncol=10);colnames(mat)<-rep("A",ncol(mat))
    expect_error(FixedExpressionData("MAS5",mat),"must have row names")

    mat<-matrix(rep(1e16,100),ncol=10);colnames(mat)<-rep("A",ncol(mat));rownames(mat)<-rep("A",nrow(mat))
    expect_error(FixedExpressionData("MAS5",mat),"duplicated row names found")
    expect_error(FixedExpressionData("MAS5",mat),"duplicated column names found")

    mat<-matrix(rep(1e16,100),ncol=10);colnames(mat)<-letters[1:10];rownames(mat)<-rep("A",nrow(mat))
    expect_error(FixedExpressionData("MAS5",mat),"duplicated row names found")

    mat<-matrix(rep(1e16,100),ncol=10);colnames(mat)<-rep("A",ncol(mat));rownames(mat)<-letters[1:10]
    expect_error(FixedExpressionData("MAS5",mat),"duplicated column names found")
    
    mat<-matrix(rep(1e16,100),ncol=10);colnames(mat)<-letters[1:10];rownames(mat)<-letters[1:10];mat[3]<-NA
    expect_error(FixedExpressionData("MAS5",mat),"NA values found in expression matrix")
   
    mat<-matrix(rep(1,100),ncol=10);colnames(mat)<-letters[1:10];rownames(mat)<-letters[1:10]
    expect_warning({
        eset.fixed<-FixedExpressionData("MAS5",mat)
        expect_true(is(eset.fixed,"FixedExpressionData"))
    },"The data seems to be log2 transformed")
    
    mat<-matrix(rep(1,100)+1e6,ncol=10);colnames(mat)<-letters[1:10];rownames(mat)<-letters[1:10]
    expect_warning({
        eset.fixed<-FixedExpressionData("MAS5",mat,isLog2Transformed=TRUE)
        expect_true(is(eset.fixed,"FixedExpressionData"))
    },"Reconsider argument isLog2Transformed = TRUE")
    
    mat<-matrix(rep(1,100),ncol=10);colnames(mat)<-letters[1:10];rownames(mat)<-letters[1:10]
    expect_warning({
        FixedExpressionData("MAS5",mat,isLog2Transformed=FALSE)
        expect_true(is(eset.fixed,"FixedExpressionData"))
    },"Reconsider argument isLog2Transformed = FALSE")
  
    mat<-matrix(rep(1,100)+1e6,ncol=10);colnames(mat)<-letters[1:10];rownames(mat)<-letters[1:10]
    expect_warning(FixedExpressionData("MAS5",mat,targetValue=500),"Using the given targetValue, but it does not correspond")

    mat<-matrix(rep(1,100)+1e6,ncol=10);mat[,2]<-1e5;colnames(mat)<-letters[1:10];rownames(mat)<-letters[1:10]
    expect_error(FixedExpressionData("MAS5",mat )  ,"Cannot determine the targetValue that was used for MAS5.0" )
   
    mat<-matrix(rep(1,100)+1e6,ncol=10);mat[,2]<-1e5;colnames(mat)<-letters[1:10];rownames(mat)<-letters[1:10]
    expect_true(is(FixedExpressionData("MAS5",mat,targetValue=100 ) ,"FixedExpressionData"))

    mat<-matrix(rep(1,10)+1e6,ncol=1);colnames(mat)<-letters[1];rownames(mat)<-letters[1:10]
    expect_warning(FixedExpressionData("MAS5",mat,targetValue=100 ),"Using the given targetValue, but it does not correspond")

    mat<-matrix(rep(1,100),ncol=10);colnames(mat)<-letters[1:10];rownames(mat)<-letters[1:10]
    expect_true(is(FixedExpressionData("MAS5",mat,isLog2Transformed=TRUE),"FixedExpressionData"))
    
    mat<-matrix(rep(1,100)+1e6,ncol=10);colnames(mat)<-letters[1:10];rownames(mat)<-letters[1:10]
    expect_true(is(FixedExpressionData("MAS5",mat),"FixedExpressionData"))
    
}
)


test_that("getTargetValue with signature('FixedExpressionData')",{
    mat<-matrix(rep(1,100)+1e6,ncol=10);colnames(mat)<-letters[1:10];rownames(mat)<-letters[1:10]
    eset.fixed<-FixedExpressionData("MAS5",mat)
    expect_equal(getTargetValue(eset.fixed),1e6)

    expect_warning({
        mat<-matrix(rep(1,100)+1e6,ncol=10);colnames(mat)<-letters[1:10];rownames(mat)<-letters[1:10]
        eset.fixed<-FixedExpressionData("MAS5",mat,targetValue=500)
        expect_equal(getTargetValue(eset.fixed),500)
    })

})


test_that("setTargetValue with signature('FixedExpressionData')",{
    expect_equal({
        mat<-matrix(rep(1,100)+1e6,ncol=10);colnames(mat)<-letters[1:10];rownames(mat)<-letters[1:10]
        eset.fixed<-FixedExpressionData("MAS5",mat)
        setTargetValue(eset.fixed)<-100
        getTargetValue(eset.fixed)}, 100)

    expect_error({
        mat<-matrix(rep(1,100)+1e6,ncol=10);colnames(mat)<-letters[1:10];rownames(mat)<-letters[1:10]
        eset.fixed<-FixedExpressionData("GCRMA",mat)
        setTargetValue(eset.fixed)<-100
        }, "Function 'setTargetValue' does only apply to MAS5.0 normalized data")
})


test_that("getExpressionEnvironment with signature('FixedExpressionData')",{
   a <- expect_type({
        mat<-matrix(rep(1,100)+1e6,ncol=10);colnames(mat)<-letters[1:10];rownames(mat)<-letters[1:10]
        eset.fixed<-FixedExpressionData("MAS5",mat)
        getExpressionEnvironment(eset.fixed)
    }, "environment")

    
    expect_true("expressionMatrix"%in%ls(getExpressionEnvironment(eset.fixed)))
})

test_that("exprs function with signature('FixedExpressionData')",{
    a <- expect_true({
        mat<-matrix(rep(1,100)+1e6,ncol=10);colnames(mat)<-letters[1:10];rownames(mat)<-letters[1:10]
        eset.fixed<-FixedExpressionData("MAS5",mat)
        is(exprs(eset.fixed),"matrix")
    })
 
})



test_that("rawExprs function with signature('FixedExpressionData')",{
    a <- expect_true({
        mat<-matrix(rep(1,100)+1e6,ncol=10);colnames(mat)<-letters[1:10];rownames(mat)<-letters[1:10]
        eset.fixed<-FixedExpressionData("MAS5",mat)
        is(rawExprs(eset.fixed),"matrix")
    })
 
})


test_that("dim for FixedExpressionData",{
    expect_equal({
        mat<-matrix(rep(1,50)+1e6,ncol=5);colnames(mat)<-letters[1:5];rownames(mat)<-letters[1:10]
        eset.fixed<-FixedExpressionData("MAS5",mat)
        dim(eset.fixed)
    },c(10,5))
})


test_that("rownames for FixedExpressionData",{
    expect_equal({
        mat<-matrix(rep(1,50)+1e6,ncol=5);colnames(mat)<-letters[1:5];rownames(mat)<-letters[1:10]
        eset.fixed<-FixedExpressionData("MAS5",mat)
        rownames(eset.fixed)
    },letters[1:10])
})


test_that("colnames for FixedExpressionData",{
    expect_equal({
        mat<-matrix(rep(1,50)+1e6,ncol=5);colnames(mat)<-letters[1:5];rownames(mat)<-letters[1:10]
        eset.fixed<-FixedExpressionData("MAS5",mat)
        colnames(eset.fixed)
    },letters[1:5])
})



test_that("Extract Parts of an FixedExpressionData using '['",{

    mat<-matrix(seq(50)+1e6,ncol=5);colnames(mat)<-letters[1:5];rownames(mat)<-letters[1:10]
    eset.fixed<-FixedExpressionData("MAS5",mat,targetValue=100)
    expect_equal(exprs(eset.fixed[2:3,3:5]),mat[2:3,3:5])
    expect_equal(exprs(eset.fixed[2:3,]),mat[2:3,])
    expect_equal(exprs(eset.fixed[,3:4]),mat[,3:4])
    expect_equal(exprs(eset.fixed[2,3]),mat[2,3,drop=FALSE])

})


test_that("Extract Parts of an FixedExpressionData using '[['",{
    mat<-matrix(seq(50)+1e6,ncol=5);colnames(mat)<-letters[1:5];rownames(mat)<-letters[1:10]
    eset.fixed<-FixedExpressionData("MAS5",mat,targetValue=100)
    expect_equal(exprs(eset.fixed[[2:3,3:5]]),mat[2:3,3:5])
    expect_equal(exprs(eset.fixed[[2:3,]]),mat[2:3,])
    expect_equal(exprs(eset.fixed[[,3:4]]),mat[,3:4])
    expect_equal(exprs(eset.fixed[[2,3]]),mat[2,3,drop=FALSE])

})



test_that("Test explicitlyChangeExprs on a FixedExpressionData",{
    mat<-matrix(seq(50)+1e6,ncol=5);colnames(mat)<-letters[1:5];rownames(mat)<-letters[1:10]
    eset.fixed<-FixedExpressionData("MAS5",mat,targetValue=100)
    mat.new<-matrix(1,ncol=1);colnames(mat.new)<-"C";rownames(mat.new)<-"D"
    expect_equal({
        explicitlyChangeExprs(eset.fixed)<-mat.new
        exprs(eset.fixed)
    },mat.new)

    mat.new<-matrix(seq(50)+1e5,ncol=5);colnames(mat.new)<-letters[1:5];rownames(mat.new)<-letters[1:10]
        
    expect_equal({
        explicitlyChangeExprs(eset.fixed)<-mat.new
        exprs(eset.fixed)
    },mat.new)

})
 
