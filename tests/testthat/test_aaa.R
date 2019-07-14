
##packageVersionInternal must results object of class package_version.
test_that("returns package_version Object",{
    expect_true(is(packageVersionInternal(),"package_version"))
})

