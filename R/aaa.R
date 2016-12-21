#' @importFrom utils packageVersion packageName
packageVersionInternal<-function(){
    ##To handle a roxygen run
    packageName = packageName()
    if (packageName == "roxygen_devtest") { return(package_version("0.0.0")) }
    packageVersion(packageName())
}
