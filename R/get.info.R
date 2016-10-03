get.info<-function(classifier){
    classifier<-match.arg(classifier,get.classifiers())
    params<-.getParams(classifier)
    params$name<-classifier
    class(params)<-"classifier.info"
    params
}

print.classifier.info<-function(x,...){
    cat(rep("-",15),"\n",sep='')
    cat("Name : " , x$name,"\n")


    normalization<-NULL;
    normMethods<-names(.Eventchain)
    for (normMethod in normMethods){
        if (x$name%in%get.classifiers(normMethod)) {normalization<-normMethod;}
    }

    cat("Normalization : " , normalization,"\n")

    cat("preprocessing : \n" )
    print(unlist( .Eventchain[[normalization]][[x$name]]));

    weights<-x$weights
    weights<-weights[order(abs(weights),decreasing=TRUE)]
    lw<-length(weights)
    if (lw>5){ weights<-weights[1:5,drop=FALSE]; }
    weights<-paste(names(weights),signif(weights,3),sep=" : ")
    weights<- paste(weights,collapse="\n\t")
    if (lw>5){
        
        weights<-paste(c(weights,"\n\t...",lw-5," more\n"),collapse="")
    }
    cat("Weights:",weights)
    
    cat("\nCitation : " , x$citation,"\n")
    cat(rep("-",15),"\n\n",sep='')
}
