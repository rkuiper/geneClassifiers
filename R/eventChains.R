##The event chain is a nested list:
##The outer list names indicate the normalization method that is assumed to be used and must be one of : 'mas5', 'gcrma', 'rma', 'frma'
##The next level list names indicate the classifier name. This list contains a vector of keywords indicating what processes must be applied to the eset.
##Keywords: 
#    'targetValue' (Must be given for mas5 classifiers. Only used for mas5 classifiers. Indicating the required targetValue, must be applied before zscore transformation!)
#    'truncate' (Must be given. Trunctate values below a give value (applies to untransformed/unlogged data. Must be given after mas5 targetValue)
#    'to.log2' (Only used for mas5 classifiers. Log 2 transformation)
#    'allow.reweighted' (Must be given. Should reweigting be allowed in case of missing probe-set(s)?)
#    'to.meancentering' (per probe-set based centering to mean=0)
#    'to.unitvariance' (per probe-set based scaling to variance=1)
#    'to.referencemeanvar' (set probe-set means and variances to equal values as observed in training set)


#The outer list
.Eventchain<-list(
    "mas5"=list(),
    "gcrma"=list(),
    "frma"=list(),    
    "rma"=list()
);


.Eventchain$gcrma[["HM19"]]<-
    list(
        "truncate" = -Inf,
        "allow.reweighted"=TRUE,
        "to.referencemeanvar"=TRUE
    )


.Eventchain$mas5[["EMC92"]]<-
    list(
        "targetValue" = 500,
        "truncate" = -Inf,
        "allow.reweighted"=TRUE,
        "to.log2"=TRUE,
        "to.meancentering"=TRUE,
        "to.unitvariance"=TRUE    
    )


.Eventchain$mas5[["UAMS17"]]<-
    list(
        "targetValue" = 500,
        "truncate" = -Inf,
        "allow.reweighted"=TRUE,
        "to.log2"=TRUE,
        "to.meancentering"=TRUE,
        "to.unitvariance"=TRUE
    )

.Eventchain$mas5[["UAMS70"]]<-
    list(
        "targetValue" = 500,
        "truncate" = -Inf,
        "allow.reweighted"=TRUE,
        "to.log2"=TRUE,
        "to.referencemeanvar"=TRUE
    )

.Eventchain$mas5[["UAMS80"]]<-
    list(
        "targetValue" = 500,
        "truncate" = -Inf,
        "allow.reweighted"=TRUE,
        "to.log2"=TRUE,
        "to.referencemeanvar"=TRUE
    )


.Eventchain$mas5[["IFM15"]]<-
    list(
        "targetValue" = 500,
        "truncate" = -Inf,
        "to.log2"=TRUE,
        "to.meancentering"=TRUE,
        "to.unitvariance"=TRUE,
        "allow.reweighted"=FALSE    
    )


.Eventchain$mas5[["MILLENNIUM100"]]<-
    list(
        "targetValue" = 500,
        "truncate" = -Inf,
        "allow.reweighted"=TRUE,
        "to.log2"=TRUE,
        "to.meancentering"=TRUE,
        "to.unitvariance"=TRUE    
    )

.Eventchain$mas5[["MRCIX6"]]<-
    list(
        "targetValue" = 500,
        "truncate" = -Inf,
        "to.log2"=TRUE,
        "to.referencemeanvar"=TRUE,
        "allow.reweighted"=FALSE
    )






.doCheckEventChains<-function(){
    if (!all(names(.Eventchain)%in%c("mas5",  "gcrma" ,"frma" , "rma" ))){
        sublevelDiff<-setdiff(names(.Eventchain),c("mas5",  "gcrma" ,"frma" , "rma" ))
        stop("Eventchain structure error: Found invalid sublevel names (",paste(sublevelDiff,collapse=","),"). Only mas5, gcrma, frma and rma are recognized as valid sublevels!");
    }
    
    allowedKeyWords<-c("targetValue","truncate","to.log2","allow.reweighted","to.meancentering","to.unitvariance","to.referencemeanvar")
    
    for (type in names(.Eventchain)    ){
        for (item in names(.Eventchain$mas5)){
            if (!all(names(.Eventchain[[type]][[item]])%in%allowedKeyWords)){
                diff<-setdiff(names(.Eventchain[[type]][[item]]),allowedKeyWords);
                stop("Eventchain structure error: Unrecognized keyword(s) in ",type," ", item,": ",paste(diff,collapse=","));
            }    
        }
    }

    #'targetValue' (Must be given for mas5 classifiers. Only used for mas5 classifiers. Indicating the required targetValue, must be applied before zscore transformation!)
    for (type in names(.Eventchain)    ){
        for (item in names(.Eventchain$mas5)){
            if (type=="mas5") {
                if (!"targetValue"%in%names(.Eventchain[[type]][[item]])){stop("Eventchain structure error: keyword 'targetValue' must be given for a mas5 classifier. Error in ",type," ",item)}
                if (!is.numeric(.Eventchain[[type]][[item]][["targetValue"]])){stop("Eventchain structure error: keyword 'targetValue' must have a numeric value. Error in ",type," ",item)}
                if (.Eventchain[[type]][[item]][["targetValue"]]<1){stop("Eventchain structure error: keyvalue 'targetValue' must be a positive real number. Error in ",type," ",item)}
            } else {
                if ("targetValue"%in%names(.Eventchain[[type]][[item]])){stop("Eventchain structure error: keyword 'targetValue' must be given for a mas5 classifiers only. Error in ",type," ",item)}
            }
        }
    }
    #'truncate' (Must be given. Trunctate values below a give value (applies to untransformed/unlogged data. Must be given after mas5 targetValue)
    for (type in names(.Eventchain)    ){
        for (item in names(.Eventchain[[type]])){
            if (!"truncate"%in%names(.Eventchain[[type]][[item]])){stop("Eventchain structure error: keyword 'truncate' must be given. Error in ",type," ",item)}
            if (!is.numeric(.Eventchain[[type]][[item]][["truncate"]])){stop("Eventchain structure error: keyword 'truncate' must have a numeric value. Error in ",type," ",item)}
            if ("targetValue"%in%names(.Eventchain[[type]][[item]])){
                pos1<-which(names(.Eventchain[[type]][[item]])=="targetValue")
                pos2<-which(names(.Eventchain[[type]][[item]])=="truncate")
                if (pos2<pos1) {stop("Eventchain structure error: keyword 'truncate' must be given after 'targetValue'. Error in ",type," ",item)}
            }
        }
    }
    #'to.log2' (Only used for mas5 classifiers. Log 2 transformation)
    for (type in names(.Eventchain)    ){
        for (item in names(.Eventchain$mas5)){
            if (type=="mas5") {
                if (!"to.log2"%in%names(.Eventchain[[type]][[item]])){stop("Eventchain structure error: keyword 'to.log2' must be given for a mas5 classifier. Error in ",type," ",item)}
                if (!is.logical(.Eventchain[[type]][[item]][["to.log2"]])){stop("Eventchain structure error: keyword 'to.log2' must be a logical TRUE or FALSE. Error in ",type," ",item)}
            } else {
                if ("to.log2"%in%names(.Eventchain[[type]][[item]])){stop("Eventchain structure error: keyword 'to.log2' must be given for a mas5 classifiers only. Error in ",type," ",item)}
            }
        }
    }

    #'allow.reweighted' (should reweigting be allowed in case of missing probe-set(s)?)
    for (type in names(.Eventchain)    ){
        for (item in names(.Eventchain[[type]])){
            if ("allow.reweighted"%in%names(.Eventchain[[type]][[item]])){
                if (!is.logical(.Eventchain[[type]][[item]][["allow.reweighted"]])){stop("Eventchain structure error: keyword 'allow.reweighted' must be a logical TRUE or FALSE. Error in ",type," ",item)}
            } else {
                stop("Eventchain structure error: keyword 'allow.reweighted' must be given. Error in ",type," ",item)
            }
        }
    }
    #'to.meancentering' (per probe-set based centering to mean=0)
    for (type in names(.Eventchain)    ){
        for (item in names(.Eventchain[[type]])){
            if ("to.meancentering"%in%names(.Eventchain[[type]][[item]])){
                if (!is.logical(.Eventchain[[type]][[item]][["to.meancentering"]])){stop("Eventchain structure error: keyword 'to.meancentering' must be a logical TRUE or FALSE. Error in ",type," ",item)}
            } 
        }
    }
    #'to.unitvariance' (per probe-set based centering to mean=0)
    for (type in names(.Eventchain)    ){
        for (item in names(.Eventchain[[type]])){
            if ("to.unitvariance"%in%names(.Eventchain[[type]][[item]])){
                if (!is.logical(.Eventchain[[type]][[item]][["to.unitvariance"]])){stop("Eventchain structure error: keyword 'to.unitvariance' must be a logical TRUE or FALSE. Error in ",type," ",item)}
            } 
        }
    }
    #'to.referencemeanvar' (set probe-set means and variances to equal values as observed in training set)
    for (type in names(.Eventchain)    ){
        for (item in names(.Eventchain[[type]])){
            if ("to.referencemeanvar"%in%names(.Eventchain[[type]][[item]])){
                if (!is.logical(.Eventchain[[type]][[item]][["to.referencemeanvar"]])){stop("Eventchain structure error: keyword 'to.referencemeanvar' must be a logical TRUE or FALSE. Error in ",type," ",item)}
            } 
        }
    }
}
