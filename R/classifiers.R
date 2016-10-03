.onAttach<-function(libname, pkgname) {
    packageStartupMessage('See vignette("geneClassifiers") for help');
}
.onLoad<-function(libname, pkgname) {
    .doCheckEventChains()
}

get.classifiers<-function(normalization=c("mas5","gcrma","rma","frma")){
    normalization<-match.arg(normalization, several.ok=TRUE);
    res<-vector();
    for (norm in normalization){
        res<-c(res,names(.Eventchain[[norm]]));
    }
    res
}

applyClassifier.mas5<-function(which,eset,dataislog2=FALSE,used.MAS5TargetValue=NULL,...){
    which=match.arg(which,get.classifiers("mas5"))
    which.eventChain = .Eventchain$mas5[[which]];
    tmp<-.doEventChain(which,eset,which.eventChain,dataislog2=dataislog2,used.MAS5TargetValue=used.MAS5TargetValue,...);        
    return(tmp)
}

applyClassifier.gcrma<-function(which,eset,...){
    which=match.arg(which,get.classifiers("gcrma"))
    which.eventChain = .Eventchain$gcrma[[which]];
    result<-.doEventChain(which,eset,which.eventChain,...);
    return(result)
}
applyClassifier.frma<-function(which,eset,...){
    which=match.arg(which,get.classifiers("frma"))
    which.eventChain = .Eventchain$frma[[which]];
    result<-.doEventChain(which,eset,which.eventChain,...);
    return(result)
}


applyClassifier.rma<-function(which,eset,...){
    which=match.arg(which,get.classifiers("rma"))
    which.eventChain = .Eventchain$rma[[which]];
    result<-.doEventChain(which,eset,which.eventChain,...);
    return(result)
}

.getParams<-function(which){
    local({
        txt<-paste("data(",which,",envir=environment())",sep="")
        eval(parse(text=txt))
        eval(parse(text=paste("return(.",which,"_params)",sep="")))
    })
}

.doEventChain<-function(which,eset,eventChain,params,...){
    if ( ! inherits(x=eset,"ExpressionSet") ) {
        stop("Expect argument eset to be of class 'ExpressionSet'")
    }
    dots <- list(...)
    dataislog2 = FALSE;

    used.MAS5TargetValue = NULL
    allow.reweighted = FALSE;
    do.batchcorrection = TRUE;
    if (!is.null(dots$dataislog2)) dataislog2<-dots$dataislog2;
    if (!is.null(dots$used.MAS5TargetValue)) used.MAS5TargetValue<-dots$used.MAS5TargetValue;
    if (!is.null(dots$allow.reweighted)) allow.reweighted<-dots$allow.reweighted;
    if (!is.null(dots$do.batchcorrection)) do.batchcorrection<-dots$do.batchcorrection;    
    if (missing(params)) params<-.getParams(which)
    if (do.batchcorrection & ncol(eset)<3)  {
        stop("Only a small number of samples available (n=",ncol(eset),"). Please consider setting 'do.batchcorrection = FALSE' but be aware this could severely induce biases!");
    } 
    if (do.batchcorrection & ncol(eset)<20) {
        warning("Small number of available samples (n=",ncol(eset),"). Include more samples or consider setting 'do.batchcorrection = FALSE' but be aware this could severely induces biases!");
    } 

    weightingType = "complete"
    batchCorrection = "yes"

    ##Check if all probe-sets needed are in eset...
    invent<-.probeInventorization(eset,params);
    if (length(invent$outs)>0) {
        if (length(invent$ins)==0){ stop("No classifier probesets found in provided data for ",which); }
        
        
        if (allow.reweighted & weightingType == "complete" & "allow.reweighted"%in%names(eventChain)) {
            if (as.numeric(eventChain["allow.reweighted"])>0){
                if (inherits(x=params$data,"ExpressionSet") ) { ##Training set is given
                    betas<-.reweight(esetTrain = params$data,betas=c(params$intercept,params$weights),invent=invent)
                    params$intercept<-betas[1];
                    params$weights<-betas[-1,drop=FALSE];
                    params<-.filterParams(params,invent)

                    weightingType = "reweighted";
                    warning("Using a reweighted classifier for ",which," because ", length(invent$outs), " of " , length(invent$ins)+length(invent$outs)," covariates are missing");
                } else {
                    warning("Cannot create a reweighted classifier ",which, " because the trainingset is unknown. ",length(invent$outs), " of " , length(invent$ins)+length(invent$outs)," covariates are missing.");
                }
            }
        }

        if (weightingType == "complete") {
            if ("allow.reweighted"%in%names(eventChain)) {
                if (as.numeric(eventChain["allow.reweighted"])>0){
                    stop("Found ", length(invent$outs), " of " , length(invent$ins)+length(invent$outs)," covariates to be missing for classifier ", which,". Consider using 'allow.reweighted = TRUE'");
                }
            }
            stop("Found ", length(invent$outs), " of " , length(invent$ins)+length(invent$outs)," covariates to be missing for classifier ", which,". Reweighting is not allowed for this classifier");
                    
        }

    }


    eset<-eset[names(params$weights),] ###Filter probes

    dataismeancentered = FALSE;
    datahasunitvariance = FALSE;
    datatransformedtoreference = FALSE;

    for (eventIdx in c(1:length(eventChain))){
        keyword = names(eventChain[eventIdx])
        keyvalue = as.numeric(eventChain[eventIdx])
        if (keyword=="targetValue"){
            if (!is.null(used.MAS5TargetValue)) {
                if (dataislog2==TRUE){ ##Handle data as already logged2
                    exprs(eset)<-exprs(eset)+log(keyvalue,2)-log(used.MAS5TargetValue,2)
                } else { 
                    exprs(eset)<-exprs(eset)*keyvalue/used.MAS5TargetValue
                }
            }
        } else if (keyword=="truncate"){
            if (which%in%get.classifiers("mas5")){
                if (is.null(used.MAS5TargetValue) & datatransformedtoreference==FALSE){
                    stop("Please set 'used.MAS5TargetValue' to the right value. This is needed for classifiers that require truncation of lowest expressed genes.") 
                }    
            }
            if (dataislog2==TRUE){ ##Handle data as already logged2
                exprs(eset)[exprs(eset)<log(keyvalue,2)]<-log(keyvalue,2);
            } else { 
                exprs(eset)[exprs(eset)<keyvalue]<-keyvalue;
            }
        }
        else if (keyword=="to.log2"){
            if (keyvalue>0) {
                if (dataislog2==FALSE){ eset<-.log2(eset); }
                dataislog2<-TRUE; 
            }            
        } else if (keyword=="to.meancentering"){
            if (keyvalue>0) {
                result<-NULL
                if (do.batchcorrection) {  result <- .meancentering(eset,NULL);}
                else {
                    if (which%in%get.classifiers("mas5")){
                        if (is.null(used.MAS5TargetValue)){
                            stop("Please set 'used.MAS5TargetValue' to the right value. This is needed when skipping batch normalization.") 
                        }    
                    }
                    result <- .meancentering(eset,params);
                } 
                params$means=result$means;
                eset<-result$eset            
                dataismeancentered<-TRUE
            }
        } else if (keyword=="to.unitvariance" ){
            if (keyvalue>0) {
                result<-NULL
                if (do.batchcorrection) { result<-.unitvariance(eset,NULL); }
                else { 
                    if (which%in%get.classifiers("mas5")){
                        if (is.null(used.MAS5TargetValue)){
                            stop("Please set 'used.MAS5TargetValue' to the right value. This is needed when skipping batch normalization.") 
                        }    
                    }
                    result<-.unitvariance(eset,params);
                } 
                params$sds=result$sds;
                eset<-result$eset
                datahasunitvariance <-TRUE
            }
        }
        else if (keyword=="to.referencemeanvar") {
            if (keyvalue>0){
                if (do.batchcorrection) { eset<-.reftransform(eset,params); }
                else if (which%in%get.classifiers("mas5")){
                    if (is.null(used.MAS5TargetValue)){
                        stop("Please set 'used.MAS5TargetValue' to the right value. This is needed when skipping batch normalization.") 
                    }    
                }
        
                datatransformedtoreference=TRUE
            }
        }
        else if (keyword=="allow.reweighted") {} 
        else {
            stop("Encountered unknown keyword (",keyword,") in the eventChain of ", which);
        }
    }
    score<-.do.run(eset,params)
    new("Classifier",params=params,score=score,ins=invent$ins,batchCorrection=do.batchcorrection,weightingType=weightingType);
}
.probeInventorization<-function(eset,params){
    ins<-intersect(names(params$weights),rownames(eset))
    outs<-setdiff(names(params$weights),rownames(eset))
    return(list(ins=ins,outs=outs))
}
.reweight<-function(esetTrain,betas,invent){
    ins<-invent$ins;
    outs<-invent$outs

    esetTrain<-t(exprs(esetTrain))
    mat<-cbind(1,esetTrain[,c(ins,outs),drop=FALSE]) #complete data
    mat0<-cbind(1,esetTrain[,ins,drop=FALSE]); #known data
    mat1<-esetTrain[,outs,drop=FALSE] #unknown data
    Omega<-solve(t(mat0)%*%mat0)%*%t(mat0)%*%mat1
    fi<-c(betas[1],betas[ins])+Omega%*%betas[outs];
    ##Assume normality to adjust new variance and bias toward the old values ensure compatability with classification thresholds
    sdNew<-sd(mat0%*%fi)
    sdOld<-sd(mat%*%c(betas[1],betas[c(ins,outs)]))
    fi[-1]<-fi[-1]*sdOld/sdNew
    orgMu<-sum(colMeans(mat)*c(betas[1],betas[c(ins,outs)]))-betas[1]
    newMu<-sum(colMeans(mat0)*fi)-fi[1]
    fi[1]<-betas[1] + orgMu-newMu ##reset intercept to original value
    names(fi)<-c("Intercept",ins)
    betas<-fi
    return(betas)
}
.do.run<-function(eset,params){
    params$do.run(eset,params)
}
.isalreadylog2Transformed.check<-function(eset){
    quantile(exprs(eset),0.75,na.rm=TRUE)<30 | any(exprs(eset)<0)
}
.log2<-function(eset){
    if (.isalreadylog2Transformed.check(eset)) {
        warning("Performing a log2 transformation but the data already seems to be log2 transformed. If this is true pass 'ignore.log2 = TRUE' as an argument" );
    }
    exprs(eset)<-log(exprs(eset),2);
    eset
}
.filterprobes<-function(eset,params){
    probe_weights<-names(params$weights);
    eset[probe_weights,];
}
.meancentering<-function(eset,params){
    if (is.null(params)){
        params<-list(means=apply(exprs(eset),1,mean));
    }
    exprs(eset)<-exprs(eset)-params$means;
    list(eset=eset,means=params$means);
}
.reftransform<-function(eset,params){
    exprs(eset)<- (exprs(eset)-apply(exprs(eset),1,mean))*params$sds/apply(exprs(eset),1,sd)+params$means;
    eset;
}
.unitvariance<-function(eset,params){
    if (is.null(params)){
        params<-list(sds=apply(exprs(eset),1,sd));
    }
    exprs(eset)<-exprs(eset)/params$sds;
    list(eset=eset,sds=params$sds);
}
.filterParams<-function(params,invent){
    params$data     <-params$data[invent$ins,];
    params$weights  <-params$weights[invent$ins,drop=FALSE];
    params$means    <-params$means[invent$ins,drop=FALSE];
    params$sds      <-params$sds[invent$ins,drop=FALSE];
    params
}
