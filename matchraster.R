matchraster<-function(datamatrix,pointposition,dataposition=NULL,sample=1,datamatrix.names=paste(substitute(datamatrix),1:sample,sep='.')){

if(is.null(dataposition)){dataposition<-list(as.numeric(rownames(datamatrix)),as.numeric(colnames(datamatrix)))}

if(length(dataposition[[1]])<2){stop('dataposition or rownames on datamatrix must be provided')}

closestlat<-(sapply(pointposition[[1]],function(x){order(sqrt((x-as.numeric(dataposition[[1]]))^2))}))[1:sample,]
closestlong<-(sapply(pointposition[[2]],function(x){order(sqrt((x-as.numeric(dataposition[[2]]))^2))}))[1:sample,]

values<-datamatrix[closestlat+(closestlong-1)*nrow(datamatrix)]

matrix(values,nrow=sample,dimnames=list(datamatrix.names,names(pointposition[[1]])))

}

