randomizeff<-function(table,runs=100){

splist<-rep(colnames(table),colSums(table))

allentry<-list()

for(i in 1:runs){

entry<-matrix(NA,nrow=nrow(table),ncol=ncol(table))
colnames(entry)=colnames(table)

ind1<-cumsum(c(0,rowSums(table)))+1
ind2<-cumsum(rowSums(table))

S1<-sample(splist)

for(k in 1:nrow(table)){
ver<-table(as.factor(S1)[ind1[k]:ind2[k]])
entry[k,]<-ver
}

#poncho3.1(entry)

#plot3d(entry)

allentry[[i]]<-entry

}

allentry
}

