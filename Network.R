### Community Network (By csdambros 18Jun2012)

migmat<-function(net){
M<-matrix(0,length(net),length(net))
for(i in 1:length(net)){
M[i,c(net[[i]])]<-1/length(c(net[[i]]))
}
M}

########################

netsimu<-function(N,M,v=0,nruns=100000,stepshow=1000,meta=NULL){

record<-vector('numeric',length=(nruns-(nruns%%stepshow))/stepshow)

plot(NA,xlim=c(0,nruns),ylim=c(0,1))

M1<-M/matrix(N,1,length(N))[rep(1,length(N)),]

if(is.null(meta)){
meta=matrix(1,max(N),length(N))
for(i in 1:length(N)){
meta[-c(1:N[i]),i]<-NA
}}

meta2=as.vector(meta/meta)
meta3=(meta2*1:(length(N)*max(N)))

for (g in 1:nruns){

died<-sample(meta3[is.na(meta3)==F],1)

#died=200
row<-1+(died-1)%%max(N)
col<-(died-row)/max(N)+1

if(runif(1)<=v){

meta[died]<-max(meta)+1

}else{
replacement<-sample(meta3[is.na(meta3)==F],1,prob=rep(M1[col,],N))
meta[died]<-meta[replacement]}

#meta<-meta-min(meta)+1

#meta<-(as.numeric(as.factor(meta)))

if((g%%stepshow)==0){
selfsim<-simpson(meta[,1])
points(g,selfsim,pch=21)
record[(g-(g%%stepshow))/stepshow]<-selfsim
}}
meta
}
###############################3

netsimu2<-function(N,M,v=0,nruns=100000,stepshow=1000,meta=NULL,plotnew=T){

if(is.null(meta)){
meta=matrix(1,max(N),length(N))
for(i in 1:length(N)){
meta[-c(1:N[i]),i]<-NA
}}

for(i in 1:length(N)){
meta[-c(1:N[i]),i]<-NA
}
nnodes<-dim(M)[1]
N1<-c(0,cumsum(N))
totalN<-sum(N)

M1<-(1-v)*(t(t(M)/N))

if(plotnew==T){
plot(NA,xlim=c(0,nruns),ylim=c(0,1))}

record<-vector('numeric',length=(nruns-(nruns%%stepshow))/stepshow)
#record=record*NA

meta6<-meta[is.na(meta)==F]

for(g in 1:nruns){

meta5<-meta6

for (i in 1:nnodes){

newsp<-(max(meta5)+1):(max(meta5)+totalN)

meta5[(N1[i]+1):N1[i+1]]<-sample(c(meta6,newsp),N[i],prob=c(rep(M1[i,],N),rep(v/totalN,totalN)),replace=T)

}

meta6<-meta5

if((g%%stepshow)==0){

S4<-1-(sum(table(meta6[1:N[1]])*(table(meta6[1:N[1]])-1)))/(N[1]*(N[1]-1))

if(S4==0&&v==0){break}
record[(g-(g%%stepshow))/stepshow]<-S4
points(g,S4,pch=21,bg=r,cex=.5,col=0)
#points(1:((nruns-(nruns%%stepshow))/stepshow)*stepshow,record,type='l',col=1,lwd=2)
}}
points(1:((nruns-(nruns%%stepshow))/stepshow)*stepshow,record,type='l',col=1,lwd=2)
record[is.na(record)]<-0
fmeta<-as.numeric(as.factor(meta6))
lev<-rep(c(1:nnodes),N)
fmeta<-matrix(unlist(tapply(fmeta,lev,function(x){c(x,rep(NA,max(N)-length(x)))})),max(N),nnodes)
network<-list(Nnodes=nnodes, M=M,finalmeta=fmeta,simpson11record=record)
class(network)<-'n2'
return(network)
}

##############################################

anetwork<-function(N,M,v=0,nruns=100000,stepshow=1000,F1=NULL,plotnew=T){

nnodes<-dim(M)[1]

if(is.null(F1)){
F1<-matrix(1,nnodes,nnodes)
}

record<-{}
record[1]<-(1-F1[1,1])
if(plotnew==T){
plot(0,1-F1[1,1],xlim=c(0,nruns),ylim=c(0,1))}

F3=F1

N<-matrix(N,nnodes,nnodes)
kronecker<-(1-diag(diag(M)*0+1))

for(g in 1:nruns){

F3<-(1-v)^2*(M%*%t(M%*%(kronecker*t(F1)))+
M%*%(t(M)*matrix(diag(F1),nnodes,nnodes)*(1-1/N))+
M%*%(t(M)*(1/N)))


#for(i in 1:nnodes){
#for(j in 1:nnodes){
#rec={}
#for(k in 1:nnodes){
#rec[k]=sum(M[i,k]*M[j,-k]*F1[k,-k])+
#M[i,k]*M[j,k]*F1[k,k]*(1-1/N[k])+#change the 1 for k in the N
#M[i,k]*M[j,k]*(1/N[k])#change the 1 for k in the N
#}
#F3[i,j]=((1-v)^2)*sum(rec)
#}}

if((g%%stepshow)==0){
if(mean(F1==F3)==1){break}
record[1+(g-(g%%stepshow))/stepshow]<-(1-mean(F3))
points(g,1-mean(F3),pch=21,bg=2,cex=.5,col=0)
points(g,1-mean(diag(F3)),pch=21,bg=3,cex=.5,col=0)
#points(1:((nruns-(nruns%%stepshow))/stepshow)*stepshow,record,type='l',col=1,lwd=2)
}
F1=F3
}
record[(1+((g-(g%%stepshow))/stepshow)):(1+((nruns-(nruns%%stepshow))/stepshow))]<-(1-mean(F3))
points(c(1,1:((nruns-(nruns%%stepshow))/stepshow)*stepshow),record,type='l',col=1,lwd=2)
network<-list(Nnodes=nnodes,M=M,simpson11record=record,finalF=F3)
class(network)<-'Analiticalnet'
return(network)
}


#####################################


extractMH<-function(x){

resu<-matrix(NA,nrow(x),ncol(x))

for (i in 1:nrow(x)){

resu[i,]<-2*x[i,]/(diag(x)[i]+diag(x))

}
resu
}





