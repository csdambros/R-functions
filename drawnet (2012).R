### Create network topology (By csdambros 05Jun2012)

drawnet<-function(nsites,self=T){

nsites=20

if(is.list(nsites)){
ver=nsites
nsites=length(nsites)
ver[[length(ver)+1]]=100
}else{
ver<-as.list(as.data.frame(matrix(NA,1,(nsites+1)),col.names=1:(nsites+1)))
ver[[nsites+1]]=100
}


coords<-seq(360/(nsites),360,360/(nsites))/180*pi

par(mar=c(0,0,0,0))
plot.new()
plot.window(xlim=c(-1.4,1.1),ylim=c(-1.3,1.3))

for(i in 1:nsites){
arrows(sin(rep(coords[i],length(ver[[i]]))),cos(rep(coords[i],length(ver[[i]]))),sin(coords[ver[[i]]]),cos(coords[ver[[i]]]),col=3,angle=10,length=0.2)
}

points(sin(coords),cos(coords))
text(sin(coords)*1.05,cos(coords)*1.05,names(ver)[1:nsites],cex=0.5)

points(-1,-1,cex=8,pch=21,bg=rgb(1,0.2,0.2),col=0)
points(-1,-1,cex=5.4,pch=21,col=0,bg='dark grey')
points(-0.99,-0.99,cex=5,pch=21,col=0,bg='white');
text(-0.995,-0.985,'Done',cex=0.7,col='dark grey')
text(-0.99,-0.98,'Done',cex=0.7,col='dark grey')

points(-1,1,cex=8,pch=21,bg=rgb(0.2,0.2,1),col=0)
points(-1,1,cex=5.4,pch=21,col=0,bg='dark grey')
points(-0.99,1.01,cex=5,pch=21,col=0,bg='white');
text(-0.995,1.015,'RBN',cex=0.7,col='dark grey')
text(-0.99,1.02,'RBN',cex=0.7)


row=identify(c(sin(coords),-1),c(cos(coords),1),n=1,plot=F)
points(sin(coords[row]),cos(coords[row]),pch=21,bg=2)

if(row==(nsites+1)){

points(-1,1,cex=8,pch=21,bg=rgb(0.2,0.2,1),col=0);points(-1,1,cex=5.4,pch=21,col=0,bg='dark grey');points(-0.995,1.005,cex=5,pch=21,col=0,bg='white');text(-1,1.01,'RBN',cex=0.7)

inputs<-readline(prompt='How many average input nodes?')

connex<-rpois(nsites,as.numeric(inputs))

for(i in 1:nsites){
randbind<-sample((1:nsites)[-i],connex[i])
ver[[i]]<-c(ver[[i]],randbind)
arrows(sin(rep(coords[i],connex[i])),cos(rep(coords[i],connex[i])),sin(coords[randbind]),cos(coords[randbind]),col=4,angle=10,length=0.2)
}

points(-1,1,cex=5.4,pch=21,col=0,bg='dark grey')
points(-0.99,1.01,cex=5,pch=21,col=0,bg='white');
text(-0.995,1.015,'RBN',cex=0.7,col='dark grey')
text(-0.99,1.02,'RBN',cex=0.7)

}else{
text(-0.99,1.02,'RBN',cex=0.7,col='dark grey')
text(-0.99,-0.98,'Done',cex=0.7)

while(1){

ver[[row]]=c(ver[[row]],identify(c(sin(coords),-0.99),c(cos(coords),-0.98),n=1,plot=F,tolerance=5))

current<-ver[[row]]
value<-current[length(current)]

if(value==(nsites+1)){
points(-1,-1,cex=8,pch=21,bg=rgb(1,0.2,0.2),col=0);points(-1,-1,cex=5.4,pch=21,col=0,bg='dark grey');points(-0.995,-0.995,cex=5,pch=21,col=0,bg='white');text(-1,-0.99,'Done',cex=0.7)
break}

if(value!=row){
arrows(sin(coords[row]),cos(coords[row]),sin(coords[value]),cos(coords[value]),col=2,angle=10,length=0.2)
}

if(value==row){

points(sin(coords[row]),cos(coords[row]),pch=21,bg='white')
row=identify(c(sin(coords),-1),c(cos(coords),-1),n=1,plot=F,tolerance=5)
points(sin(coords[row]),cos(coords[row]),pch=21,bg=2)
}

current<-ver[[row]]
value<-current[length(current)]
}
points(-1,-1,cex=8,pch=21,bg=rgb(1,0.2,0.2),col=0);points(-1,-1,cex=5.4,pch=21,col=0,bg='dark grey');points(-0.99,-0.99,cex=5,pch=21,col=0,bg='white');text(-0.995,-0.985,'Done',cex=0.7,col='dark grey')
text(-0.99,-0.98,'Done',cex=0.7)
}
ver2<-lapply(1:nsites,function(x)as.numeric(levels(factor(ver[[x]][if(self==F){ver[[x]]!=x&ver[[x]]!=(nsites+1)}else
{ver[[x]]!=(nsites+1)}]))))
names(ver2)=names(ver)[1:nsites]
ver2
}

drawnet(3)





