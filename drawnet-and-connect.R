
drawnet<-function(M){

###########################################

nsites=nrow(M)

coords<-seq(2/(nsites),2,2/(nsites))*pi

###########################################
par(mar=c(0,0,0,0))
plot.new()
plot.window(xlim=c(-1.4,1.1),ylim=c(-1.3,1.3))

points(sin(coords),cos(coords))
text(sin(coords)*1.05,cos(coords)*1.05,1:nsites,cex=0.5)
points(sin(coords)*.98,cos(coords)*.98,cex=.4,bg=2,pch=21)

###########################################

non.empties<-(1:nsites^2)[M!=0 ]

select.2<-((non.empties-1)%%nsites)+1
select.1<-1+(non.empties-(((non.empties-1)%%nsites)+1))/nsites

arrows(sin(coords[select.1])*0.98,cos(coords[select.1])*0.98,sin(coords[select.2])*0.98,cos(coords[select.2])*0.98,
length = 0,angle=25,col='dark grey',lwd=1)

text.x<-sin(coords[select.2])+(((sin(coords[select.1])-sin(coords[select.2])))/5)
text.y<-cos(coords[select.2])+(((cos(coords[select.1])-cos(coords[select.2])))/5)

text(text.x*.98+.01,text.y*.98,
M[non.empties],
cex=.7,pos=1,offset = 0.2)

invisible(M)

}



drawnet(matrix(.001,10,10))


##########################################################################

connect<-function(M,m=0.1,type='add'){

diag(M)<-1-rowSums(M)+diag(M)

nsites=nrow(M)

coords<-seq(2/(nsites),2,2/(nsites))*pi

c('+','-')[match(type,c('add','sub'))]

addsub<-c('+','-')[match(type,c('add','sub'))]

app<-match.fun(addsub)

#####################################

drawnet(M)

points(-1,-1,cex=8,pch=21,bg=rgb(1,0.2,0.2),col=0)
points(-1,-1,cex=5.4,pch=21,col=0,bg='dark grey')
points(-0.99,-0.99,cex=5,pch=21,col=0,bg='white');
text(-0.995,-0.985,'Done',cex=0.7,col='dark grey')
text(-0.99,-0.98,'Done',cex=0.7,col='black')

repeat{

select.2=nsites+2

select.1<-identify(c(sin(coords),-1),c(cos(coords),-1),n=1,plot=F)
if(select.1<nsites){points(sin(coords[select.1]),cos(coords[select.1]),bg=2,pch=21)}
if(select.1==(nsites+1)|select.2==(nsites+1)){break}


while(select.2!=select.1){

points(sin(coords[select.1]),cos(coords[select.1]),bg=2,pch=21)
select.2<-identify(c(sin(coords),-1),c(cos(coords),-1),n=1,plot=F)

if(select.1==(nsites+1)|select.2==(nsites+1)){break}

M[select.2,select.1]=app(M[select.2,select.1],m)
diag(M)<-1-rowSums(M)+diag(M)

points(sin(coords[select.2]),cos(coords[select.2]),bg=3,pch=21)

drawnet(M)
points(-1,-1,cex=8,pch=21,bg=rgb(1,0.2,0.2),col=0)
points(-1,-1,cex=5.4,pch=21,col=0,bg='dark grey')
points(-0.99,-0.99,cex=5,pch=21,col=0,bg='white');
text(-0.995,-0.985,'Done',cex=0.7,col='dark grey')
text(-0.99,-0.98,'Done',cex=0.7,col='black')

}
if(select.1==(nsites+1)|select.2==(nsites+1)){break}
}

M

}

#M<-matrix(0,30,30,dimnames=list(LETTERS[1:30],LETTERS[1:30]))

#drawnet(M)

#drawnet(cbind(M[,1]+0.1,M))

#connect(M,type='add')

#M2<-matrix(0.001,10,10)


#drawnet(M2)
#M2.1<-connect(M2)

#drawnet(M2.1)

#rowSums(M2.1)



#source('anetwork.R')

