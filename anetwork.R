##############################################
#Neutral simulation of a metacommunity using a network approach (semi analytical approach). See Economo and Keitt (2008)
#Includes two supplentary functions (not required for 'anetwork')
#Developed by Dambros CS in Jun 2012

#This function was implemented mostly using matrix algebra to run faster
#Download the script 'neutralworkedexample.R' to see examples of how to use the function ant integrate with real data

#################################
#Use: anetwork (N, M, V, ...)

#Parameters:

	# N   	The number of individuals in each local community, by default N is a single number and all local communities have the same number of individuals. A vector of local community sizes can also be provided
	# M		Migration matrix among all local communities in the metacommunity where the rows represent the areas (i) receiving migrants with rate mij from the areas on the columns (j).
	# V		The speciation rate, single numeric value. The same for all local communities.
	# ...		Further parameters that can be used (see function)

#################################
#Function

anetwork<-function(N,M,v=0,nruns=100000,stepshow=1000,F1=NULL,plotnew=T){

	nnodes<-dim(M)[1]

	if(is.null(F1)){
		F1<-matrix(1,nnodes,nnodes)}

	record<-{}
	record[1]<-(1-F1[1,1])

	if(plotnew==T){
		plot(0,1-mean(F1),ylab='Diversity ( 1 - fij )',xlab='Generation',xlim=c(0,nruns),ylim=c(0,1))}

	F3=F1

	N<-matrix(N,nnodes,nnodes)
	kronecker<-(1-diag(diag(M)*0+1))

	for(g in 1:nruns){

#The following 3 lines represent the full calculation presented by Economo and Keitt (2008):

		F3<-(1-v)^2*(M%*%t(M%*%(kronecker*t(F1)))+
		M%*%(t(M)*matrix(diag(F1),nnodes,nnodes)*(1-1/N))+
		M%*%(t(M)*(1/N)))

###############
#the followning code generates the same results not using matrix algegra (take longer time to run)

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
			points(g,1-mean(diag(F3)),pch=21,bg=3,cex=.5,col=0)}

		F1=F3
		}

	record[(1+((g-(g%%stepshow))/stepshow)):(1+((nruns-(nruns%%stepshow))/stepshow))]<-(1-mean(F3))
	points(c(1,1:((nruns-(nruns%%stepshow))/stepshow)*stepshow),record,type='l',col=1,lwd=2)
	network<-list(Nnodes=nnodes,M=M,simpson11record=record,finalF=F3)
	class(network)<-'Analiticalnet'
	return(network)
}

#End of the anetwork funtion



#####################################
#Auxiliary function to extract the Morisita-Horn diversity index among samples from the anetwork function output (x)

extractMH<-function(x){

	resu<-matrix(NA,nrow(x),ncol(x))

	for (i in 1:nrow(x)){

		resu[i,]<-2*x[i,]/(diag(x)[i]+diag(x))

	}
	resu
}

#End of the extractMH function

#####################################
#Auxiliary function to calculate the difference of the predicted diversity values given N,M,m and V from empirical data (or from other simulation)

#Output: differences in alpha diversities squared + differences in beta-diversity squared

diff.anetwork<-function(x,N,M,V,F1=1-x,...){


	if(sum(M)==nrow(M)&V<1&V>0){

		model<-anetwork(N,M,V,F1=F1,...)

		result<-sum(sqrt(((1-x)-model$finalF)^2))

		result

	}else{100}

}





#End of the function diff.anetwork

###########################################
#Auxiliary function to draw the network scheme given M (nice to visualize the migrations) 

drawnet<-function(M,cex=1){

###########################################

nsites=nrow(M)

coords<-seq(2/(nsites),2,2/(nsites))*pi

###########################################
par(mar=c(0,0,0,0))
plot.new()
plot.window(xlim=c(-1.4,1.1),ylim=c(-1.3,1.3))

points(sin(coords),cos(coords))
text(sin(coords)*1.05,cos(coords)*1.05,rownames(M),cex=.7)
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
cex=cex,pos=1,offset = 0.2)

invisible(M)

}



#drawnet(matrix(.001,10,10))


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




#End