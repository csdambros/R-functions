##Developed by csdambros in 07May2012

#_______________WRAP TORUS

wrap=function(a,b,dim){((a-1+b+dim)%%dim)+1}

#_______________NEIGHBORS


#__Moore neighborhood function

neighbors=function(map,select,dim){

af<-wrap(select,1,dim)
be<-wrap(select,-1,dim)

neighborhood=matrix(NA,nrow(select),8)

for(i in 1:nrow(select)){
neighborhood[i,]<-map[c(be[i,1],be[i,1],be[i,1],select[i,1],select[i,1],af[i,1],af[i,1],af[i,1])+
(c(af[i,2],select[i,2],be[i,2],be[i,2],af[i,2],be[i,2],select[i,2],af[i,2])-1)*dim]
}
neighborhood}

#__Neighbors in distance (including Moore neighborhood, dist=1)

neighborsdist=function(map,select,dist=1,dim){

neighborhood=list()

for(i in 1:nrow(select)){

disti=dist[(i-1)%%(length(dist))+1]

lins=wrap(select[i,1],(-(disti):(disti)),dim)
cols=wrap(select[i,2],c(-disti,disti),dim)
colsleft=wrap(select[i,2],(-(disti-1):(disti-1)),dim)
linsleft=wrap(select[i,1],c(-disti,disti),dim)

neighborhood[[i]]<-map[c((lins+(cols[1]-1)*dim),
(lins+(cols[2]-1)*dim),
(linsleft[1]+(colsleft-1)*dim),
(linsleft[2]+(colsleft-1)*dim))]

}
neighborhood}


#neighborsdist(community,rbind(lin,col),1,dim)


#__DISPERSION PROBABILITY FUNCTIONS


CK=function(x,C,a,b){C*(a^x+(b^2)/(x^2+b^2))}

hist(sample(100,100,prob=CK(seq(0.01,1,.01),C=0.3,a=0,b=0.0001),replace=T))

seq(0.01,1,(1-0.01)/(100-1))


###########################################################################
#__________________________________________________________________________




#__PROPER FUNCTION TO EVOLVE THE COMMUNITY

evolveCommunity=function(community,deathrate=1,speciationrate=0,
migrationrate=0,ngenerations=100,dispersionFUN=CK,FUNparameters=list(C=0.3,a=0,b=0.0001),
metacommunity=community,totalmetacommunity=metacommunity){

dispersionFUN=match.fun(dispersionFUN)
dim=ncol(community)
phylogeny=paste('(',paste(levels(factor(community)),'s',sep='',collapse=','),')',sep='')

#FUNparameters=list(C=0.3,a=0,b=0.0001)

#match.arg(names(formals(CK)),names(FUNparameters),several.ok=T)

#FUNparameters$x=seq(0.01,1,(1-0.01)/(dim-1))*dim
FUNparameters$x=1:dim


prob=do.call(dispersionFUN,FUNparameters)
par(mfrow=c(2,2))
#plot(1:dim,prob,type='o')
hist(sample(dim,1000,prob=prob,replace=T),0:dim)
image(community)

if (deathrate>=1){

	ndie=deathrate
}else{
	
	ndie=round(deathrate*(dim^2))
}


for (i in 1:ngenerations){

select1=sample(dim^2,ndie)
distcolonizer=sample(dim,ndie,prob=prob,replace=T)
#distcolonizer=1
lin=((select1-1)%%dim)+1
col=1+(select1-lin)/dim

#Determine the neighbors of each vacant spot in each distance
neighborhood=neighborsdist(community,cbind(lin,col),distcolonizer,dim)

#Decide which one of the neighbors will be selected (sequential number) 
sampleneighbor=ceiling(runif(ndie)*(distcolonizer*8))

#A smart way to avoid the 'for' loop here (select the sampled neighbor). It is amazingly fast!!!
#Based on sampleneighbor, select the individual from the neighborhood (identify the species)

colonizer=unlist(neighborhood)[sampleneighbor+cumsum(c(0,distcolonizer*8))[1:ndie]]

speciatemigrate=runif(ndie)

#old=paste((community[select1])[speciatemigrate<speciationrate],'s',sep='')
#new=paste((1:sum(speciatemigrate<speciationrate))+(max(c(community,metacommunity))),'s',sep='')

#facold=levels(factor(old))

#for(j in facold){
#newsps=paste('(',paste(c(j,new[old==j]),collapse=','),')',sep='')
#phylogeny=gsub(j,newsps,phylogeny)
#}

#WARNING: Metacommunity and totalmetacommunity temporarily unimplemented

colonizer[speciatemigrate<speciationrate]=1:sum(speciatemigrate<speciationrate)+(max(totalmetacommunity))
colonizer[speciatemigrate>(1-migrationrate)]=sample(metacommunity,sum(speciatemigrate>(1-migrationrate)),replace=T)

community[select1]=colonizer

totalmetacommunity=max(max(totalmetacommunity),max(community))

}
image(community)
hist(community,0:max(community))

evolve=list(
#phylo=phylogeny,
community=community)
#return(evolve)
#c(sampleneighbor,prob,distcolonizer,colonizer,neighborhood,select1)
}

#########################################
#########################################











