####################################################
#Functions to rarefy the individuals in community given the species counts and estimate S without rarefying (analytical approximation)
#CSDambros March-05-2013

#################################
#Uses: rarefy.N(x,n,FUN=function(x)length(x),rep=100,replace=TRUE)
#      estimate.S(x)

#Parameters:

	# x   		A vector with the abundance of each species in a community
	# n   		The number of individuals to be randomly draw
	# FUN 		A calculation to be made on a list in the same format as x. For example, to count the number of species use FUN=function(x)length(x)
	# rep 		The number of replicates (default is 100)
	# replace 	Should the individuals be replaced in each draw? (default is TRUE)

#################################
#Functions

rarefy.N<-function(x,FUN,n,rep=100,replace=T){

FUN<-match.fun(FUN)

count.list<-rep(names(x),x)

result<-{}

for(r in 1:rep){

new<-sample(count.list,n,replace=replace)

result[r]<-FUN(table(new))

}
result
}



estimate.S<-function(x){

E<-NA

for(n in 1:sum(x)){

E[n]<-sum(1-((1-(n/sum(x)))^x))
}
E
}


#end
