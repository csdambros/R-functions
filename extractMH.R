#####################################
#Auxiliary function to extract the Morisita-Horn diversity index among samples from the anetwork function output (x)

extractMH<-function(x){

resu<-matrix(NA,nrow(x),ncol(x))

for (i in 1:nrow(x)){

resu[i,]<-2*x[i,]/(diag(x)[i]+diag(x))

}
resu
}

#end of the extractMH function
###############################################################