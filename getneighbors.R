#### Get the neighbors in a round torus for a given cell and radius

getneighbors<-function(select1,comm,radius=1){

if(is.null(dim(comm))==FALSE){dim<-dim(comm)[1]}else{dim<-comm}

select=t(matrix(select1,length(select1),radius*2+1))

x1<-((select-1+(-radius:radius))%%dim)+1
y1<-(((select-(select%%dim))/dim+(-radius:radius))%%dim)+1

extx<-rep(t(x1),radius*2+1)
exty<-rep(t(y1),radius*2+1)

putorder<-rep(select1*dim^2,(radius*2+1)^2)

if(is.null(dim(comm))==FALSE){

		matrix(comm[extx[order(putorder)]+(exty[order(exty+putorder)]-1)*dim],
		ncol=length(select1))[,order(select1)]}

	else{

		matrix(extx[order(putorder)]+(exty[order(exty+putorder)]-1)*dim,
		ncol=length(select1))[,order(select1)]}


}













