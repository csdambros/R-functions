####################################################
#Function to plot histograms of species distributions along environmental gradients 
#CSDambros November-03-2012 (latest version)
#Version 3.2

#################################
#Use: Poncho (table,gradient=NULL,col=1,col.places=2,col.species=3,
#		col.gradient="gray",highlight.species=NA,highlight.places=NA,
#		xyratio=1.3333333333,cex.species=NULL,lwd=1,
#		cex.lab=NULL,cex.title=3,cex.gradient=2,
#		file=NULL,title=NULL,xlab="Ordered places",
#		ylab="Abundance",show.grad=NULL,show.places=NULL,
#		xfactor=0,lty=1,add=F,spgradient=NULL)

#alternativelly try:

#      Poncho3.2 (table)


#Parameters:

	# table			a numeric matrix, data frame or object containing objects in rows and atributes in columns.
	# gradient			an optional vector containing a atribute to order the table
	# col				the initial color(s) to fill or shade the rectangle(s) with. The default 1 means "black"
	# col.places		the color to fill or shade the species rectangles with. Only changes if highlight.species is not NULL
	# col.species		the color to fill or shade the places rectangles with. Only changes if highlight.places is not NULL
	# col.gradient  		the color to fill or shade the gradient rectangles with. Only changes if show.grad is not NULL
	# highlight.species	a vector containing the name of species that should be highlighted. The names should be the same as in table colnames 
	# highlight.places	a vector containing the name of places that should be highlighted.The names should be the same as in table rownames 
	# xyratio			xy ratio. The default value is 1.33. Changing this factor make graph more or less rectangular
	# cex.species		The magnification to be used for the name of atributes relative to the current setting of cex
	# lwd				The line width, a positive number, default is 1
	# cex.lab			The magnification to be used for x and y labels relative to the current setting of cex
	# cex.title			The magnification to be used for main title relative to the current setting of cex
	# cex.gradient		The magnification to be used for gradient name relative to the current setting of cex
	# file			The name of the output file, up to 511 characters. Should have a .pdf at the end name. If is NULL the image is plotted and not saved as pdf file
	# title			The main title (on top)
	# xlab			Label for the x axis
	# ylab			Label for the y axis
	# gradlab			label for the gradient name
	# show.grad			If not NULL, do not show the gradient, default is TRUE
	# show.places		If not NULL, show ordered places names as in table rownames 
	# xfactor			Add space on the graph to suport longer species names
	# spgradient		If not NULL, provide the order of species position

#################################
#Function

poncho<-function(
table,col=1,col.places=2,col.species=3,col.gradient="gray",
highlight.species=NA,highlight.places=NA,
xyratio=1.3333333333,cex.species=NULL,lwd=1,
cex.lab=NULL,cex.title=3,cex.gradient=2,gradient=NULL,
file=NULL,title=NULL,xlab="Ordered places",
ylab="Abundance",show.grad=TRUE,show.places=NULL,
xfactor=0,gradlab="Gradient",lty=1,lty.lines=1,border=1,decrease=FALSE,add=F,spgradient=NULL)
{

if(is.null(file)){
if(is.null(cex.species)){cex.species=.5}
if(is.null(cex.lab)){cex.lab=2}
if(is.null(cex.gradient)){cex.gradient=1.5}
}

if(is.null(cex.gradient)){cex.gradient=2}
if(is.null(cex.lab)){cex.lab=3}

####################################

if(is.null(file)==FALSE){

pdf(file=file,,paper="A4",width=0,height=0)

}

par(mar=c(2,3,.8,3))

#################################

if(is.null(gradient)){gradient<-prcomp(table)$x[,1]}

multiplier=ifelse(decrease==F,1,-1)


if(is.null(spgradient)==F){

tabela<-table[order(multiplier*gradient),order(spgradient)]
if(length(highlight.species)==length(col.species)){
col.species<-col.species[order(spgradient)]
}
}else{

tabela<-table[order(multiplier*gradient),order(colSums(table*multiplier*gradient)/colSums(table))]

if(length(highlight.species)==length(col.species)){
col.species<-col.species[order(colSums(table*multiplier*gradient)/colSums(table))]
}


}


#################################
x<-90
y<-x*xyratio

dimx<-x/length(0:(nrow(tabela)-1))
dimy<-y/length(1:ncol(tabela))
showg=0
if(is.null(show.grad)==F){showg=22}

if(add==F){
plot.new()
plot.window(xlim=c(0,120+xfactor),ylim=c(0,120+showg))
}

col<-rep(col,nrow(tabela))
col[match(highlight.places,rownames(tabela))]<-col.places

for(i in 1:ncol(tabela)){

col.sp<-match(highlight.species,colnames(tabela)[i])

if(sum(col.sp,na.rm=T)!=0){col.2=col.species}else{col.2=col}

if(length(highlight.species)==length(col.species)){col.2=col.species[i]}

rect((0:(nrow(tabela)-1))*dimx,(i-1)*dimy,(1:nrow(tabela))*dimx-dimx*.1,(i-1)*dimy+(tabela[,i]/max(tabela))*dimy*.8,col=col.2,lty=lty,lwd=lwd,border=border)
lines(c(-1,(nrow(tabela))*dimx),c((i-1)*dimy,(i-1)*dimy),lty=lty.lines,lwd=lwd)
}
text(rep(x+2,ncol(tabela)),cumsum(rep(120/ncol(tabela),ncol(tabela)))-120/ncol(tabela)/2,gsub("_"," ",colnames(tabela)),cex=if(is.null(cex.species)){ifelse(30/ncol(tabela)>1,1,30/ncol(tabela))}else{cex.species},adj=0,font=3)

if(lty.lines!=0){
rect(-1,0,-1,(ncol(tabela)-1)*dimy+dimy*.8,lwd=lwd)
}

text(x+4,133.5,gradlab,cex=cex.gradient,adj=0)

par(mgp=c(1,4,5))
title(main=title,cex.main=cex.title)

par(mgp=c(0,4,5))
title(ylab=ylab,cex.lab=cex.lab)

par(mgp=c(1,8,10),mar=c(3.5,3,.8,8))
title(xlab=xlab,cex.lab=cex.lab)


par(mgp=c(4,8,10),mar=c(0,0,0,0))
if(is.null(show.grad)==F){
gradient<-sort(gradient,decreasing=decrease)
gradient2<-gradient/abs(max(c(gradient,0))-min(c(gradient,0)))
val=126-min(c(gradient2*15,0))

rect((0:(nrow(tabela)-1))*dimx,val,(1:nrow(tabela))*dimx-dimx*.1,val+gradient2*15,lwd=lwd,col=col.gradient)

rect(x+2,126,x+2,126+15,lwd=lwd)

rect(x+2,126,x+3,126,lwd=lwd)

rect(x+2,126+15,x+3,126+15,lwd=lwd)

text(x+3.5,126,round(min(c(0,gradient)),2),adj=0,cex=.7)
text(x+3.5,126+15,round(max(c(0,gradient)),2),adj=0,cex=.7)

}

par(mar=c(5, 4, 4, 2) + 0.1,mgp=c(0,-.5,0))
if(is.null(show.places)==F){
axis(1,at=((0:(nrow(tabela)-1))*dimx+(1:nrow(tabela))*dimx-dimx*.1)/2,labels=rownames(tabela),las=2,cex.axis=if(is.null(cex.species)){ifelse(30/ncol(tabela)>1,1,30/ncol(tabela))}else{cex.species},tick=F)
}
#title(main="B",adj=0,cex.main=2)
par(mar=c(5, 4, 4, 2) + 0.1,mgp=c(3, 1, 0))
if(is.null(file)==F){

dev.off()
paste("You have saved the image as",file)
}else{"If resolution is low, try to use the argument file=imagename.pdf to save the image in a file"}

}


#end

##Use
#poncho(table)


