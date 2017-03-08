####################################################
#Function to plot histograms of species distributions along environmental gradients 
#CSDambros July-06-2014 (latest version)
#Version 4.0

#################################
#Use: 

#		Poncho (x,env=NULL,phy=NULL,col=1,col.gradient="gray",
#		,cex.species=0.6,cex.lab=NULL,cex.gradient=2,
#		file=NULL,xlab="Ordered places",ylab.top="Abundance",
#		ylab.bottom="Environment",lty=3)

#alternativelly try:

#      poncho (x)


#Parameters:

# x			a numeric matrix, data frame or object containing objects in rows and atributes in columns.
# env		an optional vector containing a atribute to order the table
# phy		Optional; Phylogenetic tree to order the species and plot.
# col		The initial color(s) to fill or shade the rectangle(s). The default 1 means "black". To color the species, provide a vector with the color (or factor) for each species. To color by the sites, provide a vector with the color for each site. To individually color each cell, provide a matrix with the same dimensions as x (main table). The function will recognize the length of the color vector and automatically attribute to sites, species or both.
# col.gradient  		the color to fill or shade the gradient rectangles with. Only changes if show.grad is not NULL
# cex.species		The magnification to be used for the name of atributes relative to the current setting of cex
# cex.lab			The magnification to be used for x and y labels relative to the current setting of cex
# cex.title			The magnification to be used for main title relative to the current setting of cex
# cex.gradient		The magnification to be used for gradient name relative to the current setting of cex
# file			The name of the output file, up to 511 characters. Should have a .pdf at the end name. If is NULL the image is plotted and not saved as pdf file
# xlab			Label for the x axis
# ylab.top			Label for the y axis
# ylab.bottom			label for the gradient name

#################################
#Function

poncho<-function(x,env=gradient,col=col.places,border="grey",col.gradient="grey",cex.species=0.5,cex.lab=2,cex.gradient=1,file=NULL,xlab="Ordered sites",ylab.top=ylab,ylab="Abundance",ylab.bottom=gradlab,lty.lines=3,decrease=FALSE,asis=FALSE,phy=NULL,symbol=1,lty=3,col.lty="darkgrey",col.places=1,col.species=1,gradlab="Environment",gradient=NULL){
  
  symbol<-((symbol-1)%%2)+1
  col.orig<-col
  col<-as.integer(as.factor(col))
  
  if(is.null(env)){env<-prcomp(x)$x[,1]}
  if(length(col)==ncol(x)){z<-(t(t(x>0)*as.integer(col)))}else{
    z<-(x>0)*as.integer(col)}
  
  z<-z[order(env),order(colSums(x*env)/colSums(x))]
  
  if(!is.null(sp.gradient)){
    
    z<-(x>0)*as.integer(col)
    z<-z[order(env),order(sp.gradient)]
    
  }
  
  
  if(!is.null(phy)){
    
    
    phy<-rotateConstr(phy,colnames(z))
    
    z<-z[,match(phy$tip.label[phy$edge[,2][phy$edge[,2]<=length(phy$tip.label)]],colnames(z))];
    colnames(z)<-phy$tip.label[phy$edge[,2][phy$edge[,2]<=length(phy$tip.label)]]
  }
  
  
  
  
  space.x=(0.70)/(nrow(z)+1)
  space.y=1/(ncol(z))
  
  x1=((z>0)*seq(.1,.8-space.x,length=nrow(z)))[(!is.na(z))&z>0]
  y1=(t(t(z>0)*(seq(0,1-space.y,length=ncol(z)))))[(!is.na(z))&z>0]
  
  # Plot the figure
  
  op<-par(no.readonly=T)
  
  par(mar=c(0,0,0,0),oma=c(0,0,0,0),xpd=NA)
  layout(matrix(c(1,1,1,1,2),5,1))
  
  xlim=c(0,1)
  ylim=c(0,1)
  
  if(!is.null(phy)){
    require(ape)
    layout(matrix(c(1,1,1,1,0,2,2,2,2,3,2,2,2,2,3),5,3))
    plot(phy, show.tip.label=F,plot=TRUE,root.edge = TRUE,no.margin = TRUE)
    ylab.top=""
    xlim=c(0.15,1.1);ylim=c(space.y/2,1-space.y/2)
    
  }
  
  plot.new()
  plot.window(xlim=xlim,ylim=ylim)
  
  segments(0.09,seq(space.y/2,1-space.y/2,length=ncol(z)),0.81,seq(space.y/2,1-space.y/2,length=ncol(z)),lty=lty,col=col.lty)
  
  segments(c(0,.09,.09),c(-.02,-0.01,-0.01),c(.82,.82,.09),c(-.02,-.01,1),lty=c(3,1,1))
  text(.02,0.5,ylab.top,srt=90,cex=2)
  
  if(symbol==1){rect(x1,y1,x1+space.x,y1+space.y,col=col.orig[z[(!is.na(z))&z>0]],border = border)}
  if(symbol==2){points(x1+space.x/2,y1+space.y/2,pch=21,col=border,bg=col.orig[z[(!is.na(z))&z>0]])}
  
  text(.82,seq(space.y/2,1-space.y/2,length=ncol(z)),gsub("_"," ",colnames(z)),adj=0,cex=cex.species,col=col.species,font=3)
  
  # Plot gradient
  
  plot.new()
  plot.window(xlim=xlim,ylim=c(-.5,1.2))
  
  gradient<-sort(env,decreasing=F)
  gradient2<-gradient/abs(max(c(gradient,0))-min(c(gradient,0)))
  val=-min(c(gradient2,0))
  
  segments(c(.085,.085,.085,.1,.1,.8,.1),c(0,0,1,-.05,-.05,-.05,val),c(.075,.085,.075,.8,.1,.8,.8),c(0,1,1,-.05,-.1,-.1,val))
  
  rect(seq(.1,0.8-space.x,length=nrow(z)),val,seq(.1,.8-space.x,length=nrow(z))+space.x,val+gradient2,col=col.gradient)
  
  text(.065,c(0,1),c(floor(min(gradient)),ceiling(max(gradient))),adj=1,cex=.8)
  text(.8-(.8-.1)/2,-.35,xlab,cex=2)
  text(.01,0.5,ylab.bottom,srt=90)
  
  # Print to pdf device
  
  if(!is.null(file)){
    dev.copy(pdf,file=paste(file,rep(".pdf",1-length(grep(".pdf",file))),sep=""),width=5+(!is.null(phy))*1.5,height=8.27)
    dev.off()
  }else{cat('use file = "filename.pdf" to save as a pdf')}
  
  par(op)
  
}

# USE

# poncho (table)
# poncho (table, env = environment, phy = phylogeny)
