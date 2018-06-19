
poncho<-function(x,env=gradient,col=x>0,xborder="grey",col.gradient="grey",cex.species=0.5,cex.lab=2,cex.gradient=1,file=NULL,xlab="Ordered sites",ylab.top="Abundance",ylab.bottom="Environment",lty.lines=3,decrease=FALSE,phy=NULL,symbol=1,lty=3,col.lty="darkgrey",sp.gradient=NULL,gradient=NULL){
  
  symbol<-((symbol-1)%%2)+1
  
  if(is.null(env)){env<-prcomp(x)$x[,1]}
  if(length(col)==ncol(x)){
    #xcol<-(t(t(x>0)*as.integer(col)))
    xcol<-matrix(rep(col,each=nrow(x)),nrow=nrow(x),ncol(x))
    print(dim(xcol))
    
    # (t(t(x>0)*as.integer(col)))
    
  }
  else{
    #xcol<-(x>0)*as.integer(col)
    xcol<-matrix(col,nrow(x),ncol(x))
    
    print(dim(xcol))
    }
  
  if(!is.null(sp.gradient)){
    
    col.order<-order(sp.gradient)
    x<-x[order(env),col.order]
    xcol<-xcol[order(env),col.order]
    
  }else{
    
    col.order<-order(colSums(x*env)/colSums(x))
    
    x<-x[order(env),col.order]
    xcol<-xcol[order(env),col.order]
    
  }
  
  if(!is.null(phy)){
    require(ape)
    
    phy<-rotateConstr(phy,colnames(x))
    
    col.order<-match(phy$tip.label[phy$edge[,2][phy$edge[,2]<=length(phy$tip.label)]],colnames(x))
    
    x<-x[,col.order]
    xcol<-xcol[,col.order]
    
    #colnames(x)<-phy$tip.label[phy$edge[,2][phy$edge[,2]<=length(phy$tip.label)]]
    
  }
  
  space.x=(0.70)/(nrow(x)+1)
  space.y=1/(ncol(x))
  
  x1=((x>0)*seq(.1,.8-space.x,length=nrow(x)))
  y1=(t(t(x>0)*(seq(0,1-space.y,length=ncol(x)))))
  # Plot the figure
  
  op<-par(no.readonly=T)
  
  par(mar=c(0,0,0,0),oma=c(0,0,0,0),xpd=NA)
  layout(matrix(c(1,1,1,1,2),5,1))
  
  xlim=c(0,1)
  ylim=c(0,1)
  
  if(!is.null(phy)){
    layout(matrix(c(1,1,1,1,0,2,2,2,2,3,2,2,2,2,3),5,3))
    plot(phy, show.tip.label=F,plot=TRUE,root.edge = TRUE,no.margin = TRUE)
    ylab.top=""
    xlim=c(0.15,1.1);ylim=c(space.y/2,1-space.y/2)
    
  }
  
  plot.new()
  plot.window(xlim=xlim,ylim=ylim)
  
  #segments(0.09,seq(space.y/2,1-space.y/2,length=ncol(x)),0.81,seq(space.y/2,1-space.y/2,length=ncol(x)),lty=lty,col=col.lty)
  segments(0.09,seq(0,1-space.y,length=ncol(x)),0.81,seq(0,1-space.y,length=ncol(x)),lty=lty,col=col.lty)
  
  
  segments(c(0,.09,.09),c(-.02,-0.01,-0.01),c(.82,.82,.09),c(-.02,-.01,1),lty=c(3,1,1))
  text(.02,0.5,ylab.top,srt=90,cex=2)
  
  if(symbol==1){rect(x1[x>0],y1[x>0],x1[x>0]+space.x,y1[x>0]+space.y*(x/max(x))[x>0],col=xcol[x>0],border = 1)}
  #  if(symbol==1){rect(x1,y1,x1+space.x,y1+space.y,col=xcol,border = 1)}
  
  if(symbol==2){points(x1[x>0]+space.x/2,y1[x>0]+space.y/2,pch=21,col=border,bg=xcol[x>0])}
  #if(symbol==2){points(x1+space.x/2,y1+space.y/2,pch=21,col=border,bg=col.orig[z[(!is.na(z))&z>0]])}
  
  
  text(.82,seq(space.y/2,1-space.y/2,length=ncol(x)),gsub("_"," ",colnames(x)),adj=0,cex=cex.species,font=3)
  
  # Plot gradient
  
  plot.new()
  plot.window(xlim=xlim,ylim=c(-.5,1.2))
  
  gradient<-sort(env,decreasing=F)
  gradient2<-gradient/abs(max(c(gradient,0))-min(c(gradient,0)))
  val=-min(c(gradient2,0))
  
  segments(c(.085,.085,.085,.1,.1,.8,.1),c(0,0,1,-.05,-.05,-.05,val),c(.075,.085,.075,.8,.1,.8,.8),c(0,1,1,-.05,-.1,-.1,val))
  
  rect(seq(.1,0.8-space.x,length=nrow(x)),val,seq(.1,.8-space.x,length=nrow(x))+space.x,val+gradient2,col=col.gradient)
  
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
