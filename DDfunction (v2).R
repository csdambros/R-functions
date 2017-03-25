DDfunction=function(geodist,comdist,formula=formula(y~a*x+b),start=list(a=1,b=1),modelname='user specified',cex=0.8,pch=21,col=0,xlab='Geographical distance (degrees)',ylab='Jaccard similarity',bg='dark grey',...){

#logarithm model: -a-b*exp(geodist*-c))

alldist=data.frame(x=matrix(geodist),y=matrix(comdist))

comdist2=nls(formula,start=start, data=alldist)

linear<-lm(y~x,data=alldist)
additive<-gam(y~s(x,k=50),data=rbind(alldist))


plot(geodist,comdist,cex=cex,pch=pch,col=col,xlab=xlab,ylab=ylab,bg=bg,...)

drawinpoints=seq(min(alldist$x),max(alldist$x),length.out=300)

points(drawinpoints,drawinpoints*linear$coefficients[2]+linear$coefficients[1],type='l',lwd=2)
#points(drawinpoints,-coefficients(comdist2)['a']-(coefficients(comdist2)['b']*exp(drawinpoints*-coefficients(comdist2)['c'])),type='l',col=2,lwd=2)
points(drawinpoints,predict(comdist2,newdata=list(x=drawinpoints)),type='l',col=2,lwd=2)
points(drawinpoints,predict(additive,newdata=list(x=drawinpoints)),type='l',col='dark green',lwd=2)

SStotlin<-sum(((alldist$y)-mean(alldist$y))^2)
SSerrlin<-sum(((alldist$y)- (alldist$x*linear$coefficients[2]+linear$coefficients[1]))^2)
Rsqrtlin=1-SSerrlin/SStotlin

SStot<-sum(((alldist$y)-mean(alldist$y))^2)
SSerr<-sum(((alldist$y)- predict(comdist2,newdata=list(x=(alldist$x))))^2)
Rsqrt=1-SSerr/SStot

SStotspl<-sum(((alldist$y)-mean(alldist$y))^2)
SSerrspl<-sum(((alldist$y)- predict(additive,newdata=list(x=(alldist$x))))^2)
Rsqrtspl=1-SSerrspl/SStotspl

legend('topright',title='R-squares',c(
paste('line:',round(Rsqrtlin,4)),
paste(modelname,round(Rsqrt,4)),
paste('spline:',round(Rsqrtspl,4))

),bty='n')

}




