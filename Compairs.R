####################################################
#Function to plot a panel with the regression of groups of predictor x response variables (save your time when dealing with multiple varibles)
#CSDambros October-05-2011. Latest version: November-10-2012

#################################
#Use: compairs(x,y)

#Parameters:

	# x   	a table of independent variables (each variable in a column)
	# y		a table of dependent variables (each variable in a column)
	# ...		any other graphic parameter as described in plot and par functions is accepted


#Output:

# Composite graph showing all the simple linear regressions among predictor and response variables
# Green lines are significant at level of 0.01
# Red lines are significant at level of 0.05

#################################
#Function

compairs<-function(x,y,FUN=lm,show="graph",xaxt="n",yaxt="n",backcol='grey',fg=1,highlight.05=2,highlight.01=1,...){

	op<-par(no.readonly = TRUE)

	base<-match.fun(FUN)

	dependent<-data.frame(y)
	independent<-data.frame(x)

	par(mfrow=c((ncol(dependent)+1),(ncol(independent)+1)))

	par(mar=c(0.1,0,0.1,1))

	for(i in 1:ncol(dependent)){

		plot(1,1,ylim=c(min(dependent[,i]),max(dependent[,i])),xlim=c(0,5),type="n",xaxt=xaxt,yaxt=yaxt,fg=0,...)
		rect(5-1,min(dependent[,i]),5,max(dependent[,i]))
		rect(5-1,dependent[,i],5,dependent[,i])
		text(2.5,(min(dependent[,i])+max(dependent[,i]))/2,colnames(dependent)[i],srt=90)

		for(k in 1:ncol(independent)){

			if(show=="cor"){
				plot(1,1,xlim=c(0,5),type="n",xaxt=xaxt,yaxt=yaxt,fg=fg,...)
				text(3,1,round(cor(dependent[,i],independent[,k]),3))}
			else{
				if(show=="graph"){
					if(var(independent[,k],na.rm=T)==0|is.na(var(independent[,k],na.rm=T))){}
					else{
						z<-lm(dependent[,i]~independent[,k])
						p<-summary(z)[[4]][2,4]
						a<-summary(z)[[4]][2,1]
	
						plot(dependent[,i]~independent[,k],xaxt=xaxt,yaxt=yaxt,fg=fg,pch=22,
						bg=ifelse(p<0.05,1,backcol),col=0)#(ifelse(sqrt(newbg^2)<1,1,newbg))[i,k],...)

						abline(z,col=ifelse(p<0.01,highlight.01,ifelse(p<.05,highlight.05,backcol)))
					}
				}
				else{
					plot(1,1,xlim=c(0,5),type="n",xaxt=xaxt,yaxt=yaxt,fg=fg,...)
					if(var(independent[,k],na.rm=T)==0|is.na(var(independent[,k],na.rm=T))){}
					else{
						text(2.5,1,round(summary(lm(dependent[,i]~independent[,k]))[[4]][2,4],4))
					}
				}
			}
		}
	}

	par(mar=c(0,0,1,1))

	plot(1,1,xlim=c(0,5),type="n",xaxt=xaxt,yaxt=yaxt,fg=0,...)

	for(i in 1:ncol(independent)){
		plot(independent[,i],1:nrow(independent),xlim=c(min(independent[,i]),max(independent[,i])),ylim=c(0,5),type="n",xaxt=xaxt,yaxt=yaxt,fg=0,...)

		rect(min(independent[,i]),4,max(independent[,i]),5)
		rect(independent[,i],4,independent[,i],5)
		text((min(independent[,i])+max(independent[,i]))/2,2.5,colnames(independent)[i])
	}

	par(op)

}


#compairs(x,y)

#end





