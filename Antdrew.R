Antdrew<-function(gene2){

int.exo<-c('E','I') # Specify if is an intron or exon
int.exo<-rep(int.exo,ncol(gene2))[1:ncol(gene2)]# repeat the sequence E-I for all boxes (or more)

n=nrow(gene2) #number of rows (sequences) on the table

x<-max(rowSums(gene2,na.rm=T)) #set the maximum value for the x axis
y<-100 # set the maximum value for the y axis

w=(5*y)/(6*n-1) #calculate the width of each box

s=w/5 # calculate the space between boxes (in each row)

x2<-t(apply(gene2,1,cumsum)) # Calculate the final (right) location of the x values for each box
x1<-cbind(rep(0,n),x2[,-ncol(gene2)]) # Calculate the left limit of each box

#not used #x1<-c(0,cumsum(gene[-length(gene)])) #Calculate the x values for the bottom left corner
#not used#x2<-c(cumsum(gene)) # Calculate the x values for the top right corner

Y1<-(w+s)*(0:(n-1)) #calculate the initial location of the boxes on y axis
Y2<-Y1+w #calculate the final (top) location of the boxes on y axis

plot.new() #Open a new screen to plot
plot.window(xlim=c(-0.1*x,x),ylim=c(0,y)) #set the limits for the plotting region

#The following 2 lines are just used for a calculation

y1<-rep(Y1,ncol(gene2)) # Repeat the location of the boxes (y1) for each columng of the table
y2<-rep(Y2,ncol(gene2))# Repeat the location of the boxes (y2) for each columng of the table

see<-t(matrix(int.exo,nrow=ncol(gene2),ncol=nrow(gene2)))#Repeat the intro-exon sequence for the whole table (depending on the number of lines)

y1.1<-ifelse(see=='I',colMeans(rbind(y1,y2)),y1) # replace the values on the y1 table to the mean of y1 and y2 for the places corresponding to introns
y2.1<-ifelse(see=='I',colMeans(rbind(y1,y2)),y2) # replace the values on the y2 table to the mean of y1 and y2 for the places corresponding to introns

colors<-ifelse(is.na(gene2),0,0) #Create colors for all the boxes (by default is always 0, no color)
#colors[c(3,8,17),]<-1

rect(x1,y1.1,x2,y2.1,col=colors,lwd=2)#draw rectangles on the specified x and ys

text(-0.02*x,colMeans(rbind(Y1,Y2)),rownames(gene2),cex=.5,pos=2,font=3)# put the species names

}

##########################################
#Use the function

Antdrew(Ant.Data[nrow(Ant.Data):1,])








