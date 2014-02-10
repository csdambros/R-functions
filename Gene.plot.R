#For Andrew - Creating a function to plot gene sequences

plot.new()# Create an empty screen
plot.window(xlim=c(0,100),ylim=c(0,100))# Create the plot window (specify the x and y limits)

gene<-c(10,5,30,10,5,5,20) # Input for lengths
int.exo<-c('E','I','E','I','E','I','E') # Specify if is an intron or exon

x1<-c(0,cumsum(gene[-length(gene)])) #Calculate the x values for the bottom left corner
x2<-c(cumsum(gene)) # Calculate the x values for the top right corner

y1<-ifelse(int.exo=='E',20,22.5)#Calculate the y values for the bottom left corner
y2<-ifelse(int.exo=='E',25,22.5)#Calculate the y values for the top right corner

rect(x1,y1,x2,y2,col=2:5)# Draw rectangles on the specified values of x and y. col determines the sequence of colors



