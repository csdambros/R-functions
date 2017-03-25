####################################################
#Function to calculate the PIE and Simpson's diversity indexes 
#(PIE is the same as the Simpson's index without replacement)
#CSDambros March-04-2013


#################################
#Use: PIE(x)
#     simpson(x,y=NULL)	

#Parameters:

	# x   A vector with the abundance of each species in a community  
	# y   An optional vector with the abundance of each species in a second community. If y is provided, then the index calculate the simpson diversity index between samples
  # see simpson 2 for another implementation
#################################
#Function PIE

PIE<-function(x){

(sum(x)/(sum(x)-1))*(1-(sum((x/sum(x))^2)))

}

# simpson<-function(x,y=NULL){
#   if(is.null(y)){
#     1-sum((x/sum(x))^2)
#   }
#   else{
#     1-sum((x/sum(x))*(y/sum(y)))
#   }
# }

simpson<-function(x,y=x){1-sum((x/sum(x))*(y/sum(y)))}

# Communities as vectors with the species names repeated (eg. c("sp2","sp2","sp3","sp4")). The communities dont need to have the same species pool neither the same length (zero abundances are not represented)

simpson2<-function(x,y=x){
  
  if(is.character(x)){
  
    
    1-((((1-simpson(table(c(x,y))))*(length(c(x,y)))^2)-(1-simpson(table(x)))*(length(x)^2)-(1-simpson(table(y)))*(length(y)^2))/(2*length(x)*length(y)))
    
  }else{
    
    print ("Use simpson for numeric vectors")
    
  }
      
}








