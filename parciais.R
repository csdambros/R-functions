####################################################
#Function to extract partial scores on a multiple regression
#CSDambros February-02-2010


#################################
#Use: parciais(table)

#Parameters:

	# table   A table with the first column as the response variable and the remaining (at least 1) as predictor variables (see function cbind to concatenate columns)  


#output:

	# table with pairs of partials. The first 2 columns are the partials from the response and the first predictor variable with respect to the other variables. The third and fourth columns are the partials for the response and second predictor variable and so on...


#################################
#Function

parciais<-function(table){

aha<-table[,1]
for(i in 1:(ncol(table)-1)){
primparc<-lm(table[,i+1]~table[,-c(1,i+1)])$residuals
segparc<-lm(table[,1]~table[,-c(1,i+1)])$residuals
aha<-cbind(aha,primparc,segparc)
}
colnames(aha)<-c("tabela",rep(c("X parcial","Y parcial"),ncol(table)-1))
aha[,-1]
}


#end
#parciais(table)