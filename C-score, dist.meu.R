###função criada por Cristian de Sales Dambros em 20/11/2009


###Retorna o índice C-score (Roberts & Stone 1990) para a dada matriz


cscore<-
function (x) 
{
        Nsites <- nrow(x)
        S <- colSums(x)
        R <- ncol(x)
	  P<-((R * (R - 1))/2)
        ab <- 0
        for (i in 1:R) {
            for (j in (1:R)[-i]) {
if(j<i){next}else{
                Q <- sum(x[, i] * x[, j])
               ab <-ab+ ((S[i] - Q) * (S[j] - Q))/P}
            }}

ab

}


######################

dist.COR<-
function (x) 
{
x<-x[,colSums(x)!=nrow(x)]
        R <- ncol(x)
        ab <- 0
        for (i in 1:R) {
            for (j in (1:R)[-i]) {

if(j<i){next}else{
                

ab<-ab+cor(x[,i],x[,j],method="pearson")
          }}}

ab

}

