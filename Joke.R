
resp<-numeric(0)
for(i in 40:130){
resp<-c(rep(i,100),resp)
}


cexpre<-(1:10)/5
pchpre<-c(21,2:10)
colpre<-c(0,1:9)



plot.new()
plot(rep(1:100,91),resp,type="n",ylim=c(0,130))
points((1:10)*10,rep(1,10),cex=cexpre,bg=1,pch=21)
points((1:10)*10,rep(10,10),cex=2,bg=0,pch=pchpre)
points((1:10)*10,rep(20,10),cex=2,bg=colpre,pch=21)
rect(5,30,9,30);points(10,30);rect(16,30,20,30);points(20,30,pch=4)
#points((3:10)*10,rep(30,8));text((3:10)*10-3,rep(30,8),(1:8)*3)



resu<-(0)
resu2<-(0)
resux<-(0)
resuy<-(0)
plo=0

x<-c((1:10)*10,(1:10)*10,(1:10)*10,c(10,20),rep(1:100,91))
y<-c(rep(0,10),rep(10,10),rep(20,10),c(30,30),resp)


while(1!=2){

saved<-identify(x,y,n=1,plot=F)

resux<-x[saved[length(saved)]]
resuy<-y[saved[length(saved)]]

resu<-c(resu,resux)
resu2<-c(resu2,resuy)

cex=resu[resu2<5][length(resu[resu2<5])][length(resu[resu2<5][length(resu[resu2<5])])]/10
pch=resu[resu2<15][resu2[resu2<15]>5][length(resu[resu2<15][resu2[resu2<15]>5])]/10
col=resu[resu2<25][resu2[resu2<25]>15][length(resu[resu2<25][resu2[resu2<25]>15])]/10
plo<-resu[resu2==30][length(resu[resu2==30])]


if(resuy>35){
points(resux,resuy,cex=cexpre[cex],col=1,pch=pchpre[pch],bg=colpre[col])

if(sum(plo)==10){

points(resu[length(resu):(length(resu)-1)][resu2[length(resu2):(length(resu2)-1)]>30],resu2[length(resu2):(length(resu2)-1)][resu2[length(resu2):(length(resu2)-1)]>30],cex=cexpre[cex],col=colpre[col],pch=pchpre[pch],bg=colpre[col],type="o")

}
}
}

