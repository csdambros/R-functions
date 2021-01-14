

cat("Caro aluno de Ecologia Geral,\n\nVoce acaba de instalar um conjunto de funcoes desenvolvidas especialmente para sua aula pratica.\nEstas funcoes lhe ajudarao a determinar com precisao a\ncobertura de fungos ou bacterias presentes na sua\namostra.\n\nTenha uma boa aula!\nProf. Cristian")

#require(rgdal)

plight<-function(file,thresh=0.5,showRGB=FALSE){
  
  require(raster)
  
  picture<-stack(file)
  values<-getValues(picture)
  
  pixWhite<-rowSums(values>thresh*255)==3
  percentage<-mean(pixWhite)
  cover<-sum(pixWhite)
  
  if(showRGB){
    
    picture[[1]][]<-1-pixWhite
    plot(picture[[1]])
    
  }
  
  
  return(list("cover(pixels)"=cover,"Notcover(pixels)"=length(pixWhite)-cover,percentage=percentage))
}

pdark<-function(file,thresh=0.5,showRGB=FALSE){
  
  require(raster)
  
  picture<-stack(file)
  values<-getValues(picture)
  
  pixWhite<-rowSums(values>thresh*255)==3
  percentage<-mean(pixWhite)
  cover<-sum(pixWhite)
  
  if(showRGB){
    
    picture[[1]][]<-1-pixWhite
    plot(picture[[1]])
    
  }
  
  
  return(list("cover(pixels)"=length(pixWhite)-cover,"Notcover(pixels)"=cover,percentage=percentage))
}


pred<-function(file,thresh=0.5,showRGB=FALSE){
  
  require(raster)
  
  picture<-raster(file,band=1)
  values<-getValues(picture)
  
  pixWhite<-values>(thresh*255)
  percentage<-mean(pixWhite)
  cover<-sum(pixWhite)
  
  if(showRGB){
    
    picture[]<-1-pixWhite
    plot(picture)
    
  }
  
  
  return(list("cover(pixels)"=length(pixWhite)-cover,"Notcover(pixels)"=cover,percentage=percentage))
}


pgreen<-function(file,thresh=0.5,showRGB=FALSE){
  
  require(raster)
  
  picture<-raster(file,band=2)
  values<-getValues(picture)
  
  pixWhite<-values>(thresh*255)
  percentage<-mean(pixWhite)
  cover<-sum(pixWhite)
  
  if(showRGB){
    
    picture[]<-1-pixWhite
    plot(picture)
    
  }
  
  
  return(list("cover(pixels)"=length(pixWhite)-cover,"Notcover(pixels)"=cover,percentage=percentage))
}

pblue<-function(file,thresh=0.5,showRGB=FALSE){
  
  require(raster)
  
  picture<-raster(file,band=1)
  values<-getValues(picture)
  
  pixWhite<-values>(thresh*255)
  percentage<-mean(pixWhite)
  cover<-sum(pixWhite)
  
  if(showRGB){
    
    picture[]<-1-pixWhite
    plot(picture)
    
  }
  
  
  return(list("cover(pixels)"=length(pixWhite)-cover,"Notcover(pixels)"=cover,percentage=percentage))
}


rasterGrid<-function(file,nx=20,ny=20,cex=0.3){

  require(raster)

im1<-stack(file)

op<-par(no.readonly = TRUE) 
  
extent(im1)<-c(0,100,0,100)
plotRGB(im1)
#grid(20,20,col = 1)
abline(v=seq(0,100,length.out = nx+1),lty=3)
abline(h=seq(0,100,length.out = ny+1),lty=3)
text(seq(2.5,100-2.5,length.out = nx),rep(seq(2.5,100-2.5,length.out = ny),each=nx),1:(nx*ny),cex=cex)
  
par(op)  
}


require(raster)

model.train<-function(folder=".",pattern=".png|.jpg"){
  require(raster)
  require(tree)
  
#  folder="../Aula10.MoldenBread/VideoLapses/"
#  extension=".png"
  pics<-list.files(folder,pattern = pattern)
  picsref<-list.files(folder,pattern = pattern,full.names = TRUE)
  
  cat(paste0("(",1:length(pics),") ",pics,"\n"))
  
  ref <- readline("Digite o número da foto de referência (treino):")
  cat("Arquivo",pics[as.integer(ref)],"selecionado para treino (ensinar o R)\n")
  
  picture<-stack(picsref[as.integer(ref)])

    {
plotRGB(picture,maxpixels=10000)

moldValues<-list()

extmx<-{}
extmy<-{}

i=1


stopmofo<-"S"

while(stopmofo=="S"){

cat("Selecione área com a primeira classe na foto")

  x<-"S"

  mold<-select(picture)
  
while(!x%in%c("S","N")){
x <- toupper(readline("Manter seleção (S/N)?"))
if(!x%in%c("S","N")){cat("Digite S para Sim ou N para Não")}
}
  if(x=="S"){
  moldValues[[i]]<-getValues(mold)
  i<-i+1
  

ext<-extent(mold)

extmx<-cbind(extmx,ext[c(1,2)])
extmy<-cbind(extmy,ext[c(3,4)])

rect(extmx[1,],extmy[1,],extmx[2,],extmy[2,],col=2,density = 20)

}

stopmofo<-""

while(!stopmofo%in%c("S","N")){
  stopmofo <- toupper(readline("Você gostaria de identificar mais regiões com a primeira classe (S/N)?"))
  if(!stopmofo%in%c("S","N")){cat("Digite S para Sim ou N para Não")}

  }

}

###
breadValues<-list()

extbx<-{}
extby<-{}

stopmofo<-"S"

while(stopmofo=="S"){
  
  cat("Selecione área com a segunda classe na foto")
  bread<-select(picture)
  
  
  x<-"S"
  
  while(!x%in%c("S","N")){
    x <- toupper(readline("Manter seleção (S/N)?"))
    if(!x%in%c("S","N")){cat("Digite S para Sim ou N para Não")}
  }
  
  if(x=="S"){
    breadValues[[i]]<-getValues(bread)
    i<-i+1
    
  ext<-extent(bread)
  
  extbx<-cbind(extbx,ext[c(1,2)])
  extby<-cbind(extby,ext[c(3,4)])
  
  rect(extbx[1,],extby[1,],extbx[2,],extby[2,],col=3,density = 20)
  
  }
  
  stopmofo<-""
  
  while(!stopmofo%in%c("S","N")){
    stopmofo <- toupper(readline("Você gostaria de identificar mais regiões com a segunda classe (S/N)?"))
    if(!stopmofo%in%c("S","N")){cat("Digite S para Sim ou N para Não")}
    
  }
  
}

#####

moldValues2<-do.call(rbind,moldValues)
breadValues2<-do.call(rbind,breadValues)
}

cat("Construindo modelo...")

resp<-rep(c("c1","c2"),c(nrow(moldValues2),nrow(breadValues2)))

m1<-tree(resp~.,data=data.frame(resp,rbind(moldValues2,breadValues2)))

cat("concluído\n")

return(m1)

}

model.classify<-function(m1,
                        folder=".",
                        pattern=".png|.jpg",
                        plot=FALSE){

require(raster)
require(tree)


#  folder="../Aula10.MoldenBread/VideoLapses/"
#  extension=".png"
pics<-list.files(folder,pattern = pattern)
picsref<-list.files(folder,pattern = pattern,full.names = TRUE)

cat("Mensurando cobertura nas fotos\n")

op<-par(no.readonly = TRUE)

if(plot){
  par(mfrow=c(length(pics),2))
}

resu<-matrix(NA,length(pics),5)
rownames(resu)<-pics
colnames(resu)<-c("npixels","c1","c2","c1.cover","c2.cover")


for (i in 1:length(pics)){

  picture<-stack(picsref[i])  
  allpic<-getValues(picture)

  colnames(allpic)<-labels(m1$terms)
  
  predicted<-predict(m1,newdata=data.frame(allpic))

  n<-nrow(predicted)
  c1<-sum(predicted[,"c1"]>=0.5)
  c2<-sum(predicted[,"c2"]>0.5)

  resu[i,]<-c(n,c1,c2,c1/n,c2/n)
  print(i)
  
  if(plot){
  
  picturecopy<-raster(picture)
  picturecopy<-setValues(picturecopy,predicted[,"c1"]>=0.5)

  plotRGB(picture,maxpixels=10000)
  plot(picturecopy,maxpixels=10000,axes=FALSE,box=FALSE)

  }
}
 
  par(op)
  return(resu)
   
}




require(raster)


select.pic<-function(file="."){

  pics<-basename(file)

  cat(paste0("(",1:length(pics),") ",pics),sep = "\n")
  
  ref <- readline("Type picture reference number (training):")
  cat("File",file[as.integer(ref)],"selected\n")

  return(file[as.integer(ref)])
  
}


select.region<-function(raster,plot.new=TRUE,col=2,density=20,maxpixels=10000000,...){

  if(plot.new){
    
    plotRGB(brick(raster),maxpixels=maxpixels)
    
  }
  
  
values<-list()

extmx<-{}
extmy<-{}

i=1

continue<-"Y"

while(continue=="Y"){
  
  cat("Select a region")
  x<-"Y"
  
  sel<-select(raster)
  
  while(!x%in%c("Y","N")){
    x <- toupper(readline("Keep selection (Y/N)?"))
    if(!x%in%c("Y","N")){cat("Type Y for Yes or N for No")}
  }
  
  if(x=="Y"){
    values[[i]]<-getValues(sel)
    i<-i+1
    
    ext<-extent(sel)
    
    extmx<-cbind(extmx,ext[c(1,2)])
    extmy<-cbind(extmy,ext[c(3,4)])
    
    rect(extmx[1,],extmy[1,],extmx[2,],extmy[2,],col=col,density = density,...)
  }
  
  continue<-""
  
  while(!continue%in%c("Y","N")){
    continue <- toupper(readline("Continue selecting other regions in this class (Y/N)?"))
    if(!continue%in%c("Y","N")){cat("Type Y for Yes or N for No")}
    
  }
  
}

return(values)

}


model.train<-function(file=".",nclasses=2,classnames=1:nclasses,maxpixels=500000){

  require(raster)
  require(tree)
  
  if(length(classnames)!=nclasses){
    
    stop("classnames length must equal nclasses")
    
  }
  
  picsref<-file
  
  if(!class(file)%in%c("raster","RasterStack","RasterBrick","RasterLayer")){
    
  if(dir.exists(file[1])){
    
    file<-list.files(file,full.names = TRUE)

  }
    
    if(length(file)>1){
    
    picsref<-select.pic(file = file)
    }
    
  }

picture<-stack(picsref)

plotRGB(picture,maxpixels=maxpixels)

resu<-list()

for(k in 1:nclasses){
  cat("Class",classnames[k],"\n")
  resu[[k]]<-select.region(raster=picture,col=k,plot.new=FALSE)
}

resu2<-lapply(resu,do.call,what = rbind)
resu3<-do.call(rbind,resu2)

n<-unlist(lapply(resu2,nrow))

resp<-factor(rep(classnames,n))

data<-data.frame(resp,resu3)

cat("Building model...")
m1<-tree(resp~.,data=data)
cat("Done\n")

return(m1)

}

model.classify<-function(m1,file=".",plot=FALSE){

require(raster)
require(tree)

  picsref<-file
  
  if(!class(file)%in%c("raster","RasterStack","RasterBrick","RasterLayer")){
    
  if(dir.exists(file[1])){
    
    picsref<-list.files(file,full.names = TRUE)

  }
    
    pics<-basename(picsref)
  
  }else{
      
    pics<-filename(file)    
    }
  
cat("Classifying regions in",length(pics),"files\n")

op<-par(no.readonly = TRUE)

 if(plot){
   par(mfrow=c(1,2))
 }

resu<-matrix(NA,length(pics),1+nlevels(m1$y))
rownames(resu)<-pics
colnames(resu)<-c("npixels",levels(m1$y))


for (i in 1:length(pics)){

  file2<-file
  
  if(!class(file)%in%c("raster","RasterStack","RasterBrick","RasterLayer")){
    
    file2<-stack(picsref[i])
  }
  
  allpic<-getValues(file2)
  
  colnames(allpic)<-labels(m1$terms)
  
  predicted<-predict(m1,newdata=data.frame(allpic),type="class")
    
  n<-length(predicted)
  predtable<-table(predicted)

  resu[i,]<-c(n,predtable)
  print(i)
  
  if(plot){
  
  picturecopy<-raster(file2)
  picturecopy<-setValues(picturecopy,predicted)

  plotRGB(file2)
  plotRGB(file2)
  plot(picturecopy,axes=FALSE,box=FALSE,add=TRUE,col=adjustcolor(c(NA,2:nlevels(predicted))),legend=FALSE)
  #plot(picturecopy,axes=FALSE,box=FALSE,add=TRUE,legend=FALSE)

  }
}
 
  par(op)
  return(resu)
   
}


# files<-list.files("CanopyCoverData",full.names = TRUE)
# 
# ver<-brick(files[1])
# ver2<-aggregate(ver,fact=6)
# 
# m1<-model.train(ver2)
# 
# model.classify(m1,ver2,plot=TRUE)
# 



# m1<-model.train(file = "../Aula10.MoldenBread",maxpixels = 10000)
# model.classify(m1,file = "../Aula10.MoldenBread",plot = TRUE)
# 
# ver<-brick("../Aula10.MoldenBread/canopy.jpg")
# ver2<-aggregate(ver,fact=4)
# 
# m1<-model.train(ver2,maxpixels = 10000000,nclasses = 3,classnames = c("sol","ceu","arvore"))
# model.classify(m1,ver2,plot = TRUE)
# 
# 
# getValues(ver2)
# 
# 
# 
# ver<-brick("../Aula10.MoldenBread/canopy.jpg")
# 
# 
# 
# 
# pgreen("../Aula10.MoldenBread/canopy.jpg",showRGB = TRUE)
# 
# ver<-stack("../Aula10.MoldenBread/canopy.jpg")
# ver2<-getValues(ver)
# 
# 
# 





