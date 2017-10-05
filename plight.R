

cat("Caro aluno de Ecologia Geral,\n\nVocê acaba de instalar um conjunto de funções\ndesenvolvidas especialmente para sua aula prática.\nEstas funções lhe ajudarão a determinar com precisão a\ncobertura de fungos ou bactérias presentes na sua\namostra.\n\nTenha uma boa aula!\nProf. Cristian")


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

im1<-stack("../Aula10.MoldenBread/CristianDambros2.jpg")
extent(im1)<-c(0,100,0,100)
plotRGB(im1)
#grid(20,20,col = 1)
abline(v=seq(0,100,length.out = nx+1),lty=3)
abline(h=seq(0,100,length.out = ny+1),lty=3)
text(seq(2.5,100-2.5,length.out = nx),rep(seq(2.5,100-2.5,length.out = ny),each=nx),1:(nx*ny),cex=cex)
}


