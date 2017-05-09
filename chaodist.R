# Functions to calculate dissimilarity between sets of communities using the metrics corrected for undersampling
# Based in Chao et al. (2005)
# Currently only Jaccard implemented (sorensen will probably be implemented in future versions)
# Uses only matrix operations (no for loops), it is slightly faster than the function available in the 'fossil' package
# Easy to use with multiple sampling sites (output in the same format as 'vegdist' in the 'vegan' package)
# Creator: CSDambros 06-Apr-2016 


chaodist<-function(comm,method="jaccard"){
  
  if(method!="jaccard"){warning("Only jaccard currently implemented, switching to jaccard")}
  
  mn<-rowSums(comm)
  pi<-comm/mn
  mndiv<-(mn-1)/mn
  
  commPA<-comm>0
  
  
  f1<-tcrossprod(comm==1,commPA)  #faster than f1<-(comm==1)%*%t(commPA)
  f2<-tcrossprod(comm==2,commPA)  #faster than f2<-(comm==2)%*%t(commPA)
  
  f2[f2==0]<-1 # 1 if no doubletons
  
  
  P1<-tcrossprod(pi,commPA) #faster than P1<-pi%*%t(commPA)
  P2<-t(t(f1/(2*f2))*mndiv)
  P3<-tcrossprod(pi,1-commPA) #faster than P3<-pi%*%t(1-commPA)
  
  
  U<-P1+P2*P3
  
  U[U>1]<-1
  
  UV<-t(U)*U
  
  Jmat<-UV/((U+t(U))-(UV))
  Jmat[is.nan(Jmat)]<-0
  
  return(1-as.dist(Jmat))
  
}

#USE
#
#chaodist(comm)
#comm is a matrix with objects (usually sampling sites) in rows and attributes as columns (usually species)
