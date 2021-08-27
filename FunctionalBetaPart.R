# Recreates the functional beta part function from the package betapart in R (Baselga et al. 2017) in order to run much faster.
# The script also implements methods to plot communities in 2d and 3d functional spaces
# Author: CSDambros

#' full
#' 
#' @description Converts vector object into simetric distance matrix. Modified from the full function in the ecodist package (Goslee and Urban 2007)
#' @author Cristian Dambros 
#' @param v Vector with distance values.
#' @param names Names of objects from which distances were calculated (usually rownames in original matrix).
#'
#' @return Returns square matrix with lower and upper triangles.
#' @export
#'
#' @examples

full<-function (v,names=NULL) 
{
  n <- (1 + sqrt(1 + 8 * length(v)))/2
  if (abs(n - round(n)) > 1e-07) 
    stop("Matrix not square.")
  n <- round(n)
  full <- matrix(0, n, n,dimnames = list(names,names))
  full[lower.tri(full)] <- v
  full2 <- t(full)
  diag(full2) <- 0
  full + full2
}



H_chset<-function(set1){
  require(rcdd)
  require(geometry)
  set1<-as.matrix(set1)
  set1rep <- d2q(cbind(0, cbind(1, set1)))
  polytope1 <- redundant(set1rep, representation = "V")$output
  H_chset1 <- scdd(polytope1, representation = "V")$output
  return(H_chset1)
}

intersectHchset <- function(H_chset1, H_chset2) {
  H_inter <- rbind(H_chset1, H_chset2)
  V_inter <- scdd(H_inter, representation = "H")$output
  vert_1n2 <- q2d(V_inter[, -c(1, 2)])
  #  coord_vert_inter <- rep(NA, ncol(set1))
  vol_inter <- 0
  if (is.matrix(vert_1n2) == T) 
    if (nrow(vert_1n2) > ncol(vert_1n2)) {
      #coord_vert_inter <- vert_1n2
      vol_inter <- convhulln(vert_1n2, "FA")$vol
    }
  res <- list(vol_inter = vol_inter)
  return(res)
}

#' Title
#'
#' @param fbc 
#' @param index.family 
#'
#' @return
#' @export
#'
#' @examples

functional.computations<-function (fbc,index.family="sorensen"){
  
  index.family <- match.arg(index.family, c("jaccard", "sorensen"))
  switch(index.family, sorensen = {
    funct.beta.sim <- fbc$min.not.shared/(fbc$min.not.shared + fbc$shared)
    funct.beta.sne <- ((fbc$max.not.shared - fbc$min.not.shared)/((2 * fbc$shared) + fbc$sum.not.shared)) * (fbc$shared/(fbc$min.not.shared + fbc$shared))
    funct.beta.sor <- fbc$sum.not.shared/(2 * fbc$shared + fbc$sum.not.shared)
    
    functional.pairwise <- data.frame(funct.beta.sim = funct.beta.sim, 
                                      funct.beta.sne = funct.beta.sne, funct.beta.sor = funct.beta.sor)
  }, jaccard = {
    funct.beta.jtu <- (2 * fbc$min.not.shared)/((2 * fbc$min.not.shared) + fbc$shared)
    funct.beta.jne <- ((fbc$max.not.shared - fbc$min.not.shared)/(fbc$shared + fbc$sum.not.shared)) * (fbc$shared/((2 * fbc$min.not.shared) + fbc$shared))
    funct.beta.jac <- fbc$sum.not.shared/(fbc$shared + fbc$sum.not.shared)
    functional.pairwise <- data.frame(funct.beta.jtu = funct.beta.jtu, 
                                      funct.beta.jne = funct.beta.jne, funct.beta.jac = funct.beta.jac)
  })
  
  return(functional.pairwise)
  
}


#' Title
#'
#' @param x 
#' @param traits 
#' @param multi 
#' @param warning.time 
#' @param return.details 
#' @param prefix 
#' @param ncores 
#' @param parallel 
#' @param useMPI 
#' @param parallel.breaks 
#'
#' @return
#' @export
#'
#' @examples

functional.betapart.core5<-function (x, traits, multi = TRUE, warning.time = TRUE, return.details = FALSE,prefix=NULL,ncores=NULL,parallel=FALSE,useMPI=FALSE){
  
  tinit<-proc.time()
  
  
  require(rcdd)
  require(geometry)
  
  # x = comm.test
  # traits=traits.test
  # multi = FALSE
  # warning.time = FALSE
  # return.details = FALSE
  # prefix=NULL
  # ncores=NULL
  # parallel=FALSE
  # useMPI=FALSE
  # parallel.breaks=3
  
  
  
  if(is.null(prefix)){
    
    prefix<-tempdir()
    
  }  
  
  file1<-paste0(prefix,"/step.fbc.txt")
  file2<-paste0(prefix,"/vert.txt")
  
  if (!is.matrix(x)) {
    x <- as.matrix(x)
  }
  if (!is.numeric(x)) 
    stop("The data in 'x' is not numeric.", call. = TRUE)
  xvals <- unique(as.vector(x))
  if (any(!is.element(xvals, c(0, 1)))) 
    stop("The 'x' table contains values other than 0 and 1: data should be presence/absence.", 
         call. = TRUE)
  if (!is.numeric(traits)) 
    stop("The data in 'traits' is not numeric.", call. = TRUE)
  if (any(is.na(traits))) 
    stop("NA are not allowed in 'traits'", call. = TRUE)
  if (ncol(x) != nrow(traits)) 
    stop("Number of species in 'x' and 'traits' must be identical", 
         call. = TRUE)
  D <- ncol(traits)
  Si <- apply(x, 1, sum)
  if (any(Si <= D)) 
    stop(paste("'community ", row.names(x)[which(Si <= D)], 
               " must contain at least ", D + 1, " species", sep = ""))
  
  
  N <- nrow(x)
  if (N < 2) 
    stop("Computing dissimilairty requires at least 2 communities", 
         call. = TRUE)
  nb.step <- 2
  if (multi == T) 
    nb.step <- N
  step.fbc <- as.data.frame(matrix("", nb.step, 1, dimnames = list(c("           FRi", 
                                                                     paste("intersection", 2:nb.step, sep = "_")), c("iteration"))))
  step.fbc[, 1] <- as.character(step.fbc[, 1])
  step.fbc[1, 1] <- paste("0/", N, sep = "")
  for (k in 2:nb.step) step.fbc[k, 1] <- paste("0/", choose(N, 
                                                            k), sep = "")
  FRi <- rep(NA, N)
  names(FRi) <- row.names(x)
  coord_vert_i <- list()
  
  #1. First for loop
  cat("Entered loop1\n")
  
  # This is fast, does not require speed up
  pb <- txtProgressBar(1,N,style = 3)
  
  for (i in 1:N) {
    #cat(paste0("loop ",i,"\n"))
    
    tr_i <- traits[which(x[i, ] == 1), ]
    vert0 <- convhulln(tr_i, paste0("Fx TO '",file2,"'"))
    vert1 <- scan(file2, quiet = T)
    verti <- (vert1 + 1)[-1]
    coord_vert_i[[i]] <- tr_i[verti, ]
    FRi[i] <- convhulln(tr_i[verti, ], "FA")$vol
    step.fbc["           FRi", 1] <- paste(i, "/", N, sep = "")
    step.fbc[, 1] <- as.character(step.fbc[, 1])
    write.table(step.fbc, file = file1, row.names = T, 
                col.names = F, sep = "\t")
    
    setTxtProgressBar(pb, i)
    
    
  }
  
  close(pb)
  
  sumFRi <- sum(FRi)
  
  #2. Second for loop
  cat("Entered loop2\n")
  
  ncomb<-(N*(N-1))/2
  
  ver<-apply(x,1,function(x){data.frame(traits[x>0,])})
  
  ###################
  
  if(isTRUE(parallel)){
    require(snow)
    require(parallel)
    if(is.null(ncores)){
      ncores<-detectCores()-1
    }
    #    if(exists("cl")){snow::stopCluster(cl)}
    
    if(isTRUE(useMPI)){
      cl<-snow::makeMPIcluster(ncores)
      cat("Started MPI cluster with",ncores,"nodes\n")
      #mpi.setup.rngstream()
    }else{
      cl<-snow::makeCluster(ncores,type="SOCK")
      cat("Started SOCK cluster with",ncores,"nodes\n")
    }
    cat(paste0("Done\nProcessing ",N, " rows..."))
    
    ver2<-snow::clusterApplyLB(cl,ver,H_chset)
    
    cat(paste0("Done\n"))
    
    snow::stopCluster(cl)
  }else{
    ver2<-lapply(ver,H_chset)
  }
  
  cat(paste0("Comparing ", ncomb, " communities\n"))
  
  
  fbc<-matrix(NA,ncomb,6,dimnames = list(NULL,c("shared","sum.not.shared","min.not.shared","max.not.shared","i","j")))
  
  pb <- txtProgressBar(1,ncomb,char = "*",style = 3)
  
  count<-function(i,j){(((i-1)*(i-2))/2)+j}
  
  for(i in 2:length(ver2)){
    for(j in 1:(i-1)){
      
      count2<-count(i,j)
      
      interij<-intersectHchset(ver2[[i]],ver2[[j]])
      
      sharedij <- interij$vol_inter
      not.sharedij <- FRi[i] - sharedij
      not.sharedji <- FRi[j] - sharedij
      
      sum.not.sharedij <- not.sharedij+not.sharedji
      max.not.sharedij <- pmax(not.sharedij,not.sharedji)
      min.not.sharedij <- pmin(not.sharedij,not.sharedji)
      
      fbc[count2,]<-c(shared=sharedij,sum.not.shared=sum.not.sharedij,min.not.shared=min.not.sharedij,max.not.shared=max.not.sharedij,i,j)
      
      setTxtProgressBar(pb, count(i,j))
      
      #cat(i,"/",j,":combination ",count2, " of ", ncomb," (", round(100*(count2/ncomb),2), "%) \n",sep = "")
      
    }
  }
  close(pb)
  
  fbc<-data.frame(fbc)
  
  fbc<-lapply(fbc[order(fbc[,6]),],full,names=rownames(x))
  
  class(fbc) <- "functional.betapart"
  
  print(proc.time()-tinit)
  
  return(fbc)
  
}

#' Title
#'
#' @param x 
#' @param traits 
#' @param index.family 
#' @param prefix 
#' @param ncores 
#' @param parallel 
#' @param useMPI 
#' @param parallel.breaks 
#'
#' @return
#' @export
#'
#' @examples
functional.beta.pair5<-function (x, traits, index.family = "sorensen",prefix=NULL,ncores=NULL,parallel=FALSE,useMPI=FALSE,parallel.breaks=NULL) 
{
  index.family <- match.arg(index.family, c("jaccard", "sorensen"))
  fbc <- x
  if (!inherits(x, "functional.betapart")) {
    fbc <- functional.betapart.core5(x, traits, multi = FALSE, 
                                     warning.time = FALSE, return.details = FALSE,prefix=prefix,ncores=ncores,parallel=parallel,useMPI=useMPI)
  }
  switch(index.family, sorensen = {
    funct.beta.sim <- fbc$min.not.shared/(fbc$min.not.shared + 
                                            fbc$shared)
    funct.beta.sne <- ((fbc$max.not.shared - fbc$min.not.shared)/((2 * 
                                                                     fbc$shared) + fbc$sum.not.shared)) * (fbc$shared/(fbc$min.not.shared + 
                                                                                                                         fbc$shared))
    funct.beta.sor <- fbc$sum.not.shared/(2 * fbc$shared + 
                                            fbc$sum.not.shared)
    functional.pairwise <- list(funct.beta.sim = as.dist(funct.beta.sim), 
                                funct.beta.sne = as.dist(funct.beta.sne), funct.beta.sor = as.dist(funct.beta.sor))
  }, jaccard = {
    funct.beta.jtu <- (2 * fbc$min.not.shared)/((2 * fbc$min.not.shared) + 
                                                  fbc$shared)
    funct.beta.jne <- ((fbc$max.not.shared - fbc$min.not.shared)/(fbc$shared + 
                                                                    fbc$sum.not.shared)) * (fbc$shared/((2 * fbc$min.not.shared) + 
                                                                                                          fbc$shared))
    funct.beta.jac <- fbc$sum.not.shared/(fbc$shared + fbc$sum.not.shared)
    functional.pairwise <- list(funct.beta.jtu = as.dist(funct.beta.jtu), 
                                funct.beta.jne = as.dist(funct.beta.jne), funct.beta.jac = as.dist(funct.beta.jac))
  })
  return(functional.pairwise)
}


functional.beta.pair5uni<-function (x, traits, index.family = "sorensen",prefix=NULL,ncores=NULL,parallel=FALSE,useMPI=FALSE,parallel.breaks=NULL,constant=FALSE) {
  
  n<-ncol(x)
  x<-cbind(x,x)
  if(constant){
    traits<-traits+runif(traits,0,diff(range(traits/10000)))
  }
  traits<-cbind(traits,0)
  traits<-rbind(traits,traits)
  
  traits[1:n,ncol(traits)]<-1
  resu<-functional.beta.pair5(x, traits, index.family,prefix,ncores,parallel ,useMPI,parallel.breaks)
  return(resu)
}



#' fspacepoly
#'
#' @description draw polygons showing funcitonal space for 2d traits
#' @param comm community matrix with species as columns and sites as rows
#' @param traits functional trait data. Species are rows
#'
fspacepoly<-function(comm,traits,new=TRUE,axes=TRUE,bg="white",alpha=0.5,col=1:10,border=1,plot.sites=FALSE,xlab=NULL,ylab=NULL,...){
  
  col<-rep_len(col,nrow(comm))
  border<-rep_len(border,nrow(comm))
  
  if(new){
    
    if(is.null(xlab)){xlab=colnames(traits)[1]}
    if(is.null(ylab)){ylab=colnames(traits)[2]}
    
    
    plot(range(traits[,1]),range(traits[,2]),type="n", axes=axes,xlab=xlab,ylab=ylab,...)
    
    
  }
  
  
  for(i in 1:nrow(comm)){
    A<-traits[comm[i,]>0,]
    P<-A[chull(A),]
    polygon(P,col=adjustcolor(col[i],alpha),border = adjustcolor(border[i],0.5))
    if(plot.sites){
    text(colMeans(P)[1],colMeans(P)[2],rownames(comm)[i],col=border[i])
    }
  }
  
}

#' fspacepoly3d
#'
#' @description draw polygons showing funcitonal space for 3d traits
#' @param comm community matrix with species as columns and sites as rows
#' @param traits functional trait data. Species are rows
#'
fspacepoly3d<-function(comm,traits,new=TRUE,axes=TRUE,bg="white",alpha=0.5,col=1:10,border=1,plot.sites=FALSE,...){
  require(geometry)
  require(rgl)
  
  col<-rep_len(col,nrow(comm))
  border<-rep_len(border,nrow(comm))
  
  if(new){
    rgl::clear3d()
    bg3d(bg)
    plot3d(traits,...)
  }
  for(i in 1:nrow(comm)){
    A<-traits[comm[i,]>0,]
    tr<-t(convhulln(A,options = "Tv"))
    rgl.triangles(A[tr,1],A[tr,2],A[tr,3],col=col[i],alpha=alpha,...)
    if(plot.sites){
      text3d(colMeans(A[tr,])[1],colMeans(A[tr,])[2],colMeans(A[tr,])[3],rownames(comm)[i],color=border[i])
    }
    
    #particles3d(A[tr,1],A[tr,2],A[tr,3],col=i,...)
    
  }
  
}


#' fspacepoly3d
#'
#' @description draw polygons showing funcitonal space for 3d traits
#' @param comm community matrix with species as columns and sites as rows
#' @param traits functional trait data. Species are rows
#'
fspacepoly4d<-function(comm,traits,pcoa=NULL,new=TRUE,alpha=0.5,bg="white",col=1:10,border=1,k=3,plot.species=FALSE,plot.sites=FALSE,arrows=TRUE,axes=TRUE,...){
  require(geometry)
  require(rgl)
  require(FD)
  require(vegan)
  
  col<-rep_len(col,nrow(comm))
  border<-rep_len(border,nrow(comm))
  
  
  if(k== 1){stop("Unidimensional not implemented")}
  if(k>3){warning("k higher than 3 not allowed, reduced to 3");k=3}
  
  if(is.null(pcoa)){
  cat("Calculating PCoA using",k,"dimensions\n")
  gdist<-gowdis(traits)
  
  
  pcoa<-cmdscale(gdist,k=k)
  colnames(pcoa)<-paste("PCoA", 1:k)
  }
#  cors<-(cor(pcoa,traits))*apply(abs(apply(pcoa,2,range))*1,2,min)
  cors<-cor(traits,pcoa)
  


  if(k == 2){
    
    #fspacepoly(comm,pcoa,axes = axes,col=col,border=border,plot.sites=plot.sites)
    fspacepoly(comm,pcoa,new=new,alpha=alpha,col=col,border=border,plot.sites=plot.sites,...)
    
    if(arrows){
      
      redfactor<-0.8
      
      xaxis<-pcoa[,1]*redfactor
      yaxis<-pcoa[,2]*redfactor
      
      xcor<-cors[,1]
      ycor<-cors[,2]
      
      scale<-diff(range(xaxis))/diff(range(yaxis))
      
      expand<-min(abs(c(min(yaxis)/ycor[ycor<0],
                        max(yaxis)/ycor[ycor>0],
                        min(xaxis)/xcor[xcor<0]/scale,
                        max(xaxis)/xcor[xcor>0]/scale)))
      
      
      #arrows(0,0,xcor*scale,ycor,length=0.1)
      arrows(0,0,xcor*scale*expand,ycor*expand,length=0.1)
      
      text(ordiArrowTextXY(cbind(xcor*scale*expand,ycor*expand), rescale = FALSE,labels = colnames(traits)),labels = colnames(traits))
      
    

      #tks<-apply(abs(apply(pcoa,2,range))*0.8,2,max)
      # text(ordiArrowTextXY(t(cors*tks), rescale = FALSE,labels = colnames(traits)),labels = colnames(traits))

    
      #text(cors[1,]*tks[1]/2,cors[2,]*tks[2]/2,colnames(traits),pos = 2)
      # arrows(rep(0,ncol(traits)),rep(0,ncol(traits)),cors[1,]*tks[1],cors[2,]*tks[2],length = 0.1)
      
      
      if(axes){
        
        zaxesvals<-pretty(c(xcor,ycor),n = 3)
        
        axis(3,zaxesvals*scale*expand,labels = zaxesvals)
        axis(4,zaxesvals*expand,labels = zaxesvals)
        
        # axis(3,seq(-tks[1],tks[1],length.out = 5),labels = seq(-1,1,length.out = 5))
        # axis(4,seq(-tks[2],tks[2],length.out = 5),labels = seq(-1,1,length.out = 5))
        # 
        # #p1<-pretty(cors[,1],5)
        #p2<-pretty(cors[,2],5)
        
#        axis(3,seq(min(cors.scl[,1]),max(cors.scl[,1]),along.with = p1),p1)
#       axis(4,seq(min(cors.scl[,2]),max(cors.scl[,2]),along.with = p2),p2)        
        
      }
      
    }
    
    if(plot.species){
      points(pcoa[,1],pcoa[,2],pch=21,bg=1,cex=0.5)
      text(pcoa[,1],pcoa[,2],cex=0.5,colnames(comm))
    }
    
  }
  
  if(k>2){
    
    #fspacepoly3d(comm,pcoa,new=TRUE,col=col,border=border,plot.sites=plot.sites)
    fspacepoly3d(comm,pcoa,new=new,alpha=alpha,col=col,border=border,plot.sites=plot.sites,...)
    
    if(axes){
      
      #axes3d(c('x--', 'y--', 'z--'),color= "black",alpha=1) 
      #axes3d(c('x--', 'y--', 'z--'),color= "blue") 
      
      redfactor<-0.8
      
      xaxis<-pcoa[,1]*redfactor
      yaxis<-pcoa[,2]*redfactor
      zaxis<-pcoa[,3]*redfactor
      
      xcor<-cors[,1]
      ycor<-cors[,2]
      zcor<-cors[,3]
      
      scalex<-diff(range(xaxis))/diff(range(xaxis,yaxis,zaxis))
      scaley<-diff(range(yaxis))/diff(range(xaxis,yaxis,zaxis))
      scalez<-diff(range(zaxis))/diff(range(xaxis,yaxis,zaxis))
      
      expand<-min(abs(c(min(xaxis)/xcor[xcor<0]/scalex,
                        max(xaxis)/xcor[xcor>0]/scalex,
                        min(yaxis)/ycor[ycor<0]/scaley,
                        max(yaxis)/ycor[ycor>0]/scaley,
                        min(zaxis)/zcor[zcor<0]/scalez,
                        max(zaxis)/zcor[zcor>0]/scalez)))
    
      zaxesvals<-pretty(c(xcor,ycor,zcor),n = 3)
      
      atx<-zaxesvals*scalex*expand
      aty<-zaxesvals*scaley*expand
      atz<-zaxesvals*scalez*expand
      
            axes3d('x++',at=atx[atx>min(pcoa[,1])&atx<max(pcoa[,1])],labels = zaxesvals,color="blue",alpha=1)
            axes3d('y++',at=aty[aty>min(pcoa[,2])&aty<max(pcoa[,2])],labels = zaxesvals,color="blue",alpha=1)
            axes3d('z-+',at=atz[atz>min(pcoa[,3])&atz<max(pcoa[,3])],labels = zaxesvals,color="blue",alpha=1)

#            title3d(xlab='PCoA 1',ylab='PCoA 2',zlab='PCoA 3',line = c(3,3,3),color="black",alpha=1)
    }
    

    if(arrows){
      texts3d(cors[,1]*scalex*expand*1.1,
              cors[,2]*scaley*expand*1.1,
              cors[,3]*scalez*expand*1.1,
              colnames(traits),color="black",alpha=1)
      
      for(i in 1:ncol(traits)){
        arrow3d(c(0,0,0),cors[i,]*c(scalex,scaley,scalez)*expand,color="blue",alpha=1)
      }
      

            arrow3d(c(0,0,0),cors[4,]*c(scalex,scaley,scalez)*expand,color="blue",alpha=1)

    }
    
    if(plot.species){
      points3d(pcoa[,1],pcoa[,2],pcoa[,3],color="black")
      text3d(pcoa[,1],pcoa[,2],pcoa[,3],colnames(comm))
    }
    
  }
  
}

