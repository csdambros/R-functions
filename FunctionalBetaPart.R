# Recreates the functional beta part function from the package betapart in R (Baselga et al. 2017) in order to run much faster.
# The script also implements methods to plot communities in 2d and 3d functional spaces


#' full
#' 
#' @description Converts vector object into simetric distance matrix. Modified from the full function in the ecodist package (Goslee and Urban 2007)
#' 
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


functional.beta.pair5uni<-function (x, traits, index.family = "sorensen",prefix=NULL,ncores=NULL,parallel=FALSE,useMPI=FALSE,parallel.breaks=NULL) 
{
  
  n<-ncol(x)
  x<-cbind(x,x)
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
fspacepoly<-function(comm,traits,new=TRUE){
  
  if(new){
    
    plot(range(traits[,1]),range(traits[,2]),type="n", xlab=colnames(traits)[1],ylab=colnames(traits)[2])
    
  }
  
  
  for(i in 1:nrow(comm)){
    A<-traits[comm[i,]>0,]
    P<-A[chull(A),]
    polygon(P,col=adjustcolor(i,0.2),border = adjustcolor(i,0.5))
    text(colMeans(P)[1],colMeans(P)[2],rownames(comm)[i],col=i)
  }
  
}

#' fspacepoly3d
#'
#' @description draw polygons showing funcitonal space for 3d traits
#' @param comm community matrix with species as columns and sites as rows
#' @param traits functional trait data. Species are rows
#'
fspacepoly3d<-function(comm,traits,new=TRUE){
  require(geometry)
  require(rgl)
  
  for(i in 1:nrow(comm)){
    A<-traits[comm[i,]>0,]
    tr<-t(convhulln(A,options = "Tv"))
    rgl.triangles(A[tr,1],A[tr,2],A[tr,3],col=i,alpha=0.5)
  }
  
}


