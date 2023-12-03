library(Matrix)
library(irlba)
library(AUC)
library(entropy)
library(igraph)
library(stringr)

reg.SP <- function(A,K,tau=1,lap=FALSE,nstart=30,iter.max=100){
  avg.d <- mean(colSums(A))
  A.tau <- A + tau*avg.d/nrow(A)
  if(!lap){SVD <- irlba(A.tau,nu=K,nv=K)}else{
    d.tau <- colSums(A.tau)
    L.tau <- diag(1/sqrt(d.tau))%*%A.tau%*%diag(1/sqrt(d.tau))
    #SVD <- svd(L.tau,nu=K,nv=K)
    SVD <- irlba(L.tau,nu=K,nv=K)
  }
  km <- kmeans(SVD$v[,1:K],centers=K,nstart=nstart,iter.max=iter.max)#,algorithm="Lloyd")
  return(list(cluster=km$cluster,loss=km$tot.withinss))
}

reg.SSP <- function(A,K,tau=1,lap=FALSE,nstart=30,iter.max=100){
  avg.d <- mean(colSums(A))
  A.tau <- A + tau*avg.d/nrow(A)
  if(!lap){SVD <- irlba(A.tau,nu=K,nv=K)
  V <- SVD$v[,1:K]
  V.norm <- apply(V,1,function(x)sqrt(sum(x^2)))
  V.normalized <- diag(1/V.norm)%*%V
  }else{
    d.tau <- colSums(A.tau)
    L.tau <- diag(1/sqrt(d.tau))%*%A.tau%*%diag(1/sqrt(d.tau))
    #SVD <- svd(L.tau,nu=K,nv=K)
    SVD <- irlba(L.tau,nu=K,nv=K)
    V <- SVD$v[,1:K]
    V.norm <- apply(V,1,function(x)sqrt(sum(x^2)))
    V.normalized <- diag(1/V.norm)%*%V
  }
  km <- kmeans(V.normalized,centers=K,nstart=nstart,iter.max=iter.max)#,algorithm="Lloyd")
  return(list(cluster=km$cluster,loss=km$tot.withinss))
}

holdout.evaluation.fast.all <- function(holdout.index,A,max.K,tau=0,dc.est=1,p.sample=1,kappa=NULL){
  n <- nrow(A)
  edge.index <- which(upper.tri(A))
  edge.n <- length(edge.index)
  A.new <- matrix(0,n,n)
  A.new[upper.tri(A.new)] <- A[edge.index]
  A.new[edge.index[holdout.index]] <- NA
  A.new <- A.new + t(A.new)
  degrees <- colSums(A.new,na.rm=TRUE)
  no.edge <- 0
  no.edge <- sum(degrees==0)
  
  Omega <- which(is.na(A.new))
  non.miss <- which(!is.na(A.new))
  #A.new[non.miss] <- A.new[non.miss] + 0.5
  SVD.result <- iter.SVD.core.fast.all(A.new,max.K,p.sample=p.sample)
  dc.block.sq.err <-  dc.loglike <- roc.auc <- bin.dev <- block.sq.err <- impute.sq.err <- loglike <- rep(0,max.K)
  
  for(k in 1:max.K){
    # print(k)
    #print(fast)
    tmp.est <- SVD.result[[k]]
    A.approx <- tmp.est$A.thr
    impute.sq.err[k] <- sum((A.approx[Omega]-A[Omega])^2)
    response <- A[edge.index[holdout.index]]#A[Omega]
    predictors <- A.approx[edge.index[holdout.index]]#A.approx[Omega]
    # print("AUC claculation")
    #print(system.time(tmp.roc <- pROC::roc(response=response,predictor=predictors)))
    #print(length(unique(predictors)))
    aa <- roc(predictions=predictors,labels=factor(response))
    #tmp.roc.smooth <- smooth(tmp.roc,method="binormal")
    roc.auc[k] <- auc(aa)#as.numeric(tmp.roc$auc)
    #print(tmp.roc$auc)
    #print(auc(aa))
    #roc.auc[k] <- as.numeric(tmp.roc.smooth$auc)
    trunc.predictors <- predictors
    trunc.predictors[predictors>(1-1e-6)] <- 1-1e-6
    trunc.predictors[predictors<1e-6] <- 1e-6
    bin.dev[k] <- sum((response-trunc.predictors)^2)#-sum(response*log(trunc.predictors)) - sum((1-response)*log(1-trunc.predictors))
    if(k==1){
      pb <- (sum(A.new,na.rm=TRUE)+1)/(sum(!is.na(A.new)) -sum(!is.na(diag(A.new)))+1)
      if(pb < 1e-6) pb <- 1e-6
      if(pb > 1-1e-6) pb <- 1-1e-6
      A.Omega <- A[Omega]
      block.sq.err[k] <- sum((pb-A[Omega])^2)
      loglike[k] <- -sum(A.Omega*log(pb)) - sum((1-A.Omega)*log(1-pb))
      
    }
    
    #U.approx <- eigen(A.approx)$vectors[,1:k]
    # print("SBM calculation")
    #print(k)
    #print(dim(tmp.est$SVD$v))
    ptm <- proc.time()
    if(k==1) {U.approx <- matrix(tmp.est$SVD$v,ncol=k)}else{
      U.approx <- tmp.est$SVD$v[,1:k]
      if(tau>0){
        A.approx <- A.approx + tau*mean(colSums(A.approx))/n
        d.approx <- colSums(A.approx)
        L.approx <- diag(1/sqrt(d.approx))%*%A.approx%*%diag(1/sqrt(d.approx))
        A.approx.svd <- irlba(L.approx,nu=k,nv=k)
        U.approx <- A.approx.svd$v[,1:k]
      }
    }
    
    km <- kmeans(U.approx,centers=k,nstart=30,iter.max=30)
    B <- matrix(0,k,k)
    Theta <- matrix(0,n,k)
    for(i in 1:k){
      for(j in i:k){
        N.i <- which(km$cluster==i)
        N.j <- which(km$cluster==j)
        if(i!=j){
          B[i,j] <- B[j,i] <- (sum(A.new[N.i,N.j],na.rm=TRUE)+1)/(sum(!is.na(A.new[N.i,N.j]))+1)
        } else{
          #print(max(N.i))
          #print(max(N.j))
          #print(dim(A.new))
          B[i,j] <- B[j,i] <- (sum(A.new[N.i,N.j],na.rm=TRUE)+1)/(sum(!is.na(A.new[N.i,N.j])) -sum(!is.na(diag(A.new[N.i,N.j])))+1)
        }
        
      }
      Theta[N.i,i] <- 1
    }
    P.hat <- Theta%*%B%*%t(Theta)
    diag(P.hat) <- 0
    block.sq.err[k] <- sum((P.hat[Omega]-A[Omega])^2)
    P.hat.Omega <- P.hat[Omega]
    A.Omega <- A[Omega]
    P.hat.Omega[P.hat.Omega < 1e-6] <- 1e-6
    P.hat.Omega[P.hat.Omega > (1-1e-6)] <- 1-1e-6
    loglike[k] <- -sum(A.Omega*log(P.hat.Omega)) - sum((1-A.Omega)*log(1-P.hat.Omega))
    # print(proc.time() - ptm)
    #### Degree correct model
    V <- U.approx
    # print("DCSBM calculation")
    ptm <- proc.time()
    #V.norms <- apply(V,1,function(x) sqrt(sum(x^2)))
    if(k==1) {V.norms <- as.numeric(abs(V))}else{
      V.norms <- apply(V,1,function(x) sqrt(sum(x^2)))
    }
    
    iso.index <- which(V.norms==0)
    Psi <- V.norms
    Psi <- Psi / max(V.norms)
    inv.V.norms <- 1/V.norms
    inv.V.norms[iso.index] <- 1
    
    V.normalized <- diag(as.numeric(inv.V.norms))%*%V
    
    #V.norms <- apply(V,1,function(x) sqrt(sum(x^2)))
    #Psi <- V.norms
    #Psi <- Psi / max(V.norms)
    #V.normalized <- diag(1/V.norms)%*%V
    #Psi.outer <- outer(Psi,Psi)
    if(k==1){
      if(dc.est>1){
        B <- sum(A.new,na.rm=TRUE)+0.01
        
        partial.d <- colSums(A.new,na.rm=TRUE)
        partial.gd <- B
        phi <- rep(0,n)
        B.g <- partial.gd
        phi <- as.numeric(partial.d/B.g)
        B <- B/p.sample
        P.hat <- t(t(matrix(B,n,n)*phi)*phi)
        #P.hat <- diag(phi)%*%matrix(B,n,n)%*%diag(phi)
        diag(P.hat) <- 0
      }
      dc.block.sq.err[k] <- sum((pb-A[Omega])^2)
      P.hat.Omega <- P.hat[Omega]
      A.Omega <- A[Omega]
      P.hat.Omega[P.hat.Omega < 1e-6] <- 1e-6
      P.hat.Omega[P.hat.Omega > (1-1e-6)] <- 1-1e-6
      
      dc.loglike[k] <- -sum(A.Omega*log(P.hat.Omega)) - sum((1-A.Omega)*log(1-P.hat.Omega))
      
      
    }else{
      km <- kmeans(V.normalized,centers=k,nstart=30,iter.max=30)
      if(dc.est>1){
        B <- matrix(0,k,k)
        Theta <- matrix(0,n,k)
        for(i in 1:k){
          for(j in 1:k){
            N.i <- which(km$cluster==i)
            N.j <- which(km$cluster==j)
            B[i,j] <- sum(A.new[N.i,N.j],na.rm=TRUE)+0.01
          }
          Theta[N.i,i] <- 1
        }
        Theta <- Matrix(Theta,sparse=TRUE)
        partial.d <- colSums(A.new,na.rm=TRUE)
        partial.gd <- colSums(B)
        phi <- rep(0,n)
        B.g <- Theta%*%partial.gd
        phi <- as.numeric(partial.d/B.g)
        B <- B/p.sample
        tmp.int.mat <- Theta*phi
        P.hat <-as.matrix(tmp.int.mat%*%B%*%t(tmp.int.mat))
        #P.hat <- diag(phi)%*%Theta%*%B%*%t(Theta)%*%diag(phi)
        diag(P.hat) <- 0
      }
      dc.block.sq.err[k] <- sum((P.hat[Omega]-A[Omega])^2)
      P.hat.Omega <- P.hat[Omega]
      A.Omega <- A[Omega]
      P.hat.Omega[P.hat.Omega < 1e-6] <- 1e-6
      P.hat.Omega[P.hat.Omega > (1-1e-6)] <- 1-1e-6
      dc.loglike[k] <- -sum(A.Omega*log(P.hat.Omega)) - sum((1-A.Omega)*log(1-P.hat.Omega))
    }
    # print(proc.time() - ptm)
    
    
    
  }
  return(list(impute.sq.err=impute.sq.err,block.sq.err=block.sq.err,loglike=loglike,roc.auc=roc.auc,no.edge=no.edge,dc.block.sq.err=dc.block.sq.err,dc.loglike=dc.loglike,bin.dev=bin.dev))
}


iter.SVD.core.fast.all <- function(A,Kmax,tol=1e-5,max.iter=100,sparse=TRUE,init=NULL,verbose=FALSE,tau=0,p.sample=1){
  if(sparse) A <- Matrix(A,sparse=TRUE)
  avg.p <- mean(as.numeric(A),na.rm=TRUE)
  cap <- 1#kappa*avg.p
  A[which(is.na(A))] <- 0
  A <- A/p.sample
  #svd.new <- svd(A,nu=K,nv=K)
  #print("begin SVD")
  svd.new <- irlba(A,nu=Kmax,nv=Kmax)
  #print("end SVD")
  result <- list()
  for(K in 1:Kmax){
    # print(K)
    if(K==1){
      A.new <- svd.new$d[1]*matrix(svd.new$u[,1],ncol=1)%*%t(matrix(svd.new$v[,1],ncol=1))
    }else{
      A.new <- A.new + svd.new$d[K]*matrix(svd.new$u[,K],ncol=1)%*%t(matrix(svd.new$v[,K],ncol=1))
    }
    A.new.thr <- A.new
    A.new.thr[A.new < 0+tau] <- 0+tau
    A.new.thr[A.new >cap] <- cap
    
    tmp.SVD <- list(u=svd.new$u[,1:K],v=svd.new$v[,1:K],d=svd.new$d[1:K])
    result[[K]] <- list(iter=NA,SVD=tmp.SVD,A=A.new,err.seq=NA,A.thr=A.new.thr)
  }
  return(result)
  
}

ECV.block <- function(A,max.K,cv=NULL,B=3,holdout.p=0.1,tau=0,dc.est=2,kappa=NULL){
  n <- nrow(A)
  edge.index <- which(upper.tri(A))
  edge.n <- length(edge.index)
  holdout.index.list <- list()
  if(is.null(cv)){
    holdout.n <- floor(holdout.p*edge.n)
    
    for(j in 1:B){
      holdout.index.list[[j]] <- sample(x=edge.n,size=holdout.n)
    }
  }else{
    sample.index <- sample.int(edge.n)
    max.fold.num <- ceiling(edge.n/cv)
    fold.index <- rep(1:cv,each=max.fold.num)[edge.n]
    cv.index <- fold.index[sample.index]
    B <- cv
    for(j in 1:B){
      holdout.index.list[[j]] <- which(cv.index==j)
    }
  }
  #print(fast)
  result <- lapply(holdout.index.list,holdout.evaluation.fast.all,A=A,max.K=max.K,tau=tau,dc.est=dc.est,p.sample=1-holdout.p,kappa=kappa)
  dc.block.err.mat <- dc.loglike.mat <- bin.dev.mat <- roc.auc.mat <- impute.err.mat <- block.err.mat <- loglike.mat <- matrix(0,nrow=B,ncol=max.K)
  no.edge.seq <- rep(0,B)
  Omega.list <- A.list <- Imputed.A.list <- list()
  for(b in 1:B){
    impute.err.mat[b,] <- result[[b]]$impute.sq.err
    block.err.mat[b,] <- result[[b]]$block.sq.err
    loglike.mat[b,] <- result[[b]]$loglike
    roc.auc.mat[b,] <- result[[b]]$roc.auc
    bin.dev.mat[b,] <- result[[b]]$bin.dev
    no.edge.seq[b] <- result[[b]]$no.edge
    dc.block.err.mat[b,] <- result[[b]]$dc.block.sq.err
    dc.loglike.mat[b,] <- result[[b]]$dc.loglike
    
  }
  
  
  output <- list(impute.err=colMeans(impute.err.mat),l2=colMeans(block.err.mat),dev=colSums(loglike.mat),auc=colMeans(roc.auc.mat),dc.l2=colMeans(dc.block.err.mat),dc.dev=colSums(dc.loglike.mat),sse=colMeans(impute.err.mat),auc.mat=roc.auc.mat,dev.mat=loglike.mat,l2.mat=block.err.mat,SSE.mat=impute.err.mat,dc.dev.mat=dc.loglike.mat,dc.l2.mat=dc.block.err.mat)
  
  if(min(output$dev)>min(output$dc.dev)){
    dev.model <- paste("DCSBM",which.min(output$dc.dev),sep="-")
  }else{
    dev.model <- paste("SBM",which.min(output$dev),sep="-")
  }
  if(min(output$l2)>min(output$dc.l2)){
    l2.model <- paste("DCSBM",which.min(output$dc.l2),sep="-")
  }else{
    l2.model <- paste("SBM",which.min(output$l2),sep="-")
  }
  output$l2.model <- l2.model
  output$dev.model <- dev.model
  
  return(output)
}



A<-read.csv("A_full.csv",header=FALSE)
A<-as.matrix(A)

cl <- makeCluster(6)
registerDoParallel(cl)

### when the dataset is cora
max.K = 14

start_time <- Sys.time()
ecv_para<-foreach(p=c(0.05,0.1,0.15,0.2,0.25),.packages=c("igraph","Matrix","irlba","AUC","entropy","stringr"),.combine="rbind") %dopar%{
  ecv<-ECV.block(A, max.K, cv = NULL, B = 3, holdout.p = p, tau = 0, dc.est = 2, kappa = NULL)
  
  ncom_est<-as.numeric(str_sub(ecv$l2.model,-1))
  ncom_est
}
ncom<-round(estimate_mode(ecv_para))

end_time <- Sys.time()
time=end_time - start_time

print(c("Estimation Time of full:",time))
print(c("Estimation of full:",ncom))


