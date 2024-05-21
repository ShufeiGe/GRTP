
## case2.
mcmc_infer <- function(partition.best,V.all,cts_z,n.iter,burnin.iter,save.all.par=FALSE){
  require(MCMCpack)
  require(MASS)
  Jt <- length(partition.best$Polytopes) 
  d <- dim(V.all)[2]
  
  Aj.all <- list()
  Bj.all <- list()
  Aj.inv.all <- list()
  AjinvBjt.all <- list()
  beta_j.all <- list()
  for (j in 1:Jt){
    idj <- partition.best$Polytopes[[j]]$V.ID
    idj2 <- idj[!is.na(cts_z[idj])]
    V.idj2.nonNA <- V.all[idj2,]
    cts_zj.nonNA <- cts_z[idj2]
    Aj.all[[j]] <- diag(1/a_h,nrow = d,ncol=d)+t(V.idj2.nonNA)%*%(V.idj2.nonNA)
    Bj.all[[j]] <- as.matrix(apply((V.idj2.nonNA*matrix(rep(cts_zj.nonNA,d),ncol=d,byrow = FALSE)),2,sum))
    
    Aj.M.eigen <- eigen(Aj.all[[j]],symmetric=TRUE)
    Aj.inv.all[[j]] <- Aj.M.eigen$vectors%*%diag(1/Aj.M.eigen$values)%*%t(Aj.M.eigen$vectors)
    AjinvBjt.all[[j]] <- Aj.inv.all[[j]]%*%Bj.all[[j]]
    beta_j.all[[j]] <- matrix(NA,nrow=n.iter,ncol=d)
  }
  

  

  sigma2 <- rep(NA,n.iter)
  sigma2.temp <- 1
  alpha.star <-  alpha_h+(Jt*d+length(cts_z[id.train]))/2   
  beta.star.part1 <- beta_h+sum(cts_z^2,na.rm=TRUE)/2
  
  for(i in 1:n.iter){
    beta.star.part2 <- 0
    for(j in 1:Jt){
      muj <-  AjinvBjt.all[[j]] 
      sigmaj <- sigma2.temp*Aj.inv.all[[j]]
      beta_j.temp <- mvrnorm(n=1,mu=muj,Sigma=sigmaj)
      beta_j.all[[j]][i,] <- beta_j.temp
      beta.star.part2 <- beta.star.part2 + t(beta_j.temp)%*%Aj.all[[j]]%*%beta_j.temp/2-t(Bj.all[[j]])%*%beta_j.temp/2
    }
    beta.star <- beta.star.part1+beta.star.part2
    sigma2[i] <-  MCMCpack::rinvgamma(1,shape=alpha.star,scale=beta.star)
    sigma2.temp <- sigma2[i] 
      
  }
 
  sd.sig <- sd(sigma2[(burnin.iter+1):n.iter])
  sd.beta <- sapply(beta_j.all,function(x){apply(as.matrix(x[(burnin.iter+1):n.iter,]),2,sd)})
  mean.sig <- mean(sigma2[(burnin.iter+1):n.iter])
  mean.beta <- sapply(beta_j.all,function(x){apply(as.matrix(x[(burnin.iter+1):n.iter,]),2,mean)})
                   
  ID.temp <- list()
  cts_z.predict <- rep(NA,length(cts_z))
  if(Jt==1){
    ID.temp[[1]] <- sapply(partition.best$Polytopes, "[[", 1)
  }else{
    ID.temp <- sapply(partition.best$Polytopes, "[[", 1)
  }
  
  for(id.j in 1:Jt){
    cts_z.predict[ID.temp[[id.j]]] <-  V.all[ID.temp[[id.j]],]%*%mean.beta[,id.j]
  }
  pred.t <- matrix(cts_z.predict,ncol=1)
  
  
 if(save.all.par){
   out <- list(all=list(sigma2=sigma2,beta=beta_j.all),
               pred.z=pred.t,
               sd.sig=sd.sig,
               sd.beta=sd.beta,
               mean.sig=mean.sig,
               mean.beta=mean.beta)
 }else{
   out <-  list(pred.z=pred.t,
                sd.sig=sd.sig,
                sd.beta=sd.beta,
                mean.sig=mean.sig,
                mean.beta=mean.beta)
 }
 
  return(out)
}



#library(MCMCpack)
#partition.best <- Partition.t[[which.max(Wt)]]
#n.iter <- 5000
#burnin.iter < -1000
#test.inf <- mcmc1(partition.best,V.all,cts_z,n.iter,burnin.iter)
#sigma2 <- test.inf$all$sigma2
#mu_j <- test.inf$all$mu_j
#n.polytp <- length(partition.best$Polytopes)
#par(mfrow=c(3,3))
#plot(c(sigma2),type="l")
#(j_v <- sample(1:n.polytp,size=8))
#for(j in j_v){
#  plot(mu_j[,j],type="l",ylab=paste0("polytope ",j))
#}
# 