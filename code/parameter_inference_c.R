## case1.
mcmc_infer <- function(partition.best,V.all,cts_z,n.iter,burnin.iter,save.all.par=FALSE){
  require(MCMCpack)
  aj_v <- sapply(partition.best$Polytopes, "[[", 5)
  bj_v <- sapply(partition.best$Polytopes, "[[", 6)
  Jt <- length(aj_v)
  mu_j <- matrix(NA,nrow=n.iter,ncol=Jt)
  sigma2 <- rep(NA,n.iter)
  sigma2.temp <- 0.1
  
  alpha.star <-  alpha_h+(Jt+length(cts_z[id.train]))/2   
  beta.star.part1 <- beta_h+Jt*mu_h^2/2/a_h+sum(cts_z^2,na.rm=TRUE)/2
  
  for(i in 1:n.iter){
    mu_j_t <- rep(NA,Jt)
    for(j in 1:Jt){
      mu_j_t[j]<- rnorm(1,mean=bj_v[j]/aj_v[j],sd=sqrt(a_h*sigma2.temp/aj_v[j]))
    }
    mu_j[i,] <- mu_j_t
    beta.star.part2 <-  mu_j_t^2%*%aj_v/2/a_h- mu_j_t%*%bj_v/a_h
    
    beta.star <- beta.star.part1+beta.star.part2
    
    sigma2.temp<-  MCMCpack::rinvgamma(1,shape=alpha.star,scale=beta.star)
    sigma2[i]  <- sigma2.temp
  }
  
  mu_j.M.out <- as.matrix(mu_j[(burnin.iter+1):n.iter,])
  mean.sig <- mean(sigma2[(burnin.iter+1):n.iter])
  sd.sig <- sd(sigma2[(burnin.iter+1):n.iter])
  mean.mu <- apply(mu_j.M.out,2,mean)
  sd.mu <- apply(mu_j.M.out,2,sd)
  
  ID.temp <- list()
  cts_z.predict <- rep(NA,length(cts_z))
  if(Jt==1){
    ID.temp[[1]] <- sapply(partition.best$Polytopes, "[[", 1)
  }else{
    ID.temp <- sapply(partition.best$Polytopes, "[[", 1)
  }

  for(id.j in 1:length(ID.temp)){
    cts_z.predict[ID.temp[[id.j]]] <-  mean.mu[id.j]
  }
  
  pred.t <- matrix(cts_z.predict,ncol=1)
  
  
  if(save.all.par){
    out <- list(all=list(sigma2=sigma2,mu_j=mu_j),
                pred.z =pred.t,
                sd.sig=sd.sig,
                sd.mu=sd.mu,
                mean.sig=mean.sig,
                mean.mu=mean.mu)
  }else{
    out <- list(pred.z =pred.t,
                sd.sig=sd.sig,
                sd.mu=sd.mu,
                mean.sig=mean.sig,
                mean.mu=mean.mu)
  }
  return(out)
}
