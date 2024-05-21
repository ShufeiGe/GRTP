# Initialize the partition set & return the pair-wised distance matrix between points. 
# The pair-wised distance matrix will be saved to be served as a lookup table to avoid duplicated computation.
# mode.type : ==1 for constant case; other values for the linear case.
partition_initial <- function(V.all,cts_z,l.max,sample.type,cut.type,accept.rate.min,mode.type){
  n <- length(cts_z)
  V.temp.ID <- c(1:n)
  V.temp <- V.all[V.temp.ID,]
  d <- dim(V.all)[2]
  Polytopes <- list()
  dist.all <- as.matrix(dist(V.all))
  
  countbygroup <- sum(!is.na(cts_z)) #number of non-NA counts in the polytope
  
  dist.max <- max(dist.all)
  
  if(mode.type==1){
    aj <- 1+a_h*length(id.train)
    bj <- mu_h+sum(a_h*(cts_z[id.train]))
    Polytopes[[1]] <- list(V.ID=V.temp.ID,  Lambda=dist.max/2, CountByGroup=countbygroup,dist.max=dist.max,aj=aj,bj=bj)
    
    alpha_t <- alpha_h+length(id.train)/2
    beta_t <- beta_h + mu_h^2/2/a_h +  0.5*sum(cts_z[id.train]^2)-0.5/a_h*bj^2/aj
    logaj <- log(1+a_h*length(id.train))
    logPt <- -(alpha_t)*log(beta_t)-0.5*logaj
    
    #Partition.Inital <- list(Polytopes=Polytopes,tau=rexp(1,rate=Polytopes[[1]]$Lambda),beta_t = beta_t,logPt=logPt)
  }else{
    V.all.nonNA <- V.all[!is.na(cts_z),]
    Aj.M.all <- diag(1/a_h,nrow = d,ncol=d)+t(V.all.nonNA)%*%(V.all.nonNA)
    cts_z.nonNA <- cts_z[!is.na(cts_z)]
    Bj.M.all <-  as.matrix(apply((V.all.nonNA*matrix(rep(cts_z.nonNA,d),ncol=d,byrow = FALSE)),2,sum))
    Aj.M.all.eigen <- eigen(Aj.M.all,symmetric=TRUE)
    Aj.inv.all <- Aj.M.all.eigen$vectors%*%diag(1/Aj.M.all.eigen$values)%*%t(Aj.M.all.eigen$vectors)
    
    aj  <-  prod(Aj.M.all.eigen$values) # det(Aj.M.left)
    bj <-  t(Bj.M.all)%*%Aj.inv.all%*%Bj.M.all
    
    Polytopes[[1]] <- list(V.ID=V.temp.ID,  Lambda=dist.max/2, CountByGroup=countbygroup,dist.max=dist.max,aj=aj,bj=bj)
    
    beta_t <- beta_h+ sum(cts_z[id.train]^2)/2-bj/2
    logPt  <- -3/2*log(a_h)-(alpha_h+length(id.train)/2)*log(beta_t)-log(aj)/2  
    #Partition.Inital <- list(Polytopes=Polytopes,tau=rexp(1,rate=Polytopes[[1]]$Lambda),beta_t = beta_t,logPt=logPt)
  }
  
  Partition.Inital <- list(Polytopes=Polytopes,tau=NULL,beta_t = beta_t,logPt=logPt)
  output <- list(Partition.Inital=Partition.Inital,dist.all=dist.all)
  return(output)
}


getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}

Cut_plane_standard <- function(V,V.ID,cts_z,sample.type,cut.type,accept.rate.min,dist.max=NA,w.normal=NULL){
  d <- dim(V)[2]
  P.u <- rep(NA,d)
  theta <- rep(NA,(d-1))
  mu <- c(NA)
  len.projection <- c(NA)
  normal.v <- rep(NA,d)
  sample.cut.count <- 1
  V.left.ID <-NULL
  V.right.ID <- NULL
  if(is.null(w.normal)){
    w.normal <- rep(1/d,d)
  }
  
  if(sum(!is.na(cts_z))>0){
    V.unique <- unique(V[!is.na(cts_z),])
    d1 <- dim(V.unique)[1]
    if(is.na(dist.max)){
      if(d1==1){
        dist.max <- 0
      }else{
        dist.max <- max(dist(V.unique))
      }
    }
    #skip the cut if  vertices are same
    skip.index <-  dist.max==0
  }else{
    skip.index =1 # skip if the responses of all vertices are NULL
  }
  
  
  if(!skip.index){
    if(sample.type==1){
      # use the largest distance between vertices as upper bound instead of the diameter of the
      # smallest d-sphere, since
      # length(projection.line.segment(Vertices))<= max{dist between Vertices}
      #                                         <= the diameter of the smallest sphere
      # dist.max would be the smallest boundary of the rejection sampling
      run.index <- 1
      if(cut.type==1){
        while(run.index){
          
          
          normal.v <- rnorm(d,sd=0.1)
          normal.v <- w.normal*normal.v
          normal.v <- normal.v/sqrt(sum(normal.v^2))
          
          t.scale <-  V%*%matrix(normal.v,nrow=d,ncol=1)/sum(normal.v^2)
          t.scale.train <- t.scale[!is.na(cts_z)]
          t.scale.ends <- range(t.scale.train,na.rm = TRUE)
          V.projection.ends <-rbind(t.scale.ends[1]*normal.v,t.scale.ends[2]*normal.v)
          len.projection <- diff(t.scale.ends)*sqrt(sum(normal.v^2))
          
          mu <- runif(1,min=0,max=dist.max)
          mu <- mu/len.projection
          if(mu<1){
            t.scale.mu <- (1-mu)*t.scale.ends[1]+mu*t.scale.ends[2]
            #if(any(t.scale<t.scale.mu) && any(t.scale>t.scale.mu)){
            #if(any(t.scale.train<=t.scale.mu) && any(t.scale.train>t.scale.mu)){
            if(!(sum(t.scale.train<=t.scale.mu)<d) & !(sum(t.scale.train>t.scale.mu)<d)){
              run.index <- 0
              P.u <- (1-mu)*V.projection.ends[1,]+mu*V.projection.ends[2,]
              V.left.ID <- which(t.scale<=t.scale.mu)
              V.right.ID <- which(t.scale>t.scale.mu)
            }
          }else{
            sample.cut.count <- sample.cut.count+1
            run.index <- (1/sample.cut.count>accept.rate.min)*1
            skip.index <- 1-run.index
          }
        }
      }else{
        while(run.index){
          
          normal.v<- c(rmultinom(1,size=1,prob=w.normal))
          
          
          t.scale <-  V%*%matrix(normal.v,nrow=d,ncol=1)/sum(normal.v^2)
          t.scale.train <- t.scale[!is.na(cts_z)]
          t.scale.ends <- range(t.scale.train,na.rm = TRUE)
          V.projection.ends <-rbind(t.scale.ends[1]*normal.v,t.scale.ends[2]*normal.v)
          len.projection <- diff(t.scale.ends)*sqrt(sum(normal.v^2))
          
          
          mu <- runif(1,min=0,max=dist.max)
          mu <- mu/len.projection
          
          if(mu<1){
            t.scale.mu<-(1-mu)*t.scale.ends[1]+mu*t.scale.ends[2]
            if(!(sum(t.scale.train<=t.scale.mu)<d) & !(sum(t.scale.train>t.scale.mu)<d)){
              #if(any(t.scale.train<=t.scale.mu) && any(t.scale.train>t.scale.mu)){
              run.index <- 0
              P.u <- (1-mu)*V.projection.ends[1,]+mu*V.projection.ends[2,]
              V.left.ID <- which(t.scale<=t.scale.mu)
              V.right.ID <- which(t.scale>t.scale.mu)
            }
            
          }else{
            sample.cut.count <- sample.cut.count+1
            run.index <- (1/sample.cut.count>accept.rate.min)*1
            skip.index <- 1-run.index
          }
        }
      }
    }else{
      run.index <- 1
      while(run.index){
        
        if(cut.type==1){
          normal.v <- rnorm(d,sd=0.1)
          normal.v <- w.normal*normal.v
          normal.v <- normal.v/sqrt(sum(normal.v^2))
          
          
        }else{
          normal.v<- c(rmultinom(1,size=1,prob=w.normal))
          
        }
        
        t.scale <-  V%*%matrix(normal.v,nrow=d,ncol=1)/sum(normal.v^2)
        t.scale.train <- t.scale[!is.na(cts_z)]
        t.scale.ends <- range(t.scale.train,na.rm = TRUE)
        V.projection.ends <-rbind(t.scale.ends[1]*normal.v,t.scale.ends[2]*normal.v)
        len.projection <- diff(t.scale.ends)*sqrt(sum(normal.v^2))
        
        mu <- runif(1,0,1)
        t.scale.mu<-(1-mu)*t.scale.ends[1]+mu*t.scale.ends[2]
        
        # if(any(t.scale.train<=t.scale.mu) && any(t.scale.train>t.scale.mu)){
        if(!(sum(t.scale.train<=t.scale.mu)<d) & !(sum(t.scale.train>t.scale.mu)<d)){
          run.index <- 0
          P.u <- (1-mu)*V.projection.ends[1,]+mu*V.projection.ends[2,]
          V.left.ID <- which(t.scale<=t.scale.mu)
          V.right.ID <- which(t.scale>t.scale.mu)
        }else{
          sample.cut.count <- sample.cut.count+1
          run.index <- (sample.cut.count<11)*1
          skip.index <- 1-run.index
        }
        
      }
      
    }
  }
  return(list(N_vec= normal.v,Pu=P.u,
              Len.projection=len.projection,dist.max=dist.max,
              sample.cut.count=sample.cut.count,
              skip.index=skip.index,
              V.left.ID=V.left.ID,
              V.right.ID=V.right.ID
  ))
}


Generative_Process <- function(partition,V.all,cts_z,tau,group.level,group.len,sample.type,cut.type,accept.rate.min,w.normal,mode.type){
  l <- length(partition$Polytopes)
  d <- dim(V.all)[2]
  tau.v <- partition$tau
  result <- c(0,0,0)
  Cut <- list()
  cut.idx <- 0
  skip.all <- 0  #only set it to 1 when  (1)  cond1 is not true OR (2) cost exceeds the budget.
  new_par_t <- list(beta_t=NULL,logPt_delta=NULL)
  
  
  
  if(tau.v[l]>=tau){
    skip.all <- 1
  }else{
    Lambdas <- sapply(partition$Polytopes, "[[", 2)
    
    ##shrink the candidate space, only choose from polytopes whose
    #  (i) Lambda>0
    # Or (ii) #polytopes < 2p  (make sure each terminal node contains at least p polytopes 
    # to ensure the model is identifiable.)
    cond1 <- (Lambdas>0)
    cond2 <-  sapply(partition$Polytopes, "[[", 3)
    cond2 <-  !(cond2 < 2*d)
    cond12 <- sum(cond1&cond2)
    if(sum(cond12)==0){
      skip.all <- 1
    }else{
      sample.space <- which(cond1&cond2)
      if(length(sample.space)==1){
        j <- sample.space
      }else{
        j <- sample(x=c(sample.space),size=1,prob = Lambdas[sample.space])
      }
      
      V.temp.ID <- partition$Polytopes[[j]]$V.ID
      V.temp <- V.all[V.temp.ID,]
      cts_z.temp <- cts_z[V.temp.ID]
      dist.max <- partition$Polytopes[[j]]$dist.max
      
      
      Cut <- Cut_plane_standard(V.temp,V.tempID,cts_z.temp,sample.type,cut.type,accept.rate.min,dist.max,w.normal)
      
      if(Cut$skip.index!=1){
        index.left <- Cut$V.left.ID
        index.right <- Cut$V.right.ID
        
        if(length(index.left)==1){
          V.temp.left <- t(as.matrix(V.temp[index.left,]))
        }else{
          V.temp.left <- V.temp[index.left,]
        }
        
        if(length(index.right)==1){
          V.temp.right <- t(as.matrix(V.temp[index.right,]))
        }else{
          V.temp.right <- V.temp[index.right,]
        }
        
        V.temp.ID.left <- V.temp.ID[index.left]
        V.temp.ID.right <- V.temp.ID[index.right]
        
        cut.idx <- 1
        l <- l+1
        
        
        aj_tminus1 <- partition$Polytopes[[j]]$aj
        bj_tminus1 <- partition$Polytopes[[j]]$bj
        
        dist.max.left <- max(dist.all[V.temp.ID.left,V.temp.ID.left])
        dist.max.right <-  max(dist.all[V.temp.ID.right,V.temp.ID.right])
        Lambda.left <- dist.max.left/2
        Lambda.right <- dist.max.right/2
        
        
        
        
        cts_z.temp.left <- cts_z[V.temp.ID.left]
        cts_z.temp.right <- cts_z[V.temp.ID.right]
        
        countbygroup.left <- sum(!is.na(cts_z.temp.left))
        countbygroup.right <- sum(!is.na(cts_z.temp.right))
        
        if(mode.type==1){
          aj.left <- 1+ a_h*countbygroup.left
          aj.right <- 1+a_h*countbygroup.right
          bj.left <-  mu_h+sum(a_h*cts_z.temp.left,na.rm = TRUE)
          bj.right <- mu_h+sum(a_h*cts_z.temp.right,na.rm= TRUE)
          
          beta_t <- partition$beta_t+mu_h^2/2/a_h + (bj_tminus1^2/aj_tminus1-bj.left^2/aj.left-bj.right^2/aj.right)/2/a_h
          logPt_delta <- alpha_t*(log(partition$beta_t)-log(beta_t))+(log(aj_tminus1)-log(aj.left)-log(aj.right))/2 
        }else{
          V.left.nonNA <- (V.all[V.temp.ID.left,])[!is.na(cts_z.temp.left),]
          V.right.nonNA <- (V.all[V.temp.ID.right,])[!is.na(cts_z.temp.right),]
          cts_z.temp.left.nonNA <- cts_z.temp.left[!is.na(cts_z.temp.left)]
          cts_z.temp.right.nonNA <- cts_z.temp.right[!is.na(cts_z.temp.right)]
          
          Aj.M.left <- diag(1/a_h,nrow = d,ncol=d)+t(V.left.nonNA)%*%(V.left.nonNA)
          Aj.M.right <- diag(1/a_h,nrow = d,ncol=d)+t(V.right.nonNA)%*%(V.right.nonNA)
          
          Bj.M.left <-  as.matrix(apply((V.left.nonNA*matrix(rep(cts_z.temp.left.nonNA,d),ncol=d,byrow = FALSE)),2,sum))
          Bj.M.right <-  as.matrix(apply((V.right.nonNA*matrix(rep(cts_z.temp.right.nonNA,d),ncol=d,byrow = FALSE)),2,sum))
          
          Aj.M.l.eigen <- eigen(Aj.M.left,symmetric=TRUE)
          Aj.M.r.eigen <- eigen(Aj.M.right,symmetric=TRUE)
          
          # Aj.l <- Aj.M.l.eigen$vectors%*%diag(Aj.M.l.eigen$values)%*%t(Aj.M.l.eigen$vectors)
          Aj.inv.l <- Aj.M.l.eigen$vectors%*%diag(1/Aj.M.l.eigen$values)%*%t(Aj.M.l.eigen$vectors)
          Aj.inv.r <- Aj.M.r.eigen$vectors%*%diag(1/Aj.M.r.eigen$values)%*%t(Aj.M.r.eigen$vectors)
 
          aj.left <-  prod(Aj.M.l.eigen$values) # det(Aj.M.left)
          aj.right <- prod(Aj.M.r.eigen$values) 
          bj.left <-  t(Bj.M.left)%*%Aj.inv.l%*%Bj.M.left
          bj.right <- t(Bj.M.right)%*%Aj.inv.r%*%Bj.M.right
          
          beta_t <- partition$beta_t+(bj_tminus1-bj.left-bj.right)/2
          logPt_delta <-log(a_h)/2+ alpha_t*(log(partition$beta_t)-log(beta_t))+(log(aj_tminus1)-log(aj.left)-log(aj.right))/2 #logP_{t+1} - logP_{t}
        }
        
       
        
        
        
        partition$Polytopes[[j]] <- list(V.ID=V.temp.ID.left,Lambda=Lambda.left,CountByGroup=countbygroup.left,
                                         dist.max=dist.max.left,aj=aj.left,bj=bj.left)
        partition$Polytopes[[l]] <- list(V.ID=V.temp.ID.right,Lambda=Lambda.right,CountByGroup=countbygroup.right,
                                         dist.max=dist.max.right,aj=aj.right,bj=bj.right)
        
        new_par_t <- list(beta_t=beta_t,logPt_delta=logPt_delta)
        
        Lambdas <- sapply(partition$Polytopes, "[[", 2)
        tau.v.plus1 <- rexp(n=1,sum(Lambdas))+tau.v[l-1]
        tau.v <- c(tau.v,tau.v.plus1)
        partition$tau <- tau.v
        partition$beta_t <- beta_t
        partition$logPt <- partition$logPt +  logPt_delta
        
        
      }
    }
    
  }
  return(list(partition=partition,
              cut.indx=cut.idx,
              Cut=Cut,skip.all=skip.all,new_par_t=new_par_t))
}

