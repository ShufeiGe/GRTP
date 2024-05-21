#if (opt$Mondrian) {
#  cut.type <- 2
#} else {
#  cut.type <- 1
#}

ALPHA <- as.numeric(opt$alpha)
if (ALPHA < 0) {
  print_error('Argument to --alpha must be positive')
}

split.seed <- as.numeric(opt$seed)
if (split.seed == -1) {
  split.seed <- as.integer((as.double(Sys.time())*1000+Sys.getpid()) %% 2^31)
  set.seed(split.seed)
}else{
  set.seed(split.seed)
}
tau <- as.numeric(opt$tau)
sample.type <- 1
accept.rate.min <- 0.05
N <- as.numeric(opt$particles)
rep.max <- as.numeric(opt$ntrees)

if (opt$weight) {
  if (!file.exists(opt$weight)) {
    print_error('Could not open weight file')
  }
  w.normal <- scan(opt$weight, numeric(), quote = "")
  w.normal <- w.normal/sum(w.normal)
} else {
  w.normal <- NULL
}


suppressMessages(library(purrr)    )




data <- read.table(fname, sep = ' ', row.names = NULL, header = TRUE)
d     <- dim(data)[2]-1
V.all <-  as.matrix(data[,1:d])


if(l.max==Inf){
  l.max <- dim(data)[1]
}

#data-preprocess : standardization
min <- matrix(apply(V.all,2,min),nrow=dim(V.all)[1],ncol=dim(V.all)[2],byrow=TRUE)
max <- matrix(apply(V.all,2,max),nrow=dim(V.all)[1],ncol=dim(V.all)[2],byrow=TRUE)
V.all <- (V.all-min)/(max-min)
rm("min","max")


cts_z <- data[,(d+1)]
id.test <- which(is.na(cts_z))

if(length(id.test)==0){
  id.train <- c(1:length(cts_z))
}else{
  id.train <- c(1:length(cts_z))[-id.test]
}

source(paste0(pwd.code,"/generative_functions.R"))
if(mode.type==1){
  source(paste0(pwd.code,"/parameter_inference_c.R"))
}else{
  source(paste0(pwd.code,"/parameter_inference_l.R"))
}
#----inference -----start from here--------



#------------------Partition versus no. cuts (No. of cuts =l.max or tau=Inf)------------------
t0 <- proc.time()
table.out <- c()

cat(sprintf("\nEstimated runtime: ?\n"))
pb <- txtProgressBar(min = 0, max = N*rep.max, style = 3)

set.seed(split.seed)

part_init <- partition_initial(V.all,cts_z,l.max,sample.type,cut.type,accept.rate.min,mode.type)

Partition.Inital <- part_init$Partition.Inital
dist.all <- part_init$dist.all
rm("part_init")

Partition.t0 <- list()
Lambda0 <- Partition.Inital$Polytopes[[1]]$Lambda



alpha_t <- alpha_h+length(id.train)/2
Wt0 <- rep(1/N,N)   # Weights vector for N particles at time t=0
logWt0 <- log(Wt0)   # log Weights vector for N particles at time t=0
logl.tminus0 <- Partition.Inital$logPt

for(i in 1:N){
  tau.temp <- rexp(1,rate=Lambda0)
  Partition.t0[[i]] <- Partition.Inital
  Partition.t0[[i]]$tau <- tau.temp
}

cts_z.predict.all.mu <- vector("list", rep.max)#NULL
cts_z.predict.all.sig <- vector("list", rep.max)#NULL
cts_z.predict.all.sig.sd <- vector("list", rep.max)#NULL
zzt.tau.tauli.all <- vector("list", rep.max)#NULL

len.taul.v <- length(tau.pred.t.v)
t1 <- proc.time()-t0

for(rep in 1:rep.max){
  # create progress bar
  setTxtProgressBar(pb, ((rep-1)*N))
  
  cts_z.predict.all.mu[[rep]] <-  matrix(NA,nrow=(dim(data)[1]),ncol=len.taul.v)
  cts_z.predict.all.sig[[rep]] <-  rep(NA,len.taul.v)
  cts_z.predict.all.sig.sd[[rep]] <-   rep(NA,len.taul.v)
  zzt.tau.tauli.all[[rep]] <-  matrix(NA,nrow=3,ncol=len.taul.v)
  zzt.tau.tauli.all.temp <- matrix(NA,nrow=3,ncol=len.taul.v)
  
  Partition.t <- Partition.t0
  Wt <- Wt0        # Weights vector for N particles at time t
  logWt <- logWt0  # log Weights vector for N particles at time t
  logl.tminus1 <- logl.tminus0 #log likelihood for N particles at time t-1
  logl.t <- logl.tminus1 #log likelihood for N particles at time t
  t2 <- proc.time() 
  for(tau.temp.id in 1:len.taul.v){
    
    tau.li <- min(map_dbl(Partition.t,~rev(.$tau)[1]))
    end.cond <- FALSE
    
    skip.index.N <- rep(0,N)
    tau <- tau.pred.t.v[tau.temp.id]
    t <- 1  # here t is the number of cuts; we use vector tau to denote the cost at each cut.
    while(tau.li<tau){
      Cut <- vector(mode="list",length=N)
      #resampling-------------start
      if(t>1){
        jC <- sample(x=1:N,size=N,prob =Wt,replace = TRUE )
        Partition.t.new <- list()
        for(i in 1:N){
          Partition.t.new[[i]] <- Partition.t[[jC[i]]]
        }
        Wt[1:N] <- 1/N
        
        logWt <- log(Wt)
        logl.t <- logl.t[jC]
        logl.tminus1 <- logl.tminus1[jC]
        skip.index.N <- skip.index.N[jC]
        Partition.t <- Partition.t.new
        rm("Partition.t.new")
      }
      #resampling-------------end
 
      for (i in 1:N){
        if(!skip.index.N[i]){
          partition.temp <- Partition.t[[i]]
          Out.temp <- Generative_Process(partition.temp,V.all,cts_z,tau,group.level,group.len,sample.type,cut.type,accept.rate.min,w.normal,mode.type)
          Cut[[i]] <- Out.temp$Cut
          if(Out.temp$cut.indx){
            logl.t[i] <- Out.temp$partition$logPt 
            logWt[i] <- logWt[i]+Out.temp$new_par_t$logPt_delta 
            logl.tminus1[i] <- logl.t[i]
            Partition.t[[i]] <- Out.temp$partition
          }else{
            skip.index.N[i] <- Out.temp$skip.all
          }
        }
      }
      
      
      Wt <- exp((logWt-max(logWt)))/sum(exp((logWt-max(logWt))))
      
      t <- max(length(Partition.t[[which.max(Wt)]]$tau)-1,1)
      
      tau.li <- min(map_dbl(Partition.t,~rev(.$tau)[1]))
      end.cond <- min(skip.index.N)
      tau.temp.v <- c(tau=tau,tau.li=tau.li)
      # to avoid endless while-loop
      #if(end.cond | (t+1)>l.max | tau.li>tau){
      if(end.cond | (t+1)>l.max){
        tau.li <- Inf
      }
      
      pred.out.idx <- (tau.li>tau)

      # prediction --- from here 
      # if(tau.li==Inf){
      #if(TRUE){
        if(pred.out.idx){
        partition.best <- Partition.t[[which.max(Wt)]]
        test.inf <- mcmc_infer(partition.best,V.all,cts_z,n.iter,burnin.iter,save.all.par=FALSE)
        names(test.inf)
        zz <- round((proc.time()-t2+t1)[3],4)
        zzt.tau.tauli.all.temp[,tau.temp.id]<- c(zz,tau.temp.v)
        
        cts_z.predict.all.sig[[rep]][tau.temp.id] <- test.inf$mean.sig 
        cts_z.predict.all.sig.sd[[rep]][tau.temp.id] <-test.inf$sd.sig
        pred.t <-test.inf$pred.z
        colnames(pred.t) <- t
        cts_z.predict.all.mu[[rep]][,tau.temp.id] <- pred.t
        rm("test.inf")
      }
      

      
    } #end while{}
  
    

  }#end for (tau)
  
  
  zzt.tau.tauli.all[[rep]] <- zzt.tau.tauli.all.temp
  

  if(rep==1){
    zz = round((proc.time()-t0)[3]/60/60*rep.max,4)
    if (zz < 1) {
      zz = round(zz * 60)
      if (zz <= 1) {
        zz = "1 minute"
      } else {
        zz = sprintf("%d minutes", zz)
      }
    } else if (zz > 24) {
      zz = round(zz/24)
      if (zz <= 1) {
        zz = "1 day"
      } else {
        zz = sprintf("%d days", zz)
      }
    } else {
      zz = round(zz)
      if (zz <= 1) {
        zz = "1 hour"
      } else {
        zz = sprintf("%d hours", zz)
      }
    }
    
    #cat(sprintf("\r\033[K\033[1A\r\033[K\r%s\n", paste("Estimated runtime: ", zz, sep="" )))
    cat( paste("\n Estimated runtime: ", zz,"\n",sep="" ))
  }
  
  ls1 <- ls()
  ls2 <- ls1[ls1!="dist.all"]
  (out.f <- paste0("./result/Model",mode.type,"cut",cut.type,"D",data.fd,data.id,".Rdata"))
  save(list=ls2,file=out.f) 
}

rm("dist.all")



