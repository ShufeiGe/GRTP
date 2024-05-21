rm(list=ls())
require(stringr)
data.fd <- "d01"  #data.folder 
data.id <- "01"   #data id
mode.type <- 1 #mode type - 1 for constant case; other values for the linear case.
cut.type <-  1 #cut type - 2 for Mondrain; 1 for uniform
opt <- list()
if (cut.type ==2 ) {
  opt$Mondrian<- TRUE
} else {
  opt$Mondrian<- FALSE
}


#setwd("your work path")
pwd.data <- paste0(getwd(),"/data/")
pwd.code <- paste0(getwd(),"/code")  

#------------ input data file
(fname = paste0("./data/",data.fd,"/X/",data.id,".txt"))



opt$alpha  <-  1/1000
opt$seed  <- 123
opt$ntrees <- 1
opt$weight <- FALSE 
opt$particles <- 100
l.max <- Inf #Inf  # no. of cuts 
tau.pred.t.v <- 10 # pre-specified budget (e.g. 5,10,20,25; any positive values); 
                   # could be a scalar or a vector if you want to compare the performance of the model with different budgets

#hyper-parameters
a_h <- 1
mu_h <- 0
alpha_h <- 0.5
beta_h <- 0.5

#no. of iterations & no. of burn-in in the MCMC alg.
n.iter <- 5000
burnin.iter <- 1000


source(paste0(pwd.code,"/tess22_frame_tau.R"))
 
(out.f2 <- paste0("./result/Model",mode.type,"cut",cut.type,"D",data.fd,data.id,"_summary.Rdata"))
(ggname = paste0(pwd.data,data.fd,"/Z/",data.id,".txt"))
fx.true <-  read.table(ggname, sep = ' ', row.names = NULL, header = TRUE)[,1]

(yname = paste0(pwd.data,data.fd,"/Y/",data.id,".txt"))
z.obs <-  read.table(yname, sep = ' ', row.names = NULL, header = TRUE)[,1]
z.obs[id.train] <- cts_z[id.train]

 
qr.test <- NULL
mae.test <- NULL
qr.train <- NULL
mae.train <- NULL
qr2.test <- NULL
mae2.test <- NULL
qr2.train <- NULL
mae2.train <- NULL

for(jj in 1:len.taul.v){
  temp.mu <- sapply(cts_z.predict.all.mu,function(x){x[,jj]})
  cts_z.predict.temp <- as.matrix(t(na.omit(t(temp.mu))))
  cts_z.predict.avg <- apply(cts_z.predict.temp,1,mean)
  err <- c(cts_z.predict.avg-fx.true)
  err2 <- c(cts_z.predict.avg-z.obs)
  qr.test <- c(qr.test,mean((err[id.test])^2))
  mae.test <- c(mae.test,mean(abs(err[id.test]))) #Mean absolute error (MAE) 
  qr.train <- c(qr.train,mean((err[id.train])^2))
  mae.train <- c(mae.train,mean(abs(err[id.train])))#Mean absolute error (MAE) 
  
  qr2.test <- c(qr2.test,mean((err2[id.test])^2))
  mae2.test <- c(mae2.test,mean(abs(err2[id.test]))) #Mean absolute error (MAE) 
  qr2.train <- c(qr2.train,mean((err2[id.train])^2))
  mae2.train <- c(mae2.train,mean(abs(err2[id.train])))#Mean absolute error (MAE) 
  
}

(metrics = rbind(c(qr.train,qr.test), c(mae.train,mae.test),
                c(qr2.train,qr2.test),
                c(mae2.train,mae2.test)))

base::save(opt,
           cut.type,
           mode.type,
           zzt.tau.tauli.all.temp, #each col. - (running time; pre-specified budget; real budget if implement +1 split)
           metrics, 
     #cts_z.predict.all.mu,
     #cts_z.predict.all.sig,
     #zzt.tau.tauli.all,z.obs,
     #cts_z.predict.all.sig.sd,
     #id.test,id.train, #test id; train id
     #fx.true,cts_z,
     #partition.best, #best partition 
     file=out.f2)
