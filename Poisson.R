rm(list = ls())
# source("DataGenerator.R")
Rcpp::sourceCpp("data_generator.cpp")
source("Algorithm.R")
library(glmnet)
library(MASS)
library(ncvreg)
library(LaplacesDemon)

set.seed(2023)

family = "poisson"
n <- 200
p <- 3000
#k <- floor(2*log(n)*n^(1/3)/3) 
k <- floor(n/(3*log(n)))
m <- 3000
m0 <- m-1000
S <- 1

#ture value of parameters
rho <- -0.15

beta <- M <- list()
length(beta) <- length(M) <- 3

M[[1]] <- sort(sample(1:p,8))
U <- rbinom(8,1,0.8)
U[which(U<1)] <- -1
beta[[1]] <- rep(0,p)
beta[[1]][M[[1]]] <- U*(log(n)/sqrt(n)+abs(rnorm(8))/8)

M[[2]] <- c(1,3,5,7,9)
beta[[2]] <- rep(0,p)
beta[[2]][1] <- 0.6
beta[[2]][3] <- -0.6
beta[[2]][5] <- 0.5
beta[[2]][7] <- -0.5
beta[[2]][9] <- 0.5

M[[3]] <- c(1:4)
beta[[3]] <- rep(0,p)
beta[[3]][1:4] <- 0.675

beta_cpp <- cbind(beta[[1]],beta[[2]],beta[[3]])

theta <- list()
length(theta) <- 3
theta[[1]] <- c(beta[[1]],rho)
theta[[2]] <- c(beta[[2]],rho)
theta[[3]] <- c(beta[[3]],rho)

#outcome of programme
names1 <- list()
length(names1) <- 4
names1[[1]] <- paste("Setting_",c(1:3))
names1[[2]] <- c("NW-SIS","SMLE","PMLE_NPG","PMLE_AMPG","PMLE_ANPG")
names1[[3]] <- c("RC","PSR","FDR","AMS","iter","TIME")
names1[[4]] <- paste("Replication_",c(1:S))
result1 <- array(0,dim=c(3,5,6,S),dimnames=names1)

names2 <- list()
length(names2) <- 4
names2[[1]] <- paste("Setting_",c(1:3))
names2[[3]] <- c("NW-SIS","SMLE","PMLE_NPG","PMLE_AMPG","PMLE_ANPG")
names2[[3]] <- c("MSE_beta","MSE_beta_post","MSE_rho","MSE_rho_post","Perr","Perr_post")
names2[[4]] <- paste("Replication_",c(1:S))
result2 <- array(0,dim=c(3,5,6,S),dimnames=names2)

#main loop
for(i in 1:S){
  time_start <- proc.time()
  
  #Data generation
  # W <- networkmat(n)
  # Data <- data_generation(n,p,m,m0,W,rho,beta,family)
  W <- Network_Cpp(n)
  Data <- dt_gen_rcpp(n,p,m,m0,W,rho,beta_cpp,family)
  
  for(j in 1:3){
    
    #i=1; j=1
    
    X <- Data$X[[j]]
    Y <- Data$Y[[j]]
    
    ##NW-SIS####################################################################
    ptm1 <- proc.time()
    fit11_11 <- NW_SIS(X,Y,W,k)
    ptm2 <- proc.time()
    result1[j,1,,i] <- c(round(measure(fit11_11$ix,M[[j]]),3), NA, round(ptm2[1]-ptm1[1],3))
    #fit11_12 <- Post_estimation(X, Y, W, p, fit11_11$ix, family)
    fit11_12 <- PMLE_POST(X, Y, W, p, fit11_11$ix,family,tol=1e-4, maxiter=200)
    result2[j,1,,i] <- c(NA,sqrt(sum((fit11_12$theta[1:p]-beta[[j]])^2)),NA,sqrt(sum((fit11_12$theta[p+1]-rho)^2)),
                         NA,fit11_12$res)
    
    ##SMLE####################################################################
    ptm1 <- proc.time()
    fit <- glmnet(X,Y,family,intercept = F)
    beta_initial <- fit$beta[,length(fit$lambda)]
    fit12_11 <- SMLE(X,Y,family,beta_initial, k,tol=1e-4,maxiter=200)
    ptm2 <- proc.time()
    result1[j,2,,i] <- c(round(measure(fit12_11$ix,M[[j]]),3),fit12_11$iter, round(ptm2[1]-ptm1[1],3))
    fit12_12 <- SMLE_POST(X,Y,p,fit12_11$ix,family,tol=1e-4, maxiter=200)
    result2[j,2,,i] <- c(sqrt(sum((fit12_11$beta-beta[[j]])^2)),sqrt(sum((fit12_12$beta-beta[[j]])^2)),NA,NA,
                         fit12_11$res,fit12_12$res)
    
    ##PMLE_NPG####################################################################
    theta_initial <- theta_ini(X,Y,W,family)
    ptm1 <- proc.time()
    fit13_11 <- PMLE_NPG(X, Y, W, family, theta_initial, k, tol=1e-4, maxiter=200)
    ptm2 <- proc.time()
    result1[j,3,,i] <- c(round(measure(fit13_11$ix,M[[j]]),3), fit13_11$iter, round(ptm2[1]-ptm1[1],3))
    #fit13_12 <- Post_estimation(X, Y, W, p, fit13_11$ix, family)
    fit13_12 <- PMLE_POST(X, Y, W, p, fit13_11$ix,family,tol=1e-4, maxiter=200)
    result2[j,3,,i] <- c(sqrt(sum((fit13_11$theta[1:p]-beta[[j]])^2)),sqrt(sum((fit13_12$theta[1:p]-beta[[j]])^2)),
                         sqrt(sum((fit13_11$theta[p+1]-rho)^2)),sqrt(sum((fit13_12$theta[p+1]-rho)^2)),fit13_11$res,fit13_12$res)
    
    ##PMLE_MAPG####################################################################
    ptm1 <- proc.time()
    fit14_11 <- PMLE_MAPG(X, Y, W, family, theta_initial, k, tol=1e-4, maxiter=200)
    ptm2 <- proc.time()
    result1[j,4,,i] <- c(round(measure(fit14_11$ix,M[[j]]),3), fit14_11$iter, round(ptm2[1]-ptm1[1],3))
    #fit14_12 <- Post_estimation(X, Y, W, p, fit14_11$ix, family)
    fit14_12 <- PMLE_POST(X, Y, W, p, fit14_11$ix,family,tol=1e-4, maxiter=200)
    result2[j,4,,i] <- c(sqrt(sum((fit14_11$theta[1:p]-beta[[j]])^2)),sqrt(sum((fit14_12$theta[1:p]-beta[[j]])^2)),
                         sqrt(sum((fit14_11$theta[p+1]-rho)^2)),sqrt(sum((fit14_12$theta[p+1]-rho)^2)),
                         fit14_11$res,fit14_12$res)
    
    ##PMLE_NAPG####################################################################
    ptm1 <- proc.time()
    fit15_11 <- PMLE_NAPG(X, Y, W, family, theta_initial, k, tol=1e-4, maxiter=200)
    ptm2 <- proc.time()
    result1[j,5,,i] <- c(round(measure(fit15_11$ix,M[[j]]),3), fit15_11$iter, round(ptm2[1]-ptm1[1],3))
    #fit15_12 <- Post_estimation(X, Y, W, p, fit15_11$ix, family)
    fit15_12 <- PMLE_POST(X, Y, W, p, fit15_11$ix,family,tol=1e-4, maxiter=200)
    result2[j,5,,i] <- c(sqrt(sum((fit15_11$theta[1:p]-beta[[j]])^2)),sqrt(sum((fit15_12$theta[1:p]-beta[[j]])^2)),
                         sqrt(sum((fit15_11$theta[p+1]-rho)^2)),sqrt(sum((fit15_12$theta[p+1]-rho)^2)),
                         fit15_11$res,fit15_12$res)
    
    print(paste("this is the",i,"replication",j,"setting"))
  }
  time_stop <- proc.time()
  timer <- time_stop[1]-time_start[1]
  print(timer)
}


res1 <- apply(result1,1:3,mean)
table1 <- round(rbind(res1[1,,],res1[2,,],res1[3,,]),3)
write.table(table1,paste("Poisson-","Screening-",n,"-",p,".txt",sep=""))
write.csv(table1,paste("Poisson-","Screening-",n,"-",p,".csv",sep=""))

res2 <- apply(result2,1:3,mean)
table2 <- round(rbind(res2[1,,],res2[2,,],res2[3,,]),3)
write.table(table2,paste("Poisson-","Estimation-",n,"-",p,".txt",sep=""))
write.csv(table2,paste("Poisson-","Estimation-",n,"-",p,".csv",sep=""))

table1
table2
