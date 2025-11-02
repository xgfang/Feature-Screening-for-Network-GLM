Rcpp::sourceCpp('mvrnorm_rcpp.cpp')

### Create a network structure 
networkmat <- function(n)
{
  ass_num <- sample(1:50,n,replace=TRUE)
  A <- matrix(rbinom(n*n,1,0.6),n,n)
  diag(A) <- 0
  A <- outer(ass_num,ass_num,"==")*A
  W <- A/rowSums(A)
  W[is.na(W)]=0
  return(W)
}

### function to get mode 
getmode <- function(v) 
{
  uniqv <- unique(v)
  mode <- uniqv[which.max(tabulate(match(v, uniqv)))]
  return(mode)
}

### function to generate data which we required
data_generation <- function(n,p,m,m0,W,rho,beta,family)
{
  ### S1
  X1 <- matrix(rnorm(n*p),n,p)
  
  ### S2
  sigma2 <- matrix(0.5,p,p)
  diag(sigma2) <- 1
  X2 <- mvrnormArma(n,rep(0,p),sigma2)
  
  ### S3
  sigma3 <- matrix(0.3,p,p)
  sigma3[(row(sigma3)<=4)&(col(sigma3)<=4)]<- 0.15
  sigma3[row(sigma3)==col(sigma3)]<- 1
  X3 <- mvrnormArma(n, rep(0,p),sigma3)
  
  Y_t1 <- rep(0,n)
  Y_t2 <- rep(0,n)
  Y_t3 <- rep(0,n)
  theta1 <- c()
  theta2 <- c()
  theta3 <- c()
  Y1 <- matrix(0,m,n)
  Y2 <- matrix(0,m,n)
  Y3 <- matrix(0,m,n)
  
  if(family == "gaussian")
  {
    YY1 <- X1%*% beta[[1]]
    YY2 <- X2%*% beta[[2]]
    YY3 <- X3%*% beta[[3]]
    for(i in 1:m){
      for(j in 1:n){
        theta1[j] <- YY1[j]+rho*sum(Y_t1*W[j,])
        Y_t1[j] <- rnorm(1,theta1[j],1)
        theta2[j] <- YY2[j]+rho*sum(Y_t2*W[j,])
        Y_t2[j] <- rnorm(1,theta2[j],1)
        theta3[j] <- YY3[j]+rho*sum(Y_t3*W[j,])
        Y_t3[j] <- rnorm(1,theta3[j],1)
        if(j==n){
          Y1[i,] <- Y_t1
          Y2[i,] <- Y_t2
          Y3[i,] <- Y_t3
        }
      }
      if(i%%100==0){print(i)}
    }
    Y1 <- apply(Y1[(m0+1):m,],2,mean)
    Y2 <- apply(Y2[(m0+1):m,],2,mean)
    Y3 <- apply(Y3[(m0+1):m,],2,mean)
  }
  if(family == "binomial")
  {
    YY1 <- X1%*%beta[[1]]
    YY2 <- X2%*%beta[[2]]
    YY3 <- X3%*%beta[[3]]
    
    for(i in 1:m){
      for(j in 1:n){
        theta1[j] <- YY1[j]+rho*sum(Y_t1*W[j,])
        Y_t1[j] <- rbinom(1,1,invlogit(theta1[j]))
        theta2[j] <- YY2[j]+rho*sum(Y_t2*W[j,])
        Y_t2[j] <- rbinom(1,1,invlogit(theta2[j]))
        theta3[j] <- YY3[j]+rho*sum(Y_t3*W[j,])
        Y_t3[j] <- rbinom(1,1,invlogit(theta3[j]))
        if(j==n){
          Y1[i,] <- Y_t1
          Y2[i,] <- Y_t2
          Y3[i,] <- Y_t3
        }
      }
      if(i%%100==0){print(i)}
    }
    Y1 <- apply(Y1[(m0+1):m,],2,getmode)
    Y2 <- apply(Y2[(m0+1):m,],2,getmode)
    Y3 <- apply(Y3[(m0+1):m,],2,getmode)
  }
  if( family == "poisson")
  {
    YY1 <- X1%*%beta[[1]]
    YY2 <- X2%*%beta[[2]]
    YY3 <- X3%*%beta[[3]]
    
    for (i in 1:m){
      for(j in 1:n){
        theta1[j] <- YY1[j]+rho*sum(Y_t1*W[j,])
        Y_t1[j] <- rpois(1,exp(theta1[j]))
        theta2[j] <- YY2[j]+rho*sum(Y_t2*W[j,])
        Y_t2[j] <- rpois(1,exp(theta2[j]))
        theta3[j] <- YY3[j]+rho*sum(Y_t3*W[j,])
        Y_t3[j] <- rpois(1,exp(theta3[j]))
        if(j==n){
          Y1[i,] <- Y_t1
          Y2[i,] <- Y_t2
          Y3[i,] <- Y_t3
        }
      }
      if(i%%100==0){print(i)}
    }
    Y1 <- apply(Y1[(m0+1):m,],2,getmode)
    Y2 <- apply(Y2[(m0+1):m,],2,getmode)
    Y3 <- apply(Y3[(m0+1):m,],2,getmode)
  }
  X <- list(X1,X2,X3)
  Y <- list(Y1,Y2,Y3)
  return(list(X=X,Y=Y))
}

tr<-function(x)
{
  return(sum(diag(x)))
}

adjust.Rsquare<-function(Y, Yhat, s)
{
  r2 = 1 - sum((Y - Yhat)^2)/sum((Y-mean(Y))^2)
  #return(r2)
  n = length(Y)
  adj_r2 = 1 - (1-r2)*(n-1)/(n-s-2)
  return(adj_r2)
}

lm.sparse<-function(Y, X)
{
  beta = ginv(as.matrix(crossprod(X)))%*%t(X)%*%(Y)
  Yhat = as.numeric(X%*%beta)
  return(Yhat)
}


projX<-function(X)
{
  return(X%*%ginv(as.matrix(crossprod(X)))%*%t(X))
}