Rcpp::sourceCpp('eigenvalues.cpp')

#Screening
NW_SIS <-function(X,Y,W,k)
{
  Z <- cbind(Y, W%*%Y)
  PZ <- Z%*%solve(crossprod(Z))%*%t(Z)
  R2 <- colSums(X*(PZ%*%X))/colSums(X*X)
  return(list(ix=sort(order(R2,decreasing = T)[1:k])))
}

SMLE <- function(X,Y,family,beta_initial,k,tol,maxiter)
{
  p = ncol(X)
  n <- nrow(X)
  beta_old <- beta_initial
  iter <- 1
  Umin <- 1
  Umax <- 1e8
  c <- 0.01
  tau <- 2
  U0 <- 1
  Terrul <- 1
  while((iter<=maxiter)&(Terrul>tol))
  {
    old <- LikScore(X,Y,beta_old,family)
    u <- U0
    beta_new <- beta_old + (1/u)*old$score
    index <- order(abs(beta_new),decreasing = TRUE)[1:k]
    beta_new[order(abs(beta_new),decreasing = TRUE)[(k+1):p]] = 0
    new <- Lik(X,Y,beta_new,family)
    count <- 1
    while(as.logical(new$lik<old$lik + 0.5*c*u*(sum((beta_old-beta_new)^2))) & u<Umax)
    {
      u <- min(tau*u,Umax);
      beta_new <- beta_old + (1/u)*old$score
      index <- order(abs(beta_new),decreasing = TRUE)[1:k]
      beta_new[order(abs(beta_new),decreasing = TRUE)[(k+1):p]] = 0
      new <- Lik(X,Y,beta_new,family)
      count <- count+1
    }
    new <- LikScore(X,Y,beta_new,family)
    diff <- sum((beta_old-beta_new)^2);
    U0 <- min(max(-t(beta_new-beta_old)%*%(new$score-old$score)/(diff),Umin),Umax);
    Terrul <- sqrt(sum((beta_new-beta_old))^2)/sqrt(sum(beta_new^2))
    beta_old <- beta_new;
    iter <- iter+1;
  }
  if(family == "gaussian"){
    Y_hat <- X%*%beta_new
    res <- sum((Y-Y_hat)^2)/n
  }else if(family == "binomial"){
    Y_hat <- 1*(1/(1+exp(-X%*%beta_new))>0.5)
    res <- sum(Y!=Y_hat)/n
  }else if (family == "poisson"){
    mu <- exp(X%*%beta_new)
    #res <- mean(Y*log(Y/mu)-Y+mu)
    res <- mean((Y-mu)^2/mu)
  }
  return(list(ix=sort(index),beta=beta_new,iter=iter,res=res))
}

PMLE_NPG <- function(X,Y,W,family,theta_initial,k,tol,maxiter)
{
  p <- ncol(X)
  Z <- as.matrix(cbind(X,W%*%Y))
  theta_old <- theta_initial
  iter <- 1
  Umin <- 1
  Umax <- 1e+08
  tau <- 2 #omega
  c <- 0.01
  M <- 4
  U0 <- 1
  Fl <- rep(Inf, M)
  Terrul <- 1
  ll <- c()
  while((iter<=maxiter)&(Terrul>tol))
  {
    old <- LikScore(Z,Y,theta_old,family)
    u <- U0
    if(iter==1){Fl[1] <- old$lik}
    theta_new <- prox(theta_old,u,old$score,k)
    new <- Lik(Z,Y,theta_new,family)
    count <- 1
    while(as.logical(new$lik-min(Fl)<0.5*c*u*(sum((theta_old-theta_new)^2)))&u<Umax)
    {
      u <- min(tau*u,Umax)
      theta_new <- prox(theta_old,u,old$score,k)
      new <- Lik(Z,Y,theta_new,family)
      count <- count+1
    }
    new <- LikScore(Z,Y,theta_new,family)
    diff <- sum((theta_old-theta_new)^2);
    U0 <- min(max(-t(theta_new-theta_old)%*%(new$score-old$score)/(diff),Umin),Umax);
    Fl[iter%%M+1] <-  as.numeric(new$lik)
    ll <- c(ll,as.numeric(new$lik))
    Terrul <- sqrt(sum((theta_new-theta_old))^2)/sqrt(sum(theta_new^2))
    theta_old <- theta_new
    cat("iter", iter,"\n",
        "Terrul",Terrul,"\n")
    iter <- iter+1
  }
  index <- which(theta_new[1:p]!=0)
  if(family == "gaussian"){
    Y_hat <- Z%*%theta_new
    res <- sum((Y-Y_hat)^2)/n
  }else if(family == "binomial"){
    Y_hat <- 1*(1/(1+exp(-Z%*%theta_new))>0.5)
    res <- sum(Y!=Y_hat)/n
  }else if (family == "poisson"){
    mu <- exp(Z%*%theta_new)
    #res <- mean(Y*log(Y/mu)-Y+mu)
    res <- mean((Y-mu)^2/mu)
  }
  return(list(ix=sort(index),theta=theta_new,iter=iter,ll=ll,res=res))
}

PMLE_MAPG <- function(X,Y,W,family,theta_initial,k,tol,maxiter)
{
  p <- ncol(X)
  Z <- as.matrix(cbind(X,W%*%Y))
  theta1 <- theta0 <- gamma0 <- theta_initial
  score <- LikScore(Z,Y,theta0,family)$score
  alpha1 <- delta1 <- prox(gamma0,200,score,k)
  omega <- 2
  c <- 0.01
  Umin <- 1
  Umax <- 1e8
  k0 <- 0
  k1 <- 1
  ll <- c()
  iter <- 1
  Terrul <- 1
  while((iter<=maxiter)&(Terrul>tol))
  {
    gamma1 <- theta1 + (k0/k1)*(alpha1-theta1) + ((k0-1)/k1)*(theta1-theta0)
    u_gamma1 <- t(alpha1-gamma0)%*%(LikScore(Z,Y,gamma0,family)$score-LikScore(Z,Y,alpha1,family)$score)
    u_gamma2 <- sum((alpha1-gamma0)^2)
    u_gamma <- as.numeric(u_gamma1)/u_gamma2
    u_gamma <- min(max(u_gamma,Umin),Umax)
    u_theta1 <- t(delta1-theta0)%*%(LikScore(Z,Y,theta0,family)$score-LikScore(Z,Y,delta1,family)$score)
    u_theta2 <- sum((delta1-theta0)^2)
    u_theta <- as.numeric(u_theta1)/u_theta2
    u_theta <- min(max(u_theta,Umin),Umax)
    old1 <- LikScore(Z,Y,gamma1,family)
    alpha2 <- prox(gamma1,u_gamma,old1$score,k)
    new1 <- Lik(Z,Y,alpha2,family)
    count1 <- 1
    while(as.logical(new1$lik < old1$lik + 0.5*c*u_gamma*(sum((alpha2-gamma1)^2))) & u_gamma < Umax)
    {
      u_gamma <- min(omega*u_gamma,Umax)
      alpha2 <- prox(gamma1,u_gamma,old1$score,k)
      new1 <- Lik(Z,Y,alpha2,family)
      count1 <- count1+1
    }
    #cat("count1",count1,"\n")
    old2 <- LikScore(Z,Y,theta1,family)
    delta2 <- prox(theta1,u_theta,old2$score,k)
    new2 <- Lik(Z,Y,delta2,family)
    count2 <- 1
    while(as.logical(new2$lik < old2$lik + 0.5*c*u_theta*(sum((delta2-theta1)^2))) & u_theta < Umax)
    {
      u_theta <- min(omega*u_theta,Umax)
      delta2 <- prox(theta1,u_theta,old2$score,k)
      new2 <- Lik(Z,Y,delta2,family)
      count2 <- count2+1
    }
    #cat("count2",count2,"\n")
    if(new1$lik >= new2$lik){
      theta2 <- alpha2
    }else{
      theta2 <- delta2 
    }
    k0 <- k1
    k1 <- (sqrt(4*k1^2+1)+1)/2
    ll <- c(ll,as.numeric(Lik(Z,Y,theta2,family)$lik))
    Terrul <- sqrt(sum((theta2-theta1))^2)/sqrt(sum(theta2^2))
    theta0 <- theta1
    theta1 <- theta2
    gamma0 <- gamma1
    alpha1 <- alpha2
    delta1 <- delta2
    cat("iter", iter,"\n",
        "Terrul",Terrul,"\n")
    iter <- iter+1
  }
  index <- which(theta2[1:p]!=0)
  if(family == "gaussian"){
    Y_hat <- Z%*%theta2
    res <- sum((Y-Y_hat)^2)/n
  }else if(family == "binomial"){
    Y_hat <- 1*(1/(1+exp(-Z%*%theta2))>0.5)
    res <- sum(Y!=Y_hat)/n
  }else if (family == "poisson"){
    mu <- exp(Z%*%theta2)
    #res <- mean(Y*log(Y/mu)-Y+mu)
    res <- mean((Y-mu)^2/mu)
  }
  return(list(ix=sort(index),theta=theta2,iter=iter,ll=ll,res=res))
}

PMLE_NAPG <- function(X,Y,W,family,theta_initial,k,tol,maxiter)
{
  p <- ncol(X)
  Z <- as.matrix(cbind(X,W%*%Y))
  alpha1 <- theta1 <- theta0 <- gamma0 <- theta_initial
  eta <- 0.8
  iter <- 1
  tau <- 2 #omega
  Umin <- 1
  Umax <- 1e8
  c <- 0.01
  c1 <- -Lik(Z,Y,theta1,family)$lik
  q1 <- 1
  k0 <- 0
  k1 <- 1
  Terrul <- 1
  ll <- c()
  while((iter<=maxiter)&(Terrul>tol))
  {
    gamma1 <- theta1 + (k0/k1)*(alpha1-theta1) + ((k0-1)/k1)*(theta1-theta0)
    if(iter==1){
      u_gamma <- 200
    }else{
      u_gamma1 <- t(gamma1-gamma0)%*%(LikScore(Z,Y,gamma0,family)$score-LikScore(Z,Y,gamma1,family)$score)
      u_gamma2 <- sum((gamma1-gamma0)^2)
      u_gamma <- as.numeric(u_gamma1)/u_gamma2
      u_gamma <- min(max(u_gamma,Umin),Umax)
    }
    old1 <- LikScore(Z,Y,gamma1,family)
    alpha2 <- prox(gamma1,u_gamma,old1$score,k)
    new1 <- Lik(Z,Y,alpha2,family)
    count1 <- 1
    while(as.logical(new1$lik < old1$lik + 0.5*c*u_gamma*(sum((alpha2-gamma1)^2)) || 
          new1$lik < -c1 + 0.5*c*u_gamma*(sum((alpha2-gamma1)^2))) & u_gamma < Umax)
    {
      u_gamma <- min(tau*u_gamma,Umax)
      alpha2 <- prox(gamma1,u_gamma,old1$score,k)
      new1 <- Lik(Z,Y,alpha2,family)
      count1 <- count1+1
    }
    #cat("count1",count1,"\n")
    if(new1$lik >= -c1 + 0.5*c*u_gamma*(sum((alpha2-gamma1)^2))){
      theta2 <- alpha2
    }else{
      if(iter==1){
        u_theta <- 200
      }else{
        u_theta1 <- t(theta1-gamma0)%*%(LikScore(Z,Y,gamma0,family)$score-LikScore(Z,Y,theta1,family)$score)
        u_theta2 <- sum((theta1-gamma0)^2)
        u_theta <- as.numeric(u_theta1)/u_theta2
        u_theta <- min(max(u_theta,Umin),Umax)
      }
      old2 <- LikScore(Z,Y,theta1,family)
      delta2 <- prox(theta1,u_theta,old2$score,k)
      new2 <- Lik(Z,Y,delta2,family)
      count2 <- 1
      while(as.logical(new2$lik < -c1 + 0.5*c*u_theta*(sum((delta2-theta1)^2))) & u_theta < Umax)
      {
        u_theta <- min(tau*u_theta,Umax)
        delta2 <- prox(theta1,u_theta,old2$score,k)
        new2 <- Lik(Z,Y,delta2,family)
        count2 <- count2+1
      }
      #cat("count2",count2,"\n")
      if(new1$lik >= new2$lik){
        theta2 <- alpha2
      }else{
        theta2 <- delta2 
      }
    }
    k0 <- k1
    k1 <- (sqrt(4*k1^2+1)+1)/2
    q2 <- eta*q1+1
    l2 <- Lik(Z,Y,theta2,family)$lik
    c1 <- (eta*q1*c1-l2)/q2
    q1 <- q2
    ll <- c(ll,as.numeric(l2))
    Terrul <- sqrt(sum((theta2-theta1))^2)/sqrt(sum(theta2^2))
    theta0 <- theta1
    theta1 <- theta2
    gamma0 <- gamma1
    alpha1 <- alpha2
    cat("iter", iter,"\n",
        "Terrul",Terrul,"\n")
    iter <- iter+1
  }
  index <- which(theta2[1:p]!=0)
  if(family == "gaussian"){
    Y_hat <- Z%*%theta2
    res <- sum((Y-Y_hat)^2)/n
  }else if(family == "binomial"){
    Y_hat <- 1*(1/(1+exp(-Z%*%theta2))>0.5)
    res <- sum(Y!=Y_hat)/n
  }else if (family == "poisson"){
    mu <- exp(Z%*%theta2)
    #res <- mean(Y*log(Y/mu)-Y+mu)
    res <- mean((Y-mu)^2/mu)
  }
  return(list(ix=sort(index),theta=theta2,iter=iter,ll=ll,res=res))
}


#Post estimation
SMLE_POST <- function(X,Y, p, fit_pre,family,tol,maxiter)
{
  X <- X[,fit_pre]
  fit <- glmnet(X,Y,family,intercept = F)
  beta_old <- fit$beta[,length(fit$lambda)]
  iter <- 1
  Umin <- 1
  Umax <- 1e8
  c <- 0.01
  tau <- 2
  U0 <- 1
  Terrul <- 1
  while((iter<=maxiter)&(Terrul>tol))
  {
    old <- LikScore(X,Y,beta_old,family)
    u <- U0
    beta_new <- beta_old + (1/u)*old$score
    new <- Lik(X,Y,beta_new,family)
    count <- 1
    while(as.logical(new$lik<old$lik + 0.5*c*u*(sum((beta_old-beta_new)^2))) & u<Umax)
    {
      u <- min(tau*u,Umax);
      beta_new <- beta_old + (1/u)*old$score
      new <- Lik(X,Y,beta_new,family)
      count <- count+1
    }
    new <- LikScore(X,Y,beta_new,family)
    diff <- sum((beta_old-beta_new)^2);
    U0 <- min(max(-t(beta_new-beta_old)%*%(new$score-old$score)/(diff),Umin),Umax);
    Terrul <- sqrt(sum((beta_new-beta_old))^2)/sqrt(sum(beta_new^2))
    beta_old <- beta_new;
    iter <- iter+1;
  }
  beta_hat <- rep(0,p)
  beta_hat[fit_pre] <- beta_new
  if(family == "gaussian"){
    Y_hat <- X%*%beta_new
    res <- sum((Y-Y_hat)^2)/n
  }else if(family == "binomial"){
    Y_hat <- 1*(1/(1+exp(-X%*%beta_new))>0.5)
    res <- sum(Y!=Y_hat)/n
  }else if (family == "poisson"){
    mu <- exp(X%*%beta_new)
    #res <- mean(Y*log(Y/mu)-Y+mu)
    res <- mean((Y-mu)^2/mu)
  }
  return(list(beta=beta_hat,iter=iter,res=res))
}

PMLE_POST <- function(X, Y, W, p, fit_pre,family,tol,maxiter)
{
  n <- dim(X)[1]
  X <- X[,fit_pre]
  Z <- cbind(X,W%*%Y)
  fit <- glmnet(Z,Y,family,intercept = F)
  theta_old <- fit$beta[,length(fit$lambda)]
  iter <- 1
  Umin <- 1
  Umax <- 1e+08
  tau <- 2 #omega
  c <- 0.01
  M <- 4
  U0 <- 1
  Fl <- rep(Inf, M)
  Terrul <- 1
  ll <- c()
  while((iter<=maxiter)&(Terrul>tol))
  {
    old <- LikScore(Z,Y,theta_old,family)
    u <- U0
    if(iter==1){Fl[1] <- old$lik}
    theta_new <- theta_old + (1/u)*old$score
    new <- Lik(Z,Y,theta_new,family)
    count <- 1
    while(as.logical(new$lik-min(Fl)<0.5*c*u*(sum((theta_old-theta_new)^2)))&u<Umax)
    {
      u <- min(tau*u,Umax)
      theta_new <- theta_old + (1/u)*old$score
      new <- Lik(Z,Y,theta_new,family)
      count <- count+1
    }
    new <- LikScore(Z,Y,theta_new,family)
    diff <- sum((theta_old-theta_new)^2);
    U0 <- min(max(-t(theta_new-theta_old)%*%(new$score-old$score)/(diff),Umin),Umax);
    Fl[iter%%M+1] <-  as.numeric(new$lik)
    ll <- c(ll,as.numeric(new$lik))
    Terrul <- sqrt(sum((theta_new-theta_old))^2)/sqrt(sum(theta_new^2))
    theta_old <- theta_new
    cat("iter", iter,"\n",
        "Terrul",Terrul,"\n")
    iter <- iter+1
  }
  theta_hat <- rep(0,p+1)
  theta_hat[c(fit_pre,p+1)] <- theta_new
  if(family == "gaussian"){
    Y_hat <- Z%*%theta_new
    res <- sum((Y-Y_hat)^2)/n
  }else if(family == "binomial"){
    Y_hat <- 1*(1/(1+exp(-Z%*%theta_new))>0.5)
    res <- sum(Y!=Y_hat)/n
  }else if (family == "poisson"){
    mu <- exp(Z%*%theta_new)
    #res <- mean(Y*log(Y/mu)-Y+mu)
    res <- mean((Y-mu)^2/mu)
  }
  return(list(theta=theta_hat,iter=iter,ll=ll,res=res))
}

Post_estimation_SMLE <- function(X, Y, p, fit_pre, family)
{
  n <- dim(X)[1]
  Z <- X[,fit_pre]
  if(family=="gaussian"){
    fit <- glm.fit(Z,Y,family = gaussian(),intercept = FALSE)
    beta_hat <- rep(0,p)
    beta_hat[fit_pre] <- fit$coefficients
    M_hat <- which(beta_hat!=0)
    M_hat <- M_hat[-length(M_hat)]
  }else if(family=="binomial"){
    fit <- glm.fit(Z,Y,family = binomial(),intercept = FALSE,control=list(maxit=200))
    beta_hat <- rep(0,p)
    beta_hat[fit_pre] <- fit$coefficients
    M_hat <- which(beta_hat!=0)
    M_hat <- M_hat[-length(M_hat)]
  }else{
    fit <- glm.fit(Z,Y,family = poisson(),intercept = FALSE,control=list(maxit=200))
    beta_hat <- rep(0,p)
    beta_hat[fit_pre] <- fit$coefficients
    M_hat <- which(beta_hat!=0)
    M_hat <- M_hat[-length(M_hat)]
  }
  return(list(ix=sort(M_hat),beta=beta_hat))
}

Post_estimation_PMLE <- function(X, Y, W, p, fit_pre, family)
{
  n <- dim(X)[1]
  Z <- cbind(X[,fit_pre],W%*%Y)
  if(family=="gaussian"){
    fit <- glm.fit(Z,Y,family = gaussian(),intercept = FALSE)
    theta_hat <- rep(0,p+1)
    theta_hat[c(fit_pre,p+1)] <- fit$coefficients
    M_hat <- which(theta_hat!=0)
    M_hat <- M_hat[-length(M_hat)]
  }else if(family=="binomial"){
    fit <- glm.fit(Z,Y,family = binomial(),intercept = FALSE,control=list(maxit=200))
    theta_hat <- rep(0,p+1)
    theta_hat[c(fit_pre,p+1)] <- fit$coefficients
    M_hat <- which(theta_hat!=0)
    M_hat <- M_hat[-length(M_hat)]
  }else{
    fit <- glm.fit(Z,Y,family = poisson(),intercept = FALSE,control=list(maxit=200))
    theta_hat <- rep(0,p+1)
    theta_hat[c(fit_pre,p+1)] <- fit$coefficients
    M_hat <- which(theta_hat!=0)
    M_hat <- M_hat[-length(M_hat)]
  }
  return(list(ix=sort(M_hat),theta=theta_hat))
}


#LASSO initial value of theta
theta_ini <- function(X,Y,W,family)
{ 
  Z <- cbind(X,W%*%Y)
  fit <- glmnet(Z,Y,family,intercept = F)
  theta <- fit$beta[,length(fit$lambda)]
  return(theta)
}

LikScore <- function(X,Y,beta,family)
{
  if(family == "gaussian")
  {
    likelihood <- t(Y)%*%X%*%beta-sum((X%*%beta)^2/2)
    score <- t(X)%*%(Y-(X%*%beta))
  }
  
  else if (family == "binomial")
  {
    likelihood <- t(Y)%*%X%*%beta-sum(log(1+exp(X%*%beta)))
    score <- t(X)%*%(Y-invlogit(X%*%beta))
  }
  
  else if (family == "poisson")
  {
    likelihood <- t(Y)%*%X%*%beta-sum(exp(X%*%beta))
    score <- t(X)%*%(Y-exp(X%*%beta))
  }
  
  else (print("Wrong input family"))
  return(list(lik = likelihood, score = score))
}

Lik <- function(X,Y,beta,family)
{
  if(family == "gaussian")
  {
    likelihood <- t(Y)%*%X%*%beta-sum((X%*%beta)^2/2)
  }
  
  else if (family == "binomial")
  {
    likelihood <- t(Y)%*%X%*%beta-sum(log(1+exp(X%*%beta)))
  }
  
  else if (family == "poisson")
  {
    likelihood <- t(Y)%*%X%*%beta-sum(exp(X%*%beta))
  }
  
  else (return(message("Wrong input family")))
  return(list(lik = likelihood))
}

prox <- function(gamma,u_gamma,score,k){
  p <- length(gamma)-1
  alpha <- gamma + (1/u_gamma)*score
  beta <- alpha[1:p]
  beta[order(abs(beta),decreasing = TRUE)[(k+1):p]] = 0
  index <- which(beta!=0)
  alpha[1:p] <- beta
  return(alpha)
}

measure <- function(estS,S)
{
  estS <- sort(estS)
  if(all(is.element(S,estS))){
    RC <- 1
  }else{
    RC <- 0
  }
  PSR <- length(intersect(estS,S))/length(S)
  FDR <- length(intersect(c(1:p)[-S],estS))/length(estS)
  if(length(estS)==length(S) && all(estS == S)){
    CSR <- 1
  }else{
    CSR <- 0
  }
  AMS <- sum(estS!=0)
  result <- c(RC, PSR, FDR, AMS)
  names(result) <- c("RC","PSR","FDR","AMS")
  return(result)
}
