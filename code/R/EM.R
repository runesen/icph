# All EM-related code including calculation of sqrt(n)-consistent estimator of theta as starting value
# EQUAL VARIANCES

#library(mixtools)
#library(optimx)
# source("fisher.info.R")


mle.em <- function(x=xf, y=yf, K, theta0 = theta0[[1]], intercept, arbvar, model = "IID", init, maxrestarts = 5){
  
  d <- ncol(x) + intercept
  N <- nrow(x)
  xx <- x; if(intercept) xx <- cbind(rep(1,N), xx)
  converged <- FALSE
  posdef    <- FALSE
  restarts  <- -1
  
  if(model == "HMM") stop("EM algorithm for HMM not yet implemented")
  if(model == "IID"){
    #while(!(posdef && converged)){
    if(init == "regmix"){
      inits <- mixtools::regmix.init(y, xx, k=2, addintercept = FALSE, arbvar=arbvar)
      beta0 <- inits$beta
      sigma0 <- inits$s
    }
    if(init == "random"){
      beta0 <- matrix(runif(d*K,-2,2),d,K)
      sigma0 <- runif(K,.01,1)
    }
    if(init == "true"){
      beta0 <- theta0$beta
      sigma0 <- theta0$sigma
    }
    if(ncol(x) > 0){
      em <- mixtools::regmixEM(y=y, x=x, k=K, beta=beta0, sigma=sigma0, arbmean = TRUE, arbvar = arbvar, addintercept = intercept)
      beta <- em$beta
    }else{
      em <- mixtools::normalmixEM(y, k=K, mu=c(beta0), sigma=sigma0, maxit = 1e4, arbmean = TRUE, arbvar = arbvar)
      beta <- matrix(em$mu, 1, K)
    }
    sigma  <- em$sigma; if(!arbvar) sigma <- rep(sigma,K)
    lambda  <- em$lambda
    h       <- apply(em$posterior, 1, which.max)
  }
  
  fit <- fit.model(x, y, h, beta, intercept)
  fisher <- fisher.info(x, y, beta, sigma, lambda, intercept, arbvar)
  jacobian <- diag(nrow(fisher))
  if(arbvar)  diag(jacobian)[(d*K+1):(d*K+K)] <- 1/(2*sigma)
  if(!arbvar) diag(jacobian)[d*K+1] <- 1/(2*sigma[1])
  covmat <- fisher2covmat(fisher, jacobian, d, K, arbvar, model, intercept)
  
  loglik <- loglik.fun(xx, y, beta, sigma, lambda)
  loglik0 <- NA
  if(!is.null(theta0)){
    mu0 <- theta0$mu
    beta0 <- theta0$beta
    sigma0 <- theta0$sigma
    gamma0 <- theta0$gamma
    #if(!arbvar) sigma0 <- sigma0[1]
    if(model=="IID") gamma0 <- gamma0[1,]
    if(is.null(mu0)) mu0 <- rep(0,K)
    if(intercept) beta0 <- rbind(mu0, beta0)
    xx <- x; if(intercept) xx <- cbind(rep(1,N),xx)
    loglik0 <- loglik.fun(xx, y, beta0, sigma0, gamma0)
  }
  
  # For now, since Markov procedure not yet implemented
  gamma <- lambda
  if(model == "IID"){
    gammaa0 <- lapply(1:K, function(k) theta0$gamma[1,k])
    gammaa <- lapply(1:K, function(k) gamma[k])
  }
  if(model == "HMM"){
    gammaa0 <- lapply(1:K, function(k) theta0$gamma[k,-k])
    gammaa <- lapply(1:K, function(k) gamma[k,-k])
  }
  
  if(intercept){
    theta <- lapply(1:K, function(i) list(intercept = beta[1,i], beta = beta[-1,i], sigma = sigma[i], gamma = gammaa[[i]]))
    theta0 <- lapply(1:K, function(i) list(intercept = theta0$mu[i], beta = theta0$beta[,i],
                                           sigma = theta0$sigma[i], gamma = gammaa0[[i]]))
  }
  if(!intercept){
    theta <- lapply(1:K, function(i) list(beta = beta[,i], sigma = sigma[i], gamma = gammaa[[i]]))
    theta0 <- lapply(1:K, function(i) list(beta = theta0$beta[,i], sigma = theta0$sigma[i], gamma = gammaa0[[i]]))
  }
  
  return(structure(list(N = N, K = K, d = d, x = x, y = y, 
                        theta = theta, theta0 = theta0, covmat = covmat, 
                        h = h, fitted = fit$fitted, residuals = fit$residuals,
                        loglik = loglik, loglik0 = loglik0, parameter = colnames(covmat),
                        arbvar = arbvar, restarts = restarts,
                        component = rep(1:K, each = d+2),
                        intercept = intercept, 
                        estimation = "EM", model = model),
                   class = "modelfit"))
}


fit.model <- function(x, y, h, beta, intercept){
  
  N <- length(y)
  K <- ncol(beta)
  if(intercept) x <- cbind(rep(1,N), x)
  fitted <- residuals <- numeric(N)
  
  for(i in 1:K){
    fitted[h==i]     <- x[h==i,,drop=F]%*%beta[,i]
    residuals[h==i]  <- y[h==i] - fitted[h==i]
  }
  
  return(list(fitted = fitted,
              residuals = residuals))
}

loglik.fun <- function(x=xx, y, beta=beta0, sigma=sigma0, lambda=lambda0){
  
  N <- length(y)
  K <- ncol(beta)
  
  sum <- 0
  for(i in 1:N){
    sum <- sum + log(sum(sapply(1:K, function(j) lambda[j]*dnorm(y[i], x[i,]%*%beta[,j], sigma[j]))))
  }
  return(sum)
}