# All EM-related code including calculation of sqrt(n)-consistent estimator of theta as starting value
# EQUAL VARIANCES

#library(mixtools)
#library(optimx)
# source("fisher.info.R")


mle.em <- function(x=xf, y=yf, K, theta0 = theta0[[1]], intercept, arbvar, model = "IID", init, maxrestarts = 5, rep=5){
  
  if(model == "HMM") stop("EM algorithm for HMM not yet implemented")
  
  d <- ncol(x) + intercept
  N <- nrow(x)
  xx <- x; if(intercept) xx <- cbind(rep(1,N), xx)
  converged <- FALSE
  posdef    <- FALSE
  restarts  <- -1
  mods <- vector("list")
  llks <- rep(NA,rep)
  
  for(r in 1:rep){
    if(model == "IID"){
      #while(!(posdef && converged)){
      if(init == "regmix"){
        inits <- mixtools::regmix.init(y, xx, k=2, addintercept = FALSE, arbvar=arbvar)
        beta0 <- inits$beta
        #print(beta0)
        sigma0 <- inits$s
        gamma0 <- matrix(inits$lambda, nrow=K, ncol=K, byrow=TRUE)
      }
      if(init == "random"){
        beta0 <- matrix(runif(d*K,-2,2),d,K)
        sigma0 <- runif(1+(K-1)*arbvar,.01,1)
        gamma0 <- matrix(runif(K^2,0,1/K),K,K) 
        diag(gamma0) <- sapply(1:K, function(k) 1-sum(gamma0[k,-k]))
      }
      if(init == "true"){
        beta0 <- rbind(theta0$mu, theta0$beta)
        #print(beta0)
        sigma0 <- theta0$sigma
        gamma0 <- theta0$gamma
      }
      
      if(model == "IID") gamma0 <- gamma0[1,]
      
      if(ncol(x) > 0){
        em <- mixtools::regmixEM(y=y, x=xx, k=K, beta=beta0, sigma=sigma0, lambda = gamma0, arbmean = TRUE, arbvar = arbvar, addintercept = FALSE)
        beta <- em$beta
      }else{
        em <- mixtools::normalmixEM(y, k=K, mu=c(beta0), sigma=sigma0, lambda = gamma0, maxit = 1e4, arbmean = TRUE, arbvar = arbvar)
        beta <- matrix(em$mu, 1, K)
      }
      mods[[r]] <- em
      sigma  <- em$sigma; if(!arbvar) sigma <- rep(sigma,K)
      lambda  <- em$lambda
      llks[r] <- loglik.fun(xx, y, beta, sigma, lambda)
    }
    ind <- which.min(llks[1:rep])
    mod <- mods[[ind]]
  }
  
  if(ncol(x) > 0) beta <- mod$beta
  if(ncol(x) == 0) beta <- matrix(mod$mu, 1, K)
  sigma  <- mod$sigma; if(!arbvar) sigma <- rep(sigma,K)
  lambda  <- mod$lambda
  h <- apply(mod$posterior, 1, which.max)
  
  
  fit <- fit.model(xx, y, h, beta)
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
    # if(!arbvar) sigma0 <- sigma0[1]
    if(model=="IID") gamma0 <- gamma0[1,]
    if(is.null(mu0)) mu0 <- rep(0,K)
    if(intercept) beta0 <- rbind(mu0, beta0)
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