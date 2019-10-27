# source("fisher.info.R")

mle.nlm <- function(x=xf, y=yf, K=3, theta0=theta0[[i]], intercept=TRUE, arbvar = FALSE, model = "IID", init="regmix", penalty, rep=5){
  
  N <- length(y)
  d <- ncol(x)+intercept
  xx <- x; if(intercept) xx <- cbind(rep(1,N), xx)
  
  lambda_seq <- unique(10^(-(0:5))*penalty)
  mod<-mods<-vector("list")
  llks<-rep(NA,rep)
  bics <- rep(NA,length(lambda_seq))
  
  # Fit some models and choose the best
  for(l in 1:length(lambda_seq)){
    for(r in 1:rep){
      if(init == "regmix"){
        inits <- mixtools::regmix.init(y, xx, k=K, addintercept = FALSE, arbvar=arbvar)
        beta0 <- inits$beta
        sigma0 <- inits$s + 0.01 # to comply with constraints in NLM method
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
        sigma0 <- theta0$sigma
        gamma0 <- theta0$gamma
      }
      
      # to avoid overparametrization
      if(model == "IID") gamma0 <- gamma0[1,]
      if(!arbvar) sigma0 <- sigma0[1]
      
      lambda <- lambda_seq[l]
      mods[[r]] <- suppressWarnings(mle(xx,y,beta0,sigma0,gamma0,arbvar,model,lambda,pl=0))
      llks[r]<-mods[[r]]$mllk
    }
    ind<-which.min(llks[1:rep])
    mod[[l]]<-mods[[ind]]
    bics[l] <- BIC(mod[[l]]$workingparas, y, xx, d, K, arbvar, model)
  }
  
  ind <- which.min(bics)
  mod <- mod[[ind]]
  bic <- bics[ind]
  parvect <- mod$workingparas
  
  # Estimates of natural parameters
  beta <- mod$beta
  sigma <- mod$sigma
  gamma <- mod$gamma
  
  # make sure sigma is always of length K, and gamma is K x K.
  if(!arbvar) sigma <- rep(sigma, K)
  if(model == "IID") gamma <- matrix(gamma, K, K, byrow=T)
  
  h        <- viterbi(xx, y, K, beta, sigma, gamma, arbvar, model)
  fit      <- fit.model(xx, y, h, beta)
  jacobian <- dpw2pn(parvect, d, K, arbvar, model)
  fisher   <- mod$hessian
  covmat   <- fisher2covmat(fisher, jacobian, d, K, arbvar, model, intercept)
  
  loglik <- -mllk(parvect, y, xx, d, K, arbvar, model, lambda=0)
  loglik0 <- NA
  if(!is.null(theta0)){
    mu0 <- theta0$mu
    beta0 <- theta0$beta
    sigma0 <- theta0$sigma
    gamma0 <- theta0$gamma
    if(!arbvar) sigma0 <- sigma0[1]
    if(model=="IID") gamma0 <- gamma0[1,]
    if(is.null(mu0)) mu0 <- rep(0,K)
    if(intercept) beta0 <- rbind(mu0, beta0)
    parvect0 <- pn2pw(beta0, sigma0, gamma0, arbvar, model)
    loglik0 <- -mllk(parvect0, y, xx, d, K, arbvar, model, lambda=0)
  }
  
  if(model == "IID"){
    gammaa0 <- lapply(1:K, function(k) theta0$gamma[1,k])
    gammaa <- lapply(1:K, function(k) gamma[k,k])
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
                        arbvar = arbvar, model = model, bic = bic,
                        component = rep(1:K, each = d+2),
                        intercept = intercept, #whichInfVar=w, 
                        estimation = "NLM"),
                   class = "modelfit"))
}

pn2pw <- function(beta=beta0, sigma=sigma0, gamma=gamma0, arbvar, model){
  
  K <- ncol(beta)
  tbeta <- c(beta)
  tsigma <- log(sigma-0.01*arbvar) # if arbvar=T, we impose a lower threshold on the sigmas here (otherwise the likelihood is unbounded)
  if(model == "HMM"){
    foo <- log(gamma/diag(gamma)) 
    tgamma <- as.vector(foo[!diag(K)])
  }
  if(model == "IID"){
    tgamma <- log(gamma[-1]/gamma[1])
  }
  parvect <- c(tbeta,tsigma,tgamma)
  return(parvect)
}

# depending on model and arbvar, this function expects different dimensions of parvect...
pw2pn <- function(parvect=parvect0, d, K, arbvar, model){
  
  beta <- matrix(parvect[1:(d*K)], d, K)
  if(arbvar){
    sigma <- exp(parvect[(d*K+1):(d*K+K)])+0.01
    if(model == "HMM"){
      gamma <- diag(K)
      gamma[!gamma] <- exp(parvect[(d*K+K+1):length(parvect)])
      gamma <- gamma/apply(gamma,1,sum)
    }
    if(model == "IID"){
      gamma <- rep(1,K)
      gamma[2:K] <- exp(parvect[(d*K+K+1):length(parvect)])
      gamma <- gamma / sum(gamma)
    }
  }
  if(!arbvar){
    sigma <- exp(parvect[d*K+1])
    if(model == "HMM"){
      gamma <- diag(K)
      gamma[!gamma] <- exp(parvect[(d*K+1+1):length(parvect)])
      gamma <- gamma/apply(gamma,1,sum)
    }
    if(model == "IID"){
      gamma <- rep(1,K)
      gamma[2:K] <- exp(parvect[(d*K+1+1):length(parvect)])
      gamma <- gamma / sum(gamma)
    }
  }
  return(list(beta=beta,sigma=sigma,gamma=gamma))
}

dpw2pn <- function(parvect, d, K, arbvar, model, silent = TRUE){
  
  gradient <- diag(length(parvect))
  if(arbvar){
    diag(gradient)[(d*K+1):(d*K+K)] <- exp(parvect[(d*K+1):(d*K+K)])
  }
  if(!arbvar){
    gradient[d*K+1, d*K+1] <- exp(parvect[d*K+1])
  }
  rowsInGamma <- 1+(K-1)*(model=="HMM")
  indexInPar <- length(parvect)-(K-1)*rowsInGamma
  for(r in 1:rowsInGamma){
    par <- parvect[(indexInPar+(r-1)*(K-1)+1):(indexInPar+r*(K-1))]
    for(i in 1:(K-1)){
      for(j in 1:(K-1)){
        gradient[indexInPar+(r-1)*(K-1)+i, indexInPar+(r-1)*(K-1)+j] <- ((i==j)*exp(par[i])*(1+sum(exp(par)))-exp(par[i]+par[j]))/(1+sum(exp(par)))^2
      }
    }
  }
  cutOff <- 0
  for(i in 1:nrow(gradient)){
    for(j in 1:nrow(gradient)){
      if(gradient[i,j]!=0 & abs(gradient[i,j]) < 1e-2){
        gradient[i,j] <- sign(gradient[i,j])*pmax(abs(gradient[i,j]),1e-2)
        cutOff <- cutOff + 1
      }
    }
  }
  if(cutOff>0 & !silent) warning("Some non-zero entries in the Jacobian are less than 10^(-3) and will be regularized")
  return(gradient)
}

if(0){
  mllk <- function(parvect=parvect0, y, x=xx, d, K, arbvar, model, lambda){
    
    lpn <- pw2pn(parvect, d, K, arbvar, model)
    beta <- lpn$beta
    sigma <- lpn$sigma
    gamma <- lpn$gamma
    if(!arbvar) sigma <- rep(sigma,K)
    if(model=="IID") gamma <- matrix(gamma, K, K, byrow=T)
    
    N <- length(y)
    logprobs <- matrix(rep(0,K*N),nrow=N)
    for (j in 1:K){
      logprobs[,j] <- dnorm(y,x%*%beta[,j],sigma[j], log = TRUE) # emission probabilities
    }
    
    logf      <- matrix(NA, nrow = N, ncol = K)
    logf[1,]  <- c(logprobs[1,1], rep(-Inf,K-1))
    
    for (n in 2:N) {
      for (k in 1:K) {
        logsum <- -Inf
        for (j in 1:K) {
          logsummand <- logf[n-1, j] + log(gamma[j,k])
          if (logsum - logsummand < 700) {
            logsum <- logsummand + log(1 + exp(logsum - logsummand)) # log(exp(logsummand)+exp(logsum)) = log(summand + sum)
          }
        }
        logf[n, k] <- logprobs[n,k] + logsum
      }
    }
    
    loglik <- logf[N, 1]
    
    for (k in 2:K) {
      temp <- logf[N, k]
      if (loglik - temp < 700) {
        loglik <- temp + log(1 + exp(loglik - temp))
      }
    }
    
    # loglik <- log(sum(exp(logf[N,])))
    penalty <- penalty.scad(gamma,lambda, N,d, arbvar)
    print(paste("loglik: ", loglik))
    print(paste("penalty: ", penalty))
    if(lambda != 0) loglik <- loglik - penalty
    return(-loglik)
  }
}

#if(0){
mllk <- function(parvect=parvect0, y, x=xx, d, K, arbvar, model, lambda){
  
  lpn <- pw2pn(parvect, d, K, arbvar, model)
  beta <- lpn$beta
  sigma <- lpn$sigma
  gamma <- lpn$gamma
  if(!arbvar) sigma <- rep(sigma,K)
  if(model=="IID") gamma <- matrix(gamma, K, K, byrow=T)
  
  N <- length(y)
  allprobs <- matrix(rep(1,K*N),nrow=N)
  for (j in 1:K){
    allprobs[,j] <- dnorm(y,x%*%beta[,j],sigma[j]) # emission probabilities
  }
  #foo <- lpn$delta  
  foo <- rep(1/K,K)
  lscale <- 0
  for (i in 1:N){
    #if(sum(allprobs[i,])!=0){
    foo <- foo%*%gamma*allprobs[i,] 
    sumfoo <- sum(foo)
    lscale <- lscale+log(sumfoo)
    foo <- foo/sumfoo # scaling to avoid numerical underflow
    #}
  }
  
  if(lambda != 0) lscale <- lscale - penalty.scad(gamma,lambda, N,d, arbvar)
  return(-lscale)
}
#}

BIC <- function(parvect=parvect0, y, x=xx, d, K, arbvar, model){
  
  lpn <- pw2pn(parvect, d, K, arbvar, model)
  beta <- lpn$beta
  sigma <- lpn$sigma
  gamma <- lpn$gamma
  if(!arbvar) sigma <- rep(sigma,K)
  if(model=="IID") gamma <- matrix(gamma, K, K, byrow=T)
  
  Mhat <- K
  for(k in 1:K){
    if(all(gamma[,k] < .001)){
      gamma[,k] <- 0
      Mhat <- Mhat - 1
    }
  }
  
  N <- length(y)
  allprobs <- matrix(rep(1,K*N),nrow=N)
  for (j in 1:K){
    allprobs[,j] <- dnorm(y,x%*%beta[,j],sigma[j]) # emission probabilities
  } ## Here, use log=TRUE in dnorm to avoid underflow. Change if loop below suitably 
  #foo <- lpn$delta  
  foo <- rep(1/K,K)
  lscale <- 0
  for (i in 1:N){
    if(sum(allprobs[i,])!=0){
      foo <- foo%*%gamma*allprobs[i,] 
      sumfoo <- sum(foo)
      lscale <- lscale+log(sumfoo)
      foo <- foo/sumfoo # scaling to avoid numerical underflow
    }
  }
  
  bic <- lscale - .5*Mhat*(1+d+arbvar)*log(N)
  return(-bic)
}

# function that runs the numerical minimization of 'mllk' (i.e. tries to find the MLE)
mle <- function(x, y,beta0,sigma0,gamma0,arbvar,model,lambda,pl){
  d <- nrow(beta0)
  K <- ncol(beta0)
  parvect <- pn2pw(beta0,sigma0,gamma0,arbvar,model)
  mod <- nlm(mllk,parvect,y,x,d,K,arbvar,model,lambda,print.level=pl,iterlim=1000,stepmax=5,hessian=TRUE)
  pn <- pw2pn(mod$estimate, d, K,arbvar,model)
  return(list(beta=pn$beta,sigma=pn$sigma,gamma=pn$gamma,workingparas=mod$estimate,hessian=mod$hessian,mllk=mod$minimum))
}

viterbi <- function(x, y, K, beta, sigma, gamma, arbvar, model){
  
  if(!arbvar) sigma <- rep(sigma, K)
  if(model == "IID") gamma <- matrix(gamma, K, K, byrow=T)
  N <- length(y)
  logprobs <- matrix(rep(1,K*N),nrow=N)
  for (j in 1:K){
    logprobs[,j] <- dnorm(y,x%*%beta[,j],sigma[j], log=TRUE)
  }
  logmaxProbs <- matrix(NA, N, K)
  logmaxProbs[1,] <- c(logprobs[1,1], rep(-Inf, K-1))
  
  for (n in 2:N) {
    for (k in 1:K) {
      maxi <- NULL
      for (j in 1:K) {
        temp <- logmaxProbs[n-1, j] + log(gamma[j,k])
        maxi <- max(maxi, temp) # after looping over all j, this is the max of all
      }
      logmaxProbs[n,k] <- maxi + logprobs[n,k]
    }
  }
  viterbiPath <- rep(NA, N)
  viterbiPath[N] <- which.max(logmaxProbs[N, ])
  
  for (n in (N - 1):1){
    viterbiPath[n] <- which.max(logmaxProbs[n,] + log(c(gamma[, viterbiPath[n+1]])))
  } 
  return(viterbiPath)
}


penalty.scad <- function(gamma, lambda, N, d, arbvar,eps=10^(-6)){
  
  return(lambda*N*(d+arbvar)*sum(apply(gamma,2,function(c) log(eps+SCAD(sqrt(mean(c^2)), lambda))-log(eps))))
}

SCAD <- function(g, lambda, a = 3.7){
  c1 <- g
  c2 <- (2*a*lambda*g-g^2)/(2*(a-1)*lambda)
  c3 <- a*lambda
  scad <- c1*(g<=lambda) + c2*(g > lambda && g <= a*lambda) + c3*(g>a*lambda)
  return(scad)
}
