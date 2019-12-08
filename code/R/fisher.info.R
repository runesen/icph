## Fisher information for for either arbvar = T or arbvar = F

# If arbvar = FALSE, sigma will still be supplied as K dim vector (equal entries), so that functions can run as usual
fisher.info <- function(x, y, beta, sigma=sqrt(sigma2), lambda, intercept, arbvar = TRUE){
  
  N <- length(y)
  if(intercept) x <- cbind(rep(1, N),x)
  
  K <- length(lambda)
  d <- ncol(x)
  
  Npar <- K+d*K+arbvar*K+(1-arbvar)*1
  
  tmp <- matrix(0, Npar, Npar)
  
  for(n in 1:N){
    for(i in 1:K){
      summand <- ddloglik(x, y, beta, sigma, lambda, i, n, arbvar) * posterior(x, y, beta, sigma, lambda, i, n) +
        matrix(dloglik(x, y, beta, sigma, lambda, i, n, arbvar),Npar,1) %*% dposterior(x, y, beta, sigma, lambda, i, n, arbvar)
      tmp <- tmp + summand
    }
  }
  
  return(-tmp)
}

posterior <- function(x, y, beta, sigma, lambda, i, n){
  dnorm(y[n], x[n,]%*%beta[,i], sigma[i])*lambda[i] / sum(dnorm(y[n], x[n,]%*%beta, sigma)*lambda)
}

dposterior <- function(x=X, y=Y, beta, sigma, lambda, i , n, arbvar){
  
  K <- length(lambda)
  d <- ncol(x)
  Npar <- K+d*K+arbvar*K+(1-arbvar)*1
  
  deriv <- numeric(Npar)
  
  for(j in (1:K)){ # diff wrt beta_j
    term1 <- c(dnorm(y[n], x[n,]%*%beta[,i], sigma[i]))*lambda[i] * (1/(sigma[i]^2)*x[n,]*c(y[n]-x[n,]%*%beta[,i])) * sum(dnorm(y[n], x[n,]%*%beta, sigma)*lambda)
    term2 <- c(dnorm(y[n], x[n,]%*%beta[,i], sigma[i]))*lambda[i] * c(dnorm(y[n], x[n,]%*%beta[,j], sigma[j]))*lambda[j] * (1/(sigma[j]^2)*x[n,]*c(y[n]-x[n,]%*%beta[,j]))
    deriv[((j-1)*d+1):(j*d)] <- ((i==j)*term1 - term2) / sum(dnorm(y[n], x[n,]%*%beta, sigma)*lambda)^2
  }
  
  if(arbvar){
    for(j in 1:K){ # diff wrt sigma2_j
      term1 <- c(dnorm(y[n], x[n,]%*%beta[,i], sigma[i]))*lambda[i] * (-1/(2*sigma[i]^2) + 1/(2*sigma[i]^4)*c(y[n]-x[n,]%*%beta[,i])^2) * sum(dnorm(y[n], x[n,]%*%beta, sigma)*lambda)
      term2 <- c(dnorm(y[n], x[n,]%*%beta[,i], sigma[i]))*lambda[i] * c(dnorm(y[n], x[n,]%*%beta[,j], sigma[j]))*lambda[j] * (-1/(2*sigma[j]^2) + 1/(2*sigma[j]^4)*c(y[n]-x[n,]%*%beta[,j])^2)
      deriv[d*K+j] <- ((i==j)*term1 - term2) / sum(dnorm(y[n], x[n,]%*%beta, sigma)*lambda)^2
    }
  }
  
  if(!(arbvar)){ # diff wrt sigma2 (SAME)
    term1 <- c(dnorm(y[n], x[n,]%*%beta[,i], sigma[i]))*lambda[i] * (-1/(2*sigma[i]^2) + 1/(2*sigma[i]^4)*c(y[n]-x[n,]%*%beta[,i])^2) * sum(dnorm(y[n], x[n,]%*%beta, sigma)*lambda)
    term2 <- c(dnorm(y[n], x[n,]%*%beta[,i], sigma[i]))*lambda[i] * sum(dnorm(y[n], x[n,]%*%beta, sigma)*lambda * (-1/(2*sigma^2) + 1/(2*sigma^4)*c(y[n]-x[n,]%*%beta)^2))
    deriv[d*K+1] <- (term1 - term2) / sum(dnorm(y[n], x[n,]%*%beta, sigma)*lambda)^2
  }
  
  for(j in 1:K){ # diff wrt lambda
    term1 <- c(dnorm(y[n], x[n,]%*%beta[,i], sigma[i])) * sum(dnorm(y[n], x[n,]%*%beta, sigma)*lambda)
    term2 <- c(dnorm(y[n], x[n,]%*%beta[,i], sigma[i]))*lambda[i] * c(dnorm(y[n], x[n,]%*%beta[,j], sigma[j]))
    deriv[d*K+arbvar*K+(1-arbvar)*1+j] <- ((i==j)*term1 - term2) / sum(dnorm(y[n], x[n,]%*%beta, sigma)*lambda)^2
  }
  
  return(deriv)
}

loglik <- function(x, y, beta, sigma, lambda, i, n){
  -.5 * log(sigma^2) - 1/(2*sigma^2) * c(y[n] - x[n,]%*%beta[,i])^2 + log(lambda[i])
}

dloglik <- function(x, y, beta, sigma, lambda, i, n, arbvar){
  
  K <- length(lambda)
  d <- ncol(x)
  Npar <- K+d*K+arbvar*K+(1-arbvar)*1
  
  deriv <- numeric(Npar)
  
  deriv[((i-1)*d+1):(i*d)] <- (1/sigma[i]^2)*t(x[n,])*c(y[n]-x[n,]%*%beta[,i]) # wrt beta_i
  deriv[d*K+arbvar*i+(1-arbvar)*1] <- -1/(2*sigma[i]^2) + 1/(2*sigma[i]^4)*c(y[n]-x[n,]%*%beta[,i])^2 # wrt sigma2_i
  deriv[d*K+arbvar*K+(1-arbvar)*1+i] <- 1/lambda[i] # wrt lambda_i
  
  return(deriv)
}

ddloglik <- function(x, y, beta, sigma, lambda, i, n, arbvar){
  
  K <- length(lambda)
  d <- ncol(x)
  Npar <- K+d*K+arbvar*K+(1-arbvar)*1
  
  dderiv <- matrix(0, Npar, Npar)
  
  dderiv[((i-1)*d+1):(i*d),((i-1)*d+1):(i*d)] <- -1/(sigma[i]^2) * t(x[n,,drop=F])%*%x[n,,drop=F] # beta_i beta_i
  dderiv[d*K+arbvar*i+(1-arbvar)*1,   d*K+arbvar*i+(1-arbvar)*1] <- 1/(2*sigma[i]^4) - 1/(sigma[i]^6) * c(y[n]-x[n,]%*%beta[,i])^2 # sigma2_i sigma2_i
  dderiv[d*K+arbvar*K+(1-arbvar)*1+i, d*K+arbvar*K+(1-arbvar)*1+i] <- -1/(lambda[i]^2) # lambda_i lambda_i
  dderiv[d*K+arbvar*i+(1-arbvar)*1, ((i-1)*d+1):(i*d)] <- dderiv[((i-1)*d+1):(i*d), d*K+arbvar*i+(1-arbvar)*1] <- -1/(sigma[i]^4)*t(x[n,])*c(y[n]-x[n,]%*%beta[,i])  # sigma2_i beta_i & beta_i sigma2_i
  
  return(dderiv)
}

fisher2covmat <- function(fisher, jacobian, d, K, arbvar, model, intercept){
  
  fisher <- (fisher + t(fisher))/2
  e <- eigen(fisher)
  v <- e$values
  v <- makePositive(v)
  U <- e$vectors
  
  covmat <- U%*%diag(1/v)%*%t(U)
  covmat <- t(jacobian) %*% covmat %*% jacobian
  covmat <- insert.entries(covmat, d, K, arbvar, model)
  
  if(model=="IID") order <- c(sapply(1:K, function(i) c((i-1)*d+1:d, d*K+i, (d+1)*K+i)))
  if(model=="HMM") order <- c(sapply(1:K, function(i) c((i-1)*d+1:d, d*K+i, (d+1)*K+(1+(i-1)*(K-1)):(i*(K-1)))))
  
  covmat <- covmat[order,]
  covmat <- covmat[,order]
  
  pars <- c(rep("beta", d-intercept), "sigma", rep("gamma", 1+(K-2)*(model=="HMM")))
  
  if(intercept) pars <- c("intercept", pars)
  pars <- rep(pars, K)
  rownames(covmat) <- colnames(covmat) <- pars
  
  return(covmat)
}

makePositive <- function(v,silent=TRUE){
  w <- which(v < 10^(-14))
  if(length(w)>0 & !silent)  warning("Some eigenvalues are below 10^(-14) and will therefore be regularized")
  for(ww in w){
    if(v[ww] < 0){
      v[ww] <- v[ww] - 2 * v[ww] + 10^(-14)
    } else {
      v[ww] <- v[ww] + 10^(-14)
    }
  }
  return(v)
}


insert.entries <- function(covmat, d, K, arbvar, model){
  
  if(!arbvar){
    newcovmat <- matrix(0, nrow(covmat)+K-1, nrow(covmat)+K-1)
    indexGamma <- setdiff(1:nrow(covmat),1:(d*K+1))
    indexGammaNew <- setdiff(1:nrow(newcovmat),1:(d*K+K))
    nGamma <- length(indexGamma)
    
    newcovmat[1:(d*K),1:(d*K)] <- covmat[1:(d*K),1:(d*K)]
    newcovmat[1:(d*K),(d*K+1):(d*K+K)] <- covmat[1:(d*K),d*K+1]
    newcovmat[1:(d*K),indexGammaNew] <- covmat[1:(d*K),indexGamma]
    
    newcovmat[(d*K+1):(d*K+K),1:(d*K)] <- matrix(covmat[d*K+1,1:(d*K)],K,d*K,byrow=T)
    newcovmat[(d*K+1):(d*K+K),(d*K+1):(d*K+K)] <- covmat[d*K+1,d*K+1]
    newcovmat[(d*K+1):(d*K+K),indexGammaNew] <- matrix(covmat[d*K+1,indexGamma],K,nGamma,byrow=T)
    
    newcovmat[indexGammaNew,1:(d*K)] <- covmat[indexGamma,1:(d*K)]
    newcovmat[indexGammaNew,(d*K+1):(d*K+K)] <- covmat[indexGamma,d*K+1]
    newcovmat[indexGammaNew,indexGammaNew] <- covmat[indexGamma,indexGamma]
    
    covmat <- newcovmat
  }
  if(model == "IID"){ # Create artificial extra entry for K'th mixing component
    covmat <- rbind(covmat, covmat[nrow(covmat),])
    covmat <- cbind(covmat, covmat[,ncol(covmat)])
  }
  
  return(covmat)
}


Q <- function(x=rep(c(1,2,.1,.1,.5,.5),2), xx=dat$x, yy=dat$y, K=2){
  
  d <- ncol(xx)
  
  x1 <- x[1:(length(x)/2)]
  x2 <- x[(length(x)/2+1):length(x)]
  
  beta <- matrix(x1[1:(d*K)], d,K)
  beta.t <- matrix(x2[1:(d*K)], d,K)
  
  sigma2 <- x1[(d*K+1)]
  sigma2.t <- x2[(d*K+1)]
  
  lambda <- x1[(d*K+1+1):(d*K+1+K)]
  lambda.t <- x2[(d*K+1+1):(d*K+1+K)]
  
  tmp <- 0
  for(n in 1:N){
    for(i in 1:K){
      summand <- loglik(xx, yy, beta, sqrt(sigma2), lambda, i, n) * posterior(xx, yy, beta.t, sqrt(sigma2.t), lambda.t, i, n)
      tmp <- tmp + summand
    }
  }
  
  return(tmp)
}

