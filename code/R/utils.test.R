
fit.model <- function(x, y, h, beta){
  
  N <- length(y)
  K <- ncol(beta)
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

decoupledTest <- function(x, y, e, intercept, plot){
  
  n <- nrow(x)
  changes <- sapply(1:n, function(i) abs(e[i]-e[i+1]) != 0)
  grid <- c(0, which(changes)-1, n)
  d <- ncol(x)
  if(d==0) S <- numeric(0)
  if(d>0)  S <- 1:d
  
  if(0){
    
    if(d==0){
      p <- fitplot.intercept.nomix(y, e)
    }else if(d==1){
      p <- fitplot.1pred.nomix(x, y, e, intercept)
    }else{
      p <- fitplot.nomix(x, y, e, intercept)
    }
    #print(p)
  }
  
  res <- seqICP::seqICP.s(Y=y, X=x, S=S, test="decoupled",
                          par.test = list(grid=grid),
                          model="iid")
  
  return(res = res)
}


output <- function(MFs){
  
  ne <- length(MFs)
  K <- MFs[[1]]$K
  d <- MFs[[1]]$d
  intercept <- MFs[[1]]$intercept
  Npars <- length(unlist(MFs[[1]]$theta[[1]]))
  # dd <- ifelse(K==2, d+2, d+3)
  
  estimates <- data.frame(parameter = c(sapply(1:ne, function(i) unlist(MFs[[i]]$parameter))),
                          component = rep(rep(1:K, each = Npars),ne),
                          environment = rep(1:ne, each = Npars*K),
                          estimate = c(sapply(1:ne, function(i) unlist(MFs[[i]]$theta))),
                          std.err = c(sapply(1:ne, function(i) sqrt(diag(MFs[[i]]$covmat)))))
  
  estimates <- estimates[order(estimates$parameter, estimates$component),]
  
  loglik <- matrix(c(sapply(1:ne, function(i) MFs[[i]]$loglik),
                     sapply(1:ne, function(i) MFs[[i]]$loglik0)), 2, ne, byrow=T)
  rownames(loglik) <- c("estimate", "true")
  
  return(list(estimates = estimates, loglik = loglik))
}
