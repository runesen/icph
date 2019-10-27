# source("fisher.info.R")
# Takes model fit and true parameter as argument and returns p-value for hypothesis H: theta=theta_0
check.cover <- function(MF=MFs[[1]], pars=testpars){
  
  K <- MF$K
  perms <- gtools::permutations(K,K)
  alphaseq <- seq(.01,1,.01)
  # intersect <- TRUE
  cover <- TRUE
  i <- 0
  
  while(cover == TRUE){
    i <- i+1
    alpha <- alphaseq[i]
    
    tmp.cover <- logical(factorial(K))
    tmp.mse <- numeric(factorial(K))
    for(k in 1:factorial(K)){
      MF     <- permute.labels(MF, perms[k,])
      tmp <- check.cover.alpha(MF, alpha, pars)
      tmp.cover[k] <- tmp$cover
      tmp.mse[k] <- tmp$mse
    }
    cover <- any(tmp.cover)
  }
  
  p.value <- alpha
  mse <- min(tmp.mse)
  return(list(p.value = p.value, mse = mse))
}

# Checks if CR at level alpha covers the true parameter
check.cover.alpha <- function(MF, alpha, pars){
  
  # if(alpha==1) return(FALSE)
  # if(alpha==0) return(TRUE)
  
  arbvar <- MF$arbvar
  
  w <- which(MF$parameter %in% pars)
  
  if(!arbvar & "sigma" %in% pars){
    ww <- which(MF$parameter == "sigma")
    w <- w[!(w %in% ww[2:length(ww)])]
  }
  
  c0 <- unlist(MF$theta0)[w]
  c  <- unlist(MF$theta)[w]
  C  <- MF$covmat[w,w]
  
  m <- length(c)
  # q <- sqrt(qchisq(1-alpha, df = m))
  q <- qchisq(1-alpha, df = m)
  
  e <- eigen(C)
  d <- e$values
  d <- makePositive(d,silent=TRUE)
  d <- sqrt(d)
  U <- e$vectors
  
  # c0 <- 1/q * diag(1/d) %*% t(U) %*%(c-c0)
  # return(sum(c0^2)<1)
  
  c0    <- diag(1/d) %*% t(U) %*%(c-c0)
  mse   <- t(c0)%*%c0
  cover <- mse < q
  
  return(list(mse = mse, cover = cover))
}
