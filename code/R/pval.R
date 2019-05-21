# Intersecting confidence regions and calculating pvalues
# source("fisher.info.R")

intersect.fun <- function(MF1, MF2, testpars, alpha=l/ne, arbvar){
  
  w <- which(MF1$parameter %in% testpars)
  
  if(!arbvar & "sigma" %in% testpars){
    ww <- which(MF1$parameter == "sigma")
    w <- w[!(w %in% ww[2:length(ww)])]
  }
  
  mle1 <- unlist(MF1$theta)[w]
  mle2 <- unlist(MF2$theta)[w]
  
  covmat1 <- MF1$covmat[w, w]
  covmat2 <- MF2$covmat[w, w]
  
  int <- intersect.ellipses(mle1, mle2, covmat1, covmat2, alpha)
  
  return(int)
}

# checks if two ellipsoid confidence regions with centers ci, covariance matrices Ci and coverage 1-alpha intersect
## BRAINSTORM:  Would it be possible to express distance in terms of alpha (appears only in c)? 
##              Then we would not have to call this function over a grid of alpha values..
intersect.ellipses <- function(c1=mle1, c2=mle2, C1=covmat1, C2=covmat2, alpha){
  
  if(alpha==1) return(FALSE)
  if(alpha==0) return(TRUE)
  
  m <- length(c1)
  q <- sqrt(qchisq(1-alpha, df = m))
  
  if(m==1){
    intersect <- c(abs(c1-c2) < (sqrt(C1)+sqrt(C2))*q)
  }else{
    
    #print(C1)
    e1 <- eigen(C1)
    d1 <- e1$values
    # d1 <- pmax(d1,min(1e-14, min(d1[d1>0])))
    # d1 <- sqrt(e1$values)
    d1 <- makePositive(d1)
    d1 <- sqrt(d1)
    U1 <- e1$vectors
    
    c <- 1/q * diag(1/d1) %*% t(U1) %*%(c2-c1)
    C <- diag(1/d1) %*% t(U1) %*% C2 %*% U1 %*% diag(1/d1)
    
    #print(C)
    e <- eigen(C)
    U <- e$vectors
    d <- e$values
    d <- makePositive(d,silent=TRUE)
    d <- sqrt(d)
    #print(e$values)
    #d <- sqrt(e$values)
    
    y <- -t(U)%*%c
    y <- abs(y) # 0 expressed in coordinate system of l and rotated to the first quadrant
    
    if(sum((y/d)^2) <= 1){ # y inside the ellipse
      intersect <- TRUE
    }else{ # Newton-Rhapson iterations
      
      # f goes monotone, quadratically to -1, so sure and fast convergence
      f <- function(t) sum((d*y/(t+d^2))^2)-1
      df <- function(t) -2*sum((y*d)^2/(t+d^2)^3)
      
      t0 <- 0
      ft0 <- f(t0)
      
      while(ft0>1e-4){
        t0 <- t0 - f(t0)/df(t0)
        ft0 <- f(t0)
      }
      
      x0 <- y*d^2/(d^2+t0) # projection of y onto (c,C)
      dist <- sqrt(sum((y-x0)^2))
      intersect <- (dist < 1)
    }
  }
  return(intersect)
}

# em is EM object, perm is a permutation of the vector (1,...,K)
permute.labels <- function(MF=MFs[[1]], perm){
  
  Npar <- length(unlist(MF$theta[[1]]))
  arbvar <- MF$arbvar
  theta.new <- lapply(perm, function(i) MF$theta[[i]])
  MF$theta <- theta.new
  MF$h <- perm[MF$h]
  new.order <- c(sapply(1:MF$K, function(i) (perm[i]-1)*Npar + 1:Npar))
  MF$covmat <- MF$covmat[new.order,]
  MF$covmat <- MF$covmat[,new.order]
  
  return(MF)
}


pval.fun <- function(MFs, pars, level, output.pvalue){
  
  ne <- length(MFs)
  K <- MFs[[1]]$K
  d <- MFs[[1]]$d
  N <- MFs[[1]]$N
  arbvar <- MFs[[1]]$arbvar
  
  pval <- 1
  reject <- FALSE
  perm <- lapply(1:ne, function(i) list(1:K))
  
  # If only interecpt is tested for invariance, the problem reduces to testing the equality of (unconditional) normals. 
  if(d==1 & MFs[[1]]$intercept){
    pval <- numeric(ne)
    for(i in 1:ne){
      y <- MFs[[i]]$y
      yy <- unlist(sapply((1:ne)[-i], function(j) MFs[[j]]$y))
      pval[i] <- ks.test(y, yy, alternative = "two.sided")$p.value
    }
    pval <- min(1, ne*min(pval))
    reject <- (pval < level)
    if(!output.pvalue) pval <- ifelse(reject, paste("<", level), paste(">=", level))
  }
  else{
    # If !output.pvalue, only check whether confidence regions intersect at level "level"
    if(output.pvalue){
      lseq <- 10^(seq(-4,0,length=100))
    } else {
      lseq <- level
    }
    count <- 0
    
    ## Loop over coverage of confidence regions
    for(l in lseq){
      
      ## All permutations of estimate in first environment
      comb <- lapply(1:factorial(K), function(i) list(gtools::permutations(K,K)[i,])) 
      comb <- list(list(1:K))
      keepIndex <- 1
      ## Loop over number of environments to compare
      for(nCompare in 2:ne){
        # keep only those configurations for which the first nCompare-1 environments overlap
        # If none overlap, we are done
        if(length(keepIndex) == 0) break
        comb <- lapply(keepIndex, function(i) comb[[i]])
        keepIndex <- c()
        
        tmp  <- lapply(1:factorial(K), function(i) gtools::permutations(K,K)[i,])
        ## For each overlapping configuration, add each permutation for the nCompare'st
        comb <- unlist(lapply(1:length(comb), function(i){ 
          lapply(1:length(tmp), function(j){ 
            lst <- comb[[i]]
            lst[[length(lst)+1]] <- tmp[[j]] 
            return(lst)
          }
          )}), recursive = FALSE)
        ## Loop over "active" combinations of the first nCompare environments
        for(k in 1:length(comb)){
          #if(k > length(comb)) break ## entries of comb are successively removed
          perms <- comb[[k]]
          all.intersect <- TRUE
          
          # count <- count+1
          # print(count)
          ## Check for pairwise intersection
          for(i in 1:(nCompare-1)){
            for(j in (i+1):nCompare){
              MFi <- permute.labels(MFs[[i]], unlist(perms[[i]]))
              MFj <- permute.labels(MFs[[j]], unlist(perms[[j]]))
              all.intersect <- all.intersect*intersect.fun(MFi, MFj, pars, l / ne, arbvar)
              if(!all.intersect) break
            }
            if(!all.intersect) break
          }
          ## If the first nCompare environments intersect, keep this combination of permutations
          if(all.intersect) keepIndex <- c(keepIndex, k)
          ## When having found a combination for which all environments intersect, break and increase \ell
          if(nCompare==ne & all.intersect){perm  <- perms; break}
        }
      }
      ## If no permutation lead to an intersection, set p-value to \ell and break. Otherwise increase \ell.
      if(!all.intersect){pval <- l; reject <- (pval < level); break}
    }
    ## If no intersection at this level, pval has been set to level and we reject.
    if(!output.pvalue){reject <- pval == level; pval <- ifelse(reject, paste("<", level), paste(">=", level))}
  }
  
  return(list(pval = pval, perm = perm, reject = reject))
}

## Old version
if(0){
pval.fun <- function(MFs, pars, level, output.pvalue){
  
  ne <- length(MFs)
  K <- MFs[[1]]$K
  d <- MFs[[1]]$d
  N <- MFs[[1]]$N
  arbvar <- MFs[[1]]$arbvar
  
  pval <- 1
  reject <- FALSE
  perm <- lapply(1:ne, function(i) list(1:K))
  
  # If only interecpt is tested for invariance, the problem reduces to testing the equality of (unconditional) normals. 
  if(d==1 & MFs[[1]]$intercept){
    pval <- numeric(ne)
    for(i in 1:ne){
      y <- MFs[[i]]$y
      yy <- unlist(sapply((1:ne)[-i], function(j) MFs[[j]]$y))
      pval[i] <- ks.test(y, yy, alternative = "two.sided")$p.value
    }
    pval <- min(1, ne*min(pval))
    reject <- (pval < level)
    if(!output.pvalue) pval <- ifelse(reject, paste("<", level), paste(">=", level))
  }
  else{
    combinations  <- expand.grid(lapply(1:ne, function(k){
      lapply(1:factorial(K), function(i) gtools::permutations(K,K)[i,])
    }))
    
    # If !output.pvalue, only check whether confidence regions intersect at level "level"
    if(output.pvalue){
      lseq <- seq(.01,1,.01)
    } else {
      lseq <- level
    }
    count<-0
    ## Loop over coverage of confidence regions
    for(l in lseq){
      
      for(k in 1:factorial(K)^ne){
        perms <- combinations[k,]
        all.intersect <- TRUE
        
        # count<-count+1
        # print(count)
        for(i in 1:(ne-1)){
          for(j in (i+1):ne){
            MFi <- permute.labels(MFs[[i]], unlist(perms[[i]]))
            MFj <- permute.labels(MFs[[j]], unlist(perms[[j]]))
            all.intersect <- all.intersect*intersect.fun(MFi, MFj, pars, l / ne, arbvar)
            if(!all.intersect) break
          }
          if(!all.intersect) break
        }
        if(all.intersect) break
      }
      if(!all.intersect){pval <- l; reject <- (pval < level); break}
    }
    if(!output.pvalue){reject <- pval == level; pval <- ifelse(reject, paste("<", level), paste(">=", level))}
  }
  return(list(pval = pval, perm = perm, reject = reject))
}
}