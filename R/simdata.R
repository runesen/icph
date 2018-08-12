# library(ggplot2)

## simdata
simdata <- function(n, ne, d, K,
                    mu, beta, sigma, gamma, 
                    shift, noise, atomic, rm.inactive = FALSE){
  
  active <- 1:(d+1)
  if(missing(shift)){
    shift <- matrix(0, nrow = ne, ncol = d+1)
  }
  if(missing(noise)){
    noise <- matrix(1, nrow = ne, ncol = d+1)
  }
  if(missing(atomic)){
    atomic <- matrix(NA, nrow = ne, ncol = d+1)
  }
  if(missing(gamma)){
    gamma <- lapply(1:ne, function(i) matrix(1/K,K,K))
  }
  if(rm.inactive){
  active <- which((apply(beta[[1]], 1, function(x) sum(x^2)) == 0)*
                    (apply(beta[[1]], 2, function(x) sum(x^2)) == 0) == 0)
  beta <- lapply(beta, function(betaa) betaa[active, active])
  }
  
  K <- length(beta)
  d <- length(active)-1
  
  e <- rep(1:ne, each = floor(n/ne))
  e[(ne*floor(n/ne)):n] <- ne
  
  h <- numeric(n)
  # if(!is.null(lambda)){
  #   for(i in 1:ne){
  #     h[e==i] <- sample(1:K, sum(e==i), prob = lambda[i,], replace = T)
  #   }
  # }
  #if(!is.null(gamma)){
    for(i in 1:ne){
      delta<-solve(t(diag(K)-gamma[[i]]+1),rep(1,K))
      w <- which(e==i)
      m <- length(w)
      h[w[1]] <- sample(1:K,size=1,prob=delta)
      for(j in 2:m){
        h[w[j]] <- sample(1:K,size=1,prob=gamma[[i]][h[w[j-1]],])
      }
    }
  #}
  
  dat <- matrix(0, nrow = n, ncol = d+1)
  done <- NULL
  for(i in 1:(d+1)){
    notdone <- which(!(1:(d+1) %in% done))
    upnext <- which(apply(beta[[1]][, notdone, drop = F], 1, function(x) sum(x^2)) == 0)
    upnext <- upnext[upnext %in% notdone]; if(length(upnext)==0) stop("No cycles allowed")
    j <- min(upnext)
    xj <- numeric(n)
    for(k in 1:K){
      for(l in 1:ne){
        dat[h==k & e==l, j] <- mu[k,j] + shift[l,j] + 
          dat[h==k & e==l, ]%*%beta[[k]][j,] + sqrt(noise[l,j])*sigma[k,j] * rnorm(sum(h==k & e==l))
        if(!is.na(atomic[l,j])) dat[h==k & e==l, j] <- atomic[l,j]
      }
    }
    done <- c(done, j)
  }
  
  out.list <- list(y = dat[,1], x = dat[,-1,drop=F], e = e, h = h)
  
  return(out.list)
}

summarY <- function(dat){
  
  y <- dat$y
  x <- dat$x
  e <- dat$e
  h <- dat$h
  
  ne <- length(unique(e))
  K <- length(unique(h))
  d <- ncol(x)
  
  dd <- data.frame(variable = rep(c("Y", sapply(1:d, function(j) paste0("X",j))), each = ne*K),
                  environment = rep(rep(unique(e), each = K), d+1),
                  state = rep(rep(1:K, length(unique(e))), d+1))
  mat <- cbind(y, x)
  m <- v <- c()
  for(j in 1:(d+1)){
    mj <- c(sapply(unique(e), function(f){
      sapply(1:K, function(k) mean(mat[e==f & h==k, j]))
    }))
    vj <- c(sapply(unique(e), function(f){
      sapply(1:K, function(k) var(mat[e==f & h == k, j]))
    }))
    m <- c(m, mj)
    v <- c(v, vj)
  }
  
  dd$mean <- m
  dd$variance <- v
  return(dd)
}

## plotdata
plotdata <- function(dat){
  
  frame <- data.frame(y=dat$y, x=dat$x, e=factor(dat$e), h=factor(dat$h))
  
  n <- nrow(frame)
  d <- ncol(frame)-3
  
  melt.frame <- reshape2::melt(frame, id.vars = c("y", "e", "h"))
  
  gglot2::ggplot(melt.frame, aes(value, y)) + geom_point(aes(col = h, shape = e), size = 2, alpha = 0.7) + 
    facet_grid(variable~., scales = "free") + theme(text = element_text(size = 20)) + 
    theme_bw()
}
