##################################################
## Figure generating scripts for 
## "Switching Regression Models and Causal Inference
## In the Presence of Latent Variables"
## Submitted to JMLR
## Author: Rune Christiansen
## mail: krunechristiansen@math.ku.dk
##################################################
## simulating data for Figure 7
##################################################
onServer <- TRUE
library(pracma)
library(icph)

mean.posterior <- function(muX1, sX1, muX2, bX2, sX2, muY, bY, sY, lambda){
  l <- c(lambda, 1-lambda)
  
  # cube containing 100-.998^3 % of the data (which is used as integration domain)
  qX1 <- qnorm(c(.001, .999), mean=muX1, sd=sX1)
  qX2 <- c(min(qnorm(.001, mean=c(muX2+bX2*qX1[1], muX2+bX2*qX1[2]), sd=sX2)),
           max(qnorm(.999, mean=c(muX2+bX2*qX1[1], muX2+bX2*qX1[2]), sd=sX2)))
  qY <- c(min(qnorm(.001, mean=c(muY+bY%*%c(qX1[1],qX2[1]), muY+bY%*%c(qX1[1],qX2[2]),
                                 muY+bY%*%c(qX1[2],qX2[1]), muY+bY%*%c(qX1[2],qX2[2])), sd=sY)),
          max(qnorm(.999, mean=c(muY+bY%*%c(qX1[1],qX2[1]), muY+bY%*%c(qX1[1],qX2[2]),
                                 muY+bY%*%c(qX1[2],qX2[1]), muY+bY%*%c(qX1[2],qX2[2])), sd=sY)))
  
  # density over (X1,X2,Y) under regime i
  pi <- function(x1,x2,y,i) dnorm(y,mean=muY[i]+bY[1,i]*x1+bY[2,i]*x2,sd=sY)*dnorm(x2,mean=muX2+bX2*x1,sd=sX2)*dnorm(x1,mean=muX1,sd=sX1)*l[i]
  # density over (X1,X2,Y)
  p <- function(x1,x2,y) pi(x1,x2,y,1)+pi(x1,x2,y,2)
  # posterior probability of state i
  posti <- function(x1,x2,y,i) ifelse(pi(x1,x2,y,i)<1e-10,0,pi(x1,x2,y,i)/p(x1,x2,y))
  
  
  out <- 1
  for(i in 1:2){
    fun <- function(x1,x2,y) posti(x1,x2,y,i) * (pi(x1,x2,y,i) / l[i])
    out <- out*integral3(fun, xmin=qX1[1], xmax=qX1[2], ymin=qX2[1], ymax=qX2[2], zmin=qY[1], zmax=qY[2], reltol = 1e-2)
  }
  out
}

num.exp <- 100
n.vec <- (1:5)*100
alpha <- .0499
dist.vec <- seq(0,2,.5)
S.list <- list(numeric(0), 1, 2, 3, 1:2, c(1,3), 2:3, 1:3)
S.vec <- c("", "1", "2", "3", "1, 2", "1, 3", "2, 3", "1, 2, 3")
set.frame <- var.frame <- NULL

s <- 0
for(k in 1:length(n.vec)){
  n <- n.vec[k]
  for(j in 1:length(dist.vec)){
    dist <- dist.vec[j]
    for(i in 1:num.exp){
      s <- s+1
      set.seed(s)
      print(paste(s, "out of", num.exp*length(n.vec)*length(dist.vec)))
      
      mu1 <- runif(n = 1, min = -0.2, max = 0.2)
      mu2 <- runif(n = 1, min = -0.2, max = 0.2)
      mu22 <- runif(n = 1, min = .5, max = 1)
      mu3 <- runif(n = 1, min = -0.2, max = 0.2)
      mu32 <- runif(n = 1, min = -1, max = -0.5)
      muY <- runif(n = 2, min = -0.2, max = 0.2)
      
      sigma1 <- runif(n = 1, min = 0.1, max = 0.3) 
      sigma2 <- runif(n = 1, min = 0.1, max = 0.3) 
      sigma22 <- runif(n = 1, min = 1, max = 1.5) 
      sigma3 <- runif(n = 1, min = 0.1, max = 0.3) 
      sigma32 <- sigma3
      sigmaY <- rep(runif(n = 1, min = 0.1, max = 0.3),2) # equal error variances
      
      beta2 <- sample(c(-1,1),1) * runif(1, .5, 1.5)
      beta22 <- beta2
      beta3 <- sample(c(-1,1),1) * runif(1, .5, 1.5)
      beta32 <- 0
      betaY <- matrix(sample(c(-1,1),4,replace=T)*runif(4, .5, 1.5), 2, 2)
      betaY[,2] <- betaY[,1] + sign(betaY[,1])*dist # Controlling the distance between coefficients of different states
      
      lambda <- runif(n = 3, min = 0.3, max = 0.7)
      
      E <- sample(c(1,2,3), size = n, replace = TRUE)
      E <- sort(E)
      
      X1 <-  sqrt(sigma1) * rnorm(n)
      
      X2 <- numeric(n)
      X2[E==1] <- mu2 +   beta2*X1[E==1] +  sqrt(sigma2) *  rnorm(n)[E==1]
      X2[E==2] <- mu22 +  beta22*X1[E==2] + sqrt(sigma22) * rnorm(n)[E==2]
      X2[E==3] <- mu2 +   beta2*X1[E==3] +  sqrt(sigma32) * rnorm(n)[E==3]
      
      H <- numeric(n)
      H[E==1] <- sample(c(1,2), size = n, prob = c(lambda[1], 1-lambda[1]), replace = TRUE)[E==1]
      H[E==2] <- sample(c(1,2), size = n, prob = c(lambda[2], 1-lambda[2]), replace = TRUE)[E==2]
      H[E==3] <- sample(c(1,2), size = n, prob = c(lambda[3], 1-lambda[3]), replace = TRUE)[E==3]
      
      Y <- numeric(n)
      Y[H==1] <- muY[1] + cbind(X1,X2)[H==1,] %*% betaY[,1,drop=F] + sqrt(sigmaY)[1] * rnorm(sum(H==1))
      Y[H==2] <- muY[2] + cbind(X1,X2)[H==2,] %*% betaY[,2,drop=F] + sqrt(sigmaY)[2] * rnorm(sum(H==2))
      
      X3 <- numeric(n)
      X3[E==1] <- mu3 +  beta3 *  Y[E==1] + sqrt(sigma3) *  rnorm(n)[E==1]
      X3[E==2] <- mu3 +  beta3 *  Y[E==2] + sqrt(sigma3) *  rnorm(n)[E==2]
      X3[E==3] <- mu32 + beta32 * Y[E==3] + sqrt(sigma32) * rnorm(n)[E==3]
      
      piP <- tryCatch({
        if(n != 100 | i > 20) stop("") # for every dist, only calculate this for 20 different data sets (more is not needed)
        mean.posterior(muX1=0, sX1 = sqrt(sigma1), muX2 = mu2, bX2 = beta2,
                       sX2 = sqrt(sigma22), muY = muY, bY=betaY, sY=sqrt(sigmaY)[1], lambda = lambda[1])
      }, error=function(e) NA)
      piP
      
      tmp <- tryCatch({
        result.icp <- icph(D = cbind(X1, X2, X3, Y), E = E, target = 4, model = "IID", par.icp = list(silent = TRUE))
        
        set.framee <- data.frame(ind=1:8, p.value = NA, reject = NA, S = S.vec)
        w <- which(S.vec %in% result.icp$pvalues$S)
        set.framee$p.value[w] <- sapply(S.vec[w], function(S) result.icp$pvalues$p.value[S==result.icp$pvalues$S])
        set.framee$reject[w] <- set.framee$p.value[w]<alpha
        set.framee$Shat <- (set.framee$S == toString(result.icp$parent.set))
        set.framee$dist <- dist
        set.framee$seed <- s
        set.framee$n <- n
        set.framee$piP <- piP
        
        var.framee <- data.frame(dist = dist, 
                                 seed = s, 
                                 n = n,
                                 piP = piP,
                                 var = 1:3)
        
        var.framee$inShat <- var.framee$var %in% result.icp$parent.set
        if(all(set.framee$reject, na.rm=1)){
          var.framee$p.value.nonCausal <- NA
        } else{
          var.framee$p.value.nonCausal <- sapply(var.framee$var,
                                                 function(v){
                                                   ord <- order(set.framee$p.value)[!is.na(set.framee$p.value)]
                                                   # index of the set with largest pvalue among all sets not containing v
                                                   w <- min(which(sapply(1:length(ord), function(i) !(v %in% S.list[[rev(ord)[i]]]))))
                                                   # corresponding p-value
                                                   set.framee$p.value[rev(ord)[w]]
                                                 })
        }
        var.framee$p.value.bestSet <- sapply(var.framee$var,
                                             function(v){
                                               ord <- order(set.framee$p.value)[!is.na(set.framee$p.value)]
                                               w <- min(which(sapply(1:length(ord), function(i) (v %in% S.list[[rev(ord)[i]]]))))
                                               set.framee$p.value[rev(ord)[w]]
                                             })
        list(set.framee, var.framee)
      }, error = function(e){
        set.framee <- data.frame(p.value=NA, reject = NA, S=sapply(S.list, toString), Shat=NA, dist = dist, seed=s, n = n, piP = piP)
        var.framee <- data.frame(dist = dist, seed=s, n = n, piP= piP, var=1:3, inShat=NA, p.value.nonCausal=NA, p.value.bestSet=NA)
        list(set.framee, var.framee)
      })    
      set.frame <- rbind(set.frame, tmp[[1]])
      var.frame <- rbind(var.frame, tmp[[2]])
      write.table(set.frame, "fig7_data_set.txt", quote = FALSE, sep = "\t")
      write.table(var.frame, "fig7_data_var.txt", quote = FALSE)
    }
  }
}

