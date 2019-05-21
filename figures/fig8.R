setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
#devtools::install_github("runesen/icph/code")
library(icph)
library(InvariantCausalPrediction)
library(ggplot2)
library(reshape2)
library(gridExtra)
library(grid)
library(rmutil)
library(pcalg)
theme_set(theme_bw() + theme(axis.text = element_text(size = 11),
                             axis.title = element_text(size = 11)))

num.exp <- 100
# num.exp <- 100
n.vec  <- 100*(1:5)
# n.vec  <- 100
alpha <- .0499
dist <- 1.5
S.list <- list(numeric(0), 1, 2, 3, 1:2, c(1,3), 2:3, 1:3)
par.test <- list(method="NLM", variance.constraint="equality")
arbvar <- FALSE

set.frame <- var.frame <- NULL

# ICPH, no model violation
s <- 0
for(k in 1:length(n.vec)){
  n <- n.vec[k]
  for(i in 1:num.exp){
    print(paste("ICPH, none, n = ", n, ", sim = ", i))
    s <- s+1
    set.seed(s)
    
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
    
    tmp <- tryCatch({
      result.icp <- icph(D = cbind(X1, X2, X3, Y), E = E, target = 4, model = "IID",
                         par.icp = list(output.pvalues = FALSE, stopIfEmpty = FALSE, silent=TRUE),
                         par.test = par.test)
      parent.set <- result.icp$parent.set
      
      var.framee <- data.frame(n = n,
                               seed = s, 
                               var = 1:3,
                               method = "ICPH",
                               violation = "none")
      var.framee$reject.nonCausal <- sapply(var.framee$var, function(x) x %in% parent.set)
      var.framee
      
    }, error = function(e){
      data.frame(n=n, seed=s, var=1:3, method="ICPH", violation="none", reject.nonCausal=NA)
    })
    var.frame <- rbind(var.frame, tmp)
  }
}

var.frame.none <- var.frame


set.frame <- var.frame <- NULL

# PC algorithm
s <- 0
for(k in 1:length(n.vec)){
  n <- n.vec[k]
  for(i in 1:num.exp){
    print(paste("PC, none, n = ", n, ", sim = ", i))
    print(s)
    s <- s+1
    set.seed(s)
    
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
    
    M <- cbind(X1,X2,X3,E,Y)
    
    tmp <- tryCatch({
      pc.fit <- pc(suffStat = list(C = cor(M), n = n),
                      indepTest = gaussCItest, 
                      alpha=0.05, labels = colnames(M), 
                      verbose = FALSE)
      # dev.off()
      # plot(pc.fit)
      amat <- as(pc.fit, "amat")
      w <- which(amat["Y",]==1) # edge to Y
      parent.set <- w[amat[w,"Y"]==0] # no edge from Y
      
      var.framee <- data.frame(n = n,
                               seed = s, 
                               var = 1:3,
                               method = "PC",
                               violation = "none")
      var.framee$reject.nonCausal <- sapply(var.framee$var, function(x) x %in% parent.set)
      var.framee
      
    }, error = function(e){
      data.frame(n=n, seed=s, var=1:3, method="PC", violation="none", reject.nonCausal=NA)
    })    
    var.frame <- rbind(var.frame, tmp)
  }
}

var.frame.pc <- var.frame


set.frame <- var.frame <- NULL

# k-means ICP algorithm
s <- 0
for(k in 1:length(n.vec)){
  n <- n.vec[k]
  for(i in 1:num.exp){
    print(paste("ICP, none, n = ", n, ", sim = ", i))
    print(s)
    s <- s+1
    set.seed(s)
    
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
    
    M <- cbind(Y,X1,X2,X3)
    
    tmp <- tryCatch({
      
      km <- kmeans(M, centers=2)
      hhat <- km$cluster
      
      icp.fit1 <- ICP(X=cbind(X1,X2,X3)[hhat==1,], Y=Y[hhat==1], ExpInd = E[hhat==1], alpha=.025)
      icp.fit2 <- ICP(X=cbind(X1,X2,X3)[hhat==2,], Y=Y[hhat==2], ExpInd = E[hhat==2], alpha=.025)
      
      parent.set <- unique(c(which(icp.fit1$pvalues<.025), which(icp.fit2$pvalues<.025)))
      
      var.framee <- data.frame(n = n,
                               seed = s, 
                               var = 1:3,
                               method = "k-means ICP",
                               violation = "none")
      var.framee$reject.nonCausal <- sapply(var.framee$var, function(x) x %in% parent.set)
      var.framee
      
    }, error = function(e){
      data.frame(n=n, seed=s, var=1:3, method="k-means ICP", violation="none", reject.nonCausal=NA)
    })    
    var.frame <- rbind(var.frame, tmp)
  }
}

var.frame.icp <- var.frame


set.frame <- var.frame <- NULL

# different error variances
s <- 0
for(k in 1:length(n.vec)){
  n <- n.vec[k]
  for(i in 1:num.exp){
    print(paste("ICPH, variance heterogeneity, n = ", n, ", sim = ", i))
    s <- s+1
    set.seed(s)
    
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
    sigmaY <- runif(n = 2, min = 0.1, max = 0.3) # different error variances
    
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
    
    tmp <- tryCatch({
      result.icp <- icph(D = cbind(X1, X2, X3, Y), E = E, target = 4, model = "IID",
                         par.icp = list(output.pvalues = FALSE, stopIfEmpty = FALSE, silent=TRUE),
                         par.test = par.test)
      parent.set <- result.icp$parent.set
      
      var.framee <- data.frame(n = n,
                               seed = s, 
                               var = 1:3,
                               method = "ICPH",
                               violation = "heterogeneous variances")
      var.framee$reject.nonCausal <- sapply(var.framee$var, function(x) x %in% parent.set)
      var.framee
      
    }, error = function(e){
      data.frame(n=n, seed=s, var=1:3, method="ICPH", violation = "heterogeneous variances", reject.nonCausal=NA)
    })
    var.frame <- rbind(var.frame, tmp)
  }
}

set.frame.different_variances <- set.frame
var.frame.different_variances <- var.frame




set.frame <- var.frame <- NULL

# uniform noise
s <- 0
for(k in 1:length(n.vec)){
  n <- n.vec[k]
  for(i in 1:num.exp){
    print(paste("ICPH, uniform noise, n = ", n, ", sim = ", i))
    s <- s+1
    set.seed(s)
    
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
    sigmaY <- rep(runif(n = 1, min = 0.1, max = 0.3),2)
    
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
    
    Y <- numeric(n) # uniform noise with same variance as the corresponding Gaussian noise
    Y[H==1] <- muY[1] + cbind(X1,X2)[H==1,] %*% betaY[,1,drop=F] + sqrt(sigmaY)[1] * runif(sum(H==1),-sqrt(3), sqrt(3))
    Y[H==2] <- muY[2] + cbind(X1,X2)[H==2,] %*% betaY[,2,drop=F] + sqrt(sigmaY)[2] * runif(sum(H==2),-sqrt(3), sqrt(3))
    
    X3 <- numeric(n)
    X3[E==1] <- mu3 +  beta3 *  Y[E==1] + sqrt(sigma3) *  rnorm(n)[E==1]
    X3[E==2] <- mu3 +  beta3 *  Y[E==2] + sqrt(sigma3) *  rnorm(n)[E==2]
    X3[E==3] <- mu32 + beta32 * Y[E==3] + sqrt(sigma32) * rnorm(n)[E==3]
    
    tmp <- tryCatch({
      result.icp <- icph(D = cbind(X1, X2, X3, Y), E = E, target = 4, model = "IID",
                         par.icp = list(output.pvalues = FALSE, stopIfEmpty = FALSE, silent=TRUE),
                         par.test = par.test)
      parent.set <- result.icp$parent.set
      
      var.framee <- data.frame(n = n,
                               seed = s, 
                               var = 1:3,
                               method = "ICPH",
                               violation = "uniform noise")
      var.framee$reject.nonCausal <- sapply(var.framee$var, function(x) x %in% parent.set)
      var.framee
      
    }, error = function(e){
      data.frame(n=n, seed=s, var=1:3, method="ICPH", violation = "uniform noise", reject.nonCausal=NA)
    })
    var.frame <- rbind(var.frame, tmp)
  }
}

set.frame.uniform_noise <- set.frame
var.frame.uniform_noise <- var.frame




set.frame <- var.frame <- NULL

# laplace noise
s <- 0
for(k in 1:length(n.vec)){
  n <- n.vec[k]
  for(i in 1:num.exp){
    print(paste("ICPH, laplace noise, n = ", n, ", sim = ", i))
    s <- s+1
    set.seed(s)
    
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
    sigmaY <- rep(runif(n = 1, min = 0.1, max = 0.3),2)
    
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
    
    Y <- numeric(n) # uniform noise with same variance as the corresponding Gaussian noise
    Y[H==1] <- muY[1] + cbind(X1,X2)[H==1,] %*% betaY[,1,drop=F] + sqrt(sigmaY)[1] * rlaplace(sum(H==1),1/sqrt(2))
    Y[H==2] <- muY[2] + cbind(X1,X2)[H==2,] %*% betaY[,2,drop=F] + sqrt(sigmaY)[2] * rlaplace(sum(H==2),1/sqrt(2))
    
    X3 <- numeric(n)
    X3[E==1] <- mu3 +  beta3 *  Y[E==1] + sqrt(sigma3) *  rnorm(n)[E==1]
    X3[E==2] <- mu3 +  beta3 *  Y[E==2] + sqrt(sigma3) *  rnorm(n)[E==2]
    X3[E==3] <- mu32 + beta32 * Y[E==3] + sqrt(sigma32) * rnorm(n)[E==3]
    
    tmp <- tryCatch({
      result.icp <- icph(D = cbind(X1, X2, X3, Y), E = E, target = 4, model = "IID",
                         par.icp = list(output.pvalues = FALSE, stopIfEmpty = FALSE, silent=TRUE),
                         par.test = par.test)
      parent.set <- result.icp$parent.set
      
      var.framee <- data.frame(n = n,
                               seed = s, 
                               var = 1:3,
                               method = "ICPH",
                               violation = "Laplace noise")
      var.framee$reject.nonCausal <- sapply(var.framee$var, function(x) x %in% parent.set)
      var.framee
      
    }, error = function(e){
      data.frame(n=n, seed=s, var=1:3, method="ICPH", violation = "Laplace noise", reject.nonCausal=NA)
    })
    var.frame <- rbind(var.frame, tmp)
  }
}

set.frame.laplace_noise <- set.frame
var.frame.laplace_noise <- var.frame


set.frame <- var.frame <- NULL

# H -> X mean shift
s <- 0
for(k in 1:length(n.vec)){
  n <- n.vec[k]
  for(i in 1:num.exp){
    print(paste("ICPH, mean shift, n = ", n, ", sim = ", i))
    s <- s+1
    set.seed(s)
    
    mu1 <- runif(n = 2, min = -1, max = 1)
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
    sigmaY <- rep(runif(n = 1, min = 0.1, max = 0.3),2)
    
    beta1 <- sample(c(-1,1),1) * runif(1, .5, 1.5)
    beta2 <- sample(c(-1,1),1) * runif(1, .5, 1.5)
    beta22 <- beta2
    beta3 <- sample(c(-1,1),1) * runif(1, .5, 1.5)
    beta32 <- 0
    betaY <- matrix(sample(c(-1,1),4,replace=T)*runif(4, .5, 1.5), 2, 2)
    betaY[,2] <- betaY[,1] + sign(betaY[,1])*dist # Controlling the distance between coefficients of different states
    
    lambda <- runif(n = 3, min = 0.3, max = 0.7)
    
    E <- sample(c(1,2,3), size = n, replace = TRUE)
    E <- sort(E)
    
    H <- numeric(n)
    H[E==1] <- sample(c(1,2), size = n, prob = c(lambda[1], 1-lambda[1]), replace = TRUE)[E==1]
    H[E==2] <- sample(c(1,2), size = n, prob = c(lambda[2], 1-lambda[2]), replace = TRUE)[E==2]
    H[E==3] <- sample(c(1,2), size = n, prob = c(lambda[3], 1-lambda[3]), replace = TRUE)[E==3]
    
    X1 <- numeric(n)
    X1[H==1] <- mu1[1] + sqrt(sigma1) * rnorm(sum(H==1))
    X1[H==2] <- mu1[2] + sqrt(sigma1) * rnorm(sum(H==2))
    
    X2 <- numeric(n)
    X2[E==1] <- mu2 +   beta2*X1[E==1] +  sqrt(sigma2) *  rnorm(n)[E==1]
    X2[E==2] <- mu22 +  beta22*X1[E==2] + sqrt(sigma22) * rnorm(n)[E==2]
    X2[E==3] <- mu2 +   beta2*X1[E==3] +  sqrt(sigma32) * rnorm(n)[E==3]
    
    Y <- numeric(n) # uniform noise with same variance as the corresponding Gaussian noise
    Y[H==1] <- muY[1] + cbind(X1,X2)[H==1,] %*% betaY[,1,drop=F] + sqrt(sigmaY)[1] * rnorm(sum(H==1))
    Y[H==2] <- muY[2] + cbind(X1,X2)[H==2,] %*% betaY[,2,drop=F] + sqrt(sigmaY)[2] * rnorm(sum(H==2))
    
    X3 <- numeric(n)
    X3[E==1] <- mu3 +  beta3 *  Y[E==1] + sqrt(sigma3) *  rnorm(n)[E==1]
    X3[E==2] <- mu3 +  beta3 *  Y[E==2] + sqrt(sigma3) *  rnorm(n)[E==2]
    X3[E==3] <- mu32 + beta32 * Y[E==3] + sqrt(sigma32) * rnorm(n)[E==3]
    
    tmp <- tryCatch({
      result.icp <- icph(D = cbind(X1, X2, X3, Y), E = E, target = 4, model = "IID",
                         par.icp = list(output.pvalues = FALSE, stopIfEmpty = FALSE, silent=TRUE),
                         par.test = par.test)
      parent.set <- result.icp$parent.set
      
      var.framee <- data.frame(n = n,
                               seed = s, 
                               var = 1:3,
                               method = "ICPH",
                               violation = "mean shift")
      var.framee$reject.nonCausal <- sapply(var.framee$var, function(x) x %in% parent.set)
      var.framee
      
    }, error = function(e){
      data.frame(n=n, seed=s, var=1:3, method="ICPH", violation = "mean shift", reject.nonCausal=NA)
    })
    var.frame <- rbind(var.frame, tmp)
  }
}

set.frame.mean_shift <- set.frame
var.frame.mean_shift <- var.frame



set.frame <- var.frame <- NULL

# H -> X variance scaling
s <- 0
for(k in 1:length(n.vec)){
  n <- n.vec[k]
  for(i in 1:num.exp){
    print(paste("ICPH, variance shift, n = ", n, ", sim = ", i))
    s <- s+1
    set.seed(s)
    
    mu1 <- runif(n = 2, min = -1, max = 1)
    mu2 <- runif(n = 1, min = -0.2, max = 0.2)
    mu22 <- runif(n = 1, min = .5, max = 1)
    mu3 <- runif(n = 1, min = -0.2, max = 0.2)
    mu32 <- runif(n = 1, min = -1, max = -0.5)
    muY <- runif(n = 2, min = -0.2, max = 0.2)
    
    sigma1 <- runif(n = 2, min = 0.1, max = 1) 
    sigma2 <- runif(n = 1, min = 0.1, max = 0.3) 
    sigma22 <- runif(n = 1, min = 1, max = 1.5) 
    sigma3 <- runif(n = 1, min = 0.1, max = 0.3) 
    sigma32 <- sigma3
    sigmaY <- rep(runif(n = 1, min = 0.1, max = 0.3),2)
    
    beta1 <- sample(c(-1,1),1) * runif(1, .5, 1.5)
    beta2 <- sample(c(-1,1),1) * runif(1, .5, 1.5)
    beta22 <- beta2
    beta3 <- sample(c(-1,1),1) * runif(1, .5, 1.5)
    beta32 <- 0
    betaY <- matrix(sample(c(-1,1),4,replace=T)*runif(4, .5, 1.5), 2, 2)
    betaY[,2] <- betaY[,1] + sign(betaY[,1])*dist # Controlling the distance between coefficients of different states
    
    lambda <- runif(n = 3, min = 0.3, max = 0.7)
    
    E <- sample(c(1,2,3), size = n, replace = TRUE)
    E <- sort(E)
    
    H <- numeric(n)
    H[E==1] <- sample(c(1,2), size = n, prob = c(lambda[1], 1-lambda[1]), replace = TRUE)[E==1]
    H[E==2] <- sample(c(1,2), size = n, prob = c(lambda[2], 1-lambda[2]), replace = TRUE)[E==2]
    H[E==3] <- sample(c(1,2), size = n, prob = c(lambda[3], 1-lambda[3]), replace = TRUE)[E==3]
    
    X1 <- numeric(n)
    X1[H==1] <- sqrt(sigma1[1]) * rnorm(sum(H==1))
    X1[H==2] <- sqrt(sigma1[2]) * rnorm(sum(H==2))
    
    X2 <- numeric(n)
    X2[E==1] <- mu2 +   beta2*X1[E==1] +  sqrt(sigma2) *  rnorm(n)[E==1]
    X2[E==2] <- mu22 +  beta22*X1[E==2] + sqrt(sigma22) * rnorm(n)[E==2]
    X2[E==3] <- mu2 +   beta2*X1[E==3] +  sqrt(sigma32) * rnorm(n)[E==3]
    
    Y <- numeric(n) # uniform noise with same variance as the corresponding Gaussian noise
    Y[H==1] <- muY[1] + cbind(X1,X2)[H==1,] %*% betaY[,1,drop=F] + sqrt(sigmaY)[1] * rnorm(sum(H==1))
    Y[H==2] <- muY[2] + cbind(X1,X2)[H==2,] %*% betaY[,2,drop=F] + sqrt(sigmaY)[2] * rnorm(sum(H==2))
    
    X3 <- numeric(n)
    X3[E==1] <- mu3 +  beta3 *  Y[E==1] + sqrt(sigma3) *  rnorm(n)[E==1]
    X3[E==2] <- mu3 +  beta3 *  Y[E==2] + sqrt(sigma3) *  rnorm(n)[E==2]
    X3[E==3] <- mu32 + beta32 * Y[E==3] + sqrt(sigma32) * rnorm(n)[E==3]
    
    tmp <- tryCatch({
      result.icp <- icph(D = cbind(X1, X2, X3, Y), E = E, target = 4, model = "IID",
                         par.icp = list(output.pvalues = FALSE, stopIfEmpty = FALSE),
                         par.test = par.test)
      parent.set <- result.icp$parent.set
      
      var.framee <- data.frame(n = n,
                               seed = s, 
                               var = 1:3,
                               method = "ICPH",
                               violation = "variance shift")
      var.framee$reject.nonCausal <- sapply(var.framee$var, function(x) x %in% parent.set)
      var.framee
      
    }, error = function(e){
      data.frame(n=n, seed=s, var=1:3, method="ICPH", violation = "variance shift", reject.nonCausal=NA)
    })
    var.frame <- rbind(var.frame, tmp)
  }
}

set.frame.variance_scaling <- set.frame
var.frame.variance_scaling <- var.frame

###############################################


set.frame <- var.frame <- NULL

# H continuous
s <- 0
for(k in 1:length(n.vec)){
  n <- n.vec[k]
  for(i in 1:num.exp){
    print(paste("ICPH, continuous H, n = ", n, ", sim = ", i))
    s <- s+1
    set.seed(s)
    
    mu1 <- runif(n = 1, min = -1, max = 1)
    mu2 <- runif(n = 1, min = -0.2, max = 0.2)
    mu22 <- runif(n = 1, min = .5, max = 1)
    mu3 <- runif(n = 1, min = -0.2, max = 0.2)
    mu32 <- runif(n = 1, min = -1, max = -0.5)
    muY <- runif(n = 1, min = -0.2, max = 0.2)
    
    sigma1 <- runif(n = 1, min = 0.1, max = 1) 
    sigma2 <- runif(n = 1, min = 0.1, max = 0.3) 
    sigma22 <- runif(n = 1, min = 1, max = 1.5) 
    sigma3 <- runif(n = 1, min = 0.1, max = 0.3) 
    sigma32 <- sigma3
    sigmaY <- rep(runif(n = 1, min = 0.1, max = 0.3),2)
    
    beta1 <- sample(c(-1,1),1) * runif(1, .5, 1.5)
    beta2 <- sample(c(-1,1),1) * runif(1, .5, 1.5)
    beta22 <- beta2
    beta3 <- sample(c(-1,1),1) * runif(1, .5, 1.5)
    beta32 <- 0
    betaY <- matrix(sample(c(-1,1),4,replace=T)*runif(4, .5, 1.5), 2, 1)
    
    E <- sample(c(1,2,3), size = n, replace = TRUE)
    E <- sort(E)
    
    H <- rnorm(n)
    
    X1 <- sqrt(sigma1) * rnorm(n)
    
    X2 <- numeric(n)
    X2[E==1] <- mu2 +   beta2*X1[E==1] +  sqrt(sigma2) *  rnorm(n)[E==1]
    X2[E==2] <- mu22 +  beta22*X1[E==2] + sqrt(sigma22) * rnorm(n)[E==2]
    X2[E==3] <- mu2 +   beta2*X1[E==3] +  sqrt(sigma32) * rnorm(n)[E==3]
    
    Y <- H*(muY + cbind(X1,X2) %*% betaY) + sqrt(sigmaY) * rnorm(n)
    
    X3 <- numeric(n)
    X3[E==1] <- mu3 +  beta3 *  Y[E==1] + sqrt(sigma3) *  rnorm(n)[E==1]
    X3[E==2] <- mu3 +  beta3 *  Y[E==2] + sqrt(sigma3) *  rnorm(n)[E==2]
    X3[E==3] <- mu32 + beta32 * Y[E==3] + sqrt(sigma32) * rnorm(n)[E==3]
    
    tmp <- tryCatch({
      result.icp <- icph(D = cbind(X1, X2, X3, Y), E = E, target = 4, model = "IID",
                         par.icp = list(output.pvalues = FALSE, stopIfEmpty = FALSE, silent=TRUE),
                         par.test = par.test)
      parent.set <- result.icp$parent.set
      
      var.framee <- data.frame(n = n,
                               seed = s, 
                               var = 1:3,
                               method = "ICPH",
                               violation = "continuous H")
      var.framee$reject.nonCausal <- sapply(var.framee$var, function(x) x %in% parent.set)
      var.framee
      
    }, error = function(e){
      data.frame(n=n, seed=s, var=1:3, method="ICPH", violation = "continuous H", reject.nonCausal=NA)
    })
    var.frame <- rbind(var.frame, tmp)
  }
}

set.frame.H_continuous <- set.frame
var.frame.H_continuous <- var.frame


################

var.frame <- rbind(var.frame.none,
                   var.frame.different_variances,
                   var.frame.uniform_noise,
                   var.frame.laplace_noise,
                   var.frame.mean_shift,
                   var.frame.variance_scaling,
                   var.frame.H_continuous,
                   var.frame.pc,
                   var.frame.icp)

write.table(var.frame, "var.frame.sensitivity.txt", quote = FALSE)

# set.frame <- read.table("set.frame.sensitivity.txt", header = TRUE)
# var.frame <- read.table("var.frame.sensitivity.txt", header = TRUE)

############################################ Rejection Rates for non-causality ############################################

rej.frame <- data.frame(n=rep(unique(var.frame$n), each=length(unique(var.frame$var))*length(unique(var.frame$violation))*length(unique(var.frame$method))),
                        var = rep(rep(unique(var.frame$var), each = length(unique(var.frame$violation))*length(unique(var.frame$method))), length(unique(var.frame$n))),
                        violation = rep(rep(unique(var.frame$violation), each=length(unique(var.frame$method))), length(unique(var.frame$var))*length(unique(var.frame$n))),
                        method = rep(unique(var.frame$method), length(unique(var.frame$n))*length(unique(var.frame$var))*length(unique(var.frame$violation))),
                        rej = c(sapply(unique(var.frame$n), function(m){
                          c(sapply(unique(var.frame$var), function(v){
                            c(sapply(unique(var.frame$violation), function(vio){
                              c(sapply(unique(var.frame$method), function(me){
                                mean((subset(var.frame, n == m & var == v & violation == vio & method == me)$reject.nonCausal), na.rm = TRUE)
                              }))
                            }))
                          }))
                        }))
)

rej.frame$violation <- factor(rej.frame$violation)
rej.frame$method <- factor(rej.frame$method)


## ggplot default colors
tmp <- data.frame(x=1:7, y=1:7, z=factor(1:7))
p <- ggplot(tmp, aes(x,y,col=z)) + geom_point()
cols <- unique(ggplot_build(p)$data[[1]][,1])[1:7]

rej.frame$violation <- factor(rej.frame$violation, 
                              levels = levels(rej.frame$violation), 
                              labels = c("none", "het. var.", "uniform noise", "Laplace noise", "mean shift", "var. shift", "cont. H"))

rej.frame$method <- factor(rej.frame$method,
                           levels = c("ICPH", "k-means ICP", "PC"),
                           labels = c("ICPH", "k-means ICP", "PC"))

p11 <- ggplot(subset(rej.frame, var == 1), aes(n, rej, group = interaction(method,violation), col = violation, alpha = (violation=="none"))) + 
  geom_point(aes(shape=method, size = (violation=="none"))) + 
  scale_alpha_manual(values=c(.5,1), guide = FALSE) + 
  scale_size_manual(values = c(3,3), guide = FALSE) + 
  geom_line(size=1, aes(lty=method)) + 
  scale_y_continuous(limits = c(0,1)) + 
  ylab("rejecetion rate for non-causality") + xlab("sample size") + ggtitle(expression(bold(X^1))) + 
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5),
        plot.margin = unit(c(2,8,2,8), "mm")) +
  scale_color_manual("model violation", values = c("black",cols)) +
  scale_linetype_manual(values = c("solid","longdash","twodash")) + 
  guides(size = "none", alpha = "none", lty = "none", 
         shape = guide_legend(order=1, nrow=1, override.aes = list(size = 3)), 
         color = guide_legend(order=2,nrow=1)) 

p22 <- ggplot(subset(rej.frame, var == 2), aes(n, rej, group = interaction(method,violation), col = violation, alpha = (violation=="none"))) + 
  geom_point(aes(shape=method, size = (violation=="none"))) + 
  scale_alpha_manual(values=c(.5,1), guide = FALSE) + 
  scale_size_manual(values = c(3,3), guide = FALSE) + 
  geom_line(size=1, aes(lty=method)) + 
  scale_y_continuous(limits = c(0,1)) + 
  ylab("rejecetion rate for non-causality") + xlab("sample size") + ggtitle(expression(bold(X^2))) + 
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5),
        plot.margin = unit(c(2,8,2,8), "mm")) + 
  scale_color_manual("model violation", values = c("black",cols), guide = guide_legend(nrow=1)) +
  scale_linetype_manual(values = c("solid","longdash","twodash")) + 
  guides(size = "none", alpha = "none", lty = "none", shape = guide_legend(order=1), color = guide_legend(order=2)) 

p33 <- ggplot(subset(rej.frame, var == 3), aes(n, rej, group = interaction(method,violation), col = violation, alpha = (violation=="none"))) + 
  geom_point(aes(shape=method, size = (violation=="none"))) + 
  scale_alpha_manual(values=c(.5,1), guide = FALSE) + 
  scale_size_manual(values = c(3,3), guide = FALSE) + 
  geom_line(size=1, aes(lty=method)) + 
  scale_y_continuous(limits = c(0,1)) + 
  ylab("rejecetion rate for non-causality") + xlab("sample size") + ggtitle(expression(bold(X^3))) + 
  geom_hline(yintercept = .05, col = "red", lty = 2) + 
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5),
        plot.margin = unit(c(2,8,2,8), "mm")) + 
  scale_color_manual("model violation", values = c("black",cols), guide = guide_legend(nrow=1)) +
  scale_linetype_manual(values = c("solid","longdash","twodash")) + 
  guides(size = "none", alpha = "none", lty = "none", shape = guide_legend(order=1), color = guide_legend(order=2)) 

pp <- icph:::grid_arrange_shared_legend(p11, p22, p33, ncol = 3, nrow = 1)

pdf("fig8.pdf", width=11.5, height=3.5)
grid.draw(pp)
dev.off()
