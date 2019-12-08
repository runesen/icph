##################################################
## Figure generating scripts for 
## "Switching Regression Models and Causal Inference
## In the Presence of Latent Variables"
## Submitted to JMLR
## Author: Rune Christiansen
## email: krunechristiansen@math.ku.dk
##################################################
## genetating data for Figure 9
##################################################

library(rmutil)
library(icph)
library(InvariantCausalPrediction)
source("jci-pc.R")


num.exp <- 100
n.vec  <- 100*(1:5)
alpha <- .0499
dist <- 1.5
S.list <- list(numeric(0), 1, 2, 3, 1:2, c(1,3), 2:3, 1:3)
par.test <- list(method="NLM", variance.constraint="equality")
arbvar <- FALSE

var.frame <- NULL


  
############################################## FCI-PC ##############################################
s <- 0
for(k in 1:length(n.vec)){
  n <- n.vec[k]
  for(i in 1:num.exp){
    print(paste("PC-FCI, none, n = ", n, ", sim = ", i))
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
      pc.fit <- jci_pc(suffStat = list(C = cor(M), n = n),
                   indepTest = gaussCItest,
                   alpha=0.05, labels = colnames(M),
                   verbose = FALSE,
                   source.nodes = "E")
      # dev.off()
      # plot(pc.fit)
      amat <- as(pc.fit, "amat")
      w <- which(amat["Y",]==1) # edge to Y
      parent.set <- w[amat[w,"Y"]==0] # no edge from Y

      var.framee <- data.frame(n = n,
                               seed = s,
                               var = 1:3,
                               method = "JCI-PC",
                               violation = "none")
      var.framee$reject.nonCausal <- sapply(var.framee$var, function(x) x %in% parent.set)
      var.framee

    }, error = function(e){
      data.frame(n=n, seed=s, var=1:3, method="JCI-PC", violation="none", reject.nonCausal=NA)
    })
    var.frame <- rbind(var.frame, tmp)
  }
}

write.table(var.frame, "fig9_data.txt", quote = FALSE, sep="\t")



############################################## k-means ICP ##############################################
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

write.table(var.frame, "fig9_data.txt", quote = FALSE, sep="\t")

############################################## ICPH, no model violation ##############################################
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

write.table(var.frame, "fig9_data.txt", quote = FALSE, sep="\t")

############################################## ICPH, variance heterogeneity ##############################################
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

write.table(var.frame, "fig9_data.txt", quote = FALSE, sep="\t")

############################################## ICPH, uniform noise ##############################################
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

write.table(var.frame, "fig9_data.txt", quote = FALSE, sep="\t")


############################################## ICPH, Laplace noise ##############################################
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

write.table(var.frame, "fig9_data.txt", quote = FALSE, sep="\t")


var.frame <- read.table("fig9_data.txt", header = TRUE, sep="\t")
############################################## ICPH, H -> X mean shift ##############################################
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

write.table(var.frame, "fig9_data.txt", quote = FALSE, sep="\t")

############################################## ICPH, H -> X variance scaling ##############################################
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
                         par.icp = list(output.pvalues = FALSE, stopIfEmpty = FALSE, silent=TRUE),
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

write.table(var.frame, "fig9_data.txt", quote = FALSE, sep="\t")

############################################## ICPH, H continuous ##############################################
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

write.table(var.frame, "fig9_data.txt", quote = FALSE, sep="\t")
