##################################################
## Figure generating scripts for 
## "Switching Regression Models and Causal Inference
## In the Presence of Latent Variables"
## Submitted to JMLR
## Author: Rune Christiansen
## email: krunechristiansen@math.ku.dk
##################################################
## simulating data for Figure 14 (middle)
##################################################

onServer <- TRUE
library(pracma)
library(glmnet)
library(icph)

num.exp <- 100
n  <- 300
alpha <- .0499
dist <- 1.5
S.list <- list(numeric(0), 1, 2, 3, 1:2, c(1,3), 2:3, 1:3)
par.test <- list(method="NLM", variance.constraint="equality")
nb.noisevar <- 10^c(-Inf, 0:3)
nb.select <- 5

set.frame <- var.frame <- NULL

############################### ICPH ###############################
s <- 0
for(k in 1:length(nb.noisevar)){
  nb <- nb.noisevar[k]
  for(i in 1:num.exp){
    print(paste(s+1, "out of", num.exp*length(nb.noisevar)))
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
    
    X <- cbind(X1,X2,X3)
    ## including noise variables as children of X3 (so that they are actually correlated with Y)
    if(nb > 0){
      b <- runif(nb,-1,1)
      Xnoise <- matrix(X3,nrow=n)%*%matrix(b,ncol=nb) + matrix(rnorm(nb*n), nrow=n, ncol=nb)
      colnames(Xnoise) <- sapply(1:nb, function(i) paste0("N",i))
      X <- cbind(X1,X2,X3,Xnoise)
    }
    
    
    mod.glm <- glmnet(x=X, y=Y, family = "gaussian")
    betamat <- t(mod.glm$beta)
    nz <- sort(apply(betamat, 2, function(x) min(which(x!=0))))
    vars <- names(nz)[1:5]
    vars <- vars[!is.na(vars)]
    Xsel <- X[,vars]

    tmp <- tryCatch({
      result.icp <- icph(D = cbind(Y, Xsel), E = E, target = 1, model = "IID",
                         par.icp = list(output.pvalues = FALSE, stopIfEmpty = FALSE, silent=TRUE),
                         par.test = list(method = "NLM", variance.constraint = "equality", init = "regmix"))
      parent.set <- vars[result.icp$parent.set-1]
      
      var.framee <- data.frame(n = n,
                               seed = s, 
                               var = 1:3,
                               method = "ICPH",
                               nb.noisevar = nb)
      var.framee$reject.nonCausal <- sapply(c("X1", "X2", "X3"), function(x) x %in% parent.set)
      var.framee
      
    }, error = function(e){
      data.frame(n=n, seed=s, var=1:3, method="ICPH", nb.noisevar = nb, reject.nonCausal=NA)
    })
    var.frame <- rbind(var.frame, tmp)
  }
  write.table(var.frame, "fig14_2_data.txt", quote = FALSE, sep="\t")
}


