##################################################
## Figure generating scripts for 
## "Switching Regression Models and Causal Inference
## In the Presence of Latent Variables"
## Submitted to JMLR
## Author: Rune Christiansen
## email: krunechristiansen@math.ku.dk
##################################################
## simulating data for Figure 14 (left)
##################################################
onServer <- TRUE
library(pracma)
library(icph)

num.exp <- 100
n <- 500
K.vec <- 2:5
alpha <- .0499
s2 <- .1
S.list <- list(numeric(0), 1, 2, 3, 1:2, c(1,3), 2:3, 1:3)


set.frame <- var.frame <- NULL

s <- 0
for(j in 1:length(K.vec)){
  K <- K.vec[j]
  for(l in c("fixed2", "varying", "known")){
    for(i in 1:num.exp){
      s <- s+1
      set.seed(s)
      print(paste(s, "out of", num.exp*2*length(K.vec)))
      if(l == "fixed2") ell <- 2 # misspecified
      if(l == "varying") ell <- 2:5 # trying different values for K
      if(l == "known") ell <- K # trying different values for K
      
      mu1 <- runif(n = 1, min = -0.2, max = 0.2)
      mu2 <- runif(n = 1, min = -0.2, max = 0.2)
      mu22 <- runif(n = 1, min = .5, max = 1)
      mu3 <- runif(n = 1, min = -0.2, max = 0.2)
      mu32 <- runif(n = 1, min = -1, max = -0.5)
      muY <- runif(n = K, min = -0.2, max = 0.2)
      
      sigma1 <- runif(n = 1, min = 0.1, max = 0.3) 
      sigma2 <- runif(n = 1, min = 0.1, max = 0.3) 
      sigma22 <- runif(n = 1, min = 1, max = 1.5) 
      sigma3 <- runif(n = 1, min = 0.1, max = 0.3) 
      sigma32 <- sigma3
      
      beta2 <- sample(c(-1,1),1) * runif(1, .5, 1.5)
      beta22 <- beta2
      beta3 <- sample(c(-1,1),1) * runif(1, .5, 1.5)
      beta32 <- 0
      betaY <- cbind(runif(2,.5,1.5), 
                     runif(2,c(1.5,-2),c(2.5,-1)),
                     runif(2,-1.5,-.5),
                     runif(2,c(-2.5,1.5),c(-1.5,2.5)),
                     c(0,0))[,1:K] 
      
      E <- sample(c(1,2,3), size = n, replace = TRUE)
      E <- sort(E)
      
      lambda <- matrix(runif(n = 3*K, min = 0.1, max = 1/K), 3, K)
      lambda <- lambda / rowSums(lambda)
      
      X1 <-  sqrt(sigma1) * rnorm(n)
      
      X2 <- numeric(n)
      X2[E==1] <- mu2 +   beta2*X1[E==1] +  sqrt(sigma2) *  rnorm(n)[E==1]
      X2[E==2] <- mu22 +  beta22*X1[E==2] + sqrt(sigma22) * rnorm(n)[E==2]
      X2[E==3] <- mu2 +   beta2*X1[E==3] +  sqrt(sigma32) * rnorm(n)[E==3]
      
      H <- numeric(n)
      H[E==1] <- sample(1:K, size = n, prob = lambda[1,], replace = TRUE)[E==1]
      H[E==2] <- sample(1:K, size = n, prob = lambda[2,], replace = TRUE)[E==2]
      H[E==3] <- sample(1:K, size = n, prob = lambda[3,], replace = TRUE)[E==3]
      
      Y <- numeric(n)
      for(k in 1:K){
        Y[H==k] <- muY[k] + cbind(X1,X2)[H==k,] %*% betaY[,k,drop=F] + sqrt(s2) * rnorm(sum(H==k))
      }
      
      X3 <- numeric(n)
      X3[E==1] <- mu3 +  beta3 *  Y[E==1] + sqrt(sigma3) *  rnorm(n)[E==1]
      X3[E==2] <- mu3 +  beta3 *  Y[E==2] + sqrt(sigma3) *  rnorm(n)[E==2]
      X3[E==3] <- mu32 + beta32 * Y[E==3] + sqrt(sigma32) * rnorm(n)[E==3]
      
      if(0){
        plot(X2,Y,col=H)
        test.equality.sr(Y,X2,E,model="IID",plot=TRUE,number.of.states = 2)$p.value
        test.equality.sr(Y,X2,E,model="IID",plot=TRUE,number.of.states = 2:4)$p.value
      }
      
      
      tmp <- tryCatch({
        result.icp <- icph(D = cbind(X1, X2, X3, Y), E = E, target = 4, 
                           model = "IID", number.of.states = ell,
                           par.icp = list(output.pvalues = TRUE, silent = TRUE, stopIfEmpty=FALSE))
        
        set.framee <- result.icp$pvalues[,-1]
        set.framee$Shat <- (set.framee$S == toString(result.icp$parent.set))
        set.framee$K <- K
        set.framee$l <- l
        set.framee$seed <- s
        set.framee$sigmasq <- s2
        
        var.framee <- data.frame(K = K, 
                                 l = l, 
                                 sigmasq = s2,
                                 seed = s, 
                                 var = 1:3)
        
        var.framee$inShat <- var.framee$var %in% result.icp$parent.set
        if(sum(!set.framee$reject)==0){
          var.framee$reject.nonCausal <- NA
        } else{
          var.framee$reject.nonCausal <- sapply(var.framee$var,
                                                function(v){
                                                  w <- which(sapply(1:nrow(set.framee), function(i) !(v %in% S.list[[i]])))
                                                  sum(!set.framee$reject[w])==0
                                                })
        }
        list(set.framee, var.framee)
      }, error = function(e){
        set.framee <- data.frame(p.value=NA, reject = NA, S=sapply(S.list, toString),  Shat=NA, K=K, l=l, seed=s, sigmasq = s2)
        var.framee <- data.frame(K=K, l=l, sigmasq = s2, seed=s, var=1:3, inShat=NA, reject.nonCausal=NA)
        list(set.framee, var.framee)
      })    
      set.frame <- rbind(set.frame, tmp[[1]])
      var.frame <- rbind(var.frame, tmp[[2]])
      write.table(set.frame, "fig14_1_data_set.txt", quote = FALSE, sep = "\t")
      write.table(var.frame, "fig14_1_data_var.txt", quote = FALSE)
    }
  }
}
