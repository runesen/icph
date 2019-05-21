setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
#install_github("runesen/icph")
library(icph)
library(ggplot2)
library(reshape2)
library(gridExtra)
library(grid)
theme_set(theme_bw() + theme(axis.text = element_text(size = 11),
                             axis.title = element_text(size = 11)))

num.exp <- 100
n.vec  <- 100*(1:5)
method.vec <- c("NLM", "EM")
constraint.vec <- c("lower bound", "equality")
alpha <- .0499
dist.vec <- 1.5
S.list <- list(numeric(0), 1, 2, 3, 1:2, c(1,3), 2:3, 1:3)
par.test <- list()

set.frame.lst <- var.frame.lst <- list()

ss <- 0

for(method in method.vec){
  for(constraint in constraint.vec){
    ss <- ss+1
    par.test$method <- method
    par.test$variance.constraint <- constraint
    arbvar <- (constraint == "lower bound")
    
    set.frame <- var.frame <- NULL
    s <- 0
    
    for(k in 1:length(n.vec)){
      n <- n.vec[k]
      for(j in 1:length(dist.vec)){
        dist <- dist.vec[j]
        for(i in 1:num.exp){
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
          
          X1 <-  mu1 + sqrt(sigma1) * rnorm(n)
          
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
                               par.icp = list(output.pvalues = FALSE, stopIfEmpty=FALSE, silent = TRUE),
                               par.test = par.test)
            
            set.framee <- result.icp$pvalues[,-1]
            set.framee$Shat <- (set.framee$S == toString(result.icp$parent.set))
            set.framee$dist <- dist
            set.framee$n <- n
            set.framee$seed <- s
            set.framee$arbvar <- arbvar
            set.framee$method <- method
            set.framee$testNoMix <- FALSE
            
            var.framee <- data.frame(dist = dist, 
                                     n = n,
                                     seed = s, 
                                     var = 1:3,
                                     arbvar = arbvar,
                                     method = method,
                                     testNoMix = FALSE)
            
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
            set.framee <- data.frame(p.value = NA, reject = NA, S=sapply(S.list, toString), 
                                     Shat=NA, dist=dist, n=n, seed=s, arbvar=arbvar, method=method, testNoMix=FALSE)
            var.framee <- data.frame(dist = dist, n=n, seed=s, var=1:3, 
                                     arbvar=arbvar, method=method, testNoMix=FALSE, inShat=NA, reject.nonCausal=NA)
            list(set.framee, var.framee)
          })    
          set.frame <- rbind(set.frame, tmp[[1]])
          var.frame <- rbind(var.frame, tmp[[2]])
        }
      }
    }
    set.frame.lst[[ss]] <- set.frame
    var.frame.lst[[ss]] <- var.frame
  }
}


set.frame.nlm.arbvar <- set.frame.lst[[1]]
var.frame.nlm.arbvar <- var.frame.lst[[1]]
set.frame.nlm.samevar <- set.frame.lst[[2]]
var.frame.nlm.samevar <- var.frame.lst[[2]]
set.frame.em.arbvar <- set.frame.lst[[3]]
var.frame.em.arbvar <- var.frame.lst[[3]]
set.frame.em.samevar <- set.frame.lst[[4]]
var.frame.em.samevar <- var.frame.lst[[4]]

set.frame <- rbind(set.frame.nlm.arbvar,
                   set.frame.nlm.samevar,
                   set.frame.em.arbvar,
                   set.frame.em.samevar)
var.frame <- rbind(var.frame.nlm.arbvar,
                   var.frame.nlm.samevar,
                   var.frame.em.arbvar,
                   var.frame.em.samevar)

############################################ Results for empirical estimator \hat S ############################################

set.frame$S <- factor(set.frame$S, levels = c("1, 2", "2", "1", "", "3", "1, 3", "2, 3", "1, 2, 3"))

set.frame$method <- factor(set.frame$method, levels = c("EM", "NLM"),
                           labels = c("method: EM", "method: NLM"))
set.frame$arbvar <- factor(set.frame$arbvar, levels = c(TRUE, FALSE),
                           labels = c("variance constraint: lower bound", "variance constraint: equality"))

p <- ggplot(subset(set.frame, Shat), aes(x=n, order = S)) + geom_bar(aes(fill=S), width = 90) + 
  facet_grid(.~method+arbvar) + xlab("sample size") + ylab(expression(paste("number of simulations resulting in ", hat(S), " = ..."))) + 
  geom_hline(yintercept = 5, linetype = 2, col = "red") + 
  scale_fill_manual(name = "S",
                    limits = c("", "1", "2", "1, 2", "3", "1, 3", "2, 3", "1, 2, 3"),
                    labels = c("{}  ", "{1}  ","{2}  ","{1,2}  ","{3}  ","{1,3}  ", "{2,3}  ", "{1,2,3}"),
                    values = c("#999999", "#56B4E9", "#0072B2", "#009E73", "#D55E00", "#F0E442", "#E69F00", "#CC79A7")) + 
  theme(legend.position = "bottom",
        strip.text = element_text(face="bold")) + 
  guides(fill = guide_legend(nrow=1))

p



pdf("fig7.pdf", width=12*1.1, height=4*1.1)
print(p)
dev.off()