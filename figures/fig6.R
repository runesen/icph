setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
#install_github("runesen/icph")
library(icph)
library(ggplot2)
library(reshape2)
library(gridExtra)
theme_set(theme_bw() + theme(axis.text = element_text(size = 10),
                             axis.title = element_text(size = 10),
                             plot.title = element_text(size = 10, hjust=.5)))


#################################################### simulating data ####################################################

num.exp <- 100
n.vec <- (1:5)*100
alpha <- .0499
dist.vec <- seq(0,2,.5)
S.list <- list(numeric(0), 1, 2, 3, 1:2, c(1,3), 2:3, 1:3)
set.frame <- var.frame <- NULL

s <- 0
for(k in 1:length(n.vec)){
  n <- n.vec[k]
  for(j in 1:length(dist.vec)){
    dist <- dist.vec[j]
    for(i in 1:num.exp){
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
        result.icp <- icph(D = cbind(X1, X2, X3, Y), E = E, target = 4, model = "IID")
        
        set.framee <- result.icp$pvalues[,-1]
        set.framee$Shat <- (set.framee$S == toString(result.icp$parent.set))
        set.framee$dist <- dist
        set.framee$seed <- s
        set.framee$n <- n
        
        var.framee <- data.frame(dist = dist, 
                                 seed = s, 
                                 n = n,
                                 var = 1:3)
        
        var.framee$inShat <- var.framee$var %in% result.icp$parent.set
        if(sum(set.framee$p.value>alpha)==0){
          var.framee$p.value.nonCausal <- NA
        } else{
          var.framee$p.value.nonCausal <- sapply(var.framee$var,
                                                 function(v){
                                                   ord <- order(set.framee$p.value)
                                                   w <- min(which(sapply(1:length(ord), function(i) !(v %in% S.list[[rev(ord)[i]]]))))
                                                   set.framee$p.value[rev(ord)[w]]
                                                 })
        }
        var.framee$p.value.bestSet <- sapply(var.framee$var,
                                             function(v){
                                               ord <- order(set.framee$p.value)
                                               w <- min(which(sapply(1:length(ord), function(i) (v %in% S.list[[rev(ord)[i]]]))))
                                               set.framee$p.value[rev(ord)[w]]
                                             })
        list(set.framee, var.framee)
      }, error = function(e){
        set.framee <- data.frame(p.value=NA, reject = NA, S=sapply(S.list, toString), Shat=NA, dist = dist, seed=s, n = n)
        var.framee <- data.frame(dist = dist, seed=s, n = n, var=1:3, inShat=NA, p.value.nonCausal=NA, p.value.bestSet=NA)
        list(set.framee, var.framee)
      })    
      set.frame <- rbind(set.frame, tmp[[1]])
      var.frame <- rbind(var.frame, tmp[[2]])
    }
  }
}


#################################################### producing figure ####################################################

set.frame$S <- factor(set.frame$S, levels = c("1, 2", "2", "1", "", "3", "1, 3", "2, 3", "1, 2, 3"))
set.frame$inSstar <- rep(sapply(unique(set.frame$seed), function(s) subset(set.frame, seed == s & Shat)$S %in% c("", "1", "2", "1, 2")), 
                         each = length(unique(set.frame$S)))

typeIS <- data.frame(n = rep(unique(set.frame$n), each = length(unique(set.frame$dist))),
                     dist=rep(unique(set.frame$dist),length(unique(set.frame$n))),
                     nocover = c(sapply(unique(set.frame$n), function(m){
                       c(sapply(unique(set.frame$dist), function(d){
                         mean((!subset(set.frame, S == "1, 2" & dist == d & n == m)$inSstar), na.rm = TRUE)
                       }))
                     }))
)

S_trans <- my_trans("S", low=0, up=.1, mid = .1, seq(0,round(max(typeIS$nocover,na.rm=T),2), .01))

pS <- ggplot(typeIS, aes(xmin=n-50, xmax = n+50, ymin=dist-.25, ymax = dist+.25)) + 
  geom_rect(aes(fill=nocover)) + xlab("sample size") + ylab("difference in regression coefficients") + 
  ggtitle("False causal discovery rate") + 
  geom_text(aes(x=n, y=dist, label=round(nocover,2))) + 
  scale_fill_gradient2(name = expression(paste("Rejection rate for ", H[0])), trans = "S", low = "green", high = "red", mid= "yellow", midpoint = .05) + 
  theme(legend.position = "none", plot.margin = unit(c(.2,.2,.2,1), "cm"))


typeIH <- data.frame(n = rep(unique(set.frame$n), each = length(unique(set.frame$dist))*length(unique(set.frame$S))),
                     dist=rep(rep(unique(set.frame$dist), each = length(unique(set.frame$S))),length(unique(set.frame$n))),
                     S = factor(rep(unique(set.frame$S), length(unique(set.frame$dist))*length(unique(set.frame$n)))),
                     rej = c(sapply(unique(set.frame$n), function(m){
                       c(sapply(unique(set.frame$dist), function(d){
                         c(sapply(unique(set.frame$S), function(s){
                           mean((subset(set.frame, S == s & dist == d & n == m)$p.value<.05), na.rm = TRUE)
                         }))
                       }))
                     }))
)


maxrej <- max(subset(typeIH, S=="1, 2")$rej, na.rm=T)
H_trans <- my_trans("H", low=0, up=maxrej, mid = .1, breaks=seq(0,round(maxrej,2), .1))

pH <- ggplot(subset(typeIH, S == "1, 2"), aes(xmin=n-50, xmax = n+50, ymin=dist-.25, ymax = dist+.25)) + 
  geom_rect(aes(fill=rej)) + xlab("sample size") + ylab("difference in regression coefficients") + 
  geom_text(aes(x=n, y=dist, label=round(rej,2))) + 
  ggtitle(expression(paste("Rejection rates for (true) null Hypotheis ", H[paste("0,",S^"*")]))) + 
  scale_fill_gradient2(name = expression(paste("Rejection rate for ", H[0])), trans = "H", low = "green", high = "red", mid= "yellow", midpoint = maxrej/2) + 
  theme(legend.position = "none", plot.margin = unit(c(.2,1,.2,.2), "cm"))

pdf("fig6.pdf", width = 8.5, height = 4)
print(grid.arrange(pH, pS, ncol=2))
dev.off()

