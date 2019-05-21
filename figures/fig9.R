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
n <- 500
K.vec <- 2:5
alpha <- .0499
s2 <- .1

S.list <- list(numeric(0), 1, 2, 3, 1:2, c(1,3), 2:3, 1:3)

set.frame <- var.frame <- NULL

s <- 0

for(j in 1:length(K.vec)){
  K <- K.vec[j]
  for(i in 1:num.exp){
    s <- s+1
    set.seed(s)
    
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
    betaY <- cbind(runif(2,1,2), 
                   runif(2,c(1,-2),c(2,-1)),
                   runif(2,-2,-1),
                   runif(2,c(-2,1),c(-1,2)),
                   c(0,0))[,1:K] 
    
    lambda <- matrix(runif(n = 3*K, min = 0.1, max = 1/K), 3, K)
    lambda <- lambda / rowSums(lambda)
    
    E <- sample(c(1,2,3), size = n, replace = TRUE)
    E <- sort(E)
    
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
    
    tmp <- tryCatch({
      result.icp <- icph(D = cbind(X1, X2, X3, Y), E = E, target = 4, 
                         model = "IID", number.of.states = K,
                         par.icp = list(output.pvalues = FALSE))
      
      set.framee <- result.icp$pvalues[,-1]
      set.framee$Shat <- (set.framee$S == toString(result.icp$parent.set))
      set.framee$K <- K
      set.framee$seed <- s
      set.framee$sigmasq <- s2
      
      var.framee <- data.frame(K = K, 
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
      set.framee <- data.frame(p.value=NA, reject = NA, S=sapply(S.list, toString),  Shat=NA, K=K, seed=s, sigmasq = s2)
      var.framee <- data.frame(K=K, sigmasq = s2, seed=s, var=1:3, inShat=NA, reject.nonCausal=NA)
      list(set.framee, var.framee)
    })    
    set.frame <- rbind(set.frame, tmp[[1]])
    var.frame <- rbind(var.frame, tmp[[2]])
  }
}

###### producing plots ######

Slevels <- c("1, 2", "2", "1", "", "3", "1, 3", "2, 3", "1, 2, 3")
set.frame$S <- factor(set.frame$S, levels = Slevels)
K.vec <- unique(set.frame$K)

dat <- data.frame(Shat = factor(rep(Slevels, length(K.vec)), levels = Slevels),
                  K = rep(K.vec, each = length(Slevels)),
                  sim = c(sapply(K.vec, function(k){
                    sapply(Slevels, function(s) sum(subset(set.frame, Shat & K == k)$S == s))
                  })),
                  percent = c(sapply(K.vec, function(k){
                    sapply(Slevels, function(s) 100*sum(subset(set.frame, Shat & K == k)$S == s) / 
                             nrow(subset(set.frame, Shat & K == k)))
                  })))


p1 <- ggplot(dat, aes(x=K, y = sim, fill = Shat)) + geom_bar(stat = "identity") + 
  xlab("number of latent states") + ylab(expression(paste("number of simulations resulting in ", hat(S), " = ..."))) + 
  geom_hline(yintercept = 5, linetype = 2, col = "red") + 
  scale_fill_manual(name = "S",
                    limits = c("", "1", "2", "1, 2", "3", "1, 3", "2, 3", "1, 2, 3"),
                    labels = c("{}  ", "{1}  ","{2}  ","{1,2}  ","{3}  ","{1,3}  ", "{2,3}  ", "{1,2,3}"),
                    values = c("#999999", "#56B4E9", "#0072B2", "#009E73", "#D55E00", "#F0E442", "#E69F00", "#CC79A7")) + 
  theme(plot.margin = unit(c(5,10,5,5), "mm"), legend.position = "bottom") + 
  guides(fill = guide_legend(nrow = 1))

Slevels <- c("", "1", "2", "3", "1, 2", "1, 3", "2, 3", "1, 2, 3")
rej.frame <- data.frame(S = factor(rep(Slevels, length(K.vec)), levels = Slevels),
                        K = rep(K.vec, each = length(Slevels)),
                        rej = c(sapply(K.vec, function(k){
                          sapply(Slevels, function(s) mean(subset(set.frame, K == k & S == s)$reject, na.rm=T))
                        })))

p2 <- ggplot(rej.frame, aes(x=K, y = rej, group = S, col = S)) + geom_point(size = 3, alpha = .7) + geom_line(size = 1, alpha = .7) + 
  ylab(expression(paste("rejecetion rates for ", H["0,S"]))) + xlab("number of latent states") + ylim(c(0,1)) + 
  geom_hline(yintercept = .05, linetype = 2, col = "red") + 
  scale_color_manual(name = "S",
                     limits = c("", "1", "2", "1, 2", "3", "1, 3", "2, 3", "1, 2, 3"),
                     labels = c("{}  ", "{1}  ","{2}  ","{1,2}  ","{3}  ","{1,3}  ", "{2,3}  ", "{1,2,3}"),
                     values = c("#999999", "#56B4E9", "#0072B2", "#009E73", "#D55E00", "#F0E442", "#E69F00", "#CC79A7")) + 
  theme(plot.margin = unit(c(5,5,5,10), "mm"))


p <- icph:::grid_arrange_shared_legend(p1, p2)

pdf("fig9.pdf", width=8.5, height=4)
grid.draw(p)
dev.off()

