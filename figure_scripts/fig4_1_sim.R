##################################################
## Figure generating scripts for 
## "Switching Regression Models and Causal Inference
## In the Presence of Latent Variables"
## Submitted to JMLR
## Author: Rune Christiansen
## email: krunechristiansen@math.ku.dk
##################################################
## simulating data for Figure 4 (left)
##################################################
library(pracma)
library(icph)

#################################################### simulating data ####################################################


## function that calculates mean posterior probabilites for the hidden states
mean.posterior <- function(sY, mu, beta, muX, sX, lambda){
  
  l <- c(lambda, 1-lambda)
  
  qX <- qnorm(c(.001, .999), mean=muX, sd=sX)
  qY <- c(min(qnorm(.001, mean=c(mu+beta*qX[1], mu+beta*qX[2]), sd=sY)), 
          max(qnorm(.999, mean=c(mu+beta*qX[1], mu+beta*qX[2]), sd=sY)))
  
  pi <- function(x,y,i) dnorm(y,mean=mu[i]+beta[i]*x,sd=sY)*dnorm(x,mean=muX,sd=sX)*l[i]
  p <- function(x,y) pi(x,y,1)+pi(x,y,2)
  
  MP <- matrix(NA,2,2)
  
  for(i in 1:2){
    for(j in 1:2){
      
      fun <- function(x,y) (pi(x,y,i) / p(x,y)) * (pi(x,y,j) / l[j])

      MP[i,j] <- integral2(fun, xmin=qX[1], xmax=qX[2], ymin=qY[1], ymax=qY[2])$Q
    }
  }
  diag(MP)
}

#########################
## generate data
#########################
if(onServer){
n <- 100
N <- 50
Nvec <- 1:N
num.exp <- 1000
pars <- rep(c("beta", "gamma", "intercept", "sigma"), each = 2)
l <- 8
out <- NULL
s <- 0

out <- tryCatch({read.table("~/server_scripts/cis_mean_posterior.txt", header = TRUE)}, error=function(e) NULL)
s <- ifelse(is.null(out), 0, max(out$s))
if(!is.null(out)) Nvec <- setdiff(Nvec, unique(out$sim.para))

for(j in Nvec){
  
  set.seed((j-1)*1000+1)
  
  sY <- runif(1,.1,.5)
  mu <- rep(runif(1,-1,1),2)
  beta <- runif(2,-1,1)
  muX <- runif(1,-1,1)
  sX <- runif(1,.1,1)
  lambda <- runif(1,.3,.7)
  
  mp <- mean.posterior(sY, mu, beta, muX, sX, lambda)
  
  if(0){
    x <- muX + sX*rnorm(n)
    h <- sample(1:2, size = n, prob = c(lambda, 1-lambda), replace = TRUE)
    y <- numeric(n)
    y[h==1] <- mu[1] + beta[1]*x[h==1] + sY*rnorm(n)[h==1]
    y[h==2] <- mu[2] + beta[2]*x[h==2] + sY*rnorm(n)[h==2]
    plot(x,y,col=h)
    
    mp
  }
  
  
  theta0 <- list(list(mu = mu, beta = matrix(beta,1,2), sigma = rep(sY,2), gamma = matrix(c(lambda, 1-lambda),2,2)))
  
    for(i in 1:num.exp){
      set.seed((j-1)*1000+i)
      s <- s+1
      print(paste0(s, " out of ", num.exp*N))
      
      x <- muX + sX*rnorm(n)
      h <- sample(1:2, size = n, prob = c(lambda, 1-lambda), replace = TRUE)
      y <- numeric(n)
      y[h==1] <- mu[1] + beta[1]*x[h==1] + sY*rnorm(n)[h==1]
      y[h==2] <- mu[2] + beta[2]*x[h==2] + sY*rnorm(n)[h==2]
      #plot(x,y,col=h)
      
      outt <- tryCatch({
        nlm.regmix <- test.equality.sr(Y=y, X=x, E=rep(1,n), model = "IID", plot = FALSE,
                                       par.test = list(method = "NLM", init = "regmix", theta0 = theta0, optim.n.inits=5))
        em.regmix <- em.true <- nlm.regmix
        nlm.true <- test.equality.sr(Y=y, X=x, E=rep(1,n), model = "IID", plot = FALSE,
                                       par.test = list(method = "NLM", init = "true", theta0 = theta0, optim.n.inits=1))
        # em.regmix <- test.equality.sr(Y=y, X=x, E=rep(1,n), model = "IID", plot = FALSE,
        #                                par.test = list(method = "EM", init = "regmix", theta0 = theta0, optim.n.inits=5))
        # em.true <- test.equality.sr(Y=y, X=x, E=rep(1,n), model = "IID", plot = FALSE,
        #                                par.test = list(method = "EM", init = "true", theta0 = theta0, optim.n.inits=1))
        
        outt <- rbind(data.frame(sY = sY, muX1 = mu[1], muX2 = mu[2], bX1 = beta[1], bX1 = beta[2], sX = sX, lambda = lambda),
                      data.frame(sY = sY, muX1 = mu[1], muX2 = mu[2], bX1 = beta[1], bX1 = beta[2], sX = sX, lambda = lambda),
                      data.frame(sY = sY, muX1 = mu[1], muX2 = mu[2], bX1 = beta[1], bX1 = beta[2], sX = sX, lambda = lambda),
                      data.frame(sY = sY, muX1 = mu[1], muX2 = mu[2], bX1 = beta[1], bX1 = beta[2], sX = sX, lambda = lambda))
        outt$mp1 <- mp[1]
        outt$mp2 <- mp[2]
        outt$p.value <- c(nlm.regmix$p.values.cover, nlm.true$p.values.cover, em.regmix$p.values.cover, em.true$p.values.cover)
        outt$loglik <- c(nlm.regmix$loglik[1,1], nlm.true$loglik[1,1], em.regmix$loglik[1,1], em.true$loglik[1,1])
        outt$method <- c("NLN", "NLM", "EM", "EM")
        outt$init <- c("regmix", "true", "regmix", "true")
        outt$sim.para <- j
        outt$sim.data <- i
        outt$s <- s
        outt
      }, error = function(e){
        outt <- rbind(data.frame(sY = sY, muX1 = mu[1], muX2 = mu[2], bX1 = beta[1], bX1 = beta[2], sX = sX, lambda = lambda),
                      data.frame(sY = sY, muX1 = mu[1], muX2 = mu[2], bX1 = beta[1], bX1 = beta[2], sX = sX, lambda = lambda),
                      data.frame(sY = sY, muX1 = mu[1], muX2 = mu[2], bX1 = beta[1], bX1 = beta[2], sX = sX, lambda = lambda),
                      data.frame(sY = sY, muX1 = mu[1], muX2 = mu[2], bX1 = beta[1], bX1 = beta[2], sX = sX, lambda = lambda))
        outt$mp1 <- mp[1]
        outt$mp2 <- mp[2]
        outt$p.value <- rep(NA,4)
        outt$loglik <- rep(NA,4)
        outt$method <- c("NLN", "NLM", "EM", "EM")
        outt$init <- c("regmix", "true", "regmix", "true")
        outt$sim.para <- j
        outt$sim.data <- i
        outt$s <- s
        outt
      })
      out <- rbind(out, outt)
    }
  write.table(out, "fig4_1_data.txt", quote = FALSE)
}

}
