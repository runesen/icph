setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

source("../code/R/icph.R")
source("../code/R/test.equality.sr.R")
source("../code/R/EM.R")
source("../code/R/NLM.R")
source("../code/R/pval.R")
source("../code/R/fisher.info.R")
source("../code/R/utils.test.R")
source("../code/R/plotfuns.R")
source("../code/R/simdata.R")
source("../code/R/check.coverage.R")




# library(icph)
library(ggplot2)
theme_set(theme_bw())

set.seed(1)

# sample size
n <- 600
# environment indicator
E <- rep(1:3, each = n/3)

## Transition probabilities for hidden Markov Chain
# Different in each environment
Gamma <- list(matrix(c(.8, .2,
                       .5, .5), 2,2,byrow=T),
              matrix(c(.5, .5,
                       .5, .5),2,2, byrow=T),
              matrix(c(.5, .5,
                       .2, .8),2,2, byrow=T)) 

# hidden Markov Chain
H <- rep(1, n)
for(t in 2:n){
  H[t] <- sample(1:2, size=1, prob = Gamma[[E[t]]][H[t-1],])
}

# Observed data -- here S* = {2}
X3 <- H * rnorm(n)
X2 <- 1 + H + cos((1:n)*2*pi/n)*X3 + .5*rnorm(n)
Y <- (H==1)*(1 + X2 + .5*rnorm(n)) + (H==2)*(1 + 2*X2 + .7*rnorm(n))
X1 <- .5 - X2 + exp(-3*(1:n)/n)*Y + .7*rnorm(n)

d <- data.frame(Y=Y, X1=X1, X2=X2, X3=X3, H=factor(H), 
                E = factor(E, labels = c("environment 1", "environment 2", "environment 3")))

## S={1} is NOT h-invariant
ggplot(d, aes(X1, Y)) + geom_point(aes(shape = H), alpha = 0.7, size = 2.5) + 
  facet_grid(.~E)

## S={2} IS h-invariant
ggplot(d, aes(X2, Y)) + geom_point(aes(shape = H), col= "#339900", alpha = 0.7, size = 2.5) + 
  facet_grid(.~E)




## h-invariance of S={1} is rejected
res.test1 <- test.equality.sr(Y = Y, X = X1, E = E, number.of.states = 2)
res.test1$p.value

## h-invariance of S={2} is not rejected
res.test2 <- test.equality.sr(Y = Y, X = X2, E = E, number.of.states = c(2,3))
res.test2$p.value
