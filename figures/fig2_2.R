setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
#install_github("runesen/icph")
library(icph)
library(ggplot2)
library(reshape2)
library(gridExtra)
theme_set(theme_bw())

n <- 600
ne <- 3
K <- 2
E <- rep(1:3, each = n/3)

set.seed(18)

Gamma <- list(matrix(c(.98, .02,
                       .1, .9), 2,2,byrow=T),
              matrix(c(.9, .1,
                       .1, .9),2,2, byrow=T),
              matrix(c(.9, .1,
                       .02, .98),2,2, byrow=T)) 

H <- rep(1, n)
for(t in 2:n){
  H[t] <- sample(1:2, size=1, prob = Gamma[[E[t]]][H[t-1],])
}

X4 <- H * rnorm(n)
X3 <- H + cos((1:n)*2*pi/n)*X4 + .5*rnorm(n)
Y <- (H==1)*(1 + X3 + .5*rnorm(n)) + (H==2)*(1 + 2*X3 + .7*rnorm(n))
X2 <- .5 - X3 + exp(-3*(1:n)/n)*Y + .7*rnorm(n)

dd <- data.frame(H=H, Y=Y, X1=X2, X2=X3, X3=X4, E=E)

d <- data.frame(value = c(X2, X3, X4, Y), 
                variable = factor(rep(c("X1", "X2", "X3", "Y"), each=n),
                                  levels = c("Y", "X1", "X2", "X3")),
                H = factor(rep(H, 4)), 
                time = rep(1:n, 4), 
                E = factor(rep(E,4)))

Y <- dd$Y; X1 <- dd$X1; X2 <- dd$X2; X3 <- dd$X3; E <- dd$E; H <- dd$H

switches <- which(H[-1]-H[-n] !=0)

rects <- data.frame(xstart = c(1,switches), xend = c(switches,n), 
                    col = factor(rep(c(1:2), n)[1:(length(switches)+1)]))

## 200 data points rolling window

step <- 200
preds <- c("X1", "X2", "X3")
beta <- c(sapply(preds, function(p){
  
  d <- data.frame(Y = dd$Y, X = dd[,p], H = dd$H)
  
  beta <- sapply((step/2+1):(n-step/2), function(i){
    #tryCatch({
    LM <- lm(Y ~ X*H, data = d[(i-step/2):(i+step/2),])
    coefficients(LM)[2] + d$H[i]*coefficients(LM)[4]
    #}, error = function(e) NA)
  })
  beta
}
))

dbeta <- data.frame(time = (step/2+1):(n-step/2), 
                    beta = beta, 
                    variable = rep(1:3, each = (n-step)),
                    H = dd$H[(step+1):n], 
                    col = rep(c("1","2","1"), each=n-step))


t <- test.equality.sr(Y, X2, E)

plot.frame <- data.frame(y = Y, x = X2, hhat = factor(t$hhat), e = E)

mu    <- c(sapply(1:ne, function(i) c(sapply(1:2, function(j) subset(t$estimates, parameter == "intercept" & environment == i & component == j)$estimate))))
beta  <- c(sapply(1:ne, function(i) c(sapply(1:2, function(j) subset(t$estimates, parameter == "beta" & environment == i & component == j)$estimate))))

linedata <- data.frame(e=rep(1:ne,each=K), 
                        component = factor(rep(1:K, ne)), intercept=mu, slope=beta)

plot.frame$e <- factor(plot.frame$e, labels = c("environment 1", "environment 2", "environment 3"))
linedata$e <- factor(linedata$e, labels = c("environment 1", "environment 2", "environment 3"))

p <- ggplot(plot.frame, aes(x, y)) + 
  geom_abline(data = linedata, aes(intercept=intercept, slope=slope, lty = component)) + 
  geom_point(aes(shape = hhat), alpha = 0.7, col = "#339900", size = 2.5) + 
  facet_grid(.~e) + 
  xlab(expression(X^2)) + ylab("Y") + 
  scale_shape_manual(values = c(16, 2)) + 
  scale_linetype_manual(values = c("twodash", "dashed")) + 
  theme(legend.position = "none", 
        axis.text = element_text(size = 16),
        strip.text = element_text(size = 20),
        axis.title = element_text(size = 16))

pdf("fig2_2.pdf", width = 12, height = 4.3)
print(p)
dev.off()