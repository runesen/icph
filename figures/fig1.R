setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
library(ggplot2)
library(reshape2)
theme_set(theme_bw())

n <- 600
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

line.data <- data.frame(variable = c(1,1,2,2,3,3), z = c(NA, NA, 1, 2, NA, NA), w = factor(c(NA, NA, 1, 2, NA, NA)))
break.data <- data.frame(variable = c(1,1,1,2,2,2,3,3,3), z = c(-2,-1,0, 1, 1.5, 2, -1, 0, 1))

p <- ggplot() + 
  geom_rect(data = rects, aes(xmin = xstart, xmax = xend, ymin = -Inf, ymax = Inf, fill=col), alpha=.4) + 
  geom_hline(data = line.data, aes(yintercept = z, lty = w), alpha = .7) +
  geom_line(data = dbeta, aes(time, beta, col=col), size = 1) + 
  facet_grid(variable~., scales = "free", labeller = label_bquote(X ^ .(variable))) + 
  #scale_y_continuous(data = break.data, aes(breaks = z)) + 
  scale_color_manual(values = c("black", "#339900")) +
  scale_linetype_manual(values = c("twodash", "dashed")) + 
  scale_fill_manual(name = "True (unknown) latent state",
                    values = c("#333333", "#999999"), 
                    labels = c("H=1", "H=2")) + 
  ylab("rolling window estimates of regression coefficients") + 
  guides(color = "none") + 
  theme(legend.position = "none", 
        axis.text = element_text(size = 16),
        strip.text = element_text(size = 16),
        axis.title = element_text(size = 16))

pdf("fig1.pdf", width = 11, height = 6.5)
print(p)
dev.off()
