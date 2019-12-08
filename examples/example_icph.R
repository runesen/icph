library(icph)
library(ggplot2)
library(lattice)
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

D <- cbind(Y, X1, X2, X3)

dd <- data.frame(value = c(X1, X2, X3, Y), 
                variable = factor(rep(c("X1", "X2", "X3", "Y"), each=n),
                                  levels = c("Y", "X1", "X2", "X3")),
                H = factor(rep(H, 4)), 
                time = rep(1:n, 4), 
                E = factor(rep(E,4)))

switches <- which(H[-1]-H[-n] !=0)
rects <- data.frame(xstart = c(1,switches), xend = c(switches,n), 
                    col = factor(rep(c(1:2), n)[1:(length(switches)+1)]))

## Timeseries plots
ggplot() + 
  geom_rect(data = rects, aes(xmin = xstart, xmax = xend, ymin = -Inf, ymax = Inf, fill=col), alpha=.4) +
  geom_line(data = dd, aes(time, value)) + 
  facet_grid(variable~., scales = "free") + 
  scale_fill_manual(name = "", 
                    values = c("#333333", "#999999"), 
                    labels = c("H=1", "H=2")) + 
  ylab("")

## Scatter plot matrix
splom(D, groups = H, col=c(1,2), pch = c(1,2), size=.1,
      panel=panel.superpose,
      key=list(title="",
               columns=2,
               points=list(pch=c(1,2),
                           col=c(1,2)),
               text=list(c("H=1","H=2"))))

## ICPH
res.icph <- icph(D=D, E=E, target=1, model = "HMM", par.icp = list(output.pvalues = TRUE))
res.icph
