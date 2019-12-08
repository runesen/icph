##################################################
## Figure generating scripts for 
## "Switching Regression Models and Causal Inference
## In the Presence of Latent Variables"
## Submitted to JMLR
## Author: Rune Christiansen
## mail: krunechristiansen@math.ku.dk
##################################################
## producing Figure 5 
##################################################
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
#devtools::install_github("runesen/icph/code")
library(icph)
library(ggplot2)
library(reshape2)
library(gridExtra)
theme_set(theme_bw() + theme(axis.text = element_text(size = 11),
                             axis.title = element_text(size = 13)))
  
n <- 500
set.seed(408)
  
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
sigmaY <- runif(n = 2, min = 0.1, max = 0.3) 
  
beta2 <- sample(c(-1,1),1) * runif(1, .5, 1.5)
beta22 <- beta2
beta3 <- sample(c(-1,1),1) * runif(1, .5, 1.5)
beta32 <- 0
betaY <- matrix(sample(c(-1,1),4,replace=T)*runif(4, .5, 1.5), 2, 2)
  
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

dat <- data.frame(X1=X1, X2=X2, X3=X3, Y=Y, H=H, E=E)
dat$H <- factor(dat$H-1)
dat$E <- factor(dat$E)
dat.melt <- melt(dat, id.vars = c("Y", "H", "E"))

####################### With empirical regression lines #######################

tmp <- icph(cbind(Y,X1,X2,X3), E, target = 1, model = "IID")
lines.data <- data.frame(intercept = NA, slope = NA, variable = rep(c("X1", "X2", "X3"), each = 6),
                         H = factor(rep(rep(c(0,1), each=3),3)), E = factor(rep(1:3,6)))
lines.data$intercept <- subset(tmp$estimates, ind %in% c(2,3,4) & parameter == "intercept")$estimate
lines.data$slope <- subset(tmp$estimates, ind %in% c(2,3,4) & parameter == "beta")$estimate

p1hat <- ggplot(subset(dat.melt, variable == "X1"), aes(value, Y)) + geom_point(aes(col=E, pch=H), alpha = .7) + 
  geom_abline(data=subset(lines.data, variable == "X1"), aes(intercept = intercept, slope = slope, col = E, lty = H), alpha = .7) + 
  scale_color_manual(name = "environment", 
                     values = c("black", "red", "#00BFC4"),
                     labels = c(paste("t = 1,...,", sum(dat$E==1)), 
                                paste("t = ", sum(dat$E==1)+1, ",...,", sum(dat$E!=3)),
                                paste("t = ", sum(dat$E!=3)+1, ",...,", length(dat$E)))) + 
  scale_shape_manual(name = "hidden state",
                     values = c(16,2),
                     labels = c(expression(paste(H["t"], " = 0")),
                                expression(paste(H["t"], " = 1")))) + 
  scale_linetype_manual(name = "model fit",
                        values = c(1,2),
                        labels = c(expression(paste(hat(H)["t"], " = 0")),
                                   expression(paste(hat(H)["t"], " = 1")))) + 
  scale_y_continuous(breaks = c(-2,0,2,4)) + 
  xlab(expression(X^1))  + 
  guides(color = guide_legend(order=1),
         shape = guide_legend(order=2),
         linetype = guide_legend(order=3))
  #theme(legend.position = "none")

p2hat <- ggplot(subset(dat.melt, variable == "X2"), aes(value, Y)) + geom_point(aes(col=E, pch=H), alpha = .7) +
  geom_abline(data=subset(lines.data, variable == "X2"), aes(intercept = intercept, slope = slope, col = E, lty = H), alpha = .7) + 
  scale_color_manual(name = "environment", 
                     values = c("black", "red", "#00BFC4"),
                     labels = c(paste("t = 1,...,", sum(dat$E==1)), 
                                paste("t = ", sum(dat$E==1)+1, ",...,", sum(dat$E!=3)),
                                paste("t = ", sum(dat$E!=3)+1, ",...,", length(dat$E)))) + 
  scale_shape_manual(name = "hidden state",
                     values = c(16,2),
                     labels = c(expression(paste(H["t"], " = 0")),
                                expression(paste(H["t"], " = 1")))) + 
  scale_linetype_manual(name = "model fit",
                        values = c(1,2),
                        labels = c(expression(paste(hat(H)["t"], " = 0")),
                                   expression(paste(hat(H)["t"], " = 1")))) + 
  scale_y_continuous(breaks = c(-2,0,2,4)) + 
  xlab(expression(X^2))  
  #theme(legend.position = "none")

p3hat <- ggplot(subset(dat.melt, variable == "X3"), aes(value, Y)) + geom_point(aes(col=E, pch=H), alpha = .7) +
  geom_abline(data=subset(lines.data, variable == "X3"), aes(intercept = intercept, slope = slope, col = E, lty = H), alpha = .7) + 
  scale_color_manual(name = "environment", 
                     values = c("black", "red", "#00BFC4"),
                     labels = c(paste0("t = 1,...,", sum(dat$E==1)), 
                                paste0("t = ", sum(dat$E==1)+1, ",...,", sum(dat$E!=3)),
                                paste0("t = ", sum(dat$E!=3)+1, ",...,", length(dat$E)))) + 
  scale_shape_manual(name = "hidden state (unknown)",
                     values = c(16,2),
                     labels = c(expression(paste(H["t"], " = 1")),
                                expression(paste(H["t"], " = 2")))) + 
  scale_linetype_manual(name = "model fit",
                        values = c(1,2),
                        labels = c(expression(paste(hat(H)["t"], " = 1")),
                                   expression(paste(hat(H)["t"], " = 2")))) + 
  scale_y_continuous(breaks = c(-2,0,2,4)) + 
  guides(color = guide_legend(order=1),
         shape = guide_legend(order=2),
         linetype = guide_legend(order=3)) + 
  xlab(expression(X^3))


pdf("../figures/fig5.pdf", width=13, height=3.5)
grid.arrange(p1hat + theme(legend.position = "none",
                           plot.margin = unit(c(2,5,2,2), "mm")), 
             p2hat + theme(legend.position = "none",
                           plot.margin = unit(c(2,5,2,5), "mm")),
             p3hat + theme(plot.margin = unit(c(2,2,2,5), "mm")),
             ncol=3,widths=c(1,1.02,1.45))
dev.off()


