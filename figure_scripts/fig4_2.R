setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
#devtools::install_github("runesen/icph/code")
library(icph)
library(ggplot2)
library(reshape2)
library(gridExtra)
library(scales)
theme_set(theme_bw() + theme(axis.text = element_text(size = 11),
                             axis.title = element_text(size = 11)))

#################################################### simulating data ####################################################
n <- 100
b.vec <- seq(0,.7,.1)
sigma <- .3
num.exp <- 1000
lambda <- .5
pars <- rep(c("beta", "gamma", "intercept", "sigma"), each = 2)
l <- 8
out <- NULL
s <- 0

for(i in 1:num.exp){
  for(j in 1:length(b.vec)){
    s <- s+1
    set.seed(s)
    
    b <- b.vec[j]
    mu <- runif(1,-1,1)
    beta <- runif(1,-1,1)
    
    theta0 <- list(list(mu = c(mu,mu), beta = matrix(c(beta,beta+b),1,2), sigma = rep(sigma,2), gamma = matrix(c(lambda, 1-lambda),2,2)))
    
    x <- rnorm(n)
    h <- sample(1:2, size = n, prob = c(lambda, 1-lambda), replace = TRUE)
    y <- numeric(n)
    
    
    y[h==1] <- mu + beta*x[h==1] + sigma*rnorm(n)[h==1]
    y[h==2] <- mu + (beta+b)*x[h==2] + sigma*rnorm(n)[h==2]
    
    outt <- tryCatch({
      tmp1 <- test.equality.sr(Y=y, X=x, E=rep(1,n), model = "IID", plot = FALSE,
                               par.test = list(method = "NLM", init = "regmix", theta0 = theta0))
      tmp3 <- test.equality.sr(Y=y, X=x, E=rep(1,n), model = "IID", plot = FALSE,
                               par.test = list(method = "EM", init = "regmix", theta0 = theta0))
      tmp11 <- test.equality.sr(Y=y, X=x, E=rep(1,n), model = "IID", plot = FALSE,
                                par.test = list(method = "NLM", init = "true", theta0 = theta0))
      tmp31 <- test.equality.sr(Y=y, X=x, E=rep(1,n), model = "IID", plot = FALSE,
                                par.test = list(method = "EM", init = "true", theta0 = theta0))
      
      outt <- rbind(tmp1$estimates, 
                    tmp3$estimates, 
                    tmp11$estimates, 
                    tmp31$estimates
      )
      outt$beta <- b
      outt$seed <- s
      outt$method <- rep(c(rep("NLM", l), rep("EM", l)),2)
      outt$arbvar <- FALSE
      outt$init <- c(rep("regmix", 2*l), rep("true", 2*l))
      outt$p.value <- rep(c(tmp1$p.values.cover, 
                            tmp3$p.values.cover, 
                            tmp11$p.values.cover, 
                            tmp31$p.values.cover 
      ), each = l)
      outt$mse <- rep(c(tmp1$mse, 
                        tmp3$mse, 
                        tmp11$mse, 
                        tmp31$mse 
      ), each = l)
      outt$loglik <- rep(c(tmp1$loglik[1,1], 
                           tmp3$loglik[1,1], 
                           tmp11$loglik[1,1], 
                           tmp31$loglik[1,1] 
      ), each = l)
      outt
    }, error = function(e){
      outt <- data.frame(parameter = pars, component = rep(1:2, 4), environment = 1, 
                         estimate = NA, std.err = NA, beta = b, seed = s, method = rep(c(rep("NLM", l), rep("EM", l)),2),
                         arbvar = FALSE, init = c(rep("regmix", 2*l), rep("true", 2*l)), p.value = NA, mse = NA, loglik = NA)
      outt
    })
    out <- rbind(out, outt)  
  }
}


#################################################### producing figure ####################################################

pval.frame <- data.frame()
pval.vec <- seq(.05,.95,.1)

for(m in c("EM", "NLM")){
  for(i in c("regmix", "true")){
    for(b in unique(out$beta)){
      tmp <- subset(out, method == m & arbvar == FALSE & beta == b & init == i)
      value <- sapply(pval.vec, function(p) mean(tmp$p.value > p-0.05 & tmp$p.value <= p+.05, na.rm=T))
      tmpp <- data.frame(method = m, arbvar = FALSE, beta = b, init = i, p.value = pval.vec, value = value)
      pval.frame <- rbind(pval.frame, tmpp)
    }
  }
}

scale_trans <- function(x){
  w1 <- which(!(x < 0 | x > 1) &  x <= .1)
  w2 <- which(!(x < 0 | x > 1) &  x > .1)
  out <- x
  out[w1] <- 5*out[w1]
  out[w2] <- (5/9)*out[w2] + (.5-5/90)
  return(out)
}

ascale_trans <- function(x){
  w1 <- which(!(x < 0 | x > 1) &  x <= .5)
  w2 <- which(!(x < 0 | x > 1) &  x > .5)
  out <- x
  out[w1] <- out[w1]/5
  out[w2] <- (9/5)*(out[w2] - (.5-5/90))
  return(out)
}

my_trans <- function() trans_new("my", scale_trans, ascale_trans, breaks = function(x) seq(0,1,.1), domain = c(0,1))

p <- ggplot(subset(pval.frame, method == "NLM" & arbvar == FALSE), aes(xmin=beta-.05, xmax = beta+.05, ymin=p.value-.05, ymax = p.value+0.05)) + 
  geom_rect(aes(fill=value)) + xlab(expression(paste("difference in regression coefficients "))) + 
  ylab(expression(paste("p-value for (true) hypothesis ", H[0], ": ", theta, "=", theta^0))) + 
  scale_fill_gradient2(name = "Density", trans = "my", low = "black", high = "red", mid = "green", midpoint = 0.5) + 
  scale_x_continuous(breaks=seq(0,1,.2)) + 
  theme(legend.position = "none", 
        plot.margin = unit(c(5,2,5,5), "mm"),
        plot.title = element_text(hjust = 0.5)) + 
  ggtitle(expression(paste(n, " = 100")))

p.legend <- ggplot(data.frame(x=0, y=seq(0,.99,.0001)), aes(xmin=x-.05, xmax=x+.05, ymin=y,ymax=y+.01)) + 
  geom_rect(aes(fill=y)) + 
  scale_fill_gradient2(name = "Density", trans = "my", low = "black", high = "red", mid = "green", midpoint = 0.5) + 
  scale_y_continuous(breaks = seq(0,1,.1), position = "right") + 
  theme(legend.position = "none", 
        axis.ticks.x = element_blank(), 
        axis.text = element_text(size=7),
        axis.ticks.y = element_line(size=.3),
        axis.text.x = element_blank(), 
        panel.grid = element_blank(), 
        panel.border = element_blank(), 
        panel.background = element_blank(),
        plot.margin = unit(c(10,2,15,2), "mm"))

pdf("fig4_2.pdf", width = 4.7, height = 4)
print(grid.arrange(p, p.legend, ncol=2, widths=c(1,.15)))
dev.off()