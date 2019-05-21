setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
#install_github("runesen/icph")
library(icph)
library(ggplot2)
library(reshape2)
library(gridExtra)
library(scales)
theme_set(theme_bw() + theme(axis.text = element_text(size = 11),
                             axis.title = element_text(size = 11)))

#################################################### simulating data ####################################################
n.vec <- c(100, 200, 500, 1000)
b <- 0.3
sigma <- .3
num.exp <- 1000
lambda <- .5
pars <- rep(c("beta", "gamma", "intercept", "sigma"), each = 2)
l <- 8

outn <- NULL
s <- 0
for(i in 1:num.exp){
  for(j in 1:length(n.vec)){
    n <- n.vec[j]
    mu <- runif(1,-1,1)
    beta <- runif(1,-1,1)
    
    theta0 <- list(list(mu = c(mu,mu), beta = matrix(c(beta,beta+b),1,2), sigma = rep(sigma,2), gamma = matrix(c(lambda, 1-lambda),2,2)))
    
    s <- s+1
    set.seed(s)
    
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
      outt$n <- n
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
                         estimate = NA, std.err = NA, n = n, seed = s, method = rep(c(rep("NLM", l), rep("EM", l)),2),
                         arbvar = FALSE, init = c(rep("regmix", 2*l), rep("true", 2*l)), p.value = NA, mse = NA, loglik = NA)
      outt
    })
    outn <- rbind(outn, outt)  
  }
}


#################################################### producing figure ####################################################
cov.frame.n <- data.frame()

for(m in c("EM", "NLM")){
  for(i in c("regmix", "true")){
    for(nn in unique(outn$n)){
      tmp <- subset(outn, method == m & arbvar == FALSE & n == nn & init == i)
      cov <- mean(tmp$p.value >= alpha, na.rm=T)
      tmpp <- data.frame(method = m, arbvar = FALSE, n = nn, init = i, p.value = pval.vec, cov=cov)
      cov.frame.n <- rbind(cov.frame.n, tmpp)
    }
  }
}
cov.frame.n$method <- factor(cov.frame.n$method, levels = c("NLM", "EM"))
cov.frame.n$init <- factor(cov.frame.n$init, levels = c("regmix", "true"))

pdf("fig4_3.pdf", width = 4.4, height = 4)
print(
  ggplot(subset(cov.frame.n, arbvar == FALSE), aes(x=n, y=cov, group=interaction(method,init))) +
    geom_point(aes(col=method, shape = init), size = 2, alpha = .7) + 
    geom_line(aes(col=method, lty = init), size = 1, alpha = .7) + 
    geom_hline(yintercept = .95, lty=2, col="red") +
    xlab("sample size") + ylab("empirical coverage of \n 95% confidence regions") + 
    coord_cartesian(ylim = c(0,1)) + 
    scale_x_continuous(breaks=c(0,100,200,500,1000)) + 
    scale_color_manual(name = "method",
                       values = c("dodgerblue2", "limegreen")) + 
    scale_linetype_manual(name = "initialization",
                          values = c("solid", "dotted"),
                          labels = c("data driven", "true values")) + 
    guides(color = guide_legend(order=1),
           shape = "none",
           linetype = guide_legend(order=3)) + 
    theme(legend.position = c(.8,.35), 
          plot.margin = unit(c(3,5,6,5), "mm"), 
          plot.title = element_text(hjust = 0.5, vjust = -1)) + 
    ggtitle(expression(paste(Delta, beta, " = 0.3")))
)
dev.off()