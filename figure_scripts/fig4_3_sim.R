##################################################
## Figure generating scripts for 
## "Switching Regression Models and Causal Inference
## In the Presence of Latent Variables"
## Submitted to JMLR
## Author: Rune Christiansen
## email: krunechristiansen@math.ku.dk
##################################################
## simulating data for Figure 4 (right)
##################################################
onServer <- TRUE
library(pracma)
library(icph)

df <- read.table("~/server_scripts/figure_scripts/fig4_1_data.txt", header=TRUE)
df <- subset(df, sim.data==1 & init == "regmix")
df$piP <- sqrt(df$mp1*df$mp2)
piP.vec <- c(.5,.55,.6,.65,.7)
j.vec <- sapply(piP.vec, function(p) which.min(abs(p-df$piP)))

n.vec <- c(100,200,500,1000)
num.exp <- 1000
out <- NULL
s <- 0

for(j in j.vec){
  set.seed((j-1)*1000+1)
  sY <- runif(1,.1,.5)
  mu <- rep(runif(1,-1,1),2)
  beta <- runif(2,-1,1)
  muX <- runif(1,-1,1)
  sX <- runif(1,.1,1)
  lambda <- runif(1,.3,.7)
  piP <- df$piP[j]
  
  theta0 <- list(list(mu = mu, beta = matrix(beta,1,2), sigma = rep(sY,2), gamma = matrix(c(lambda, 1-lambda),2,2)))
  
  for(n in n.vec){
    for(i in 1:num.exp){
      s <- s+1
      set.seed(s)
      print(paste0(s, " out of ", length(j.vec)*length(n.vec)*num.exp))
      
      x <- muX + sX*rnorm(n)
      h <- sample(1:2, size = n, prob = c(lambda, 1-lambda), replace = TRUE)
      y <- numeric(n)
      y[h==1] <- mu[1] + beta[1]*x[h==1] + sY*rnorm(n)[h==1]
      y[h==2] <- mu[2] + beta[2]*x[h==2] + sY*rnorm(n)[h==2]
      #plot(x,y,col=h)
      
      outt <- tryCatch({
        nlm.regmix <- test.equality.sr(Y=y, X=x, E=rep(1,n), model = "IID", plot = FALSE,
                                       par.test = list(method = "NLM", init = "regmix", theta0 = theta0, optim.n.inits=5))
        
        outt <- data.frame(piP = piP, 
                           p.value = nlm.regmix$p.values.cover,
                           method = "NLM", 
                           n = n, 
                           sim.data = i,
                           seed = s)
        outt
      }, error = function(e){
        print(e)
        outt <- data.frame(piP = piP, 
                           p.value = NA, 
                           method = "NLM", 
                           n = n, 
                           sim.data = i,
                           seed = s)
        outt
      })
      out <- rbind(out, outt)
    }
    write.table(out, "fig4_3_data.txt", quote = FALSE)
  }
}
