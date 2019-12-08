##################################################
## Figure generating scripts for 
## "Switching Regression Models and Causal Inference
## In the Presence of Latent Variables"
## Submitted to JMLR
## Author: Rune Christiansen
## email: krunechristiansen@math.ku.dk
##################################################
## producing Figure 4 (left) 
## run script fig4_1_sim.R to generate fig4_1_data.txt
##################################################

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
library(icph)
library(ggplot2)
library(reshape2)
library(gridExtra)
library(scales)
theme_set(theme_bw() + theme(axis.text = element_text(size = 11),
                             axis.title = element_text(size = 11)))

df <- read.table("fig4_1_data.txt", header=TRUE)
df <- subset(df, init == "regmix")
df$piP <- round(2*sqrt(df$mp1*df$mp2),1)/2

pval.frame <- data.frame()
pval.vec <- seq(.05,.95,.1)

for(p in unique(df$piP)){
  tmp <- subset(df, piP == p)
  value <- sapply(pval.vec, function(p) mean(tmp$p.value > p-0.05 & tmp$p.value <= p+.05, na.rm=T))
  tmpp <- data.frame(piP = p,p.value = pval.vec, value = value)
  pval.frame <- rbind(pval.frame, tmpp)
}

scale_trans <- function(x){
  w1 <- which(!(x < 0 | x > 1) &  x <= .1)
  w2 <- which(!(x < 0 | x > 1) &  x > .1)
  df <- x
  df[w1] <- 5*df[w1]
  df[w2] <- (5/9)*df[w2] + (.5-5/90)
  return(df)
}

ascale_trans <- function(x){
  w1 <- which(!(x < 0 | x > 1) &  x <= .5)
  w2 <- which(!(x < 0 | x > 1) &  x > .5)
  df <- x
  df[w1] <- df[w1]/5
  df[w2] <- (9/5)*(df[w2] - (.5-5/90))
  return(df)
}

my_trans <- function() trans_new("my", scale_trans, ascale_trans, breaks = function(x) seq(0,1,.1), domain = c(0,1))

p2 <- ggplot(pval.frame, aes(xmin=piP-.025, xmax = piP+.025, ymin=p.value-.05, ymax = p.value+.05)) + 
  geom_rect(aes(fill=value)) + 
  # xlab(expression(pi(P))) + 
  xlab("GMEP") + 
  ylab(expression(paste("p-value for (true) hypothesis ", H[0], ": ", theta, "=", theta^0))) + 
  scale_fill_gradient2(name = "Density", trans = "my", low = "black", high = "red", mid = "green", midpoint = 0.5) + 
  scale_x_continuous(breaks=seq(.5,.9,.1)) + 
  theme(legend.position = "none", 
        plot.margin = unit(c(6,2,5,5), "mm"),
        plot.title = element_text(hjust = 0.5))

p2.legend <- ggplot(data.frame(x=0, y=seq(0,.99,.0001)), aes(xmin=x-.05, xmax=x+.05, ymin=y,ymax=y+.01)) + 
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

print(grid.arrange(p2, p2.legend, ncol=2, widths=c(1,.15)))

pdf("~/Google Drive/phd/causing/drafts/icph/figures/fig4_2.pdf", width = 4.9, height = 4)
print(grid.arrange(p2, p2.legend, ncol=2, widths=c(1,.15)))
dev.off()

