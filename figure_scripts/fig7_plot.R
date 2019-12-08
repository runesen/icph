##################################################
## Figure generating scripts for 
## "Switching Regression Models and Causal Inference
## In the Presence of Latent Variables"
## Submitted to JMLR
## Author: Rune Christiansen
## mail: krunechristiansen@math.ku.dk
##################################################
## producing Figure 7
## run script fig7_sim.R to generate 
## fig7_data_set.txt and fig7_data_set.txt
##################################################
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
library(ggplot2)
library(reshape2)
library(gridExtra)
library(scales)
theme_set(theme_bw() + theme(axis.text = element_text(size = 10),
                             axis.title = element_text(size = 10),
                             plot.title = element_text(size = 10, hjust=.5)))

## loading data
set.frame <- read.table("fig7_data_set.txt", header = TRUE, sep = "\t")
var.frame <- read.table("fig7_data_var.txt", header = TRUE)
set.frame$S <- factor(set.frame$S, levels = c("1, 2", "2", "1", "", "3", "1, 3", "2, 3", "1, 2, 3"))
set.frame$inSstar <- rep(sapply(unique(set.frame$seed), function(s) subset(set.frame, seed == s & Shat)$S %in% c("", "1", "2", "1, 2")), 
                         each = length(unique(set.frame$S)))

aggregate(sqrt(piP) ~ dist, var.frame, mean, na.rm=1)

## ICPH false discovery rate
typeIS <- data.frame(n = rep(unique(set.frame$n), each = length(unique(set.frame$dist))),
                     dist=rep(unique(set.frame$dist),length(unique(set.frame$n))),
                     nocover = c(sapply(unique(set.frame$n), function(m){
                       c(sapply(unique(set.frame$dist), function(d){
                         mean((!subset(set.frame, S == "1, 2" & dist == d & n == m)$inSstar), na.rm = TRUE)
                       }))
                     }))
)

## type-I error for true hypothesis H_{0S*}
typeIH <- data.frame(n = rep(unique(set.frame$n), each = length(unique(set.frame$dist))*length(unique(set.frame$S))),
                     dist=rep(rep(unique(set.frame$dist), each = length(unique(set.frame$S))),length(unique(set.frame$n))),
                     S = factor(rep(unique(set.frame$S), length(unique(set.frame$dist))*length(unique(set.frame$n)))),
                     rej = c(sapply(unique(set.frame$n), function(m){
                       c(sapply(unique(set.frame$dist), function(d){
                         c(sapply(unique(set.frame$S), function(s){
                           mean((subset(set.frame, S == s & dist == d & n == m)$p.value<.05), na.rm = TRUE)
                         }))
                       }))
                     }))
)
maxrej <- max(subset(typeIH, S=="1, 2")$rej, na.rm=T)

## transformation needed for colors
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


## plots
pS <- ggplot(typeIS, aes(xmin=n-50, xmax = n+50, ymin=dist-.25, ymax = dist+.25)) + 
  geom_rect(aes(fill=nocover/maxrej)) + xlab("sample size") + 
  ylab("diff. in regr. coefficients (incrasing GMEP)") + 
  ggtitle("False causal discovery rate") + 
  geom_text(aes(x=n, y=dist, label=round(nocover,2))) + 
  scale_fill_gradient2(name = expression(paste("Rejection rate for ", H[0])), trans = "my", low = "green", high = "red", mid= "yellow", midpoint = .5) + 
  theme(legend.position = "none", plot.margin = unit(c(.2,.2,.2,1), "cm"))
pS

pH <- ggplot(subset(typeIH, S == "1, 2"), aes(xmin=n-50, xmax = n+50, ymin=dist-.25, ymax = dist+.25)) + 
  geom_rect(aes(fill=rej/maxrej)) + xlab("sample size") + 
  ylab("diff. in regr. coefficients (incrasing GMEP)") + 
  geom_text(aes(x=n, y=dist, label=round(rej,2))) + 
  ggtitle(expression(paste("Rejection rates for (true) null Hypotheis ", H[paste("0,",S^"*")]))) + 
  scale_fill_gradient2(name = expression(paste("Rejection rate for ", H[0])), trans = "my", low = "green", high = "red", mid= "yellow", midpoint = .5) + 
  theme(legend.position = "none", plot.margin = unit(c(.2,1,.2,.2), "cm"))
pH

pdf("~/Google Drive/phd/causing/drafts/icph/figures/fig7.pdf", width = 8.5, height = 4)
print(grid.arrange(pH, pS, ncol=2))
dev.off()

