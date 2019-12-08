##################################################
## Figure generating scripts for 
## "Switching Regression Models and Causal Inference
## In the Presence of Latent Variables"
## Submitted to JMLR
## Author: Rune Christiansen
## email: krunechristiansen@math.ku.dk
##################################################
## generating Figure 9
## run script fig9_sim.R to generate the file
## fig9_data.txt
##################################################
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
library(ggplot2)
library(reshape2)
library(gridExtra)
library(grid)
theme_set(theme_bw() + theme(axis.text = element_text(size = 11),
                             axis.title = element_text(size = 11)))


var.frame <- read.table("fig9_data.txt", header = TRUE, sep="\t")

rej.frame <- data.frame(n=rep(unique(var.frame$n), each=length(unique(var.frame$var))*length(unique(var.frame$violation))*length(unique(var.frame$method))),
                        var = rep(rep(unique(var.frame$var), each = length(unique(var.frame$violation))*length(unique(var.frame$method))), length(unique(var.frame$n))),
                        violation = rep(rep(unique(var.frame$violation), each=length(unique(var.frame$method))), length(unique(var.frame$var))*length(unique(var.frame$n))),
                        method = rep(unique(var.frame$method), length(unique(var.frame$n))*length(unique(var.frame$var))*length(unique(var.frame$violation))),
                        rej = c(sapply(unique(var.frame$n), function(m){
                          c(sapply(unique(var.frame$var), function(v){
                            c(sapply(unique(var.frame$violation), function(vio){
                              c(sapply(unique(var.frame$method), function(me){
                                mean((subset(var.frame, n == m & var == v & violation == vio & method == me)$reject.nonCausal), na.rm = TRUE)
                              }))
                            }))
                          }))
                        }))
)

rej.frame$violation <- factor(rej.frame$violation)
rej.frame$method <- factor(rej.frame$method)


## ggplot default colors
tmp <- data.frame(x=1:7, y=1:7, z=factor(1:7))
p <- ggplot(tmp, aes(x,y,col=z)) + geom_point()
cols <- unique(ggplot_build(p)$data[[1]][,1])[1:7]

rej.frame$violation <- factor(rej.frame$violation, 
                              levels = c("none", "heterogeneous variances", "uniform noise", "Laplace noise", 
                                         "mean shift", "variance shift", "continuous H"),
                              labels = c("none", "het. var.", "uniform noise", "Laplace noise", "mean shift", "var. shift", "cont. H"))

rej.frame$method <- factor(rej.frame$method,
                           levels = c("ICPH", "k-means ICP", "JCI-PC"),
                           labels = c("ICPH", "k-means ICP", "JCI-PC"))

p11 <- ggplot(subset(rej.frame, var == 1), aes(n, rej, group = interaction(method,violation), col = violation, alpha = (violation=="none"))) + 
  geom_point(aes(shape=method, size = (violation=="none"))) + 
  scale_alpha_manual(values=c(.5,1), guide = FALSE) + 
  scale_size_manual(values = c(3,3), guide = FALSE) + 
  geom_line(size=1, aes(lty=method)) + 
  scale_y_continuous(limits = c(0,1)) + 
  ylab("rejecetion rate for non-causality") + xlab("sample size") + ggtitle(expression(bold(X^1))) + 
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5),
        plot.margin = unit(c(2,8,2,8), "mm")) +
  scale_color_manual("model violation", values = c("black",cols)) +
  scale_linetype_manual(values = c("solid","longdash","twodash")) + 
  guides(size = "none", alpha = "none", lty = "none", 
         shape = guide_legend(order=1, nrow=1, override.aes = list(size = 3)), 
         color = guide_legend(order=2,nrow=1)) 

p22 <- ggplot(subset(rej.frame, var == 2), aes(n, rej, group = interaction(method,violation), col = violation, alpha = (violation=="none"))) + 
  geom_point(aes(shape=method, size = (violation=="none"))) + 
  scale_alpha_manual(values=c(.5,1), guide = FALSE) + 
  scale_size_manual(values = c(3,3), guide = FALSE) + 
  geom_line(size=1, aes(lty=method)) + 
  scale_y_continuous(limits = c(0,1)) + 
  ylab("rejecetion rate for non-causality") + xlab("sample size") + ggtitle(expression(bold(X^2))) + 
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5),
        plot.margin = unit(c(2,8,2,8), "mm")) + 
  scale_color_manual("model violation", values = c("black",cols), guide = guide_legend(nrow=1)) +
  scale_linetype_manual(values = c("solid","longdash","twodash")) + 
  guides(size = "none", alpha = "none", lty = "none", shape = guide_legend(order=1), color = guide_legend(order=2)) 

p33 <- ggplot(subset(rej.frame, var == 3), aes(n, rej, group = interaction(method,violation), col = violation, alpha = (violation=="none"))) + 
  geom_point(aes(shape=method, size = (violation=="none"))) + 
  scale_alpha_manual(values=c(.5,1), guide = FALSE) + 
  scale_size_manual(values = c(3,3), guide = FALSE) + 
  geom_line(size=1, aes(lty=method)) + 
  scale_y_continuous(limits = c(0,1)) + 
  ylab("rejecetion rate for non-causality") + xlab("sample size") + ggtitle(expression(bold(X^3))) + 
  geom_hline(yintercept = .05, col = "red", lty = 2) + 
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5),
        plot.margin = unit(c(2,8,2,8), "mm")) + 
  scale_color_manual("model violation", values = c("black",cols), guide = guide_legend(nrow=1)) +
  scale_linetype_manual(values = c("solid","longdash","twodash")) + 
  guides(size = "none", alpha = "none", lty = "none", shape = guide_legend(order=1), color = guide_legend(order=2)) 

pp <- icph:::grid_arrange_shared_legend(p11, p22, p33, ncol = 3, nrow = 1)
grid.draw(pp)

pdf("~/Google Drive/phd/causing/drafts/icph/figures/fig9.pdf", width=11.5, height=3.5)
grid.draw(pp)
dev.off()
