##################################################
## Figure generating scripts for 
## "Switching Regression Models and Causal Inference
## In the Presence of Latent Variables"
## Submitted to JMLR
## Author: Rune Christiansen
## email: krunechristiansen@math.ku.dk
##################################################
## genetating Figure 14 (left)
## run fig14_1_sim.R to simulate data fig14_1_data_*.txt
##################################################
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
library(ggplot2)
library(Cairo)
theme_set(theme_bw() + theme(text = element_text(size=13)))

set.frame <- read.table("fig14_1_data_set.txt", header = TRUE, sep = "\t")
var.frame <- read.table("fig14_1_data_var.txt", header = TRUE)

set.frame$method <- set.frame$l
var.frame$method <- var.frame$l

rej.frame <- data.frame(var = rep(unique(var.frame$var), each = length(unique(var.frame$K))*length(unique(var.frame$method))),
                        K = rep(rep(unique(var.frame$K), each=length(unique(var.frame$method))), length(unique(var.frame$var))),
                        method = rep(unique(var.frame$method), length(unique(var.frame$var))*length(unique(var.frame$K))),
                        rej = c(sapply(unique(var.frame$var), function(v){
                          c(sapply(unique(var.frame$K), function(k){
                            c(sapply(unique(var.frame$method), function(me){
                              mean((subset(var.frame, var == v & K == k & method == me)$reject.nonCausal), na.rm = TRUE)
                            }))
                          }))
                        })
                        ))

rej.frame$method <- factor(rej.frame$method, 
                           levels = c("fixed2", "varying", "known"),
                           labels = c("k = 2", "k \u2264 5", expression(k^0)))
rej.frame$var <- factor(rej.frame$var)

p <- ggplot(rej.frame, aes(K, rej, group = interaction(var, method))) + 
  geom_point(aes(shape=method, col=var), size=3, alpha=.8) + 
  # scale_alpha_manual(values=c(.5,1), guide = FALSE) + 
  # scale_size_manual(values = c(3,3), guide = FALSE) +
  geom_line(size=1, aes(lty=method,col=var)) + 
  scale_y_continuous(limits = c(0,1)) + 
  geom_hline(yintercept = .05, lty = "dashed", col="red", size=1) + 
  ylab("rejecetion rate for non-causality") + 
  xlab(expression(paste("number of latent states ", k[0]))) + 
  ggtitle("non-binary latent variables") + 
  theme(plot.title = element_text(hjust = 0.5),
        plot.margin = unit(c(2,2,2,2), "mm"), 
        legend.position = "right", 
        legend.text.align = 0) +
  scale_color_manual("variable", 
                     labels = c(expression(X^1), expression(X^2), expression(X^3)),
                     values = c("#56B4E9", "#0072B2", "#D55E00")) +
  scale_shape_manual("testing for\nh-invariance \nof degree", 
                     #levels = c("fixed2", "varying", "known"), 
                     labels = c("k = 2", "k \u2264 5", expression(paste("k = ", k[0]))), 
                     values = c(16,17,15)) + 
  scale_linetype_manual(values = c("solid","longdash","twodash")) + 
  guides(size = "none", alpha = "none", lty = "none", 
         color = guide_legend("variable         ", order=1))
p

cairo_pdf("~/Google Drive/phd/causing/drafts/icph/figures/fig12_1.pdf", width=5.5, height=4)
p
dev.off()
