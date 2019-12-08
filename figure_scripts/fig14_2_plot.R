##################################################
## Figure generating scripts for 
## "Switching Regression Models and Causal Inference
## In the Presence of Latent Variables"
## Submitted to JMLR
## Author: Rune Christiansen
## email: krunechristiansen@math.ku.dk
##################################################
## genetating Figure 14 (middle)
## run fig14_2_sim.R to simulate data fig14_2_data.txt
##################################################
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
library(ggplot2)
theme_set(theme_bw() + theme(text = element_text(size=13)))

var.frame <- read.table("fig14_2_data.txt", header = TRUE)


rej.frame <- data.frame(var = rep(unique(var.frame$var), each = length(unique(var.frame$nb.noisevar))),
                        nb.noisevar = rep(unique(var.frame$nb.noisevar), length(unique(var.frame$var))),
                        rej = c(sapply(unique(var.frame$var), function(v){
                          c(sapply(unique(var.frame$nb.noisevar), function(nb){
                            mean((subset(var.frame, var == v & nb.noisevar == nb)$reject.nonCausal), na.rm = TRUE)
                          }))
                        })
                        ))

rej.frame$var <- factor(rej.frame$var)


p <- ggplot(subset(rej.frame, nb.noisevar > 0), aes(nb.noisevar, rej, group = var)) + 
  geom_point(aes(col=var), size=3, alpha=.8) + 
  scale_x_continuous(trans = "log10") + 
  geom_line(size=1, aes(col=var)) + 
  scale_y_continuous(limits = c(0,1)) + 
  geom_hline(yintercept = .05, lty = "dashed", col="red", size=1) + 
  ylab("rejecetion rate for non-causality") + 
  xlab(expression(paste("# of predictors in addition to ", X^1, ", ", X^2, " and ",  X^3))) + 
  ggtitle("large systems of variables") + 
  theme(plot.title = element_text(hjust = 0.5),
        plot.margin = unit(c(2,2,2,2), "mm"), 
        legend.position = "right") +
  scale_color_manual("variable", 
                     labels = c(expression(X^1), expression(X^2), expression(X^3)),
                     values = c("#56B4E9", "#0072B2", "#D55E00")) +
  guides(size = "none", alpha = "none", lty = "none", 
         color = guide_legend("variable         ", order=1))
p


pdf("~/Google Drive/phd/causing/drafts/icph/figures/fig14_2.pdf", width=5.5, height=4)
p
dev.off()
