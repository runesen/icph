##################################################
## Figure generating scripts for 
## "Switching Regression Models and Causal Inference
## In the Presence of Latent Variables"
## Submitted to JMLR
## Author: Rune Christiansen
## email: krunechristiansen@math.ku.dk
##################################################
## generating Figure 14 (right)
## run fig14_3_sim.R to generate the data file
## fig14_3_data.txt
##################################################
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
library(ggplot2)
theme_set(theme_bw() + theme(text = element_text(size=13)))

var.frame <- read.table("fig14_3_data.txt", header = TRUE, sep="\t")

rej.frame <- data.frame(var = rep(unique(var.frame$var), each = length(unique(var.frame$intervention))*length(unique(var.frame$method))),
                        intervention = rep(rep(unique(var.frame$intervention), each=length(unique(var.frame$method))), length(unique(var.frame$var))),
                        method = rep(unique(var.frame$method), length(unique(var.frame$var))*length(unique(var.frame$intervention))),
                        rej = c(sapply(unique(var.frame$var), function(v){
                          c(sapply(unique(var.frame$intervention), function(interv){
                            c(sapply(unique(var.frame$method), function(me){
                              mean((subset(var.frame, var == v & intervention == interv & method == me)$reject.nonCausal), na.rm = TRUE)
                            }))
                          }))
                        })
                        ))

rej.frame$intervention <- factor(rej.frame$intervention)
rej.frame$method <- factor(rej.frame$method)
rej.frame$var <- factor(rej.frame$var)

p <- ggplot(rej.frame, aes(intervention, rej, group = var)) + 
  geom_point(aes(col=var), size=3, alpha=.8) + 
  geom_hline(yintercept = .05, lty = "dashed", col="red", size=1) + 
  geom_line(size=1, aes(col=var)) + 
  scale_y_continuous(limits = c(0,1)) + 
  ylab("rejecetion rate for non-causality") + 
  xlab(expression(paste("strength "~Delta~" of intervention on Y"))) + 
  ggtitle("violations of h-invariance") + 
  theme(plot.title = element_text(hjust = 0.5),
        plot.margin = unit(c(2,2,2,2), "mm"), 
        legend.position = "right") +
  scale_color_manual("variable", 
                     labels = c(expression(X^1), expression(X^2), expression(X^3)),
                     values = c("#56B4E9", "#0072B2", "#D55E00")) +
  scale_linetype_manual(values = c("solid","longdash","twodash")) + 
  guides(size = "none", alpha = "none", lty = "none", 
         color = guide_legend("variable         ", order=1))
p

pdf("~/Google Drive/phd/causing/drafts/icph/figures/fig12_3.pdf", width=5.5, height=4)
p
dev.off()
