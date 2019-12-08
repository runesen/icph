##################################################
## Figure generating scripts for 
## "Switching Regression Models and Causal Inference
## In the Presence of Latent Variables"
## Submitted to JMLR
## Author: Rune Christiansen
## email: krunechristiansen@math.ku.dk
##################################################
## generating Figure 8
## run script fig8_sim.R to generate the files 
## fig8_data_set.txt and fig8_data_var.txt
##################################################
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
library(ggplot2)
library(reshape2)
library(gridExtra)
library(grid)
theme_set(theme_bw() + theme(axis.text = element_text(size = 11),
                             axis.title = element_text(size = 11)))

set.frame <- read.table("fig8_data_set.txt", header = TRUE, sep = "\t")
var.frame <- read.table("fig8_data_var.txt", header = TRUE)

set.frame <- subset(set.frame, method == "NLM")
var.frame <- subset(var.frame, method == "NLM")


set.frame$S <- factor(set.frame$S, levels = c("1, 2", "2", "1", "", "3", "1, 3", "2, 3", "1, 2, 3"))

# set.frame$method <- factor(set.frame$method, levels = c("EM", "NLM"),
#                            labels = c("method: EM", "method: NLM"))
set.frame$arbvar <- factor(set.frame$arbvar, levels = c(TRUE, FALSE),
                           labels = c("variance constraint: lower bound", "variance constraint: equality"))

Slevels <- c("", "1", "2", "3", "1, 2", "1, 3", "2, 3", "1, 2, 3")
rej.frame <- data.frame(arbvar = rep(unique(set.frame$arbvar), each = length(unique(set.frame$n))*length(Slevels)),
                        n = rep(rep(unique(set.frame$n), each = length(Slevels)), length(unique(set.frame$arbvar))),
                        S = factor(rep(Slevels, length(unique(set.frame$n))*length(unique(set.frame$arbvar))), levels = Slevels),
                        rej = c(sapply(unique(set.frame$arbvar), function(a){
                          c(sapply(unique(set.frame$n), function(m){
                            c(sapply(Slevels, function(s){
                              mean(subset(set.frame, arbvar == a & n == m & S == s)$reject, na.rm=T)
                            }))
                          }))
                        })))


#############

p11 <- ggplot(subset(set.frame, Shat & arbvar == "variance constraint: lower bound"), aes(x=n, order = S)) + 
  geom_bar(aes(fill=S), width = 90) + 
  #facet_grid(.~arbvar) + 
  xlab("sample size") + 
  ylab(expression(paste("number of simulations resulting in ", hat(S), " = ..."))) + 
  geom_hline(yintercept = 5, linetype = 2, col = "red") + 
  scale_fill_manual(name = "S",
                    limits = c("", "1", "2", "1, 2", "3", "1, 3", "2, 3", "1, 2, 3"),
                    labels = c("{}  ", "{1}  ","{2}  ","{1,2}  ","{3}  ","{1,3}  ", "{2,3}  ", "{1,2,3}"),
                    values = c("#999999", "#56B4E9", "#0072B2", "#009E73", "#D55E00", "#F0E442", "#E69F00", "#CC79A7")) + 
  theme(legend.position = "bottom",
        strip.text = element_text(face="bold")) + 
  guides(fill = guide_legend(nrow=1))

p12 <- ggplot(subset(set.frame, Shat & arbvar == "variance constraint: equality"), aes(x=n, order = S)) + 
  geom_bar(aes(fill=S), width = 90) + 
  #facet_grid(.~arbvar) + 
  xlab("sample size") + 
  ylab(expression(paste("number of simulations resulting in ", hat(S), " = ..."))) + 
  geom_hline(yintercept = 5, linetype = 2, col = "red") + 
  scale_fill_manual(name = "S",
                    limits = c("", "1", "2", "1, 2", "3", "1, 3", "2, 3", "1, 2, 3"),
                    labels = c("{}  ", "{1}  ","{2}  ","{1,2}  ","{3}  ","{1,3}  ", "{2,3}  ", "{1,2,3}"),
                    values = c("#999999", "#56B4E9", "#0072B2", "#009E73", "#D55E00", "#F0E442", "#E69F00", "#CC79A7")) + 
  theme(legend.position = "bottom",
        strip.text = element_text(face="bold")) + 
  guides(fill = guide_legend(nrow=1))

p21 <- ggplot(subset(rej.frame, arbvar == "variance constraint: lower bound"), aes(x=n, y = rej, group = S, col = S)) + geom_point(size = 3, alpha = .7) + geom_line(size = 1, alpha = .7) + 
  ylab(expression(paste("rejecetion rates for ", H["0,S"]))) + 
  xlab("sample size") + ylim(c(0,1)) + 
  #facet_grid(.~arbvar) + 
  geom_hline(yintercept = .05, linetype = 2, col = "red") + 
  scale_color_manual(name = "S",
                     limits = c("", "1", "2", "1, 2", "3", "1, 3", "2, 3", "1, 2, 3"),
                     labels = c("{}  ", "{1}  ","{2}  ","{1,2}  ","{3}  ","{1,3}  ", "{2,3}  ", "{1,2,3}"),
                     values = c("#999999", "#56B4E9", "#0072B2", "#009E73", "#D55E00", "#F0E442", "#E69F00", "#CC79A7")) + 
  theme(legend.position = "bottom",
        strip.text = element_text(face="bold"))

p22 <- ggplot(subset(rej.frame, arbvar == "variance constraint: equality"), aes(x=n, y = rej, group = S, col = S)) + geom_point(size = 3, alpha = .7) + geom_line(size = 1, alpha = .7) + 
  ylab(expression(paste("rejecetion rates for ", H["0,S"]))) + 
  xlab("sample size") + ylim(c(0,1)) + 
  #facet_grid(.~arbvar) + 
  geom_hline(yintercept = .05, linetype = 2, col = "red") + 
  scale_color_manual(name = "S",
                     limits = c("", "1", "2", "1, 2", "3", "1, 3", "2, 3", "1, 2, 3"),
                     labels = c("{}  ", "{1}  ","{2}  ","{1,2}  ","{3}  ","{1,3}  ", "{2,3}  ", "{1,2,3}"),
                     values = c("#999999", "#56B4E9", "#0072B2", "#009E73", "#D55E00", "#F0E442", "#E69F00", "#CC79A7")) + 
  theme(legend.position = "bottom",
        strip.text = element_text(face="bold")) 


p1 <- icph:::grid_arrange_shared_legend(p11, p21)
p2 <- icph:::grid_arrange_shared_legend(p12, p22)

grid.draw(p1)
grid.draw(p2)

pdf("~/Google Drive/phd/causing/drafts/icph/figures/fig8_1.pdf")
grid.draw(p1)
dev.off()

pdf("~/Google Drive/phd/causing/drafts/icph/figures/fig8_2.pdf")
grid.draw(p2)
dev.off()
