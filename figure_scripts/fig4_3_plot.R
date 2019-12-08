##################################################
## Figure generating scripts for 
## "Switching Regression Models and Causal Inference
## In the Presence of Latent Variables"
## Submitted to JMLR
## Author: Rune Christiansen
## email: krunechristiansen@math.ku.dk
##################################################
## generating Figure 4 (right)
## run script fig4_3_sim.R to generate the file
## fig4_3_data.txt
##################################################
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
library(ggplot2)
library(reshape2)
library(gridExtra)
library(grid)
theme_set(theme_bw() + theme(axis.text = element_text(size = 11),
                             axis.title = element_text(size = 11)))


outn <- read.table("fig4_3_data.txt", header = TRUE)
outn.aggr <- aggregate(p.value ~ piP + n, outn, function(p) mean(p>=.05))
colnames(outn.aggr)[3] <- "coverage"


p3 <- ggplot(outn.aggr, aes(n, coverage, group = piP, col = piP)) + 
  geom_point(size=3) + geom_line(size=1) + 
  scale_x_continuous(breaks=c(0,100,200,500,1000)) + 
  coord_cartesian(ylim = c(0,1)) + 
  xlab("sample size") + 
  ylab("empirical coverage of \n 95% confidence regions") + 
  geom_hline(yintercept = .95, size=1, col="red", lty=2) + 
  guides(color = guide_legend(order=1, 
                              #title = expression(pi(P)), 
                              title = "GMEP",
                              reverse = T),
         shape = "none",
         linetype = "none") + 
  theme(legend.position = c(.8,.35), 
        plot.margin = unit(c(5,5,6,5), "mm"), 
        plot.title = element_text(hjust = 0.5, vjust = -1))

p3

pdf("~/Google Drive/phd/causing/drafts/icph/figures/fig4_3.pdf", width = 4.4, height = 3.8)
print(p3)
dev.off()
