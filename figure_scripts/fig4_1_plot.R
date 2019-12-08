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
library(ggplot2)
theme_set(theme_bw() + theme(axis.text = element_text(size = 11),
                             axis.title = element_text(size = 11)))

# change path when final
df <- read.table("fig4_1_data.txt", header=TRUE)
df$init <- factor(df$init, levels = c("regmix", "true"))
df.aggr <- subset(df, sim.data == 1)
df.aggr$cov <- aggregate(p.value~init+sim.para, df, FUN=function(p.value) mean(p.value>=.05))$p.value

mean(subset(df, init == "true")$loglik > subset(df, init == "regmix")$loglik+.0001)
mean(subset(df, init == "regmix")$loglik > subset(df, init == "true")$loglik+.0001)


p1 <- ggplot(df.aggr, aes(sqrt(mp1*mp2), cov, group = init)) +
  geom_point(aes(shape=init), alpha=.7, size=4) +
  guides(color="none") + 
  scale_shape_manual(name="initialization", 
                     values=c(16,1), 
                     labels = c("data driven", "true values")) +
  scale_linetype_manual(name="initialization", 
                        values=c(1,2), 
                        labels = c("data driven", "true values")) +
  geom_smooth(aes(lty=init),se=FALSE) + 
  geom_hline(yintercept = .95, size=1, col="red", lty=2) + 
  coord_cartesian(ylim=c(0,1)) + 
  ylab("empirical coverage of \n 95% confidence regions") + 
  # xlab(expression(pi(P))) + 
  xlab("GMEP") + 
  theme(legend.position = c(.8,.35), 
        plot.margin = unit(c(5,5,5,5), "mm"), 
        plot.title = element_text(hjust = 0.5)) 

p1

pdf("~/Google Drive/phd/causing/drafts/icph/figures/fig4_1.pdf", width = 4.4, height = 3.8)
print(p1)
dev.off()

