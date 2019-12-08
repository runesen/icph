##################################################
## Figure generating scripts for 
## "Switching Regression Models and Causal Inference
## In the Presence of Latent Variables"
## Submitted to JMLR
## Author: Rune Christiansen
## mail: krunechristiansen@math.ku.dk
##################################################
## producing Figure 12
##################################################
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
setwd("~/Google Drive/phd/causing/data/SIF/NA_GPP")
#devtools::install_github("runesen/icph/code")
library(raster)
library(rgdal)
library(lattice)
library(ggplot2)
library(reshape2)
library(icph)
library(gridExtra)
theme_set(theme_bw())

aparfile<-character(12)
SIFfile<-character(12)

doy <- c("001", "032", "060", "091", "121", "162", "182", "213", "244", "274", "305", "335")
mon <- c("01", "02", "03","04","05","06","07","08","09","10","11","12")

for (i in 1:5){
  for (j in 1:12){
    aparfile[(i-1)*12+j]<-paste('Monthly_APAR_hd/apar.',(i+2009)*100+j,'.hd.SG.tif',sep='')
  }
}


SIFfile<-list.files('monthly_SIF/',pattern='tif',full.name=TRUE)

lc<-raster('max.2010.hd.tif')
pct<-raster('max_per.2010.hd.tif')
lc.copy <- lc
lc[pct<0.75] <- NA
lc[!(lc %in% c(1,12))] <- NA # only ENF and CRO

lc.frame <- as.data.frame(lc, xy=TRUE)

SIF <- APAR <- VEG <- XY <- S <- NULL
for(j in 1:12){
  
  sif <- raster(SIFfile[j])
  apar <- raster(aparfile[j])/100
  
  sif[is.na(lc)] <- NA
  apar[is.na(lc)] <- NA
  
  sif <- getValues(sif)
  apar <- getValues(apar)
  veg <- 1+floor(getValues(lc)/11)
  
  w <- which(!is.na(sif) & !is.na(apar) & sif!=0 & apar!=0)
  
  xy <- lc.frame[w,1:2]
  
  XY <- rbind(XY, xy)
  S <- c(S, w)
  SIF <- c(SIF, sif[w])
  APAR <- c(APAR, apar[w])
  VEG <- c(VEG, veg[w])
}

set.seed(2)
res <- test.equality.sr(Y=SIF, X=APAR, E=rep(1,length(SIF)), intercept=FALSE, model="IID")
betahat <- subset(res$estimates, parameter == "beta")$estimate
sigmahat <- subset(res$estimates, parameter == "sigma")$estimate

dens <- sapply(1:2, function(i){
  beta <- betahat[i]
  sigma <- sigmahat[i]
  dnorm(SIF, mean = beta * APAR, sd = sigma)
})

lambda <- matrix(subset(res$estimates, parameter == "gamma")$estimate, nrow = nrow(dens), ncol= ncol(dens), byrow=TRUE)
post.probs <- dens*lambda / matrix(rowSums(dens*lambda), nrow = nrow(dens), ncol= ncol(dens))

d <- data.frame(x = XY$x, y = XY$y, sif = SIF, apar = APAR, h = VEG, hhat = post.probs%*%c(1,2))
ord <- order(d$x, d$y)
d <- d[ord,]
dd <- aggregate(hhat~x+y+h, d, function(hhat) round(mean(hhat)))
count <- aggregate(hhat~x+y+h, d, function(hhat) length(hhat))$hhat
df <- data.frame(x=rep(dd$x, count), y=rep(dd$y, count), hhat=rep(dd$hhat, count))
df <- df[order(df$x, df$y),]
d$hhat <- factor(df$hhat)
mean(d$hhat == d$h)
d$hmatch <- factor(ifelse(d$h == d$hhat, d$hhat, rep(3,nrow(d))))
w <- which(d$hmatch == 3);w<-c()
ww <- setdiff(1:nrow(d), w)
ww <- sample(ww, size = length(ww), replace = FALSE)

p1 <- ggplot(d[c(ww,w),], aes(apar, sif)) + geom_point(aes(col = hhat, shape = hhat), alpha = .7) + 
  scale_color_manual(name = "", 
                     values = c("darkorange2", "chartreuse4"), 
                     labels = c(expression(paste(hat(H), " = 1")), expression(paste(hat(H), " = 2")))) + 
  geom_abline(intercept = c(0,0), slope = betahat, col = c("darkorange2", "chartreuse4"), lty = c(2,1))  + 
  xlab(expression(X^2)) + ylab("Y") + 
  # ggtitle("Fluorescence yield") + 
  guides(colour = guide_legend(override.aes = list(size = 5, shape = c(16:17)))) + 
  guides(shape = "none") + 
  #guides(shape = "none", 
  #       color = guide_legend(override.aes = list(size=3))) + 
  scale_shape_manual(name = "", 
                     breaks = c(1,2), 
                     labels = c(expression(paste(hat(H), " = 1")), expression(paste(hat(H), " = 2"))),
                     values = 16:17) + 
  theme(legend.position = "bottom", 
        plot.title = element_text(hjust = 0.5))
p1

lc.copy[pct<0.75] <- 0;plot(lc.copy)
lc.copy[lc.copy %in% (1:16)[-c(1,12)]] <- 0;plot(lc.copy)
lc.full <- as.data.frame(lc.copy, xy=TRUE)
colnames(lc.full)[3] <- "h"
lc.full <- subset(lc.full, h != 17)
lc.full$h[lc.full$h == 12] <- 2
lc.full$hhat <- lc.full$h

w <- which(duplicated(rbind(dd[,1:2],lc.full[,1:2])))
frame <- rbind(dd, lc.full)[-w,]
frame$h <- factor(frame$h)
frame$hhat <- factor(frame$hhat)
tmp.melt <- melt(frame, id.vars =c("x", "y"))
tmp.melt$value <- factor(tmp.melt$value)

p2 <- ggplot(subset(frame, y < 65 & x > -140), aes(xmin=x-.25, xmax=x+.25, ymin=y-.25, ymax=y+.25,fill=h, alpha=h),col="transparent") + 
  geom_rect() + 
  scale_fill_manual(name = "", 
                    breaks = c(1,2), 
                    labels = c("ENF", "CRO"),
                    values = c("#999999", "darkorange2", "chartreuse4")) + 
  guides(alpha = "none") + 
  scale_alpha_manual(values = c(.5, 1,1)) + 
  xlab("longitude") + ylab("latitude") + 
  # ggtitle("IGBP land cover classification") + 
  theme(legend.position = "bottom", 
        plot.title = element_text(hjust = 0.5))
p2

p3 <- ggplot(subset(frame, y < 65 & x > -140), aes(xmin=x-.25, xmax=x+.25, ymin=y-.25, ymax=y+.25,fill=hhat, alpha=hhat),col="transparent") + 
  geom_rect() + 
  scale_fill_manual(name = "", 
                    breaks = c(1,2), 
                    labels = c(expression(paste(hat(H), " = 1")), expression(paste(hat(H), " = 2"))),
                    values = c("#999999", "darkorange2", "chartreuse4")) + 
  guides(alpha = "none") + 
  scale_alpha_manual(values = c(.5, 1,1)) + 
  xlab("longitude") + ylab("latitude") + 
  # ggtitle("Classification based on reconstruction of H") + 
  theme(legend.position = "bottom", 
        plot.title = element_text(hjust = 0.5)) 
p3

grid.arrange(p2, p3, p1, widths = c(1.3, 1.3, .95), ncol = 3)

pdf("~/Google Drive/phd/causing/drafts/icph/figures/fig12.pdf", width = 12, height = 3.5)
grid.arrange(p2, p3, p1, widths = c(1.3, 1.3, .95), ncol = 3)
dev.off()
