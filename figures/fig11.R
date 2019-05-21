setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
library(raster)
library(rgdal)
library(lattice)
library(ggplot2)
library(reshape2)
library(icph)
library(gridExtra)
theme_set(theme_bw())

vpmfile<-character(12)
aparfile<-character(12)
SIFfile<-character(12)
tempfile<-character(12)
swfile<-character(12)

doy <- c("001", "032", "060", "091", "121", "162", "182", "213", "244", "274", "305", "335")
mon <- c("01", "02", "03","04","05","06","07","08","09","10","11","12")

for (j in 1:12){
  aparfile[j]<-paste('data/Monthly_APAR_hd/apar.',201000+j,'.hd.SG.tif',sep='')
  swfile[j]<-paste("data/SW/GLDAS_NOAH025_M.A2010", mon[j],".021.nc4",sep="")
}

vpmfile <- list.files('data/MPI_HD/',pattern='tif',full.name=TRUE)

SIFfile<-list.files('data/monthly_SIF/',pattern='tif',full.name=TRUE)

lc<-raster('data/max.2010.hd.tif')
pct<-raster('data/max_per.2010.hd.tif')
lc[pct<0.75] <- NA
lc[!(lc %in% c(1,12))] <- NA # only ENF and CRO


factx <- facty <- 3:16


pval <- n <- c()
for(i in 1:length(factx)){
  
  lcfrac <- aggregate(lc, fact = c(factx[i], facty[i]), fun=mean, na.rm=TRUE)
  lcfrac[!(lcfrac %in% c(1,12))] <- NA
  
  SIF <- APAR <- GPP <- SW <- E <- VEG <- c()
  for(j in 1:12){
    
    sif <- raster(SIFfile[j])
    gpp <- raster(vpmfile[j])/100
    apar <- raster(aparfile[j])/100
    sw <- aggregate(crop(raster(swfile[j], varname = "SWdown_f_tavg"), extent(lc)), fact = c(2,2), fun=mean, na.rm=TRUE); sw[sw==0] <- NA
    
    sif[is.na(lc)] <- NA
    gpp[is.na(lc)] <- NA
    apar[is.na(lc)] <- NA
    sw[is.na(lc)] <- NA
    
    siffrac <- aggregate(sif, fact = c(factx[i], facty[i]), fun=mean, na.rm=TRUE)
    gppfrac <- aggregate(gpp, fact = c(factx[i], facty[i]), fun=mean, na.rm=TRUE)
    aparfrac <- aggregate(apar, fact = c(factx[i], facty[i]), fun=mean, na.rm=TRUE)
    swfrac <- aggregate(sw, fact = c(factx[i], facty[i]), fun=mean, na.rm=TRUE)
    
    siffrac[is.na(lcfrac)] <- NA
    gppfrac[is.na(lcfrac)] <- NA
    aparfrac[is.na(lcfrac)] <- NA
    swfrac[is.na(lcfrac)] <- NA
    
    sif <- getValues(siffrac)
    gpp <- getValues(gppfrac)
    apar <- getValues(aparfrac)
    sw <- getValues(swfrac)
    veg <- getValues(lcfrac)
    
    w <- which(!is.na(sif) & !is.na(gpp) & !is.na(apar) & !is.na(sw) & sif!=0 & gpp!=0 & apar!=0 & sw!=0)
    
    SIF <- c(SIF, sif[w])
    GPP <- c(GPP, gpp[w])
    APAR <- c(APAR, apar[w])
    SW <- c(SW, sw[w])
    E <- c(E, rep(1+(j %in% 2:7), length(w)))
    VEG <- c(VEG, veg[w])
  }
  
  set.seed(i)
  
  res0 <- tryCatch({
    test.equality.sr(Y=SIF, X=matrix(APAR, ncol=1)[,-1], E=E, model = "IID", intercept = TRUE, plot=FALSE)
  }, error = function(e) list(p.value=NA))
  
  res1 <- tryCatch({
    test.equality.sr(Y=SIF, X=APAR, E=E, model = "IID", intercept = FALSE, plot=FALSE)
  }, error = function(e) list(p.value=NA))
  
  res2 <- tryCatch({
    test.equality.sr(Y=SIF, X=GPP, E=E, model = "IID", intercept = FALSE, plot=FALSE)
  }, error = function(e) list(p.value=NA))
  
  res3 <- tryCatch({
    test.equality.sr(Y=SIF, X=SW, E=E, model = "IID", intercept = FALSE, plot=FALSE)
  }, error = function(e) list(p.value=NA))
  
  res12 <- tryCatch({
    test.equality.sr(Y=SIF, X=cbind(APAR, GPP), E=E, model = "IID", intercept = FALSE, plot=FALSE)
  }, error = function(e) list(p.value=NA))
  
  res13 <- tryCatch({
    test.equality.sr(Y=SIF, X=cbind(APAR, SW), E=E, model = "IID", intercept = FALSE, plot=FALSE)
  }, error = function(e) list(p.value=NA))
  
  res23 <- tryCatch({
    test.equality.sr(Y=SIF, X=cbind(GPP, SW), E=E, model = "IID", intercept = FALSE, plot=FALSE)
  }, error = function(e) list(p.value=NA))
  
  res123 <- tryCatch({
    test.equality.sr(Y=SIF, X=cbind(APAR, GPP, SW), E=E, model = "IID", intercept = FALSE, plot=FALSE)
  }, error = function(e) list(p.value=NA))
  
  pval <- c(pval, res0$p.value, res1$p.value, res2$p.value, res3$p.value, 
            res12$p.value, res13$p.value, res23$p.value, res123$p.value)
  n <- c(n, rep(length(E), 8))
  
  print(paste("after i =", i))
}

frame <- data.frame(S = factor(rep(c("{}", "APAR", "GPP", "SW", "APAR, GPP", "APAR, SW", "GPP, SW", "APAR, GPP, SW"), length(factx)), 
                               levels = c("{}", "APAR", "GPP", "SW", "APAR, GPP", "APAR, SW", "GPP, SW", "APAR, GPP, SW")),
                    n = n, 
                    p.value = pmax(pval, 1e-4))

tmp <- c()
for(m in unique(frame$n)){
  
  sub <- subset(frame, n == m)
  S.list <- list(numeric(0), 1, 2, 3, 1:2, c(1,3), 2:3, 1:3)
  w <- which(sub$p.value >= 0.05)
  Shat <- Reduce(intersect, S.list[w])
  if(is.null(Shat) | length(Shat) == 0){
    wShat <- 1
  } else{
    wShat <- 1+which(sapply(2:length(S.list), function(i) (all.equal(S.list[[i]], Shat) == TRUE)))
  }
  
  Shat <- rep(FALSE, length(S.list))
  Shat[wShat] <- TRUE
  tmp <- c(tmp, Shat)
}
frame$Shat <- tmp


p <- ggplot(frame[order(frame$Shat),], aes(n, log10(p.value), group = S, col = S)) + geom_point(aes(shape = Shat, size = Shat)) + 
  geom_line(size = .7, alpha = .8) + geom_hline(yintercept = log10(.05), lty = 2) + 
  scale_x_continuous(trans = "log10", breaks = c(300, 500, 1000, 2000), limits = c(min(n), max(n))) + 
  scale_color_manual(name = "S",
                     limits = c("{}", "SW", "APAR", "APAR, SW", "GPP", "APAR, GPP", "GPP, SW", "APAR, GPP, SW"),
                     labels = c("{}", "{1}", "{2}", "{1,2}", "{3}", "{2,3}", "{1,3}", "{1,2,3}"),
                     values = c("#999999", "#0072B2", "#009E73", "#56B4E9", "#D55E00", "#F0E442", "#E69F00", "#CC79A7")) +
  xlab("sample size") + ylab(expression(paste(log[10], " p-value for h-invariance"))) + 
  scale_shape_manual(limits = c(TRUE, FALSE), values = c(17,1)) + 
  scale_size_manual(limits = c(TRUE, FALSE), values = c(3,2)) + 
  guides(shape = "none", size = "none") + 
  theme(plot.margin = unit(c(5,5,5,5), "mm"))
p

pval.noncaus <- c()
for(m in unique(frame$n)){
  
  sub <- subset(frame, n == m)
  if(all(sub$p.value < .05)){
    pv <- c(1,1,1)
  } else{
    pv <- c(max(sub$p.value[c(1,3,4,7)]), max(sub$p.value[c(1,2,4,6)]), max(sub$p.value[c(1,2,3,5)]))
  }
  pval.noncaus <- c(pval.noncaus, pv)
}

pval.frame <- data.frame(n = rep(unique(frame$n), each = 3), 
                         var = factor(rep(c("APAR", "GPP", "SW"), length(unique(frame$n)))),
                         p.value = pval.noncaus)

pval.frame$lty <- factor(1+(pval.frame$n>1000))

pval.frame$lty <- factor(1+(pval.frame$n>1300))

pp <- ggplot() + 
  geom_point(data=pval.frame, aes(n, log10(p.value), group = var, col = var)) + 
  geom_line(data=subset(pval.frame, n < 1800), aes(n, log10(p.value), group = var, col = var), size = .7, alpha = .8) + 
  geom_line(data=subset(pval.frame, n > 1300), aes(n, log10(p.value), group = var, col = var, lty = var), size = .7, alpha = .8) + 
  geom_hline(yintercept = log10(.05), lty = 2) + 
  scale_x_continuous(trans = "log10", breaks = c(300, 500, 1000, 2000), limits = c(min(n), max(n))) + 
  scale_color_manual(name = "variable",
                     limits = c( "SW","APAR", "GPP"),
                     labels = c( expression(X^1),expression(X^2), expression(X^3)),
                     values = c( "#0072B2","#009E73", "#D55E00")) + 
  xlab("sample size") + ylab(expression(paste(log[10], " p-value for non-causality"))) + 
  theme(plot.margin = unit(c(5,5,5,5), "mm")) + 
  guides(linetype="none")
pp

pdf("fig11.pdf", width = 12, height = 4)
grid.arrange(p, pp, ncol = 2, widths = c(1.06,1))
dev.off()
 