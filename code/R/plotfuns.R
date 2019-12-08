# Plot functions for diagnostics of test.fisher.cr 
# use unlist instead of c(), since if environments have different lengths, sapply sorts content in list

# library(ggplot2)
# library(grid)
# library(gridExtra)
# library(scales)
# library(reshape2)

plotfun <- function(MFs, alpha){
  
  ne <- length(MFs)
  intercept <- MFs[[1]]$intercept
  d <- MFs[[1]]$d
  K <- MFs[[1]]$K
  
  if(intercept){
    for(i in 1:ne){
      for(k in 1:K){
        MFs[[i]]$theta[[k]]$beta <- rbind(MFs[[i]]$theta[[k]]$intercept, MFs[[i]]$theta[[k]]$beta)
      }
    }
  }
  
  
  if(d==1 & intercept){
    p <- fitplot.intercept(MFs)
    p2 <- betasigmaplot(MFs, alpha/ne)
  }else if((d==2 & intercept) | (d==1 & !(intercept))){
    p <- fitplot.1pred(MFs)
    #p2 <- betaplot(MFs, alpha / ne)
    #p <- arrangeGrob(p1, p2, nrow=2)
  }else{
    p1 <- fitplot(MFs)
    p2 <- betaplot(MFs, alpha / ne)
    p <- gridExtra::arrangeGrob(p1, p2, nrow=2)
  }
  
  grid::grid.draw(p)

}

betaplot <- function(MFs, alpha){
  
  ne <- length(MFs)
  K  <- MFs[[1]]$K
  d  <- MFs[[1]]$d
  df <- (d+1)*K
  
  out.frame <- NULL
  for(i in 1:ne){
    loop.frame <- ellipses(MFs[[i]], alpha)
    loop.frame$environment <- i
    out.frame <- rbind(out.frame, loop.frame)
  }
  out.frame$environment <- factor(out.frame$environment)
  out.frame$component <- factor(out.frame$component)
  
  p <- ggplot2::ggplot(out.frame, aes(x,y)) + 
    geom_point(aes(col = component, shape = environment), alpha = 0.7, size = 2/((d-1))) + 
    geom_path(aes(col = component, group = interaction(environment, component)), alpha = 0.7, size = 1/((d-1))) + 
    facet_grid(parx~pary) + theme(text = element_text(size = 20)) + 
    ggtitle(paste("pairwise confidence regions at level alpha = ", round(alpha,4))) + 
    theme_bw()
    
  
  return(p)
}

# Only ever called when beta=beta0
betasigmaplot <- function(MFs, alpha){
  
  if(length(MFs[[1]]$d) > 1) stop("Works only for d=1")
  
  ne <- length(MFs)
  K  <- MFs[[1]]$K
  d  <- MFs[[1]]$d
  df <- (d+1)*K
  
  out.frame <- NULL
  for(i in 1:ne){
    loop.frame <- ellipses.betasigma(MFs[[i]], alpha)
    loop.frame$environment <- i
    out.frame <- rbind(out.frame, loop.frame)
  }
  out.frame$environment <- factor(out.frame$environment)
  out.frame$component <- factor(out.frame$component)
  
  p <- ggplot2::ggplot(out.frame, aes(x,y)) + 
    geom_point(aes(col = component, shape = environment), alpha = 0.7, size = 2) + 
    geom_path(aes(col = component, group = interaction(environment, component)), alpha = 0.7, size = 1) + 
    theme(text = element_text(size = 20)) + 
    xlab(paste0("beta",1-MFs[[1]]$intercept)) + ylab("sigmasquared") + 
    ggtitle(paste("pairwise confidence regions at level alpha = ", round(alpha,4))) + 
    theme_bw()
  
  
  return(p)
}

ellipses.betasigma <- function(MF=MFs[[1]], alpha){
  
  K  <- MF$K
  d  <- MF$d
  df <- (d+1)*K
  r <- sqrt(qchisq(1-alpha, df=df))
  
  w <- which(MF$parameter != "lambda")
  C <- MF$covmat[w,w]
  
  out.frame <- NULL
  for(k in 1:K){
    CC <- C[(k-1)*(d+1)+1:2, (k-1)*(d+1)+1:2]
    cc <- c(MF$theta[[k]]$beta, MF$theta[[k]]$sigma2)
    loop.frame <- data.frame(x=ellipse(cc, CC, r)$x, 
                             y=ellipse(cc, CC, r)$y,
                             parx = paste0("beta",1-MF$intercept),
                             pary = "sigma2",
                             component = k)
    out.frame <- rbind(out.frame, loop.frame)
  }
  
  return(out.frame)
}

# TODO: when d==1 & !intercept the below code does not work. Alternative: betasigma plot..
ellipses <- function(MF=MFs[[1]], alpha){
  
  K  <- MF$K
  d  <- MF$d 
  df <- (d+1)*K - !MF$arbvar
  r <- sqrt(qchisq(1-alpha, df=df))
  
  w <- which(MF$parameter %in% c("intercept", "beta"))
  C <- MF$covmat[w,w]
  
  out.frame <- NULL
  for(k in 1:K){
    for(i in 1:(d-1)){
      for(j in (i+1):d){
        CC <- C[(k-1)*d+c(i,j), (k-1)*d+c(i,j)]
        cc <- MF$theta[[k]]$beta[c(i,j)] 
        loop.frame <- data.frame(x=ellipse(cc, CC, r)$x, 
                                 y=ellipse(cc, CC, r)$y,
                                 parx = paste0("beta",i-MF$intercept),
                                 pary = paste0("beta",j-MF$intercept),
                                 component = k)
        out.frame <- rbind(out.frame, loop.frame)
      }
    }
  }
  
  return(out.frame)
}

ellipse <- function(c, C, r){
  
  e <- eigen(C)
  S <- e$vectors%*%diag(sqrt(e$values))
  len <- 30
  
  circ.seq <- seq(0,2*pi,length=len)
  circ <- rbind(cos(circ.seq), sin(circ.seq))
  el <- matrix(c,2,len) + r*S%*%circ
  return(list(x=el[1,], y=el[2,]))
}

fitplot <- function(MFs){
  
  ne <- length(MFs)
  N <- MFs[[1]]$N
  
  plot.frame <- data.frame(y = unlist(lapply(1:ne, function(i) MFs[[i]]$y)),
                           fitted = unlist(lapply(1:ne, function(i) MFs[[i]]$fitted)),
                           residuals = unlist(lapply(1:ne, function(i) MFs[[i]]$residuals)),
                           component = factor(unlist(lapply(1:ne, function(i) MFs[[i]]$h))),
                           environment = factor(unlist(lapply(1:ne, function(i) rep(i, MFs[[i]]$N)))))
  
  p <- ggplot2::ggplot(plot.frame, aes(fitted, y)) + 
    geom_point(aes(col = component, shape = environment), alpha = 0.7, size = 2/ne) + 
    facet_grid(.~environment) + theme(text = element_text(size = 20)) + 
    ggtitle("Observations against fitted values using most probable hidden states") + 
    theme_bw()
  
  return(p)
}

fitplot.1pred <- function(MFs){
  
  ne <- length(MFs)
  N <- MFs[[1]]$N
  K <- MFs[[1]]$K
  intercept <- MFs[[1]]$intercept
  
  plot.frame <- data.frame(y = unlist(lapply(1:ne, function(i) MFs[[i]]$y)),
                           x = unlist(lapply(1:ne, function(i) c(MFs[[i]]$x))),
                           component = factor(unlist(lapply(1:ne, function(i) MFs[[i]]$h))),
                           environment = factor(unlist(lapply(1:ne, function(i) rep(i, MFs[[i]]$N)))))
  
  if(intercept){
  beta0 <- c(sapply(1:ne, function(i) c(sapply(1:MFs[[i]]$K, function(j) MFs[[i]]$theta[[j]]$beta[1]))))
  beta1 <- c(sapply(1:ne, function(i) c(sapply(1:MFs[[i]]$K, function(j) MFs[[i]]$theta[[j]]$beta[2]))))
  }
  if(!intercept){
    beta0 <- rep(0,K)
    beta1 <- c(sapply(1:ne, function(i) c(sapply(1:MFs[[i]]$K, function(j) MFs[[i]]$theta[[j]]$beta[1]))))
  }
  
  linedata <- data.frame(environment=rep(1:ne,each=K), h = factor(rep(1:K, ne)), intercept=beta0, slope=beta1)
  
  p <- ggplot2::ggplot(plot.frame, aes(x, y)) + 
    geom_point(aes(shape = environment, col = component), alpha = 0.7, size = 2/ne) + 
    geom_abline(data = linedata, aes(intercept=intercept, slope=slope, col = h)) + 
    facet_grid(.~environment) + theme(text = element_text(size = 20)) + 
    ggtitle("Y vs X -- classified using most probable hidden states") + 
    theme_bw()
  
  return(p)
}

fitplot.intercept <- function(MFs){
  
  ne <- length(MFs)
  N <- MFs[[1]]$N
  K <- MFs[[1]]$K
  
  beta0 <- c(sapply(1:ne, function(i) c(sapply(1:MFs[[i]]$K, function(j) MFs[[i]]$theta[[j]]$beta))))
  
  vlinedata <- data.frame(environment=rep(1:ne,each=K), beta0=beta0)
  
  plot.frame <- data.frame(y = unlist(lapply(1:ne, function(i) MFs[[i]]$y)),
                           h = factor(unlist(lapply(1:ne, function(i) MFs[[i]]$h))),
                           environment = factor(unlist(lapply(1:ne, function(i) rep(i, MFs[[i]]$N)))))
  
  p <- ggplot2::ggplot(plot.frame, aes(x=y)) + geom_histogram(aes(fill=h)) + 
    facet_grid(.~environment) + 
    ggtitle("Marginals of target for each environment -- vertical lines indicating estimated kernel means") + 
    geom_vline(data=vlinedata, aes(xintercept = beta0), col="blue") + 
    theme_bw()
    
  return(p)
}

fitplot.intercept.nomix <- function(y, e){
  
  vlinedata <- data.frame(environment = unique(e), 
                          fitted = sapply(unique(e), function(f) mean(y[e==f])))
  
  plot.frame <- data.frame(y = y,
                           environment = e)
  
  p <- ggplot2::ggplot(plot.frame, aes(x=y)) + geom_histogram() + 
    facet_grid(.~environment) + 
    ggtitle("Marginals of target for each environment -- vertical lines indicating estimated means") + 
    geom_vline(data=vlinedata, aes(xintercept = fitted), col="blue") + 
    theme_bw()
  
  return(p)
}

fitplot.1pred.nomix <- function(x, y, e, intercept){
  
  plot.frame <- data.frame(y = y,
                           x = x,
                           environment = e)
  
  p <- ggplot2::ggplot(plot.frame, aes(x, y)) + geom_point() + 
    facet_grid(.~environment) + geom_smooth(method="lm")
  
  return(p)
}

fitplot.nomix <- function(x, y, e, intercept){
  
  if(intercept) x <- cbind(x, rep(1, nrow(x)))
  fitted <- c(sapply(unique(e), function(f) fitted(lm(y[e==f] ~ -1 + x[e==f,]))))
  
  plot.frame <- data.frame(y = y,
                           fitted = fitted,
                           environment = e)
  
  p <- ggplot2::ggplot(plot.frame, aes(y, fitted)) + geom_point() + 
    facet_grid(.~environment) + 
    ggtitle("Observations against fitted values")
  
  return(p)
}

my_trans <- function(name, up, low, mid, breaks){
  
  pos <- mid-low + .5*(up-low)
  slope1 <- .5*(up-low)/(mid-low)
  intercept1 <- low-slope1*low
  slope2 <- 1/slope1
  intercept2 <- pos-slope2*(.5*(up-low))
  
  scale_trans <- function(x){
    w1 <- which(x <= mid)
    w2 <- which(x > mid & x < pos)
    out <- x
    out[w1] <- slope1*out[w1]+intercept1
    out[w2] <- slope2*out[w2]+intercept2
    return(out)
  }
  
  ascale_trans <- function(x){
    w1 <- which(x <= (up-low)/2)
    w2 <- which(x > (up-low)/2 & x < pos)
    out <- x
    out[w1] <- (out[w1]-intercept1)/slope1
    out[w2] <- (out[w2]-intercept2)/slope2
    return(out)
  }
  
  function() trans_new(name, scale_trans, ascale_trans, breaks = function(x) breaks, domain = c(low,up))
} 

grid_arrange_shared_legend <- function(..., title = NULL, ncol = length(list(...)), nrow = 1, position = c("bottom", "right")) {
  
  plots <- list(...)
  position <- match.arg(position)
  g <- ggplot2::ggplotGrob(plots[[1]] + theme(legend.position = position))$grobs
  legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]
  lheight <- sum(legend$height)
  lwidth <- sum(legend$width)
  gl <- lapply(plots, function(x) x + theme(legend.position="none"))
  gl <- c(gl, ncol = ncol, nrow = nrow)
  
  combined <- switch(position,
                     "bottom" = arrangeGrob(do.call(arrangeGrob, gl),
                                            legend,
                                            top = title,
                                            ncol = 1,
                                            heights = unit.c(unit(.9+.1*is.null(title), "npc") - lheight, lheight)),
                     "right" = arrangeGrob(do.call(arrangeGrob, gl),
                                           legend,
                                           top = title,
                                           ncol = 2,
                                           widths = unit.c(unit(.9+.1*is.null(title), "npc") - lwidth, lwidth)))
  
}
