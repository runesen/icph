## modified pc-algorithm to include background knowledge on possible source nodes
library(pcalg)

jci_pc <- function (suffStat, indepTest, alpha, labels, p=4, fixedGaps = NULL, 
                fixedEdges = NULL, NAdelete = TRUE, m.max = Inf, 
                source.nodes = "x1", 
                u2pd = c("relaxed", "rand", "retry"), 
                skel.method = c("stable", "original", "stable.fast"), 
                conservative = FALSE, maj.rule = FALSE, 
                solve.confl = FALSE, numCores = 1, verbose = FALSE) 
{
  cl <- match.call()
  if (!missing(p)) 
    stopifnot(is.numeric(p), length(p <- as.integer(p)) == 
                1, p >= 2)
  if (missing(labels)) {
    if (missing(p)) 
      stop("need to specify 'labels' or 'p'")
    labels <- as.character(seq_len(p))
  }
  else {
    stopifnot(is.character(labels))
    if (missing(p)) {
      p <- length(labels)
    }
    else if (p != length(labels)) 
      stop("'p' is not needed when 'labels' is specified, and must match length(labels)")
    else message("No need to specify 'p', when 'labels' is given")
  }
  u2pd <- match.arg(u2pd)
  skel.method <- match.arg(skel.method)
  if (u2pd != "relaxed") {
    if (conservative || maj.rule) 
      stop("Conservative PC and majority rule PC can only be run with 'u2pd = relaxed'")
    if (solve.confl) 
      stop("Versions of PC using lists for the orientation rules (and possibly bi-directed edges)\\n can only be run with 'u2pd = relaxed'")
  }
  if (conservative && maj.rule) 
    stop("Choose either conservative PC or majority rule PC!")
  skel <- skeleton(suffStat, indepTest, alpha, labels = labels, 
                   method = skel.method, fixedGaps = fixedGaps, fixedEdges = fixedEdges, 
                   NAdelete = NAdelete, m.max = m.max, numCores = numCores, 
                   verbose = verbose)
  # dev.off()
  # plot(skel@graph)
  source.ind <- which(labels %in% source.nodes)
  for(i in 1:p){
    skel@graph@edgeL[[i]]$edges <- setdiff(skel@graph@edgeL[[i]]$edges, source.ind)
  }
  # if(special){
  #   skel@graph@edgeL$x3$edges <- setdiff(skel@graph@edgeL$x3$edges,1)
  #   skel@graph@edgeL$x2$edges <- setdiff(skel@graph@edgeL$x2$edges,1)
  #   skel@graph@edgeL$x4$edges <- setdiff(skel@graph@edgeL$x4$edges,1)
  #   skel@graph@graphData$edgemode <- "directed"
  # }
  # if(length(sourcenodes)>0) skel@graph@graphData$edgemode <- "directed"
  # # dev.off()
  # plot(skel@graph)
  
  # rm.edges <- sapply(labels[labels!=sourcenode], function(l) paste0(sourcenode, "|", l))
  # skel@graph@edgeData@data[rm.edges] <- NULL
  
  skel@call <- cl
  if (!conservative && !maj.rule) {
    switch(u2pd, rand = udag2pdag(skel), retry = udag2pdagSpecial(skel)$pcObj, 
           relaxed = udag2pdagRelaxed(skel, verbose = verbose, 
                                      solve.confl = solve.confl))
  }
  else {
    pc. <- pc.cons.intern(skel, suffStat, indepTest, alpha, 
                          version.unf = c(2, 1), maj.rule = maj.rule, verbose = verbose)
    udag2pdagRelaxed(pc.$sk, verbose = verbose, unfVect = pc.$unfTripl, 
                     solve.confl = solve.confl)
  }
}
