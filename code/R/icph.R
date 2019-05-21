#source("test.equality.sr.R", chdir = TRUE)
#library(methods)

#' Invariant Causal Prediction in the presence of latent variables
#' 
#' Infers the set of causal parents of a target variable Y among a vector of covariates X. The method allows for the presence of discrete latent variables acting on Y.
#' @author Rune Christiansen \email{krunechristiansen@@math.ku.dk}
#' @seealso \code{\link{test.equality.sr}}
#' 
#' @param D A numeric \code{n} times \code{d} data matrix. Each column corresponds to a variable, each row to an observation. 
#' @param E The environment indicator, can be numeric or a factor variable. 
#'   This variable indicates which environment or interventional setting the observations belong to, i.e., 
#'   the \code{i}'th entry of \code{E} should be \code{k} if the observation in the 
#'   \code{i}'th row of the data matrix \code{D} comes from environment \code{k}.  
#' @param target Integer between 1 and d indicating which of the variables 
#'   in the data matrix \code{D} should be taken as target variable.
#' @param number.of.states An integer greater or equal to 2 specifying the number of mixture components. 
#'   Numerical experiments show poor performance of the test for more than 3 components. 
#' @param model The dependence model for the hidden variables. If ``IID'', 
#'   all observations are assumed to be independent. 
#'   If ``HMM'', the hidden variables are assumed to follow a first order Markov chain.
#' @param intercept Logical indicating whether an intercept should 
#'   be included in the linear regression of all components.
#' @param alpha Numeric between 0 and 1 indicating the desired type I error control.
#' @param par.icp A list containing the following options. 
#'   \itemize{
#'     \item{\code{max.parents}:}{ Integer between 1 and \code{d-1}, indicating the maximal size 
#'       of parent sets that should be tested for.}
#'     \item{\code{silent}:}{ Logical. If \code{TRUE}, then information on the progression of the algorithm 
#'   is printed to the console.}
#'     \item{\code{stopIfEmpty}:}{ Logical. If \code{FALSE}, then the algorithm terminates if the 
#'   empty set is accepted.}
#'     \item{\code{output.pvalues}:}{ Logical. If \code{TRUE}, p-values for the individual hypotheses 
#'   will be calculated and returned. Since p-values are computed by repeated function evaluations 
#'   along a grid of thresholds, this option is time consuming if the data matrix \code{D} is large. 
#'   If \code{FALSE}, the individual hypotheses are only tested at level \code{alpha} 
#'   (without calculating an actual p-value)}
#'   }
#' @param par.test A list of options to be passed to the testing procedure. 
#'       See documentation of \link{test.equality.sr} for content and default values.
#'   
#' @return Object of class \code{icph} containg a list with the following items. 
#'     \item{parent.set}{ A subset of the indices from 1 to \code{d} (except \code{target}) 
#'       indicating the estimated set of observable causal parents.}
#'     \item{pvalues}{ A data frame with the results for the individual hypothesis tests. 
#'       If \code{output.pvalues = TRUE}, the p-values are returned. 
#'       Otherwise, an indicator for the rejection/accept of the individual hypothesis is returned.
#'     }
#' 
#' 
#' @examples  
#'   set.seed(1)
#'   n <- 600
#'   # Binary hidden variable that follows an order 1 Markov chain
#'   Gamma <- matrix(c(.9, .1, .1, .9),2,2, byrow=TRUE)
#'   H <- rep(1, n)
#'   for(t in 2:n){
#'     H[t] <- sample(1:2, size=1, prob = Gamma[H[t-1],])
#'   }
#'   # SCM with assignments that change between every observation
#'   # True set of observable parents: S*={3}
#'   X4 <- H * rnorm(n)
#'   X3 <- H + cos((1:n)*2*pi/n)*X4 + .5*rnorm(n)
#'   Y <- (H==1)*(1 + X3 + .5*rnorm(n)) + (H==2)*(1 + 2*X3 + .7*rnorm(n))
#'   X2 <- .5 - X3 + exp(-3*(1:n)/n)*Y + .7*rnorm(n)
#'   # Data matrix
#'   D <- cbind(Y, X2, X3, X4)
#'   # Construct environments from time-ordering of data
#'   E <- rep(1:3, each=n/3)
#'   # Run icph
#'   res.icph <- icph(D=D, E=E, target=1)
#'   res.icph
#'   
#' @import methods
#' @importFrom utils combn
#'   
#' @export 

icph <- function(D, E, target,
                 number.of.states = 2,
                 model = "HMM", 
                 intercept = TRUE,
                 alpha = 0.05,
                 par.icp = list(silent = FALSE, 
                                stopIfEmpty = TRUE, 
                                output.pvalues = TRUE), 
                 par.test = list()){

  if(!exists("max.parents", par.icp)){
    par.icp$max.parents=ncol(D)-1
  }
  if(!exists("silent", par.icp)){
    par.icp$silent = FALSE
  }
  if(!exists("stopIfEmpty", par.icp)){
    par.icp$stopIfEmpty = TRUE
  }
  if(!exists("output.pvalues", par.icp)){
    par.icp$output.pvalues=TRUE
  }
  if(!exists("testNoMix", par.test)){
    par.test$testNoMix=FALSE
  }
  if(!exists("silent", par.test)){
    par.test$silent=TRUE
  }
  
  max.parents <- par.icp$max.parents
  silent <- par.icp$silent
  stopIfEmpty <- par.icp$stopIfEmpty
  par.test$output.pvalue <- par.icp$output.pvalues
  output.estimates <- !par.test$testNoMix
  
  n <- nrow(D)
  d <- ncol(D)
  n.vars <- d-1
  

  # define variable that collects all variables except target
  allS <- (1:d)[-target]
  
  ###
  # Iterate over potential parent sets
  ###
  
  # initialize variables
  parent.set <- allS
  len <- 1
  ind <- 1
  S <- list(numeric(0))

  # reorder D matrix according to S and target
  sorted_D <- D[,c(target, S[[1]]), drop=F]
  
  # check wheter model is invariant over the empty set S={}
  if(!silent){
    print("###################################################")
    print("Fitting set:")
    print(S[[1]])
  }
  
  res <- test.equality.sr(Y = sorted_D[,1],
                          X = sorted_D[,-1],
                          E = E,
                          number.of.states = number.of.states,
                          model = model,
                          intercept = intercept,
                          alpha = alpha, 
                          plot = FALSE,
                          par.test = par.test)

  if(!silent) print(paste("p-value: ", res$p.value))
  pvalues <- data.frame(ind=1,
                        p.value=res$p.value,
                        reject = res$reject)
  if(output.estimates){
    estimates <- res$estimates
    estimates$ind <- 1
  }
  
  # compute intersection of rejected sets
  if(!res$reject){
    parent.set <- numeric(0)
  }
  empty.int <- length(parent.set)==0
  # iterate over all sets until intersection is empty or max.parents is reached
  while(len <= max.parents & !(empty.int & stopIfEmpty)){
    # list of all subsets of 1:n.vars of length len
    tmp <- combn(n.vars, len, simplify=FALSE)
    for(k in 1:length(tmp)){
      S <- append(S, list(allS[tmp[[k]]]))
    }
    while(ind < length(S) & !(empty.int & stopIfEmpty)){
      ind <- ind+1
      # reorder D matrix according to S[[ind]] and target
      sorted_D <- D[,c(target, S[[ind]]), drop=F]
      # check wheter model is invariant over the set S[[ind]]
      if(!silent){
        print("###################################################")
        print("Fitting set:")
        print(S[[ind]])
      }
      
      res <- test.equality.sr(Y = sorted_D[,1],
                              X = sorted_D[,-1],
                              E = E,
                              number.of.states = number.of.states,
                              model = model,
                              intercept = intercept,
                              alpha = alpha, 
                              plot = FALSE,
                              par.test = par.test)
      
      if(!silent) print(paste("p-value: ", res$p.value))
      # read out results
      pvalues   <- rbind(pvalues, data.frame(ind=ind, p.value=res$p.value, reject = res$reject))
      if(output.estimates){
        tmp       <- res$estimates
        tmp$ind   <- ind
        estimates <- rbind(estimates, tmp)
      }    
      
      
      # compute intersection of rejected sets
      if(!res$reject){
        parent.set <- intersect(parent.set, S[[ind]])
      }
      empty.int <- length(parent.set)==0
    }
    len <- len+1
  }

  # check if all tests where rejected
  if(sum(pvalues$reject) == length(pvalues$reject)){
    empty.int <- TRUE
    parent.set <- numeric(0)
  }
  if(!exists("model", par.test)){
    par.test$model <- "independent"
  }
  pvalues$S <- sapply(1:length(S), function(i) toString(S[[i]]))
  rownames(pvalues) <- NULL
  
  out <- list(parent.set=parent.set,
              pvalues=pvalues)
  if(output.estimates){
    rownames(estimates) <- NULL
    out$estimates <- estimates
  }
  
  out <- structure(out, class = "icph")
  return(out)
}

#' @method print icph
#' @export
print.icph <- function(out){
  cat("Estimated parent set:")
  cat("\n")
  print(out$parent.set)
  cat("\n")
  cat("p-values:")
  cat("\n")
  print(out$pvalues[,c(4,2,3)])
}