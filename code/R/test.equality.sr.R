#' A test for the equality of switching regression models
#' 
#' Given several datasets containing observations from a target variable Y
#' and covariates X, this function tests the hypothesis of whether the datasets
#' are generated under the same switching regression model for the conditional
#' distribution of Y given X.
#' 
#' @author Rune Christiansen \email{krunechristiansen@@math.ku.dk}
#' @seealso \code{\link{icph}}
#' 
#' @param Y A numeric vector with observations from the target variable for all data sets.
#' @param X A numeric matrix with the predictor variables for all data sets.
#' @param E The environment indicator, can be numeric or a factor variable. 
#'   This variable indicates which dataset the observations belong to, i.e., 
#'   the \code{i}th entry of \code{E} should be \code{k} if the observation corresponding to the 
#'   \code{i}th entry of \code{Y} and the \code{i}th row of the 
#'   data matrix \code{X} comes from data set \code{k}.  
#' @param number.of.states An integer greater or equal to 2 specifying the number of mixture components. 
#'   Numerical experiments show poor performance of the test for more than 3 components. 
#' @param model The dependence model for the hidden variables. If ``IID'', 
#'   all observations are assumed to be independent. 
#'   If ``HMM'', the hidden variables are assumed to follow a first order Markov chain.
#' @param intercept Logical indicating whether an intercept should 
#'   be included in the linear regression of all components.
#' @param alpha Numeric between 0 and 1 indicating the desired type I error control.
#' @param plot Logical indicating whether diagnostic plots should be returned.#'       
#' @param par.test A list containing the following options
#'   \itemize{
#'     \item{\code{method}:}{ The algorithm to be used for likelihood optimization. 
#'       If ``EM'', the EM-algorithm is used (via the function \code{mixreg} from the package \code{mixtools}). 
#'       If ``NLM'', the Newton-type optimizer \code{nlm} is called. 
#'       Note: the EM-algorithm is currently only implemented for model ``IID''. 
#'     }
#'     \item{\code{optim.n.inits}:}{ Number of different starting values in likelihood optimization. 
#'     }
#'     \item{\code{variance.constraint}:}{ The constraint imposed on the error variance of the mixture components. 
#'       If ``equality'', all error variances are required to be equal. 
#'       If ``lower bound'', a lower bound is imposed (\code{1e-4} for ``NLM'' and \code{1e-16} for ``EM'').
#'     }
#'     \item{\code{test.parameters}:}{ A vector of strings indicating which parameters of the 
#'       switching regression models should be tested for equality. Possible entries are 
#'       ``intercept'' (the intercept terms), ``beta'' (the regression coefficients), 
#'       ``sigma'' (the error variances) and ``gamma'' (the transition matrix of the hidden Markov chain). 
#'       For example, if \code{test.parameters = c}(``beta'', ``sigma''), 
#'       then only the regression coefficients and the error variances of the switching regression models 
#'       are tested for equality, while intercept terms and transition matrices are allowed to vary across different models. 
#'     }
#'     \item{\code{testNoMix}:}{ Logical. If \code{TRUE}, an ordinary linear regression model 
#'       is fit to every data set, and a test for equality of these methods is performed (using the package \code{seqICP}). 
#'       If this test does not reject, its output is returned and the algorithm terminates. 
#'       If this test does reject, its output is disregarded and the test for equality 
#'       of swichting regression models is performed as usual.
#'     }
#'     \item{\code{penalty}:}{ Logical. If \code{TRUE}, a penalty term will be added to the likelihood function, 
#'       penalizing for the number of mixture components. Tuning parameters are automatically selected 
#'       using a BIC criterion. This option is beneficial if not all mixture components are ``active'' 
#'       in all data sets. The method will then automatically detect this, and acknowledge that no information 
#'       on the missing mixture component is available, rather than estimating it to coincide with another component 
#'       (which usually happens in the unpenalized case).
#'     }
#'     \item{\code{init}:}{ Initialization of optimization algorithms. Can be chosen to be one of ``regmix'' 
#'       (data driven starting values automatically calculated using the function \code{regmix.init} from the package \code{mixtools}), 
#'       ``random'' (randomly sampled starting values) and ``true'' (user supplied starting values, see item \code{theta0} for the format).
#'     }
#'     \item{\code{output.pvalue}:}{ Logical. If \code{TRUE}, a p-value is calculated and returned. 
#'       Since p-values are computed by repeated function evaluations along a grid of thresholds, 
#'       this option is time consuming if the number of observations is large. If \code{FALSE}, 
#'       the hypothesis of equality is only tested at level \code{alpha} (without calculating an actual p-value).
#'     }
#'     \item{\code{silent}:}{Logical. If \code{TRUE}, then information on the progression of the algorithm 
#'       is printed to the console.}
#'     \item{\code{theta0}:}{ Either \code{NULL} or a list containing paramter values. 
#'       If parameter values are supplied, 
#'       p-values for hypotheses (one for each supplied data set) of the true values 
#'       being equal to the supplied values, are computed (note: only the parameters specified in 
#'       \code{test.parameters} are included in the hypothesis). The supplied list \code{theta0} contains a list
#'       for every level of \code{E}, each of which contains the following items: 
#'         \code{mu}: a vector of intercept terms (one for each component);
#'         \code{beta}: a matrix of regression coefficients
#'           (each column corresponds to the vector of regression coefficients for one component);
#'         \code{sigma}: a vector of error variances (one for each component);
#'         \code{gamma}: transition matrix (the \code{(i,j)}th entry denotes 
#'           the transition probability from state i to state j).
#'     }
#'   }
#'   
#' @return Object of class \code{test} containing a list with the following items.
#'     \item{estimates}{ A dataframe with estimates and standard errors of all model parameters.}
#'     \item{loglik}{ The likelihood scores for the estimated parameters (one for each data set). 
#'       If \code{!is.null(theta0)}, then also likelihood scores for the supplied parameters are returned.
#'     }
#'     \item{p.value}{ A p-value for the hypothesis of equality of the switching regression models. 
#'       If \code{output.pvalue = FALSE}, this output is of the form \code{>= alpha} or \code{< alpha}.
#'     }
#'     \item{alpha}{ Test level.}
#'     \item{output.pvalue}{ Logical indicating if p-values were calculated.}
#'     \item{p.values.cover}{ If \code{!is.null(theta0)}, this item contains p-values for the hypotheses 
#'       of the true values being equal to the supplied values. If \code{output.pvalue = FALSE}, 
#'       this output is of the form \code{>= alpha} or \code{< alpha}. 
#'     }
#'     \item{reject}{ Logical indicating whether the test rejects at level \code{alpha}.}
#'     \item{model}{ The model that was used (``IID'' or ``HMM'').}
#'     \item{intercept}{ Logical indicating if intercept terms were fitted.}
#'     \item{ne}{ The number of data sets.}
#'     \item{K}{ The number of mixture components.}
#'     \item{d}{ The number of covariates (including a possible intercept term).}
#'   
#' @examples 
#'   set.seed(1)
#'   n <- 600
#'   # construct indicator for several data sets
#'   E <- rep(1:3, each=n/3)
#'   # transition matrices
#'   Gamma1 <- matrix(c(.9, .1, .1, .9), 2,2,byrow=TRUE)
#'   Gamma2 <- matrix(c(.9, .1, .3, .7), 2,2,byrow=TRUE)
#'   Gamma3 <- matrix(c(.5, .5, .5, .5), 2,2,byrow=TRUE)
#'   # hidden variables with different distributions for each data set
#'   H1 <- H2 <- H3 <- rep(1, n/3)
#'   for(t in 2:(n/3)){
#'     H1[t] <- sample(1:2, size=1, prob = Gamma1[H1[t-1],])
#'     H2[t] <- sample(1:2, size=1, prob = Gamma2[H2[t-1],])
#'     H3[t] <- sample(1:2, size=1, prob = Gamma3[H3[t-1],])
#'   }
#'   H <- c(H1,H2,H3)
#'   
#'   X1 <- rnorm(n)
#'   Y1  <- (H==1)*(1 + X1 + .5*rnorm(n)) + (H==2)*(1 + 2*X1 + .7*rnorm(n))
#'   
#'   # The intercept, regression coefficients and error variances 
#'   # for the switching regression models of Y1 given X1 coincide
#'   test1 <- test.equality.sr(Y = Y1, X = X1, E = E)
#'   test1
#'   
#'   X2 <- rnorm(n)
#'   Y2  <- (H==1)*(1 + E*X2 + .5*rnorm(n)) + (H==2)*(E + 2*X2 + .7*rnorm(n))
#'   
#'   # For the switching regression models of Y2 given X2, they do not
#'   test2 <- test.equality.sr(Y = Y1, X = X2, E = E)
#'   test2
#' 
#' @import ggplot2
#' @importFrom  mixtools regmixEM regmix.init normalmixEM 
#' @importFrom  seqICP seqICP.s
#' @importFrom  grid grid.draw 
#' @importFrom  gridExtra arrangeGrob
#' @importFrom  scales trans_new
#' @importFrom stats dnorm fitted ks.test lm nlm qchisq rnorm runif var
#' 
#' @export

test.equality.sr <- function(Y, X, E, number.of.states = 2, 
                             model = "HMM", intercept = TRUE, 
                             alpha = 0.05, plot = TRUE,
                             par.test = list(method = "NLM", 
                                             optim.n.inits = 5,
                                             variance.constraint = "equality",
                                             test.parameters = c("intercept", "beta", "sigma"), 
                                             testNoMix = FALSE, 
                                             penalty = FALSE, 
                                             init = "regmix",
                                             output.pvalue = TRUE,
                                             silent = TRUE,
                                             theta0 = NULL)
                             ){
  
  if(!exists("method", par.test)){
    par.test$method <- "NLM"
  }
  if(!exists("optim.n.inits", par.test)){
    par.test$optim.n.inits <- 5
  }
  if(!exists("variance.constraint", par.test)){
    par.test$variance.constraint <- "equality"
  }
  if(!exists("test.parameters", par.test)){
    par.test$test.parameters <- c("intercept", "beta", "sigma")
  }
  if(!exists("testNoMix", par.test)){
    par.test$testNoMix <- FALSE
  }
  if(!exists("penalty", par.test)){
    par.test$penalty <- 0
  }
  if(!exists("init", par.test)){
    par.test$init <- "regmix"
  }
  if(!exists("output.pvalue", par.test)){
    par.test$output.pvalue <- TRUE
  }
  if(!exists("silent", par.test)){
    par.test$silent <- TRUE
  }
  if(!exists("theta0", par.test)){
    par.test$theta0 <- NULL
  }
  
  K <- sort(number.of.states)
  method <- par.test$method
  optim.n.inits <- par.test$optim.n.inits
  arbvar <- (par.test$variance.constraint == "lower bound")
  testpars <- par.test$test.parameters
  testNoMix <- par.test$testNoMix
  penalty <- par.test$penalty
  init <- par.test$init
  if(init == "true") number.of.inits <- 1
  output.pvalue <- par.test$output.pvalue
  silent <- par.test$silent
  theta0 <- par.test$theta0
  
  X <- as.matrix(X)
  
  ## Input control and warnings
  if(init == "true" & is.null(theta0)){ warning("No true parameters supplied -- using generic starting values"); init <- "regmix"}
  
  if(testNoMix){
    res <- decoupledTest(X, Y, E, intercept, plot)
    print(paste("P-value for decoupled test: ",res$p.value))
    res$reject <- res$p.value<alpha
    if(!(res$reject)){
      print("Test without Mixing does not reject -- result:")
      return(res)
    }
  }
  
  ne  <- length(unique(E))
  d   <- ncol(X) + intercept
  
  p.value <- -1
  
  # run through all supplied values for K, choose model that yields the highest p-value
  for(k in K){
    
    # Object to hold all model fits (MFs)
    MFs.new <- list()
    p.values.cover.new <- mse.new <- rep(NA,ne)
    
    ## Fit mixture models
    if(method == "EM"){
      for(i in 1:ne){
        if(!silent) print(paste0("Fitting environment ", i))
        f <- unique(E)[i]
        w <- which(E==f)
        xf <- X[w,,drop=F]
        yf <- Y[w]
        MFs.new[[i]] <- mle.em(xf, yf, k, theta0[[i]], intercept, arbvar, model, init, rep=optim.n.inits)
        if(!is.null(theta0)){ tmp <- check.cover(MFs.new[[i]], testpars); p.values.cover.new[i] <- tmp$p.value; mse.new[i] <- tmp$mse}
      }
    }
    if(method == "NLM"){
      for(i in 1:ne){
        if(!silent) print(paste0("Fitting environment ", i))
        f <- unique(E)[i]
        w <- which(E==f)
        xf <- X[w,,drop=F]
        yf <- Y[w]
        MFs.new[[i]] <- mle.nlm(x=xf, y=yf, k, theta0=theta0[[i]], intercept, arbvar, model, init, penalty, rep=optim.n.inits)
        if(!is.null(theta0)){ tmp <- check.cover(MFs.new[[i]], testpars); p.values.cover.new[i] <- tmp$p.value; mse.new[i] <- tmp$mse}
      }
    }
    
    if(ne>1){
      ## Compute p-values for all label permutations and return largest p-value and the corresponding permutation
      if(!silent) print("Computing p-value")
      tmp  <- pval.fun(MFs.new, testpars, alpha, output.pvalue)
      p.value.new <- tmp$pval
      reject.new <- tmp$reject
      
      ## Permute model labels according to largest p-value
      for(i in 1:ne){
        MFs.new[[i]] <- permute.labels(MF=MFs.new[[i]], perm=unlist(tmp$perm[[i]]))
      }
    
      # if current number of states yields larger p-value than previous one, adopt it
      if(p.value.new > p.value){ # this will always be the case for first k value, since p.value is initialized at -1
        p.value <- p.value.new
        p.values.cover <- p.values.cover.new
        mse <- mse.new
        reject <- reject.new
        MFs <- MFs.new
        if(p.value >= alpha) break
      }
    }
    if(ne==1){
      p.value <- NA
      reject <- NA
      p.values.cover <- p.values.cover.new
      mse <- mse.new
      MFs <- MFs.new
      break
    }
  }
  
  if(plot) plotfun(MFs, alpha)
  
  out <- output(MFs)
  out$p.value <- p.value
  out$alpha <- alpha
  out$output.pvalue <- output.pvalue
  out$p.values.cover <- p.values.cover
  #out$mse <- mse
  out$reject <- reject
  out$model <- model
  out$intercept <- intercept
  out$ne <- ne
  out$d <- d
  out$K <- MFs[[1]]$K
  
  ## Most probable sequence of latent states
  hhat <- numeric(length(Y))
  for(i in 1:ne){
    f <- unique(E)[i]
    hhat[E==f] <- MFs[[i]]$h
  }
  out$hhat <- hhat
  
  fitted <- numeric(length(Y))
  for(i in 1:ne){
    f <- unique(E)[i]
    fitted[E==f] <- MFs[[i]]$fitted
  }
  out$fitted <- fitted
  #out$MFs <- MFs

  out <- structure(out, class = "test")
  
  return(out)
}

#' @method print test
#' @export
print.test <- function(out){
  
    ne <- out$ne
    K <- out$K
    d <- out$d-out$intercept
    model <- out$model
    
    cat("\n\n")
    if(ne > 1){
      if(!out$reject) cat("Hypothesis of h-invariance not rejected")
      if(out$reject) cat("Hypothesis of h-invariance rejected")
      cat("\n")
      cat(paste("p-value:", ifelse(out$output.pvalue, round(out$p.value, 5), out$p.value)))
    } else{
      cat(paste("p-value:", round(out$p.values.cover, 5)))
    }
    
    cat("\n\n")
    cat("Estimates:")
    cat("\n\n")
    
    for(i in 1:ne){
      cat(paste(" %%%%%%%%%%%%%%%%%%%% environment", i, "%%%%%%%%%%%%%%%%%%%% "))
      cat("\n")
      cat(" intercept: \n") 
      cat(round(subset(out$estimates, parameter == "intercept" & environment == i)$estimate,3))
      cat("\n\n")
      cat(" regression coefficients: \n")
      print(matrix(round(subset(out$estimates, parameter == "beta" & environment == i)$estimate,3), nrow = d, ncol = K, byrow=T))
      cat("\n")
      cat(" error variances: \n")
      cat(round(subset(out$estimates, parameter == "sigma" & environment == i)$estimate^2,3))
      cat("\n\n")
      cat(" transition matrix: \n")
      if(model == "IID"){
        gamma <- matrix(rep(subset(out$estimates, parameter == "gamma" & environment == i)$estimate, each=K), K, K)
      }
      if(model == "HMM"){
        gamma <- matrix(c(sapply(1:K, function(k){
          g <- numeric(K)
          g[-k] <- subset(out$estimates, parameter == "gamma" & environment == i & component == k)$estimate
          g[k] <- 1-sum(g[-k])
          return(g)
        })), K, K, byrow=T)
      }
      print(round(gamma, 3))
      cat("\n")
    }
}

