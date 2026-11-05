# The vendored, modified grplasso (Meier, van de Geer & Buhlmann 2008; the
# fork of grplasso::grplasso.default with a scaffold and two deviation
# multipliers `mscale`/`dscale`), moved here from R/: nothing in the package calls it, and its one remaining role is a
# second solver for two tests, G1 (it reproduces fit_prnode_path()) and F6(a)
# (the four-block leak record of Table 2 of the derivation notes). It is not tested itself.
#
# Only grplasso.default() and what it calls are kept: the S4 classes of the
# model and the control (defined in this helper's environment, `where =
# environment()`, so the package namespace is never touched), grpl.model(),
# LinReg(), grpl.control(), blockstand() and residualize_columns(). The
# formula interface, lambdamax(), the S3 methods (print, plot, predict,
# fitted, `[`), LogReg() and PoissReg() left with R/; the last package copy is
# at commit 2a69b09 (R/grplasso.R, helpers.R, methods.R, family.R, control.R).
# Code below is verbatim from there except the two `where =` arguments.

setClass("grpl.model", representation = representation(
    invlink          = "function",
    link             = "function",                     
    nloglik          = "function",                   
    ngradient        = "function",
    nhessian         = "function",
    check            = "function",
    name             = "character",
    comment          = "character"
), where = environment())

grpl.model <- function(invlink, link, nloglik, ngradient, nhessian,
                       check, name = "user-specified",
                       comment = "user-specified"){
  ## Purpose: Generates models to be used for the Group Lasso algorithm.
  ## ----------------------------------------------------------------------
  ## Arguments:
  ## invlink: a function with arguments "eta" implementing the inverse link
  ##          function.
  ## link:    a function with arguments "mu" implementing the link
  ##          function.
  ## nloglik: a function with arguments "y", "mu" and "weights"
  ##          implementing the negative log-likelihood function.
  ## ngradient: a function with arguments "x", "y", "mu" and "weights" 
  ##            implementing the negative gradient of the log-likelihood
  ##            function.
  ## nhessian: a function with arguments "x", "mu" and "weights"
  ##           implementing the negative} hessian of the log-likelihood
  ##           function.
  ## name: a character name
  ## comment: a character comment
  ## ----------------------------------------------------------------------
  ## Author: Lukas Meier, Date:  1 Jun 2006, 10:12

    RET <- new("grpl.model",
               invlink   = invlink,
               link      = link,
               nloglik   = nloglik,
               ngradient = ngradient,
               nhessian  = nhessian,
               check     = check,
               name      = name,
               comment   = comment)
    RET
}


## Linear Regression
LinReg <- function(){
  grpl.model(invlink  = function(eta) eta,
             link  = function(mu) mu,
             nloglik   = function(y, eta, weights, ...)
               sum(weights * (y - eta)^2),
             ngradient = function(x, y, mu, weights, ...)
               -2 * crossprod(x, weights * (y - mu)),
             nhessian  = function(x, mu, weights, ...)
               2 * crossprod(x, weights * x),
             check     = function(y) TRUE,
             name      = "Linear Regression Model",
             comment   = "Use update.hess=\"lambda\" in grpl.control because the Hessian is constant")
}

setClass("grpl.control",
         representation = representation(
           save.x       = "logical",
           save.y       = "logical",
           update.hess  = "character",
           update.every = "numeric",
           inner.loops  = "numeric",
           line.search  = "logical",
           max.iter     = "numeric",
           tol          = "numeric",
           lower        = "numeric",
           upper        = "numeric",
           beta         = "numeric",
           sigma        = "numeric",
           trace        = "numeric"),

         prototype = list(
           save.x       = FALSE,
           save.y       = TRUE,
           update.hess  = "lambda",
           update.every = 3,
           inner.loops  = 10,
           line.search  = TRUE,
           max.iter     = 500,
           tol          = 5 * 10^-8,
           lower        = 10^-2,
           upper        = 10^9,
           beta         = 0.5,
           sigma        = 0.1,
           trace        = 1),
           
         validity = function(object){
           if(ceiling(object@update.every) != floor(object@update.every) |
              object@update.every <= 0)
             return("update.every has to be a natural number")

           if(ceiling(object@inner.loops) != floor(object@inner.loops) |
              object@inner.loops < 0)
             return("inner.loops has to be a natural number or 0")

           if(ceiling(object@max.iter) != floor(object@max.iter) |
              object@max.iter <= 0)
             return("inner.loops has to be a natural number or greater than 0")

           if(object@beta <= 0 | object@beta >= 1)
             return("beta has to be in (0, 1)")
           
           if(object@sigma <= 0 | object@sigma >= 1)
             return("sigma has to be in (0, 1)")

           if(object@tol <= 0)
             return("tol has to be positive")

           if(object@lower > object@upper)
             return("lower <= upper has to hold")

           return(TRUE)
         },
         where = environment()
)

grpl.control <- function(save.x = FALSE, save.y = TRUE,
                         update.hess = c("lambda", "always"),
                         update.every = 3, inner.loops = 10,
                         line.search = TRUE, max.iter = 500,
                         tol = 5 * 10^-8, lower = 10^-2, upper = Inf,
                         beta = 0.5, sigma = 0.1, trace = 1){
  
  ## Purpose: Options for the Group Lasso Algorithm
  ## ----------------------------------------------------------------------
  ## Arguments:
  ## save.x: a logical indicating whether the design matrix should be saved.
  ## save.y: a logical indicating whether the response should be saved.
  ## update.hess: should the hessian be updated in each
  ##              iteration ("always")? update.hess = "lambda" will update
  ##              the Hessian once for each component of the penalty
  ##              parameter "lambda" based on the parameter estimates
  ##              corresponding to the previous value of the penalty
  ##              parameter.
  ## inner.loops: how many loops should be done (at maximum) when solving
  ##              only the active set (without considering the remaining
  ##              predictors)
  ## tol: convergence tolerance; the smaller the more precise, see
  ##      details below.
               
  ## lower: lower bound for the diagonal approximation of the
  ##        corresponding block submatrix of the Hessian of the negative
  ##        log-likelihood function.
  ## upper: upper bound for the diagonal approximation of the
  ##        corresponding block submatrix of the Hessian of the negative
  ##        log-likelihood function.
  ## beta: scaling factor beta < 1 of the Armijo line search.
  ## sigma: 0 < \sigma < 1 used in the Armijo line search.
  ## trace: integer. "0" omits any output,
  ##        "1" prints the current lambda value,
  ##        "2" prints the improvement in the objective function after each
  ##        sweep through all the parameter groups and additional
  ##        information.
  ## ----------------------------------------------------------------------
  ## Author: Lukas Meier, Date:  1 Jun 2006, 10:02

  
  update.hess <- match.arg(update.hess)

  RET <- new("grpl.control",
             save.x       = save.x,
             save.y       = save.y,
             update.hess  = update.hess,
             update.every = update.every,
             inner.loops  = inner.loops,
             line.search  = line.search,
             max.iter     = max.iter,
             tol          = tol,
             lower        = lower,
             upper        = upper,
             beta         = beta,
             sigma        = sigma,
             trace        = trace)
  RET
}

blockstand <- function(x, ipen.which, inotpen.which)
{
  ## Purpose:
  ## ----------------------------------------------------------------------
  ## Arguments:
  ## ----------------------------------------------------------------------
  ## Author: Lukas Meier, Date:  4 Aug 2006, 08:50
  
  n <- nrow(x)
  x.ort <- x
  scale.pen <- list(); length(scale.pen) <- length(ipen.which)
  scale.notpen <- NULL
  
  if(length(inotpen.which) > 0){
    one <- rep(1, n)
    scale.notpen <- sqrt(drop(one %*% (x[,inotpen.which]^2)) / n)
    x.ort[,inotpen.which] <- scale(x[,inotpen.which], FALSE, scale.notpen)
  }
    
  for(j in 1:length(ipen.which)){
    ind <- ipen.which[[j]]
    decomp <- qr(x[,ind])
    if(decomp$rank < length(ind)) ## Warn if block has not full rank
      stop("Block belonging to columns ", paste(ind, collapse = ", "),
              " has not full rank! \n")
    scale.pen[[j]] <- qr.R(decomp) * 1 / sqrt(n)
    x.ort[,ind] <- qr.Q(decomp) * sqrt(n)
  }
  list(x = x.ort, scale.pen = scale.pen, scale.notpen = scale.notpen)
}

residualize_columns <- function(x, index)
{
  ## Purpose: replace each penalized column (index not NA) by its residual on
  ##          the non-penalized columns (index NA); used by grplasso.default,
  ##          lambdamax.default and the tuning code, so all of them
  ##          scale penalties the same way
  ## Returns: list(x = residualized design, coef = least-squares coefficients
  ##          of the penalized columns on the non-penalized ones)
  notpen <- which(is.na(index))
  pen    <- which(!is.na(index))
  qr.notpen <- qr(x[, notpen, drop = FALSE])
  if(qr.notpen$rank < length(notpen))
    stop("Non-penalized columns are not of full rank; cannot residualize")
  coef <- qr.coef(qr.notpen, x[, pen, drop = FALSE])
  x[, pen] <- qr.resid(qr.notpen, x[, pen, drop = FALSE])
  list(x = x, coef = coef)
}

grplasso.default <- function(x, y, index, M_groups, D_groups, weights = rep(1, length(y)),
                             offset = rep(0, length(y)), lambda, dscale = 1, mscale = 1,
                             coef.init = rep(0, ncol(x)),
                             penscale = sqrt, model = LogReg(), center = TRUE,
                             standardize = TRUE, residualize = FALSE, control = grpl.control(), ...)
{
  ## Purpose: Function to fit a solution (path) of a group lasso problem
  ## ----------------------------------------------------------------------
  ## Arguments:
  ## x: design matrix (including intercept), already rescaled and
  ##    possibly blockwise orthonormalized.
  ## y: response vector
  ## index: vector which defines the grouping of the variables. Components
  ##        sharing the same number build a group. Non-penalized
  ##        coefficients are marked with "NA".
  ## weights: vector of observation weights.
  ## offset: vector of offset values; needs to have the same length as the
  ##         response vector.
  ## lambda: vector of penalty parameters. Optimization starts with the
  ##         first component. See details below.
  ## coef.init: initial vector of parameter estimates corresponding to the
  ##            first component in the vector "lambda".
  ## penscale: rescaling function to adjust the value of the penalty
  ##           parameter to the degrees of freedom of the parameter group.
  ##           See the reference below.
  ## model: an object of class "grpl.model" implementing the negative
  ##        log-likelihood, gradient, hessian etc. See the documentation
  ##        of "grpl.model" for more details.
  ## residualize: if TRUE, each penalized column is replaced by its residual on
  ##              the non-penalized columns (after centering, before
  ##              standardizing), so a group's penalty is scaled by its
  ##              variation beyond what the non-penalized terms can absorb.
  ##              This only reparametrizes the non-penalized coefficients:
  ##              fitted values and penalized coefficients keep their usual
  ##              meaning, and the returned coefficients are back-transformed
  ##              to the original parametrization. As with standardization,
  ##              coef.init, norms.pen and ngradient refer to the internal
  ##              (transformed) scale.
  ## control: options for the fitting algorithm, see "grpl.control".
  ## ...: additional arguments to be passed to the functions defined
  ##      in "model".
  ## ----------------------------------------------------------------------
  ## Author: Lukas Meier, Date: 30 Aug 2005, 09:02

  ## Do some error checking first

  ## Check the design matrix
  if(!is.matrix(x))
    stop("x has to be a matrix")

  if(any(is.na(x)))
    stop("Missing values in x not allowed!")

  ## Check the response
  if(!is.numeric(y))
    stop("y has to be of type 'numeric'")

  if(NROW(x) != length(y))
    stop("x and y have not correct dimensions")

  if(!model@check(y))
    stop("y has wrong format")

  ## Check the other arguments
  if(length(weights) != length(y))
    stop("length(weights) not equal length(y)")

  if(any(weights < 0))
    stop("Negative weights not allowed")

  if((!isTRUE(all.equal(weights, rep(weights[1], length(y))))) &
     (center | standardize))
    warning("Weights not considered for centering/scaling at the moment...")

  if(length(offset) != length(y))
    stop("length(offset) not equal length(y)")

  if(length(coef.init) != ncol(x))
    stop("length(coef.init) not equal ncol(x)")

  if(!is.numeric(index))
    stop("index has to be of type 'numeric'!")

  if(length(index) != ncol(x))
    stop("length(index) not equal ncol(x)!")

  if(is.unsorted(rev(lambda)))
    warning("lambda values should be sorted in decreasing order")

  if(all(is.na(index)))
    stop("None of the predictors are penalized.")

  if(residualize){
    if(!any(is.na(index))){
      warning("No non-penalized predictors; setting residualize = FALSE")
      residualize <- FALSE
    }else if(!isTRUE(all.equal(weights, rep(weights[1], length(y)))))
      stop("residualize = TRUE requires equal observation weights")
  }

  check <- validObject(control) ## will stop the program if error occurs

  ## Extract the control information
  update.hess  <- control@update.hess
  update.every <- control@update.every
  inner.loops  <- control@inner.loops
  line.search  <- control@line.search
  max.iter     <- control@max.iter
  lower        <- control@lower
  upper        <- control@upper
  save.x       <- control@save.x
  save.y       <- control@save.y
  tol          <- control@tol
  trace        <- control@trace
  beta         <- control@beta
  sigma        <- control@sigma

  nrlambda <- length(lambda)
  ncolx    <- ncol(x)
  nrowx    <- nrow(x)

  if(nrlambda > 1 & update.hess == "always"){
    warning("More than one lambda value and update.hess = \"always\". You may want to use update.hess = \"lambda\"")
  }

  ## For the linear model, the Hessian is constant and has hence to be
  ## computed only *once*

  if(model@name == "Linear Regression Model"){
    if(update.hess != "lambda"){
      update.hess <- "lambda"
      if(trace >= 1)
        cat("Setting update.hess = 'lambda'\n")
    }
    if(update.every <= length(lambda)){
      update.every <- length(lambda) + 1
      if(trace >= 1)
        cat("Setting update.every = length(lambda) + 1\n")
    }
  }

  ## Which are the non-penalized parameters?
  any.notpen    <- any(is.na(index))
  inotpen.which <- which(is.na(index))
  nrnotpen      <- length(inotpen.which)

  intercept.which <- which(apply(x == 1, 2, all))
  has.intercept   <- length(intercept.which)

  if(!has.intercept & center){
    message("Couldn't find intercept. Setting center = FALSE.")
    center <- FALSE
  }

  if(length(intercept.which) > 1)
    stop("Multiple intercepts!")

  if(has.intercept)
    has.intercept.notpen <- is.na(index[intercept.which])
  else
    has.intercept.notpen <- FALSE

  others.notpen   <- nrnotpen - has.intercept.notpen
  notpen.int.only <- has.intercept.notpen & !others.notpen

  if(has.intercept & !center & standardize)
    warning("Are you sure that you don't want to perform centering in a model with intercept and standardized predictors?")

  ##if(center & others.notpen)
  if(others.notpen & !residualize)
    warning("Penalization not adjusted to non-penalized predictors.")

  ## Index vector of the penalized parameter groups
  if(any.notpen){
    ipen <- index[-inotpen.which]
    ipen.which <- split((1:ncolx)[-inotpen.which], ipen)
  }else{
    if(has.intercept)
      warning("All groups are penalized, including the intercept.")
    ipen <- index
    ipen.which <- split((1:ncolx), ipen)
  }

  nrpen    <- length(ipen.which)
  dict.pen <- sort(unique(ipen))

  ## M_groups and D_groups are matched to the penalized groups by *position*
  ## among the sorted labels, so the labels must be exactly 1:(number of
  ## groups); any other labelling would silently mis-assign the multipliers
  ## (test plan Z10, G7)
  if(!isTRUE(all.equal(as.numeric(dict.pen), as.numeric(seq_along(dict.pen)))))
    stop("the penalized group labels of index must be seq_along(unique(index)), i.e. 1:(3p) for the ",
         "three-block design, because M_groups and D_groups are matched by position; got labels ",
         paste(dict.pen, collapse = ", "))

  ## Table of degrees of freedom
  ipen.tab   <- table(ipen)[as.character(dict.pen)]
  x.old <- x

  if(center){
    if(!has.intercept) ## could be removed; already handled above
      stop("Need intercept term when using center = TRUE")

    mu.x                 <- apply(x[,-intercept.which], 2, mean)
    x[,-intercept.which] <- sweep(x[,-intercept.which], 2, mu.x)
  }

  ## Residualize the penalized columns on the non-penalized ones, so that each
  ## group's penalty only measures what the non-penalized terms cannot absorb
  if(residualize){
    ipen.cols  <- (1:ncolx)[-inotpen.which]
    resid      <- residualize_columns(x, index)
    x          <- resid$x
    resid.coef <- resid$coef
  }

  ## Standardize the design matrix -> blockwise orthonormalization
  if(standardize){
    ##warning("...Using standardized design matrix.\n")
    stand        <- blockstand(x, ipen.which, inotpen.which)
    x            <- stand$x
    scale.pen    <- stand$scale.pen
    scale.notpen <- stand$scale.notpen
  }
  ## From now on x is the *normalized* design matrix!

  ## Extract the columns into lists, works faster for large matrices
  if(any.notpen){
    x.notpen <- list(); length(x.notpen) <- nrnotpen
    for(i in 1:length(inotpen.which))
      x.notpen[[i]] <- x[,inotpen.which[[i]], drop = FALSE]
  }

  x.pen <- list(); length(x.pen) <- length(nrpen)
  for(i in 1:length(ipen.which))
    x.pen[[i]] <- x[,ipen.which[[i]], drop = FALSE]

  ## Extract the needed functions
  check     <- validObject(model)
  invlink   <- model@invlink
  nloglik   <- model@nloglik
  ngradient <- model@ngradient
  nhessian  <- model@nhessian

  coef      <- coef.init
  coef.pen  <- coef.init
  if(any.notpen)
    coef.pen  <- coef[-inotpen.which]

  norms.pen    <- c(sqrt(rowsum(coef.pen^2, group = ipen)))

  norms.pen.m  <- matrix(0, nrow = nrpen, ncol = nrlambda,
                         dimnames = list(NULL, lambda))
  norms.npen.m <- matrix(0, nrow = nrnotpen, ncol = nrlambda,
                         dimnames = list(NULL, lambda))
  nloglik.v <- fn.val.v <- numeric(nrlambda)
  coef.m    <- grad.m   <- matrix(0, nrow = ncolx, ncol = nrlambda,
                                  dimnames = list(colnames(x), lambda))
  fitted    <- linear.predictors <- matrix(0, nrow = nrowx, ncol = nrlambda,
                                           dimnames = list(rownames(x), lambda))

  converged <- rep(TRUE, nrlambda)

  ## *Initial* vector of linear predictors (eta) and transformed to the
  ## scale of the response (mu)
  eta <- offset + c(x %*% coef)
  mu <- invlink(eta)

  ## Create vectors for the Hessian approximations
  if(any.notpen){
    nH.notpen <- numeric(nrnotpen)
  }
  nH.pen <- numeric(nrpen)

  for(pos in 1:nrlambda){
    l <- lambda[pos]

    if(trace >= 2)
      cat("\nLambda:", l, "\n")

    ## Initial (or updated) Hessian Matrix of the *negative* log-likelihood
    ## function (uses parameter estimates based on the last penalty parameter
    ## value)

    if(update.hess == "lambda" & pos %% update.every == 0 | pos == 1){
      ## Non-penalized groups
      if(any.notpen){
        for(j in 1:nrnotpen){ ## changed
          Xj <- x.notpen[[j]]
          nH.notpen[j] <- min(max(nhessian(Xj, mu, weights, ...), lower), upper)
        }
      }
      ## Penalized groups
      for(j in 1:nrpen){
        ind <- ipen.which[[j]]
        Xj  <- x.pen[[j]]
        diagH <- numeric(length(ind))
        for(i in 1:length(ind)){
          diagH[i] <- nhessian(Xj[, i, drop = FALSE], mu, weights, ...)
        }
        nH.pen[j] <- min(max(diagH, lower), upper)
      }
    }

    ## Start the optimization process
    fn.val <- nloglik(y, eta, weights, ...) +
      l * sum(penscale(ipen.tab[-c(M_groups, D_groups)]) * norms.pen[-c(M_groups, D_groups)]) +
      mscale * l * sum(penscale(ipen.tab[M_groups]) * norms.pen[M_groups]) +
      dscale * l * sum(penscale(ipen.tab[D_groups]) * norms.pen[D_groups])

    ## These are needed to get into the while loop the first time
    do.all <- FALSE
    d.fn   <- d.par <- 1

    counter    <- 1 ## Count the sub-loops
    iter.count <- 0 ## Count the loops through *all* groups

    ## Stop the following while loop if the convergence criterion is fulfilled
    ## but only if we have gone through all the coordinates

    ##while(d.fn > tol | d.par > sqrt(tol) | !do.all){
    while(d.fn > tol | d.par > sqrt(tol) | !do.all){
      ## Escape loop if maximal iteration reached
      if(iter.count >= max.iter){
        converged[pos] <- FALSE
        warning(paste("Maximal number of iterations reached for lambda[", pos,
                      "]", sep = ""))
        break
      }

      ## Save the parameter vector and the function value of the previous step
      fn.val.old <- fn.val
      coef.old   <- coef

      ## Check whether we have some useful information from the previous step

      ## Go through all groups if counter == 0 or if we have exceeded the
      ## number of inner loops (inner.loops)
      if(counter == 0 | counter > inner.loops){
        do.all <- TRUE
        guessed.active <- 1:nrpen
        counter <- 1
        if(trace >= 2)
          cat("...Running through all groups\n")
      }else{## Go through the groups which were identified at the previous step
        guessed.active <- which(norms.pen != 0)
        if(length(guessed.active) == 0){
          guessed.active <- 1:nrpen
          do.all <- TRUE
          if(trace >= 2)
            cat("...Running through all groups\n")
        }else{
          do.all <- FALSE
          if(counter == 1 & trace >= 2)
            cat("...Starting inner loop\n")
          counter <- counter + 1
        }
      }
      if(do.all)
        iter.count <- iter.count + 1

      ## These are used for the line search, start at initial value 1
      ## They are currently here for security reasons
      start.notpen <- rep(1, nrnotpen)
      start.pen    <- rep(1, nrpen)

      if(any.notpen){
        ## Optimize the *non-penalized* parameters
        for(j in 1:nrnotpen){
          ind <- inotpen.which[j]
          Xj  <- x.notpen[[j]]

          ## Gradient of the negative log-likelihood function
          ngrad <- c(ngradient(Xj, y, mu, weights, ...))

          ## Update the Hessian if necessary
          if(update.hess == "always"){
            nH <- min(max(nhessian(Xj, mu, weights, ...), lower), upper)
          }else{
            nH <- nH.notpen[j]
          }

          ## Calculate the search direction
          d <- -(1 / nH) * ngrad
          ## Set to 0 if the value is very small compared to the current
          ## coefficient estimate

          d <- zapsmall(c(coef[ind], d), digits = 16)[2]


          ## If d != 0, we have to do a line search
          if(d != 0){
            scale <- min(start.notpen[j] / beta, 1) ##1
            coef.test      <- coef
            coef.test[ind] <- coef[ind] + scale * d

            Xjd       <- Xj * d
            eta.test  <- eta + Xjd * scale

            if(line.search){
              qh    <- sum(ngrad * d)

              fn.val0     <- nloglik(y, eta, weights, ...)
              fn.val.test <- nloglik(y, eta.test, weights, ...)

              qh <- zapsmall(c(qh, fn.val0), digits = 16)[1]

              ## Armijo line search. Stop if scale gets too small (10^-30).
              while(fn.val.test > fn.val0 + sigma * scale * qh & scale > 10^-30){
                ##cat("Doing line search (nonpen)\n")
                scale          <- scale * beta
                coef.test[ind] <- coef[ind] + scale * d
                eta.test       <- eta + Xjd * scale
                fn.val.test    <- nloglik(y, eta.test, weights, ...)
              }
            } ## end if(line.search)
            if(scale <= 10^-30){ ## Do nothing in that case
              #cat("Running into problems with scale\n")
              #cat("qh", qh, "d", d, "\n")
              ## coef.test <- coef
              ## eta.test  <- eta
              ## mu        <- mu
              start.notpen[j] <- 1
            }else{ ## Update the information
              coef <- coef.test
              eta  <- eta.test
              mu   <- invlink(eta)
              start.notpen[j] <- scale
            }

            ## Save the scaling factor for the next iteration (in order that
            ## we only have to do very few line searches)
            ## start.notpen[j] <- scale

            ## Update the remaining information
            ## coef <- coef.test
            ## eta  <- eta.test
            ## mu   <- invlink(eta)
          } ## end if(abs(d) > sqrt(.Machine$double.eps))
        } ## end for(j in 1:nrnotpen)
      } ## if(any.notpen)

      ## Optimize the *penalized* parameter groups
      for(j in guessed.active){
        ind  <- ipen.which[[j]]
        npar <- ipen.tab[j]

        coef.ind       <- coef[ind]
        cross.coef.ind <- crossprod(coef.ind)

        ## Design matrix of the current group
        Xj <- x.pen[[j]]

        ## Negative gradient of the current group
        ngrad <- c(ngradient(Xj, y, mu, weights, ...))

        ## Update the Hessian if necessary
        if(update.hess == "always"){
          diagH <- numeric(length(ind))
          for(i in 1:length(ind)){ ## for loop seems to be faster than sapply
            diagH[i] <- nhessian(Xj[,i,drop = FALSE], mu, weights, ...)
          }
          nH <- min(max(diagH, lower), upper)
        }else{
          nH <- nH.pen[j]
        }

        cond       <- -ngrad + nH * coef.ind
        cond.norm2 <- c(crossprod(cond)) ## was: crossprod(cond)
        ## Because of warning message:
        ## Recycling array of length 1 in vector-array arithmetic is deprecated.

        ## Check the condition whether the minimum is at the non-differentiable
        ## position (-coef.ind) via the condition on the subgradient.
        lval <- ifelse(j %in% c(M_groups, D_groups), (j %in% M_groups)*mscale*l + (j %in% D_groups)*dscale*l, l)
        border <- penscale(npar) * lval
        if(cond.norm2 > border^2){
          d <- (1 / nH) *
            (-ngrad - lval * penscale(npar) * (cond / sqrt(cond.norm2)))
          ##d <- zapsmall(c(coef.ind, d))[-(1:npar)]
        }else{
          d <- -coef.ind
        }

        ## If !all(d == 0), we have to do a line search
        if(!all(d == 0)){
          scale <- min(start.pen[j] / beta, 1)

          coef.test      <- coef
          coef.test[ind] <- coef.ind + scale * d
          Xjd            <- c(Xj %*% d)
          eta.test       <- eta + Xjd * scale

          if(line.search){
            qh <- sum(ngrad * d) +
              lval * penscale(npar) * sqrt(crossprod(coef.ind + d)) -
                lval * penscale(npar)* sqrt(cross.coef.ind)

            fn.val.test    <- nloglik(y, eta.test, weights, ...)
            fn.val0        <- nloglik(y, eta, weights, ...)

            left <- fn.val.test +
              lval  * penscale(npar) * sqrt(crossprod(coef.test[ind]))

            ## lval, not l: the group's own multiplier (test plan Z12, G10)
            right <- fn.val0 + lval * penscale(npar) * sqrt(cross.coef.ind) +
              sigma * scale * qh

            while(left > right & scale > 10^-30){
              ##cat("Doing line search (pen)\n")
              scale          <- scale * beta
              coef.test[ind] <- coef.ind + scale * d
              eta.test       <- eta + Xjd * scale
              fn.val.test    <- nloglik(y, eta.test, weights, ...)

              left <- fn.val.test +
                lval  * penscale(npar) * sqrt(crossprod(coef.test[ind]))

              right <- fn.val0 + lval * penscale(npar) * sqrt(cross.coef.ind) +
                sigma * scale * qh
            } ## end while(left > right & qh != 0)
          } ## end if(line.search)
          ## If we escaped the while loop because 'scale' is too small
          ## (= we add nothing), we just stay at the current solution to
          ## prevent tiny values
          if(scale <= 10^-30){ ## Do *nothing* in that case
            ##coef.test <- coef
            ##eta.test  <- eta
            ##mu        <- mu
            start.pen[j] <- 1
          }else{
            coef <- coef.test
            eta  <- eta.test
            mu   <- invlink(eta)
            start.pen[j] <- scale
          }
        } ## end if(!all(d == 0))
        norms.pen[j] <- sqrt(crossprod(coef[ind]))
      } ## end for(j in guessed.active)

      fn.val <- nloglik(y, eta, weights, ...) +
        l * sum(penscale(ipen.tab[-c(M_groups, D_groups)]) * norms.pen[-c(M_groups, D_groups)]) +
        mscale * l * sum(penscale(ipen.tab[M_groups]) * norms.pen[M_groups]) +
        dscale * l * sum(penscale(ipen.tab[D_groups]) * norms.pen[D_groups])

      ## Relative difference with respect to parameter vector
      ##d.par <- sqrt(crossprod(coef - coef.old)) / (1 + sqrt(crossprod(coef)))

      d.par <-  max(abs(coef - coef.old) / (1 + abs(coef)))
      ##d.par <-  max(abs(coef - coef.old) / (ifelse(abs(coef), abs(coef), 1)))

      ## Relative difference with respect to function value (penalized
      ## likelihood)
      d.fn <- abs(fn.val.old - fn.val) / (1 + abs(fn.val))

      ## Print out improvement if desired (trace >= 2)
      if(trace >= 2){
        cat("d.fn:", d.fn, " d.par:", d.par,
            " nr.var:", sum(coef != 0), "\n")
      }

      ## If we are working on a sub-set of predictors and have converged
      ## we stop the optimization and will do a loop through all
      ## predictors in the next run. Therefore we set counter = 0.

      ##if(d.fn <= tol & d.par <= sqrt(tol)){
      if(d.fn <= tol & d.par <= sqrt(tol)){
        counter <- 0 ## will force a run through all groups
        if(trace >= 2 & !do.all)
          cat("...Subproblem (active set) solved\n")
      }
    } ## end of while(d.fn > tol | d.par > sqrt(tol) | !do.all)

    if(trace == 1) {
      cat("Lambda:", l, " nr.var:", sum(coef != 0), "\n")
    }


    coef.m[,pos]            <- coef
    fn.val.v[pos]           <- fn.val
    norms.pen.m[,pos]       <- norms.pen
    nloglik.v[pos]          <- nloglik(y, eta, weights, ...)
    grad.m[,pos]            <- ngradient(x, y, mu, weights, ...)
    linear.predictors[,pos] <- eta
    fitted[,pos]            <- invlink(eta)
  } ## end for(pos in 1:nrlambda){

  ## Transform the coefficients back to the original scale if the design
  ## matrix was standardized
  if(standardize){
    if(any.notpen)
      coef.m[inotpen.which,] <- (1 / scale.notpen) * coef.m[inotpen.which,]
    ## For df > 1 we have to use a matrix inversion to go back to the
    ## original scale
    for(j in 1:length(ipen.which)){
      ind <- ipen.which[[j]]
      coef.m[ind,] <- solve(scale.pen[[j]], coef.m[ind,,drop = FALSE])
    }
  }

  ## Undo the residualization: move the components absorbed by the
  ## non-penalized columns back onto their coefficients
  if(residualize)
    coef.m[inotpen.which, ] <- coef.m[inotpen.which, , drop = FALSE] -
      resid.coef %*% coef.m[ipen.cols, , drop = FALSE]

  ## Need to adjust intercept if we have performed centering
  if(center){
    coef.m[intercept.which,] <- coef.m[intercept.which,] -
      apply(coef.m[-intercept.which,,drop = FALSE] * mu.x, 2, sum)
  }

  ## Overwrite values of x.old if we don't want to save it
  if(!save.x)
    x.old <- NULL
  if(!save.y)
    y <- NULL

  out <- list(x = x.old, ## use untransformed values
              y = y,
              coefficients = coef.m,
              norms.pen    = norms.pen.m,
              lambda       = lambda,
              mscale       = mscale,
              dscale       = dscale,
              residualize  = residualize,
              resid.coef   = if(residualize) resid.coef else NULL,
              index        = index,
              penscale     = penscale,
              model        = model,
              ngradient    = grad.m,
              nloglik      = nloglik.v,
              fitted       = fitted,
              linear.predictors = linear.predictors,
              fn.val       = fn.val.v,
              converged    = converged,
              weights      = weights,
              offset       = offset,
              control      = control,
              call         = match.call())
  structure(out, class = "grplasso")
}


