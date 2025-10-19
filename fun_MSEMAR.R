##############################################################################-#
# The following codes include functions implementing MSEMAR using two parameter#
# estimation methods: ml and reml                                              #
##############################################################################-#

#'@description Emvmeta implement the parameter estimation of MSEMAR
#'@author Luo
#'@param formula An object of class "formula", see mvmeta::mvmeta.
#'@param S The covariance for estimated outcomes in a stratified analysis. See \r
#' in mvmeta::mvmeta.
#' @param data,subset,offset,na.action,model,offset,na.action seen mvmeta::mvmeta.
#' @param method Currently, only support "ml" and "reml". 
#' @param Wmatrix A matrix denoting spatial adjacence relationship
#' @param bscov Define the stucture of Psi in MSEMAR. Currently,only support "unstr".
#' @param control seen the following function Emvmeta.control

Emvmeta <- function (formula, S, data, subset,Wmatrix = NULL, method = "reml", bscov = "unstr", 
                     model = TRUE, contrasts = NULL, offset, na.action, control = list()){
  call <- match.call()
  mcall <- match.call(expand.dots = FALSE)
  mn <- match(c("formula", "data", "subset", "weights", "na.action", 
                "offset"), names(mcall), 0L)
  mcall <- mcall[c(1L, mn)]
  mcall$drop.unused.levels <- TRUE
  mcall[[1L]] <- as.name("model.frame")
  
  if (missing(data)) 
    data <- parent.frame()
  if (!inherits(eval(substitute(formula), data), "formula")) {
    formula <- as.formula(paste(deparse(substitute(formula), 
                                        width.cutoff = 499L), "~ 1"), env = parent.frame())
    environment(formula) <- parent.frame()
    call[[mn[1]]] <- mcall[[mn[1]]] <- formula
  }
  if (missing(data)) 
    data <- environment(formula)
  
  mcall$na.action <- "na.pass"
  mf <- eval(mcall, parent.frame())
  class(mf) <- c("data.frame.mvmeta", class(mf))
  if (missing(na.action)) 
    na.action <- getOption("na.action")
  if (length(na.action)) 
    mf <- do.call(na.action, list(mf))
  if (method == "model.frame") 
    return(mf)
  if (is.empty.model(mf)) 
    stop("empty model not allowed")
  
  method <- match.arg(method, c("fixed", "ml", "reml", "mm", 
                                "vc"))
  bscov <- match.arg(bscov, c("unstr", "diag", "id", "cs", 
                              "hcs", "ar1", "prop", "cor", "fixed"))
  if (bscov != "unstr" && !method %in% c("ml", "reml")) 
    stop("structured Psi only available for methods 'ml' or 'reml'")
  
  terms <- attr(mf, "terms")
  y <- as.matrix(model.response(mf, "numeric"))
  X <- model.matrix(terms, mf, contrasts)
  offset <- as.vector(model.offset(mf))
  if (!is.null(offset)) {
    if (length(offset) != NROW(y)) 
      stop("number of offsets should equal number of observations")
  }
  S <- eval(call$S, data, parent.frame())
  S <- mvmeta:::mkS(S, y, attr(mf, "na.action"), if (missing(subset))
    NULL
    else eval(call$subset, data, parent.frame()))
  if (nrow(y) < 2L) 
    stop("less than 2 valid studies after exclusion of missing")
  
  if(is.null(control$initial_V)){
    mvcontrol <- mvmeta:::mvmeta.control()
    mvcontrol <- modifyList(mvcontrol,control)
    fit0 <- mvmeta::mvmeta.fit(X, y, S, offset, method, bscov, list())
    control$initial_V <- fit0$Psi
  }
  
  fit <- Emvmeta.fit(X, y, S, Wmatrix,offset, method, bscov, control)
  fit$model <- if (model) mf  else NULL
  fit$S <- S
  fit$na.action <- attr(mf, "na.action")
  fit$call <- call
  fit$formula <- formula
  fit$terms <- terms
  fit$contrasts <- attr(X, "contrasts")
  fit$xlevels <- .getXlevels(terms, mf)
  class(fit) <- "Emvmeta"
  fit
}


Emvmeta.fit <- function (X, y, S, Wmatrix, offset = NULL, method = "reml", bscov = "unstr", 
                         control = list()){
  control <- do.call("Emvmeta.control", control)
  
  y <- as.matrix(y)
  nay <- is.na(y)
  k <- ncol(y) 
  m <- nrow(y)
  p <- ncol(X)
  nall <- length(y) 
  nk <- colnames(y)
  if (k > 1L && is.null(nk)) 
    nk <- paste("y", seq(k), sep = "")
  nm <- rownames(y)
  np <- colnames(X)
  
  if (control$inputna) {
    augdata <- inputna(y, S, inputvar = control$inputvar)
    y <- augdata[, seq(k)]
    S <- augdata[, -seq(k)]
    nay[nay] <- FALSE
  }
  
  if (!is.null(offset)) 
    y <- y - offset
  
  Xlist <- lapply(seq(m), function(i) diag(1, k) %x% matrix(X[i,],nrow = 1))
  Xstar <- do.call("rbind",Xlist)
  ylist <- lapply(seq(m),function(i) y[i,][!nay[i,]])
  Ystar <- as.numeric((unlist(ylist)))
  Wstar <- Wmatrix %x% diag(1, k)
  
  if (dim(S)[2] == k) 
    S <- inputcov(sqrt(S), control$Scor)
  Slist <- lapply(seq(m),function(i) {
    Si <- xpndMat(S[i,])[!nay[i,],!nay[i,],drop=FALSE]
    if(any(is.na(Si))) stop("missing pattern in 'y' and S' is not consistent")
    return(Si)
  })
  D <- lapply(seq(m), function(i) t(diag(1, m)[i,] %x% Slist[[i]]))
  D <- do.call("rbind",D)
  
  fun <- paste("Emvmeta", method, sep = ".")
  fit <- do.call(fun, list(Xstar = Xstar, Ystar = Ystar, D = D, Wstar = Wstar,
                           k = k, m = m, p = p, bscov = bscov,nall = nall, control = control))
  
  if (!fit$converged) {
    warning("Not convergency")
  }
  
   iter_llvalues <- fit$iter_llvalues
   names(iter_llvalues) <- c("Lambda",rep(c("Psi","Lambda"),length(iter_llvalues)/2-1),"PsiLambda")
   fit$iter_llvalues <- iter_llvalues
  #
  fit$Xstar <- Xstar
  fit$Ystar <- Ystar
  fit$Wstar <- Wstar
  U <- solve(diag(k*m)-fit$Lambda*Wstar)
  Psi_E <- fit$Psi_E
  V <- diag(m) %x% Psi_E
  G <- U %*% V %*% t(U)
  invSIGMA <- chol2inv(chol(G+D))
  ytilde <- Ystar - fit$fitted.values
  fit$xi <- G %*% invSIGMA %*% ytilde
  fit$epsilon <- ytilde - fit$xi
  # 
  fit$fitted.values.spatial <- fit$fitted.values + fit$xi
  fit$method <- method
  fit$bscov <- bscov
  fit$offset <- offset
  fit$dim <- list(k = k, m = m, p = p)
  fit$df <- list(nobs = nall - (method == "reml") * 
                   fit$rank, df = nall - fit$df.residual, fixed = fit$rank, 
                 random = ifelse(method == "fixed", 0, nall - fit$rank - 
                                   fit$df.residual))
  fit$D <- D
  fit$lab <- list(k = nk, p = np)
  temp <- as.numeric(fit$fitted.values)
  fit$fitted.values <- matrix(temp, m, k, byrow = TRUE)
  temp <- as.numeric(fit$fitted.values.spatial)
  fit$fitted.values.spatial <- matrix(temp, m, k, byrow = TRUE)
  temp <- as.numeric(fit$epsilon)
  fit$epsilon <- matrix(temp, m, k, byrow = TRUE)
  temp <- as.numeric(fit$xi)
  fit$xi <- matrix(temp, m, k, byrow = TRUE)
  if (!is.null(offset)) {
    y <- y + offset
    fit$fitted.values <- fit$fitted.values + offset
  }
  if (method != "fixed") 
    dimnames(fit$Psi_E) <- list(nk, nk)
  if (k == 1L) {
    names(fit$coefficients) <- np
    dimnames(fit$vcov) <- list(np, np)
    fit$fitted.values <- drop(fit$fitted.values)
    fit$residuals <- drop(y - fit$fitted.values)
    names(fit$residuals) <- names(fit$fitted.values) <- nm
    names(fit$fitted.values.spatial) <- names(fit$epsilon) <- names(fit$xi) <- nm
  }
  else {
    fit$coefficients <- matrix(fit$coefficients, p, k, dimnames = list(np,nk))
    rownames(fit$vcov) <- colnames(fit$vcov) <- paste(rep(nk, 
                                                          each = p), rep(np, k), sep = ".")
    fit$residuals <- y - fit$fitted.values
    dimnames(fit$residuals) <- dimnames(fit$fitted.values) <- list(nm,nk)
    dimnames(fit$fitted.values.spatial) <- dimnames(fit$epsilon) <- dimnames(fit$xi) <- list(nm,nk)
  }
  fit
}

#'@description Define the optimal parameters in Emvmeta
#'@param opt.iter.method The optimal iterative method for estimation Psi
#'@param factr Controls the convergence of the "L-BFGS-B" method.Only for the \r
#'final optimal function in MSEMAR, i.e., the function with respect to beta, Lambda, \r
#'and Psi at the same time.
#'@param hessian Logical. Should a numerically differentiated Hessian matrix be returned?
#'@param opt.iter A integer larger 0. The iterative time for Lambda-Psi.
#'@param opt.iter.show Logical. show whether the iterative values of AIC are printed.
#'@param lltol judge if the log-likelihood is convergency in the iterative process of Lambda-Psi.
#'@param others seen mvmeta::mvmeta.control 
Emvmeta.control <- function (optim = list(), initial_V = NULL, showiter = FALSE, maxiter = 100, initPsi = NULL, 
                             Psifix = NULL, Psicor = 0, Scor = 0, inputna = FALSE, inputvar = 10^4, 
                             hessian = FALSE, vc.adj = TRUE, factr = 1e7, lltol = 1e-5,
                             set.negeigen = sqrt(.Machine$double.eps),opt.iter.method = "BFGS",
                             opt.iter = 1,opt.iter.show = F) {
  optim <- modifyList(list(fnscale = -1, maxit = maxiter, factr = factr), 
                      optim)
  if (showiter) {
    optim$trace <- 6
    optim$REPORT <- 1
  }
  list(optim = optim, initial_V = initial_V, showiter = showiter, maxiter = maxiter, 
       hessian = hessian, initPsi = initPsi, Psifix = Psifix, 
       Psicor = Psicor, Scor = Scor, inputna = inputna, inputvar = inputvar, 
       vc.adj = vc.adj, factr = factr, set.negeigen = set.negeigen,lltol= lltol,
       opt.iter.method = opt.iter.method,opt.iter = opt.iter,opt.iter.show = opt.iter.show)
}


library(mvmeta)
mlprof.fn_E <- function (par,Xstar, Ystar, Wstar, D, k, m){
  Lambda <- par[1]
  Psi_E <- par2V(par[-1],k) 
  gls <- glsfit_E(Xstar, Ystar, Wstar, D, Psi_E, Lambda, k, m, onlycoef = FALSE)
  pconst <- -0.5 * m * k * log(2 * pi) 
  pres <- -0.5 * (crossprod(gls$invtMY - gls$invtMX %*% gls$coef))
  pdet <- -sum(log(diag(gls$M)))
  as.numeric(pconst + pdet + pres)
}
mlprof.fnLambda_E <- function (par,Xstar, Ystar, Wstar, Psi_E, D, k, m){
  Lambda <- par
  gls <- glsfit_E(Xstar, Ystar, Wstar, D, Psi_E, Lambda, k, m, onlycoef = FALSE) 
  pconst <- -0.5 * m * k * log(2 * pi)
  pres <- -0.5 * (crossprod(gls$invtMY - gls$invtMX %*% gls$coef))
  pdet <- -sum(log(diag(gls$M)))
  as.numeric(pconst + pdet + pres)
}
mlprof.fnPsi_E <- function (par,Xstar, Ystar, Wstar, Lambda, D, k, m){
  Psi_E <- par2V(par,k)
  gls <- glsfit_E(Xstar, Ystar, Wstar, D, Psi_E, Lambda, k, m, onlycoef = FALSE)
  pconst <- -0.5 * m * k * log(2 * pi)
  pres <- -0.5 * (crossprod(gls$invtMY - gls$invtMX %*% gls$coef))
  pdet <- -sum(log(diag(gls$M)))
  as.numeric(pconst + pdet + pres)
}
glsfit_E <- function (Xstar, Ystar, Wstar, D, Psi_E, Lambda, k, m, onlycoef = TRUE){
  
  ycoef <- solve(diag(m*k)-Lambda*Wstar)
  Psi_Estar <- diag(m) %x% Psi_E
  SIGMA <- ycoef %*% Psi_Estar %*% t(ycoef) + D  
  M <- chol(SIGMA)
  Xnew <- Xstar
  invM <- backsolve(M,diag(ncol(M)))
  invtMX <- crossprod(invM,Xnew) 
  invtMY <- crossprod(invM,Ystar)
  coef <- as.numeric(qr.solve(invtMX, invtMY))
  
  if (onlycoef) 
    return(coef)
  list(coef = coef, SIGMA = SIGMA, M = M, invM = invM,invtMX = invtMX, invtMY = invtMY)
}

opt.iter.UV_E <- function(fnLambda,fnPsi_E,Psi_E,opt.iter,opt.iter.method,Xstar,Ystar,
                        D,k,m,Wstar,lower,upper,opt.iter.show){
  optLambda <- optimize(f = fnLambda,lower = lower,upper = upper,maximum = T, Xstar = Xstar, 
                     Ystar = Ystar, D = D, k = k, m = m, Psi_E = Psi_E, Wstar = Wstar)
  if(opt.iter.show) cat(optLambda$objective," ")
  llvalues <- optLambda$objective
  Lambda <- optLambda$maximum
  Psi_E <-as.matrix(Matrix::nearPD(Psi_E)[[1]])  
  parPsi_E <- vechMat(t(chol(Psi_E))) 
  if(opt.iter < 2) par <- c(Lambda,parPsi_E) else{
    for (i in 1:(opt.iter-1)) {
      optPsi_E <- optim(par = parPsi_E, fn = fnPsi_E, gr = NULL, Xstar = Xstar, 
                        Ystar = Ystar, D = D, k = k, m = m, Lambda = Lambda, Wstar = Wstar,
                        method = opt.iter.method, control = list(fnscale = -1), hessian = F)
      Psi_E <- par2V(optPsi_E$par,k)
      optLambda <- optimize(f = fnLambda,lower = lower,upper = upper,maximum = T, Xstar = Xstar, 
                         Ystar = Ystar, D = D, k = k, m = m, Psi_E = Psi_E, Wstar = Wstar)
      Lambda <- optLambda$maximum
      parPsi_E <- optPsi_E$par
      if(opt.iter.show) cat(optPsi_E$value,optLambda$objective," ")
      llvalues <- c(llvalues,optPsi_E$value,optLambda$objective)
    }
    par <- c(Lambda,parPsi_E)
  }
  list(par=par,llvalues=llvalues)
}


Emvmeta.ml <- function (Xstar, Ystar, D = D, Wstar = Wstar, k, m, p, bscov, 
                        nall,control, ...){
  fn <- mlprof.fn_E
  fnLambda <- mlprof.fnLambda_E
  fnPsi_E <- mlprof.fnPsi_E
  Psi_E <- control$initial_V
  yita <- 1
  lower <- -yita + 1e-6  
  upper <- yita - 1e-6  
  par <- opt.iter.UV_E(fnLambda = fnLambda,fnPsi_E = fnPsi_E,Psi_E = Psi_E,opt.iter = control$opt.iter,
                     opt.iter.method = control$opt.iter.method,Xstar = Xstar,
                     Ystar = Ystar,D = D,k = k,m = m,Wstar = Wstar,
                     lower = lower,upper = upper,opt.iter.show = control$opt.iter.show)
  llvalues <- par[[2]]
  par <- par[[1]]
  nparPsi_E <- length(par) - 1
  opt <- optim(par = par, fn = fn, gr = NULL, Xstar = Xstar,
               Ystar = Ystar, D = D, k = k, m = m, Wstar = Wstar,
               lower = c(lower,rep(-1 + 1e-6,nparPsi_E)), upper = c(upper,rep(1 - 1e-6,nparPsi_E)),
               method = "L-BFGS-B", control = control$optim, hessian = control$hessian)  
  if(control$opt.iter.show) cat(opt$value," ")
  llvalues <- c(llvalues,opt$value)
  Psi_E <- par2V(opt$par[-1],k)
  Lambda <- opt$par[1]
  
  gls <- glsfit_E(Xstar, Ystar, Wstar, D, Psi_E, Lambda, k, m, onlycoef = FALSE)
  
  qrinvtMX <- qr(gls$invtMX)
  R <- qr.R(qrinvtMX)
  Qty <- qr.qty(qrinvtMX, gls$invtMY)
  vcov <- tcrossprod(backsolve(R, diag(1, ncol(gls$invtMX))))
  
  res <- NULL
  fitted <- Xstar %*% gls$coef
  rank <- qrinvtMX$rank
  converged <- opt$convergence == 0 | abs(llvalues[length(llvalues)]-llvalues[length(llvalues)-1]) < control$lltol
  
  c(list(coefficients = gls$coef, vcov = vcov, Psi_E = Psi_E, Lambda = Lambda,residuals = res, 
         fitted.values = fitted, df.residual = nall - rank - length(par), iter_llvalues = llvalues,
         rank = rank, logLik = opt$value, converged = converged, par = opt$par), if (!is.null(opt$hessian)) 
           list(hessian = opt$hessian, se.Lambda = sqrt(-solve(opt$hessian)[1,1])), 
    list(niter = opt$counts[[2]], control = control))
}



Emvmeta.reml <- function (Xstar, Ystar, D = D, Wstar = Wstar, k, m, p, bscov, 
                          nall,control, ...){
  fn <- remlprof.fn_E
  fnLambda <- remlprof.fnLambda_E
  fnPsi_E <- remlprof.fnPsi_E
  Psi_E <- control$initial_V
  yita <- 1
  lower <- -yita + 1e-6  
  upper <- yita - 1e-6  
  par <- opt.iter.UV_E(fnLambda = fnLambda,fnPsi_E = fnPsi_E,Psi_E = Psi_E,opt.iter = control$opt.iter,
                     opt.iter.method = control$opt.iter.method,Xstar = Xstar,
                     Ystar = Ystar,D = D,k = k,m = m,Wstar = Wstar,
                     lower = lower,upper = upper,opt.iter.show = control$opt.iter.show)
  llvalues <- par[[2]]
  par <- par[[1]]
  nparPsi_E <- length(par) - 1
  opt <- optim(par = par, fn = fn, gr = NULL, Xstar = Xstar,
               Ystar = Ystar, D = D, k = k, m = m, Wstar = Wstar,
               lower = c(lower,rep(-1 + 1e-6,nparPsi_E)), upper = c(upper,rep(1 - 1e-6,nparPsi_E)),
               method = "L-BFGS-B", control = control$optim, hessian = control$hessian)  
  if(control$opt.iter.show) cat(opt$value," ")
  llvalues <- c(llvalues,opt$value)
  Psi_E <- par2V(opt$par[-1],k)
  Lambda <- opt$par[1]
  
  gls <- glsfit_E(Xstar, Ystar, Wstar, D, Psi_E, Lambda, k, m, onlycoef = FALSE)
  
  qrinvtMX <- qr(gls$invtMX)
  R <- qr.R(qrinvtMX)
  Qty <- qr.qty(qrinvtMX, gls$invtMY)
  vcov <- tcrossprod(backsolve(R, diag(1, ncol(gls$invtMX))))
  
  res <- NULL
  fitted <- Xstar %*% gls$coef
  rank <- qrinvtMX$rank
  converged <- opt$convergence == 0 | abs(llvalues[length(llvalues)]-llvalues[length(llvalues)-1]) < control$lltol
  
  c(list(coefficients = gls$coef, vcov = vcov, Psi_E = Psi_E, Lambda = Lambda,residuals = res, 
         fitted.values = fitted, df.residual = nall - rank - length(par), iter_llvalues = llvalues,
         rank = rank, logLik = opt$value, converged = converged, par = opt$par), if (!is.null(opt$hessian)) 
           list(hessian = opt$hessian, se.Lambda = sqrt(-solve(opt$hessian)[1,1])), 
    list(niter = opt$counts[[2]], control = control))
}
remlprof.fn_E <- function (par,Xstar, Ystar, Wstar, D, k, m){
  Lambda <- par[1]
  Psi_E <- par2V(par[-1],k) 
  gls <- glsfit_E(Xstar, Ystar, Wstar, D, Psi_E, Lambda, k, m, onlycoef = FALSE)
  pconst <- -0.5 * (m * k - length(gls$coef)) * log(2 * pi) 
  pres <- -0.5 * (crossprod(gls$invtMY - gls$invtMX %*% gls$coef))
  pdet1 <- -sum(log(diag(gls$M)))
  tXWXtot <- crossprod(gls$invtMX)
  pdet2 <- -sum(log(diag(chol(tXWXtot))))
  
  as.numeric(pconst + pdet1 + pdet2 + pres)
}
remlprof.fnLambda_E <- function (par,Xstar, Ystar, Wstar, Psi_E, D, k, m){
  Lambda <- par
  gls <- glsfit_E(Xstar, Ystar, Wstar, D, Psi_E, Lambda, k, m, onlycoef = FALSE)
  pconst <- -0.5 * (m * k - length(gls$coef)) * log(2 * pi)
  pres <- -0.5 * (crossprod(gls$invtMY - gls$invtMX %*% gls$coef))
  pdet1 <- -sum(log(diag(gls$M)))
  tXWXtot <- crossprod(gls$invtMX)
  pdet2 <- -sum(log(diag(chol(tXWXtot))))
  
  as.numeric(pconst + pdet1 + pdet2 + pres)
}
remlprof.fnPsi_E <- function (par,Xstar, Ystar, Wstar, Lambda, D, k, m){
  Psi_E <- par2V(par,k)
  gls <- glsfit_E(Xstar, Ystar, Wstar, D, Psi_E, Lambda, k, m, onlycoef = FALSE)
  pconst <- -0.5 * (m * k - length(gls$coef)) * log(2 * pi)
  pres <- -0.5 * (crossprod(gls$invtMY - gls$invtMX %*% gls$coef))
  pdet1 <- -sum(log(diag(gls$M)))
  tXWXtot <- crossprod(gls$invtMX)
  pdet2 <- -sum(log(diag(chol(tXWXtot))))
  
  as.numeric(pconst + pdet1 + pdet2 + pres)
}


#Output overall ERR coeffients and their var/cov matrix 
Outoerr_E <- function(object){
  
  coef <- object$coefficients ; vcov <- object$vcov
  pred <- NULL
  pred$coefsim <- coef
  pred$vcovsim <- vcov
    
  pred
  }


blup.Emvmeta <- function (object,nsim=5000) 
{
  mf <- model.frame(object)
  X <- model.matrix(object)
  y <- as.matrix(model.response(mf, "numeric"))
  S <- object$S
  offset <- object$offset
  if (!is.null(object$na.action)) {
    X <- napredict(object$na.action, X)
    y <- napredict(object$na.action, y)
    S <- napredict(object$na.action, S)
    offset <- napredict(object$na.action, offset)
  }
  if (!is.null(offset)) 
    y <- y - offset
  Lambda <- object$Lambda
  se.Lambda <- object$se.Lambda
  Wstar <- object$Wstar
  m <- nrow(y)
  k <- ncol(y)
  Xstar <- object$Xstar
  Ystar <- object$Ystar
  fittedstar <- as.numeric(t(object$fitted.values))
  
  U <- solve(diag(k*m)-Lambda*Wstar)
  Psi_E <- object$Psi_E
  V <- diag(m) %x% Psi_E
  G <- U %*% V %*% t(U)
  invSIGMA <- chol2inv(chol(G+object$D))
  
  ytilde <- Ystar - fittedstar
  xi <- G %*% invSIGMA %*% ytilde
  epsilon <- ytilde - xi
  fitted.values.spatial <- fittedstar + xi
 
  covfit <- Xstar %*% tcrossprod(object$vcov,Xstar) + G - G %*% invSIGMA %*% G
  temp <- as.numeric(fitted.values.spatial)
  fitted.values.spatial <- matrix(temp, m, k, byrow = TRUE)
  temp <- lapply(1:m, function(x) {
    index <- (x - 1) * k + 1:k
    a <- covfit[index, index]
    colnames(a) <- rownames(a) <- colnames(y)
    a
  })
  names(temp) <- rownames(y)
  blup <- vector(mode = "list", length = m)
  names(blup) <- rownames(y)
  for (i in 1:m) {
    a <- fitted.values.spatial[i, ]
    names(a) <- colnames(y)
    blup[[i]] <- list(blup = a, vcov = temp[[i]])
  }
  blup
}

par2V <- function(par,k){
  L <- diag(0, k)
  L[lower.tri(L, diag = TRUE)] <- par
  R_L <- t(L)
  t(R_L) %*% R_L
}
