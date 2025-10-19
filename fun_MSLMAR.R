##############################################################################-#
# The following codes include functions implementing MSLMAR using two parameter#
# estimation methods: ml and reml                                              #
##############################################################################-#

#'@description Lmvmeta implement the parameter estimation of MSLMAR
#'@author Luo
#'@param formula An object of class "formula", see mvmeta::mvmeta.
#'@param S The covariance for estimated outcomes in a stratified analysis. See \r
#' in mvmeta::mvmeta.
#' @param data,subset,offset,na.action,model,offset,na.action seen mvmeta::mvmeta.
#' @param method Currently, only support "ml" and "reml". 
#' @param Wmatrix A matrix denoting spatial adjacence relationship
#' @param bscov Define the stucture of Psi in MSLMAR. Currently,only support "unstr".
#' @param control seen the following function Lmvmeta.control

Lmvmeta <- function (formula, S, data, subset,Wmatrix = NULL, method = "reml", bscov = "unstr", 
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
  
  fit <- Lmvmeta.fit(X, y, S, Wmatrix,offset, method, bscov, control)
  fit$model <- if (model) mf  else NULL
  fit$S <- S
  fit$na.action <- attr(mf, "na.action")
  fit$call <- call
  fit$formula <- formula
  fit$terms <- terms
  fit$contrasts <- attr(X, "contrasts")
  fit$xlevels <- .getXlevels(terms, mf)
  class(fit) <- "Lmvmeta"
  fit
}



Lmvmeta.fit <- function (X, y, S, Wmatrix, offset = NULL, method = "reml", bscov = "unstr", 
                         control = list()){
  control <- do.call("Lmvmeta.control", control)
  
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
  
  fun <- paste("Lmvmeta", method, sep = ".")
  fit <- do.call(fun, list(Xstar = Xstar, Ystar = Ystar, D = D, Wstar = Wstar,
                           k = k, m = m, p = p, bscov = bscov,nall = nall, control = control))
  
  if (!fit$converged) {
    warning("Not convergency")
  }
  
  iter_llvalues <- fit$iter_llvalues
  names(iter_llvalues) <- c("Rho",rep(c("Psi_L","Rho"),length(iter_llvalues)/2-1),"PsiRho")
  fit$iter_llvalues <- iter_llvalues
  #
  fit$Xstar <- Xstar
  fit$Ystar <- Ystar
  fit$Wstar <- Wstar
  U <- solve(diag(k*m)-fit$Rho*Wstar)
  Psi_L <- fit$Psi_L
  V <- diag(m) %x% Psi_L
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
    dimnames(fit$Psi_L) <- list(nk, nk)
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

#'@description Define the optimal parameters in Lmvmeta
#'@param opt.iter.method The optimal iterative method for estimation Psi
#'@param factr Controls the convergence of the "L-BFGS-B" method.Only for the \r
#'final optimal function in MSLMAR, i.e., the function with respect to beta, Rho, \r
#'and Psi at the same time.
#'@param hessian Logical. Should a numerically differentiated Hessian matrix be returned?
#'@param opt.iter A integer larger 0. The iterative time for Rho-Psi.
#'@param opt.iter.show Logical. show whether the iterative values of AIC are printed.
#'@param lltol judge if the log-likelihood is convergency in the iterative process of Rho-Psi
#'@param others seen mvmeta::mvmeta.control 
Lmvmeta.control <- function (optim = list(), initial_V = NULL, showiter = FALSE, maxiter = 100, initPsi = NULL, 
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
mlprof.fn_L <- function (par,Xstar, Ystar, Wstar, D, k, m){
  Rho <- par[1]
  Psi_L <- par2V(par[-1],k) 
  gls <- glsfit_L(Xstar, Ystar, Wstar, D, Psi_L, Rho, k, m, onlycoef = FALSE)
  pconst <- -0.5 * m * k * log(2 * pi) 
  pres <- -0.5 * (crossprod(gls$invtMY - gls$invtMX %*% gls$coef))
  pdet <- -sum(log(diag(gls$M)))
  as.numeric(pconst + pdet + pres)
}
mlprof.fnRho_L <- function (par,Xstar, Ystar, Wstar, Psi_L, D, k, m){
  Rho <- par
  gls <- glsfit_L(Xstar, Ystar, Wstar, D, Psi_L, Rho, k, m, onlycoef = FALSE)
  pconst <- -0.5 * m * k * log(2 * pi)
  pres <- -0.5 * (crossprod(gls$invtMY - gls$invtMX %*% gls$coef))
  pdet <- -sum(log(diag(gls$M)))
  as.numeric(pconst + pdet + pres)
}
mlprof.fnPsi_L <- function (par,Xstar, Ystar, Wstar, Rho, D, k, m){
  Psi_L <- par2V(par,k)
  gls <- glsfit_L(Xstar, Ystar, Wstar, D, Psi_L, Rho, k, m, onlycoef = FALSE)
  pconst <- -0.5 * m * k * log(2 * pi)
  pres <- -0.5 * (crossprod(gls$invtMY - gls$invtMX %*% gls$coef))
  pdet <- -sum(log(diag(gls$M)))
  as.numeric(pconst + pdet + pres)
}
glsfit_L <- function (Xstar, Ystar, Wstar, D, Psi_L, Rho, k, m, onlycoef = TRUE){
  
  ycoef <- solve(diag(m*k)-Rho*Wstar)
  Psi_Lstar <- diag(m) %x% Psi_L
  SIGMA <- ycoef %*% Psi_Lstar %*% t(ycoef) + D
  
  M <- chol(SIGMA)
  Xnew <- solve(diag(m*k)-Rho*Wstar) %*% Xstar
  invM <- backsolve(M,diag(ncol(M)))
  invtMX <- crossprod(invM,Xnew) 
  invtMY <- crossprod(invM,Ystar)
  coef <- as.numeric(qr.solve(invtMX, invtMY))
  
  if (onlycoef) 
    return(coef)
  list(coef = coef, SIGMA = SIGMA, M = M, invM = invM,invtMX = invtMX, invtMY = invtMY)
}

opt.iter.UV_L <- function(fnRho,fnPsi_L,Psi_L,opt.iter,opt.iter.method,Xstar,Ystar,
                        D,k,m,Wstar,lower,upper,opt.iter.show){
  optRho <- optimize(f = fnRho,lower = lower,upper = upper,maximum = T, Xstar = Xstar, 
                     Ystar = Ystar, D = D, k = k, m = m, Psi_L = Psi_L, Wstar = Wstar)
  if(opt.iter.show) cat(optRho$objective," ")
  llvalues <- optRho$objective
  Rho <- optRho$maximum
  Psi_L <-as.matrix(Matrix::nearPD(Psi_L)[[1]])  
  parPsi_L <- vechMat(t(chol(Psi_L))) 
  if(opt.iter < 2) par <- c(Rho,parPsi_L) else{
    for (i in 1:(opt.iter-1)) {
      optPsi_L <- optim(par = parPsi_L, fn = fnPsi_L, gr = NULL, Xstar = Xstar, 
                    Ystar = Ystar, D = D, k = k, m = m, Rho = Rho, Wstar = Wstar,
                    method = opt.iter.method, control = list(fnscale = -1), hessian = F)
      Psi_L <- par2V(optPsi_L$par,k)
      optRho <- optimize(f = fnRho,lower = lower,upper = upper,maximum = T, Xstar = Xstar, 
                         Ystar = Ystar, D = D, k = k, m = m, Psi_L = Psi_L, Wstar = Wstar)
      Rho <- optRho$maximum
      parPsi_L <- optPsi_L$par
      if(opt.iter.show) cat(optPsi_L$value,optRho$objective," ")
      llvalues <- c(llvalues,optPsi_L$value,optRho$objective)
    }
    par <- c(Rho,parPsi_L)
  }
  list(par=par,llvalues=llvalues)
}


Lmvmeta.ml <- function (Xstar, Ystar, D = D, Wstar = Wstar, k, m, p, bscov, 
                        nall,control, ...){
  fn <- mlprof.fn_L
  fnRho <- mlprof.fnRho_L
  fnPsi_L <- mlprof.fnPsi_L
  Psi_L <- control$initial_V
  yita <- 1
  lower <- -yita + 1e-6  
  upper <- yita - 1e-6  
  par <- opt.iter.UV_L(fnRho = fnRho,fnPsi_L = fnPsi_L,Psi_L = Psi_L,opt.iter = control$opt.iter,
                     opt.iter.method = control$opt.iter.method,Xstar = Xstar,
                     Ystar = Ystar,D = D,k = k,m = m,Wstar = Wstar,
                     lower = lower,upper = upper,opt.iter.show = control$opt.iter.show)
  llvalues <- par[[2]]
  par <- par[[1]]
  nparPsi_L <- length(par) - 1
  opt <- optim(par = par, fn = fn, gr = NULL, Xstar = Xstar,
               Ystar = Ystar, D = D, k = k, m = m, Wstar = Wstar,
               lower = c(lower,rep(-1 + 1e-6,nparPsi_L)), upper = c(upper,rep(1 - 1e-6,nparPsi_L)),
               method = "L-BFGS-B", control = control$optim, hessian = control$hessian)  
  if(control$opt.iter.show) cat(opt$value," ")
  llvalues <- c(llvalues,opt$value)
  Psi_L <- par2V(opt$par[-1],k)
  Rho <- opt$par[1]
  
  gls <- glsfit_L(Xstar, Ystar, Wstar, D, Psi_L, Rho, k, m, onlycoef = FALSE)
  
  qrinvtMX <- qr(gls$invtMX)
  R <- qr.R(qrinvtMX)
  Qty <- qr.qty(qrinvtMX, gls$invtMY)
  vcov <- tcrossprod(backsolve(R, diag(1, ncol(gls$invtMX))))
  
  res <- NULL
  fitted <- solve(diag(m*k)-Rho*Wstar) %*% Xstar %*% gls$coef
  rank <- qrinvtMX$rank
  converged <- opt$convergence == 0 | abs(llvalues[length(llvalues)]-llvalues[length(llvalues)-1]) < control$lltol
  
  c(list(coefficients = gls$coef, vcov = vcov, Psi_L = Psi_L, Rho = Rho,residuals = res, 
         fitted.values = fitted, df.residual = nall - rank - length(par), iter_llvalues = llvalues,
         rank = rank, logLik = opt$value, converged = converged, par = opt$par), if (!is.null(opt$hessian)) 
           list(hessian = opt$hessian, se.Rho = sqrt(-solve(opt$hessian)[1,1])), 
    list(niter = opt$counts[[2]], control = control))
}



Lmvmeta.reml <- function (Xstar, Ystar, D = D, Wstar = Wstar, k, m, p, bscov, 
                        nall,control, ...){
  fn <- remlprof.fn_L
  fnRho <- remlprof.fnRho_L
  fnPsi_L <- remlprof.fnPsi_L
  Psi_L <- control$initial_V
  lower <- -yita + 1e-6  
  upper <- yita - 1e-6 
  par <- opt.iter.UV_L(fnRho = fnRho,fnPsi_L = fnPsi_L,Psi_L = Psi_L,opt.iter = control$opt.iter,
                     opt.iter.method = control$opt.iter.method,Xstar = Xstar,
                     Ystar = Ystar,D = D,k = k,m = m,Wstar = Wstar,
                     lower = lower,upper = upper,opt.iter.show = control$opt.iter.show)
  llvalues <- par[[2]]
  par <- par[[1]]
  nparPsi_L <- length(par) - 1
  opt <- optim(par = par, fn = fn, gr = NULL, Xstar = Xstar,
               Ystar = Ystar, D = D, k = k, m = m, Wstar = Wstar,
               lower = c(lower,rep(-1 + 1e-6,nparPsi_L)), upper = c(upper,rep(1 - 1e-6,nparPsi_L)),
               method = "L-BFGS-B", control = control$optim, hessian = control$hessian)  
  if(control$opt.iter.show) cat(opt$value," ")
  llvalues <- c(llvalues,opt$value)
  Psi_L <- par2V(opt$par[-1],k)
  Rho <- opt$par[1]
  
  gls <- glsfit_L(Xstar, Ystar, Wstar, D, Psi_L, Rho, k, m, onlycoef = FALSE)
  
  qrinvtMX <- qr(gls$invtMX)
  R <- qr.R(qrinvtMX)
  Qty <- qr.qty(qrinvtMX, gls$invtMY)
  vcov <- tcrossprod(backsolve(R, diag(1, ncol(gls$invtMX))))
  
  res <- NULL
  fitted <- solve(diag(m*k)-Rho*Wstar) %*% Xstar %*% gls$coef
  rank <- qrinvtMX$rank
  converged <- opt$convergence == 0 | abs(llvalues[length(llvalues)]-llvalues[length(llvalues)-1]) < control$lltol
  
  c(list(coefficients = gls$coef, vcov = vcov, Psi_L = Psi_L, Rho = Rho,residuals = res, 
         fitted.values = fitted, df.residual = nall - rank - length(par), iter_llvalues = llvalues,
         rank = rank, logLik = opt$value, converged = converged, par = opt$par), if (!is.null(opt$hessian)) 
           list(hessian = opt$hessian, se.Rho = sqrt(-solve(opt$hessian)[1,1])), 
    list(niter = opt$counts[[2]], control = control))
}
remlprof.fn_L <- function (par,Xstar, Ystar, Wstar, D, k, m){
  Rho <- par[1]
  Psi_L <- par2V(par[-1],k) 
  gls <- glsfit_L(Xstar, Ystar, Wstar, D, Psi_L, Rho, k, m, onlycoef = FALSE)
  pconst <- -0.5 * (m * k - length(gls$coef)) * log(2 * pi)
  pres <- -0.5 * (crossprod(gls$invtMY - gls$invtMX %*% gls$coef))
  pdet1 <- -sum(log(diag(gls$M)))
  tXWXtot <- crossprod(gls$invtMX)
  pdet2 <- -sum(log(diag(chol(tXWXtot))))
  
  as.numeric(pconst + pdet1 + pdet2 + pres)
}
remlprof.fnRho_L <- function (par,Xstar, Ystar, Wstar, Psi_L, D, k, m){
  Rho <- par
  gls <- glsfit_L(Xstar, Ystar, Wstar, D, Psi_L, Rho, k, m, onlycoef = FALSE)
  pconst <- -0.5 * (m * k - length(gls$coef)) * log(2 * pi)
  pres <- -0.5 * (crossprod(gls$invtMY - gls$invtMX %*% gls$coef))
  pdet1 <- -sum(log(diag(gls$M)))
  tXWXtot <- crossprod(gls$invtMX)
  pdet2 <- -sum(log(diag(chol(tXWXtot))))
  
  as.numeric(pconst + pdet1 + pdet2 + pres)
}
remlprof.fnPsi_L <- function (par,Xstar, Ystar, Wstar, Rho, D, k, m){
  Psi_L <- par2V(par,k)
  gls <- glsfit_L(Xstar, Ystar, Wstar, D, Psi_L, Rho, k, m, onlycoef = FALSE)
  pconst <- -0.5 * (m * k - length(gls$coef)) * log(2 * pi)
  pres <- -0.5 * (crossprod(gls$invtMY - gls$invtMX %*% gls$coef))
  pdet1 <- -sum(log(diag(gls$M)))
  tXWXtot <- crossprod(gls$invtMX)
  pdet2 <- -sum(log(diag(chol(tXWXtot))))
  
  as.numeric(pconst + pdet1 + pdet2 + pres)
}


#Output overall ERR coeffients and their var/cov matrix 
Outoerr_L <- function(object,nsim=5000){
  
  coef <- object$coefficients ; vcov <- object$vcov
  Rho <- object$Rho ; se.Rho <- object$se.Rho
  Wstar <- object$Wstar ; Xstar <- object$Xstar
  dim <- object$dim
  
  oerr_coef <- object$fitted.values[1,]
  Z <- matrix(rnorm(length(coef)*nsim),nrow=nsim,ncol=length(coef))
  coefsim <- t(c(coef) + t(Z%*%chol(vcov)))
  rhosim <- rnorm(nsim,Rho,se.Rho)
  rhosimlim <- rhosim[rhosim>-1 & rhosim<1]
  coefsimlim <- coefsim[rhosim>-1 & rhosim<1,]
  simlim <- cbind(rhosimlim,coefsimlim)
  
  oerr_coefsim <- apply(simlim,1, function(coefi) {
    weight <- solve(diag(dim$k*dim$m)-coefi[1]*Wstar)[1:dim$k,]
    oerr_coefstari <- weight %*% Xstar %*% as.vector(coefi[-1])
    oerr_coefstari
  })
  
  pred <- NULL
  pred$coefsim <- oerr_coef
  pred$vcovsim <- cov(t(oerr_coefsim))
  pred$nsim <- nsim
  
  pred
  }


blup.Lmvmeta <- function (object,nsim=5000) 
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
  Rho <- object$Rho
  se.Rho <- object$se.Rho
  Wstar <- object$Wstar
  m <- nrow(y)
  k <- ncol(y)
  Xstar <- object$Xstar
  Ystar <- object$Ystar
  fittedstar <- as.numeric(t(object$fitted.values))
  
  U <- solve(diag(k*m)-Rho*Wstar)
  Psi_L <- object$Psi_L
  V <- diag(m) %x% Psi_L
  G <- U %*% V %*% t(U)
  invSIGMA <- chol2inv(chol(G+object$D))
  
  ytilde <- Ystar - fittedstar
  xi <- G %*% invSIGMA %*% ytilde
  epsilon <- ytilde - xi
  fitted.values.spatial <- fittedstar + xi
  
  coef <- object$coefficients
  vcov <- object$vcov
  Z <- matrix(rnorm(length(coef)*nsim),nrow=nsim,ncol=length(coef))
  coefsim <- t(c(coef) + t(Z%*%chol(vcov)))
  rhosim <- rnorm(nsim,Rho,se.Rho)
  rhosimlim <- rhosim[rhosim>-1 & rhosim<1]
  coefsimlim <- coefsim[rhosim>-1 & rhosim<1,]
  simlim <- cbind(rhosimlim,coefsimlim)
  oerr_coefsim <- apply(simlim,1, function(coefi) {
    weight <- solve(diag(k*m)-coefi[1]*Wstar)
    oerr_coefstari <- weight %*% Xstar %*% as.vector(coefi[-1])
    oerr_coefstari
  })
  
  covfit <- cov(t(oerr_coefsim)) + G - G %*% invSIGMA %*% G
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
