## roxygen2 package
##############################################################
##############################################################
##############################################################
##############################################################
coxlasso<- function(fix=formula, rnd = NULL, vary.coef = NULL, 
                       data, xi, adaptive.weights = NULL, control = list()){
  
  est <- est.coxlasso(fix=fix,rnd=rnd,vary.coef=vary.coef,data=data,
                         xi=xi,adaptive.weights=adaptive.weights,control=control)
  
  est$StdDev <- est$Q
  est$call <- match.call()
  class(est) <- "coxlasso"
  est
}



est.coxlasso <- function(fix, rnd, vary.coef, data, xi, adaptive.weights = NULL, control = list() )
{
  if(grepl("\\*", fix[3]))
    stop("Usage of '*' not allowed in formula! Please specify the corresponding variables separately.")  
  
  ## Print stuff.
  if(is.null(control$flushit))
    control$flushit <- TRUE

  ia <- if(control$flushit) interactive() else FALSE

  ic.dummy<-attr(terms(fix),"intercept")  
  
  Resp <- model.response(model.frame(fix, data))
  
  if(ncol(Resp)==2)
  {
    y0 <- rep(0,nrow(Resp))  
    y <- Resp[,1]
    delta <- Resp[,2]
  }else{
    y0 <- Resp[,1]  
    y <- Resp[,2]
    delta <- Resp[,3]
  }
  
  xl <- 0
  N <- length(y)
  
  ## loop index for problems in Fisher matrix: default is NULL as hopefully not needed
  loop.ind <- NULL
  
  ## read all specifications from control function
  control <- do.call(coxlassoControl, control)
  
  ## shall the knots be placed quantile-based?
  quant.knots <- control$quant.knots
  
  ## this is the Oelker & Tutz (2017) quadratic-approximation-of-the-penalty parameter 
  c.app <- control$c.app
  
  ## fixed effects below this threshold will be treated as zero
  eps <- control$eps
  
  ## what is the right boundary of time horizon
  if(is.null(control$xr))
  {
    xr <- max(y)
  }else{
    xr <- control$xr
  }
  
  ## save fix for coxFL
  fix.fix <- fix
  
  ## build design matrix without intercept
  if(ic.dummy)
    fix <- update(fix,~ .-1)

  ## check whether fixed effects are present 
  fix.names <- attr(terms(fix),"term.labels")
  
  ## change terms with as.factor(...) in formula into "real" factors
  ind2.fix <- sapply(fix.names,factor.test)
  if(sum(ind2.fix)>0)
  {
    fix.names.old <- fix.names
    fix.names[ind2.fix] <- sapply(fix.names[ind2.fix], remove.fac)
    data[fix.names[ind2.fix]] <- lapply(data[fix.names[ind2.fix]],as.factor)
    fix <- update(fix, paste("~ .-",paste(fix.names.old[ind2.fix],collapse="-"),"+",
                               paste(fix.names[ind2.fix],collapse="+")))
  }
  
  
  fix.temp <- update(fix,~ . +1)
  
  model.matrix.temp <- model.matrix(fix.temp,data)

  ## check wether there are any NAs in the design matrix
  if(length(fix.names)>0)
    na.vec <- apply(data[,fix.names,drop=FALSE],2,is.na.fct)
  
  if(sum(na.vec)>0){
    stop(paste("There are NAs in the following covariate(s): ", 
               paste(fix.names[na.vec>0],collapse="; ")))
  }
  
  orig.names <- colnames(model.matrix.temp)[-1]
  
  ## set up index vector for groups and for unpenalized fixed effects
  if(is.null(control$index))
    control$index<-1:length(fix.names)
  
  names(control$index) <- fix.names

  ## check whether some fixed effects are already factors
  if(ncol(data[,is.element(colnames(data),fix.names),drop=FALSE])>0)
  {
    ## check whether fixed effects formula contains factors
    ind1.fix <- sapply(data[,fix.names,drop=FALSE],is.factor)
    ## check whether for some fixed effects are grouped in the index vector
    ind3.order <- order(control$index); index.reor <- control$index[ind3.order]
    index.reor[is.na(index.reor)] <- c((max(index.reor,na.rm=T)+1):(max(index.reor,na.rm=T)+sum(is.na(index.reor))))
    ind3.help <- (table(index.reor)>1)
    ind3.fix <- rep(ind3.help, times = table(index.reor))[order(ind3.order)] 
    names(ind3.fix) <- fix.names
  }else{
    ind1.fix <- 0
    ind3.fix <- 0
  }
  
  
  ind.fix <- (ind1.fix + ind3.fix)>0
  
  number.cat.fix <- sum(ind.fix)
  number.noncat.fix <- length(fix.names)-number.cat.fix
  
  ## dummy, whether categorical covariates / groups are present
  cat.present.fix <- number.cat.fix>0
  noncat.present.fix <- number.noncat.fix>0 

  real.factor.names <- character()
  lev.fac <- list()
  
  ## formula for categorical predictors (always with intercept)
  if(cat.present.fix)
  {
    ## these are factor variables; this will then be used in the predict-function
    real.factor.names <- fix.names[ind1.fix]
    lev.fac <- lapply(data[real.factor.names],levels)
    
     fix.coef.cat <- formula(paste("~ 1+",paste(fix.names[ind.fix],collapse="+")))
  }else{
    fix.coef.cat <- formula("~ 1")
  }
  
  
  ## 2nd formula for metric predictors (always with intercept)
  if(noncat.present.fix)
  {
    fix.coef <- formula(paste("~ 1+",paste(fix.names[!ind.fix],collapse="+")))
  }else{
    fix.coef <- formula("~ 1")
  }
  
  ## design matrix metric covariates
  X <- model.matrix(fix.coef, data)
  metric.names <- attr(terms(fix.coef),"term.labels")
  
  ## design matrix categorical covariates
  X2 <- model.matrix(fix.coef.cat, data)
  cat.names <- attr(terms(fix.coef.cat),"term.labels")

  ## remove intercept from both design parts
  X <- X[,-1,drop=F]
  X2 <- X2[,-1,drop=F]
  m <- ncol(X)
  m2 <- ncol(X2)
  lin <- m + m2
  
  ## adjust index for categorical covariates => count levels -1
  rep.index <- rep(1, length(control$index))
  rep.index[ind1.fix] <- apply(data[,fix.names[ind1.fix],drop=F],2,nelem)-1
  names(rep.index) <- fix.names
  rep.index <- rep.index[c(metric.names,cat.names)]
  
  ## combine both matrices (note that now covariate columns are ordered: first metric, then categorical)
  X.both <- cbind(X,X2)
  
  ####### Center & Standardization
  center <- control$center
  standardize <- control$standardize
  
  ## re-order index vector
  index.new <- control$index[c(metric.names,cat.names)]
  index.new <- rep(index.new, times = rep.index); names(index.new) <- colnames(X.both)
  
  scale.obj <- scaleDesign(X = X.both, index.vec =  index.new, center = control$center, 
                            standardize = control$standardize)
  X.both <- scale.obj$X
  
  if(m>0)
  X <- X.both[,1:m] 
  
  if(m2>0)
  X2 <- X.both[,(m+1):lin]
  
  ## initialize adaptive weights, if null, run CoxFL and use betas
  if(is.null(adaptive.weights))
    #adaptive.weights <- rep(1,m+m2)
  {
    control.iter.temp <- control$print.iter
    control$print.iter <- FALSE
    # set control starting values to zero
    control.temp <- control
    control.temp$q_start <- NULL
    control.temp$start <- NULL
    control.temp$standardize <- TRUE
    control.temp$center <- FALSE
    try(coxFL.mod <- coxFL(fix=fix.fix,rnd=rnd,vary.coef=vary.coef,data=data,control=control.temp))
      if(inherits(class(coxFL.mod),"try-error"))
        stop("Unpenalized ML adaptive weights don't exist! Please use either coxridge in advance to get suitable adaptive weights or cv.coxlasso directly.")
    adaptive.weights <- abs(coxFL.mod$coef.stand)
    names(adaptive.weights) <- NULL
    control$print.iter <- control.iter.temp
  }
    
  ## stop if there are no covariates present
  if(ncol(X.both)==0)
    stop("No terms to select! Use coxph or coxme!")
  
  ## check if random effects are present
  if(!is.null(rnd))
  {  
      ## set up random formula
      rnd.len<-length(rnd)
      
      rndformula <- as.character(rnd)
      
      trmsrnd <- terms(rnd[[1]])
      
      if(!is.factor(data[,names(rnd)[1]]))
      {
        data[,names(rnd)[1]] <- as.factor(data[,names(rnd)[1]])
        warning("Cluster variable should be specified as a factor variable!")  
      }
      
      newrndfrml <- "~ -1"
      newrndfrml <- paste(newrndfrml,  if(attr(trmsrnd, "intercept")) names(rnd)[1] else "", sep=" + ")
      
      if(length(attr(trmsrnd, "variables"))>1)
      {
        newrndfrml <- paste(newrndfrml,  
                            paste(sapply(attr(trmsrnd,"term.labels"), function(lbl){
                              paste(lbl, names(rnd)[1], sep=":")
                            }), collapse=" + "), sep="+") 
      }
      
      ## random effects design matrix part
      W_start <- model.matrix(formula(newrndfrml), data)
      
      rnlabels <- terms(formula(newrndfrml))
      random.labels <- attr(rnlabels,"term.labels")
      ran.tab <- table(data[,colnames(data)==(names(rnd)[1])])   
      
      ## number clusters
      n <- length(ran.tab)
      
      ## number of random effect types (e.g. 1 if only intercept)
      s <- ncol(W_start)/n
      
      ## if more than one random effect type, structure random effect design matrix suitably
      if(s>1)
      {
        W<-W_start[,seq(from=1,to=1+(s-1)*n,by=n)]
        for (i in 2:n)
          W<-cbind(W,W_start[,seq(from=i,to=i+(s-1)*n,by=n)])
      }else{
        W<-W_start
      }
      subject.names<-names(rnd)  
      
      ## create starting values for random effects covariance matrix
      if(is.null(control$q_start) & sum(s)>1)
        control$q_start<-diag(rep(0.1,sum(s)))
      
      Q1 <- control$q_start
      
      dimr <- n*s
  }else{
    ## create RE objects such that code still works
    W <- NULL  
    Q1 <- NULL  
    dimr <- 0
    s <- 0
  }

  ## final design matrix fixed and random effect parts
  XW <- cbind(X.both,W)
  
  ## if no vary.coef are specified, create formula with intercept (for baseline hazard)
  if(is.null(vary.coef))
    vary.coef <- formula("~ 1")
    
  ## check whether time-varying effects are present 
  vary.names <- attr(terms(vary.coef),"term.labels")
  
  ## add time varying coef if some are present
  if(length(vary.names)>0)
  vary.coef <- formula(paste("~ 1+",paste(vary.names,collapse="+")))
  
  ## design matrix metric time-varying covariates
  U <- model.matrix(vary.coef, data)
  m.vary <- ncol(U)-1 
  
  ## preparation of spline parts for the baseline hazard (and all other time-varying effects)
  smooth <- control$smooth
  start <- control$start
  exact <- control$exact
  
  ## if value of exact is not specified by the user, set it automatically such that about 100 equally sized Rieman bars are obtained
  if (is.null(exact)){
    y.max <- max(y)
    if (y.max<1) exact <- round(y.max, digits = 2) / 100  
    if (y.max > 1 & y.max < 1000) exact <- floor(y.max) / 100  
    if (y.max > 1000) exact <- round(y.max, digits = -2) / 100
  }

  nbasis <- smooth$nbasis
  # diff.ord <- 2
  diff.ord <- 1
  spline.degree <- 3
  penal <- smooth$penal
  
  ## repeat each column nbasis-times
  U.each <- t(apply(U,1,rep,each=nbasis))
  
  ## build huge double matrix with nbasis^2 * ncol(U)^2 columns 
  # nbasis times ( 1..1 age ...age ) nbasis times (age ...age, age*age ...age*age) 
   U.double.each <- matrix(NA,N,nbasis^2*ncol(U)^2)
  for(j in 1:(m.vary+1))
    U.double.each[,((j-1)*(m.vary+1)*nbasis^2+1):(j*(m.vary+1)*nbasis^2)] <- t(apply(t(apply(U*U[,j],1,rep,each=nbasis)),1,rep,nbasis))
  
  ## generate penalty matrix
  D<-diag(nbasis)
  D<-diff(D);t.D<-t(D)
  
  ## prepare setting for B-splines
  xmin<-xl-(xr-xl)/100
  xmax<-xr+(xr-xl)/100
  dx<-(xmax-xmin)/(nbasis-3)
  knots<-seq(xmin-spline.degree*dx,xmax+spline.degree*dx,by=dx)
  quant.vec <- NULL
  if(quant.knots)
  {
    quant.vec <- quantile(y,seq(0,1,by=1/(nbasis-1))[-c(1,nbasis)])
    knots[3+(1:(nbasis-2))] <- quant.vec
  }
  knots<-knots[1:nbasis]
  
  ## prepare modified parametrization of the B-splines (see Fahrmeir et al., 2004; Kneib et al., 2009)
  Pen <- t.D%*%solve(D%*%t.D)
  B.unpen.fact<-rep(1, length(knots))
  
  K.base <- rep(c(0,rep(1,nbasis-1)),m.vary+1)
  K.base <- penal * K.base
  
  if(!(diff.ord<spline.degree))
    stop("Order of differences must be lower than degree of B-spline polynomials!")
  
  ## set up classical B-spline design matrix evaluated at event times
  Phi.base <- bs.design.asy(y,xl=xl,xr=xr,spline.degree=spline.degree,nbasis=nbasis, quant = quant.vec)  
  
  ## use modified parametrization of the b-splines (see Fahrmeir et al., 2004; Kneib et al., 2009)
  Phi.base  <-cbind(Phi.base%*%B.unpen.fact,Phi.base%*%Pen)
  colnames(Phi.base) <- paste("baseline",rep(1:dim(Phi.base)[2],each=1), sep=".")
  
  ## multiply B-spline design matrix with covariates (changes nothing for baseline hazard)
  ## v_ijk := z_ijk*B(t) (Biometrics page 2 above equ. 1)
  Phi<-matrix(0,N,(m.vary+1)*nbasis)
  for(j in 1:(m.vary+1))
    Phi[,((j-1)*nbasis+1):(j*nbasis)] <- Phi.base * U[,j]

  colnames(Phi) <- paste0(rep(colnames(U),each = nbasis),rep(1:nbasis,times=m.vary+1))
  dimb <- ncol(Phi)
  
  ## 
  time.grid <- seq(0,max(y),by=exact)
  time.length <- length(time.grid)
  
  ## huge matrix specifying which fine time grid points have to be integrated for each observation
  T.mat <- matrix(NA,time.length,N)
  event.ma <- cbind(y0,y)
  T.mat[1:length(time.grid),1:N] <- apply(event.ma,1,function(z,time.grid){(time.grid<z[2] & time.grid>=z[1])},time.grid=time.grid)
  
  ## set up classical b-spline design matrix evaluated at integral grid times
  Phi.big.temp <- bs.design.asy(time.grid,xl=xl,xr=xr,spline.degree=spline.degree,nbasis=nbasis, quant = quant.vec)
  colnames(Phi.big.temp) <- paste("baseline",rep(1:dim(Phi.base)[2],each=1), sep=".")

  ## use modified parametrization of the B-splines (see Fahrmeir et al., 2004; Kneib et al., 2009)
  Phi.big.temp  <-cbind(Phi.big.temp%*%B.unpen.fact,Phi.big.temp%*%Pen)
  
  Phi.big <- t(apply(Phi.big.temp,1,rep,m.vary+1))
  
  Phi.double.big <- matrix(NA,length(time.grid),(m.vary+1)*nbasis^2)
  for(j in 1:nbasis)
    Phi.double.big[,((j-1)*nbasis*(m.vary+1)+1):(j*nbasis*(m.vary+1))] <- t(apply(Phi.big.temp*Phi.big.temp[,j],1,rep,m.vary+1))
  Phi.double.big <- (t(apply(Phi.double.big,1,rep,m.vary+1)))
  
  ## define transposed matrices
  t.U.each <- t(U.each)
  t.XW <- t(XW)

  Design <- cbind(Phi,XW)
  
  ## define products
  delta.Phi <- delta %*% Phi
  delta.fixran <- t.XW %*% delta

  if(is.null(start) & is.null(rnd))
    start <- c(log(sum(delta)/sum(y-y0)),rep(0,dimb + lin + dimr-1))

  if(!is.null(rnd) & (is.null(Q1) | is.null(start)))
  {
    if(s==1)
    {
      if(is.null(Q1))
      {
        Q1 <- 0.01  
      }
      if(is.null(start))
      {
        start <- c(log(sum(delta)/sum(y-y0)),rep(0,dimb + lin + dimr-1))
      }
    }else{
      stop("Provide starting values for q_start and start!")
    }  
  }
  
  ## permut coefs of starting vector
  permut.ind <- match(colnames(X.both),orig.names)
  start[(dimb+1):(dimb+lin)] <- start[(dimb+1):(dimb+lin)][permut.ind]
  
  Delta <- matrix(NA,control$max.iter+1, dimb + lin + dimr)
  colnames(Delta) <- c(colnames(Phi),colnames(X.both),colnames(W))

  Delta[1,] <- start
  
  ##
  
  Q<-list()
  Q.base<-numeric()
  
  if(m.vary>0)
    Q.vary <- list()
  
  Q[[1]]<-Q1
  Q.base[1]<-1/penal
  
  if(m.vary>0)
    Q.vary[[1]] <- rep(1/penal,m.vary) 
  
  alpha <- start[1:dimb]
  
  if(lin>0)
  {
    fixef <- start[(1+dimb):(dimb+lin)]
  }else{
    fixef <- NULL
  }
  
  if(dimr>0)
  {
    ranef <- start[(1+dimb+lin):(dimb+lin+dimr)]
  }else{
    ranef <- NULL
  }
  
  alphama <- matrix(alpha,ncol(U.each),nrow(U.each))
  
  ## start with estimation
  if(control$print.iter)
     {
    cat(if(ia) "\r" else NULL)
    cat("Iteration  1")
    if(.Platform$OS.type != "unix" & ia) flush.console()
  }
  
  multi.obj <- IntegrMultiCpp(alphama, c(fixef,ranef), Phi.big,
                              t.U.each, U.each, T.mat, Phi.double.big, U.double.each, XW)
  multi.obj$int.ma <- exact*multi.obj$int.ma
  multi.obj$int.array <- exact*multi.obj$int.array
  multi.obj$eta <- drop(multi.obj$eta)
  
    ## intitialize loglik, starting values
  loglik <- numeric()
  loglik.temp <- sum(delta * (Design %*% start)) - exact * sum(multi.obj$int.likeli) 
  loglik <- c(loglik,loglik.temp)
  
  int.vec.fix <-  exact * rowSums(multi.obj$help.calc)
  
  Fisher.upp <- matrix(multi.obj$int.array, dimb, dimb)
  
  P <- matrix(0,dimb+lin+dimr,dimb+lin+dimr)  

  ## simple ridge penalty on baseline & all time-varying effects
  diag(P)[1:((m.vary+1)*nbasis)] <- K.base
  
  ## here comes the LASSO on the (metric) fixed effects
   if(m>0)
     diag(P)[(dimb+1):(dimb+m)] <- xi * (fixef[1:m]^2 + c.app)^(-0.5) / adaptive.weights[1:m]

  ## Group LASSO for categorical variables
  if (m2>0)
  {
      # vec incl number of levels -1 for each cat var 
      levels.vec <- numeric()
      for(i in 1:number.cat.fix)
        levels.vec[i] <- nlevels(data[,real.factor.names[i]])-1
      cum.levels <- c(0,cumsum(levels.vec))

      # for each cat var create A.gr
      for(k in 1:number.cat.fix)
      {
      beta.vec <- fixef[(cum.levels[k] + m+1):(cum.levels[k+1]+ m)]
      
      A.gr <-  sqrt(levels.vec[k]) * (sum(beta.vec^2)+c.app)^(-0.5) * diag(length(beta.vec))
      # additional term sqrt(levels.vec[k]) not from Oelker & Tutz, but needed to adjust group lasso to normal lasso terms
      
      # write A.gr in penalty matrix
      P[(cum.levels[k] + dimb+m+1):(cum.levels[k+1]+ dimb+m), (cum.levels[k] + dimb+m+1):(cum.levels[k+1]+ dimb+m)] <- xi * A.gr /  eucl.norm(adaptive.weights[(cum.levels[k]+m+1):(cum.levels[k+1]+m)])
      }
   }  
  

  ## update RE covariance matrix (if RE are present)
  if(!is.null(rnd))
  {
      if(s==1)
      {
        diag(P)[(dimb+lin+1):(dimb+lin+dimr)] <- 1/Q1
      }else{
        Q_inv.start <- chol2inv(chol(Q1))
        for(jf in 1:n)
          P[(dimb+lin+(jf-1)*s+1):(dimb+lin+jf*s),(dimb+lin+(jf-1)*s+1):(dimb+lin+jf*s)] <- Q_inv.start
      }
  }


  ## calculate score vector
  score <- c(delta.Phi - multi.obj$t.eta %*% multi.obj$int.ma,
             delta.fixran -t.XW %*% (int.vec.fix * multi.obj$eta)) - P %*% start
  
  xi.max <- max(abs(score[(dimb+1):(dimb+lin)]))
  
  ## calculate Fisher matrix
  Fisher <- matrix(NA,dimb+lin+dimr,dimb+lin+dimr)
  
  Fisher[1:dimb,1:dimb] <- -Fisher.upp
  
  Fisher[(dimb+1):(dimb+lin+dimr),1:dimb] <-  - t.XW %*% (multi.obj$int.ma*multi.obj$eta)
  
  Fisher[1:dimb,(dimb+1):(dimb+lin+dimr)] <- 
    t(Fisher[(dimb+1):(dimb+lin+dimr),1:dimb])
  
  Fisher[(dimb+1):(dimb+lin+dimr),(dimb+1):(dimb+lin+dimr)] <- 
    -t.XW %*% (XW*multi.obj$eta*int.vec.fix)
  
  ## try to invert Fisher matrix
  InvFisher<-try(chol2inv(chol(P - Fisher)),silent=T)
  if(inherits(class(InvFisher)[1],"try-error"))
    InvFisher<-try(solve(P - Fisher),silent=T)
  if(inherits(class(InvFisher)[1],"try-error") || sum(is.na(InvFisher))>0)
    stop("Current Fisher matrix not invertible")  
  
  ## Fisher scoring update of all estimates
  Delta[2, ] <- start + InvFisher %*% score
  
  alpha <- Delta[2,1:dimb]
  if(lin>0)
  {
    fixef <- Delta[2,(dimb+1):(dimb+lin)]
  }else{
    fixef <- NULL
  }
  
  ### random effects variance
  if(!is.null(rnd))
  {
    ranef <- Delta[2,(1+dimb+lin):(dimb+lin+dimr)]

    if(s==1)
    {
      Q1 <- sum(diag(InvFisher)[(dimb+lin+1):(dimb+lin+dimr)])+sum(ranef^2)
    }else{ 
      Q1<-InvFisher[(dimb+lin+1):(dimb+lin+s),(dimb+lin+1):(dimb+lin+s)]+ranef[1:s]%*%t(ranef[1:s])
      for (i in 2:n)
        Q1<-Q1+InvFisher[(dimb+lin+(i-1)*s+1):(dimb+lin+i*s),(dimb+lin+(i-1)*s+1):(dimb+lin+i*s)]+ranef[((i-1)*s+1):(i*s)]%*%t(ranef[((i-1)*s+1):(i*s)])
    }  
    Q1<-Q1/n
  }  
  
  Q[[2]]<-Q1
  
  ### ran eff baseline
  ranef.base <- Delta[2,2:nbasis]
  Q1.base <- (sum(diag(InvFisher)[2:nbasis])+sum(ranef.base^2))/(nbasis-1)
  Q.base[2]<-Q1.base
  
  ### ran eff vary
  if(m.vary>0)
  {
    Q1.vary <- numeric()
    for(r in 1:m.vary) 
    Q1.vary[r] <- (sum(diag(InvFisher)[(r*nbasis+2): ((r+1)*nbasis)]) + sum(Delta[2, (r*nbasis+2): ((r+1)*nbasis)]^2))/(nbasis-1)
    
    Q.vary[[2]] <- Q1.vary
  }
  
  
  ########################################################################
  ########################################################################
  ########################################################################
  for(ll in 2:control$max.iter)
  {

    if(control$print.iter)
    {
      cat(if(ia) "\r" else if(ll > 1) "\n" else NULL)
      cat(paste("Iteration ",ll))
      if(.Platform$OS.type != "unix" & ia) flush.console()
    }
    
    alphama <- matrix(alpha,ncol(U.each),nrow(U.each))
    multi.obj <- IntegrMultiCpp(alphama, c(fixef,ranef), Phi.big,
                                t.U.each, U.each, T.mat, Phi.double.big, U.double.each, XW)
    multi.obj$int.ma <- exact*multi.obj$int.ma
    multi.obj$int.array <- exact*multi.obj$int.array
    multi.obj$eta <- drop(multi.obj$eta)
    
    ## loglik for starting values
    loglik.temp <- sum(delta * (Design %*% Delta[ll,])) - exact * sum(multi.obj$int.likeli) 
    loglik <- c(loglik,loglik.temp)
    
    int.vec.fix <-  exact * rowSums(multi.obj$help.calc)
    
    Fisher.upp <- matrix(multi.obj$int.array, dimb, dimb)
    
    ## update P 
    ## here comes the LASSO part
    if(m>0)
      diag(P)[(dimb+1):(dimb+m)] <- xi * (fixef[1:m]^2 + c.app)^(-0.5) / adaptive.weights[1:m]
    
    ## Group LASSO for categorical variables
    if (m2>0)
    {
      # for each cat var create A.gr
      for(k in 1:number.cat.fix)
      {
        beta.vec <- fixef[(cum.levels[k] + m+1):(cum.levels[k+1]+ m)]
        
        A.gr <-  sqrt(levels.vec[k]) * (sum(beta.vec^2)+c.app)^(-0.5) * diag(length(beta.vec))
        # additional term sqrt(levels.vec[k]) not from Oelker & Tutz, but needed to adjust group lasso to normal lasso terms
        
        # write A.gr in penalty matrix
        P[(cum.levels[k] + dimb+m+1):(cum.levels[k+1]+ dimb+m), (cum.levels[k] + dimb+m+1):(cum.levels[k+1]+ dimb+m)] <- xi * A.gr /  eucl.norm(adaptive.weights[(cum.levels[k]+m+1):(cum.levels[k+1]+m)]) 
      }
      
    }  
    
  
    ## update RE covariance matrix (if RE are present)
    if(!is.null(rnd))
    {
      if(s==1)
      {
        diag(P)[(dimb+lin+1):(dimb+lin+dimr)] <- 1/Q1
      }else{
        Q_inv <- chol2inv(chol(Q1))
        for(jf in 1:n)
          P[(dimb+lin+(jf-1)*s+1):(dimb+lin+jf*s),(dimb+lin+(jf-1)*s+1):(dimb+lin+jf*s)] <- Q_inv
      }
    }
    
    
    ## update baseline Random Intercept penalty matrix 
    diag(P)[2:nbasis] <- 1/Q1.base
    
    ## update varying coef Random Intercept penalty matrix
    if(m.vary>0)
    {
      for(r in 1:m.vary) 
       diag(P)[(r*nbasis+2): ((r+1)*nbasis)] <- 1/Q1.vary[r]  
    }
    
    
    score <- c(delta.Phi - multi.obj$t.eta %*% multi.obj$int.ma,
               delta.fixran -t.XW %*% (int.vec.fix * multi.obj$eta)) - P %*% Delta[ll,]
    
    Fisher <- matrix(NA,dimb+lin+dimr,dimb+lin+dimr)
    
    Fisher[1:dimb,1:dimb] <- -Fisher.upp
    
    Fisher[(dimb+1):(dimb+lin+dimr),1:dimb] <-  - t.XW %*% (multi.obj$int.ma*multi.obj$eta)
    
    Fisher[1:dimb,(dimb+1):(dimb+lin+dimr)] <- 
      t(Fisher[(dimb+1):(dimb+lin+dimr),1:dimb])
    
    Fisher[(dimb+1):(dimb+lin+dimr),(dimb+1):(dimb+lin+dimr)] <- 
      -t.XW %*% (XW*multi.obj$eta*int.vec.fix)
    
    ###
    InvFisher<-try(chol2inv(chol(P - Fisher)),silent=T)
    if(inherits(class(InvFisher)[1],"try-error"))
      InvFisher<-try(solve(P - Fisher),silent=T)
    if(inherits(class(InvFisher)[1],"try-error") || sum(is.na(InvFisher))>0)
      stop("Current Fisher matrix not invertible")  

    Delta[ll+1, ] <- Delta[ll, ] + InvFisher %*% score
    
    ### random effects variance
    alpha <- Delta[ll+1,1:dimb]
    if(lin>0)
    {
      fixef <- Delta[ll+1,(dimb+1):(dimb+lin)]
    }else{
      fixef <- NULL
    }
    
    ### random effects variance
    if(!is.null(rnd))
    {  
        ranef <- Delta[ll+1,(1+dimb+lin):(dimb+lin+dimr)]
        
        if(s==1)
        {
          Q1 <- sum(diag(InvFisher)[(dimb+lin+1):(dimb+lin+dimr)])+sum(ranef^2)
        }else{ 
          Q1<-InvFisher[(dimb+lin+1):(dimb+lin+s),(dimb+lin+1):(dimb+lin+s)]+ranef[1:s]%*%t(ranef[1:s])
          for (i in 2:n)
            Q1<-Q1+InvFisher[(dimb+lin+(i-1)*s+1):(dimb+lin+i*s),(dimb+lin+(i-1)*s+1):(dimb+lin+i*s)]+ranef[((i-1)*s+1):(i*s)]%*%t(ranef[((i-1)*s+1):(i*s)])
        }  
        Q1<-Q1/n
    }
    
    Q[[ll+1]]<-Q1
    
    
    ### ran eff baseline
    #ranef.base <- Delta[ll+1,3:nbasis]
    ranef.base <- Delta[ll+1,2:nbasis]
    Q1.base <- (sum(diag(InvFisher)[2:nbasis])+sum(ranef.base^2))/(nbasis-1)
    Q.base[ll+1]<-Q1.base
    
    
    ### ran eff varying coef   
    if(m.vary>0)
    {
      Q1.vary <- numeric()
      for(r in 1:m.vary) 
      Q1.vary[r] <- (sum(diag(InvFisher)[(r*nbasis+2): ((r+1)*nbasis)]) + sum(Delta[ll+1, (r*nbasis+2): ((r+1)*nbasis)]^2))/(nbasis-1)
      Q.vary[[ll+1]] <- Q1.vary
    }
    
    
    ### convergence criteria   
    
    finish <- eucl.norm(Delta[ll, ] - Delta[ll+1, ])/eucl.norm(Delta[ll, ]) < control$conv.eps
    
    if(finish)
      break
  }

  ##### calculate df's of final fit and AIC & BIC
  index.fix <- c((dimb+1):(dimb+lin))
  InvFisher.small<-try(chol2inv(chol((P - Fisher)[-index.fix,-index.fix])),silent=T)
  if(inherits(class(InvFisher.small)[1],"try-error"))
    InvFisher.small<-try(solve((P - Fisher)[-index.fix,-index.fix]),silent=T)
  if(inherits(class(InvFisher.small)[1],"try-error") || sum(is.na(InvFisher.small))>0)
    stop("Current Fisher matrix not invertible")  
  
  df.smooth <- sum(diag(-Fisher[-index.fix,-index.fix]%*%InvFisher.small))
  
  df.fixef <- sum(abs(fixef) > eps)
  
  df <- df.smooth + df.fixef

  multi.obj <- IntegrMultiCpp(alphama, c(fixef,ranef), Phi.big,
                              t.U.each, U.each, T.mat, Phi.double.big, U.double.each, XW)
  
  loglik.temp <- sum(delta * (Design %*% Delta[ll+1,])) - exact * sum(multi.obj$int.likeli) 
  loglik <- c(loglik,loglik.temp)
  
  ## once again loglik as only nu * update
  AIC <- -2*loglik.temp + 2*df
  BIC <- -2*loglik.temp + log(N)*df
  
  if(!is.null(rnd) && s==1)
    Q1<-sqrt(Q1)
  
  ada.w.metric <- NULL
  ada.w.cat <- NULL

  if(xi==0)
  {
    if(m>0)  
      ada.w.metric <- 1/abs(fixef[1:m])
    
    if(cat.present.fix)
    {
      index.cat <- unique(index.new[colnames(X2)])
      ada.w.cat <- rep(NA,length(index.cat))

      for(jv in index.cat)
        ada.w.cat[jv] <- 1/eucl.norm(fixef[index.new==jv])
    }
  }
  
  ## if time-varying effects, create final object
  time.vary <- NULL
  if(m.vary>0)
  {  
    time.vary <- as.numeric(alpha[(nbasis+1):((m.vary+1)*nbasis)])
    names(time.vary) <- paste(rep(colnames(U)[-1],each=nbasis),rep(1:nbasis,m.vary),sep=".")
  }
  
  ##if no time-varying
  if(m.vary == 0)
  {
    Q1.vary <- NULL
    Q.vary <- NULL
  }
  
  base.haz <- base.haz.cen <- as.numeric(alpha[1:nbasis])
  
  if(m>0)
    which.notsel <- abs(fixef[1:m]) < eps
  
  if(m2>0)
  {
    which.notsel.cat <- rep(FALSE,m2)
    for(k in 1:number.cat.fix)
    {
      fix.temp <- fixef[(cum.levels[k] + m+1):(cum.levels[k+1]+ m)]
      if(sqrt(sum(fix.temp^2)*levels.vec[k]) < eps)
        which.notsel.cat[(cum.levels[k] + 1):(cum.levels[k+1])] <- TRUE
    }
  }
  
  ## Transform the coefficients back to the original scale if the design
  ## matrix was standardized
  if(control$standardize){
    if(scale.obj$any.notpen)
      fixef[scale.obj$inotpen.which] <- (1 / scale.obj$scale.notpen) * fixef[scale.obj$inotpen.which]

    ## For df > 1 we have to use a matrix inversion to go back to the
    ## original scale
    for(j in 1:length(scale.obj$ipen.which)){
      ind <- scale.obj$ipen.which[[j]]
      fixef[ind] <- solve(scale.obj$scale.pen[[j]], fixef[ind,drop = FALSE])
     }
    if(s>1)
      warning("Random slopes are not standardized back!")
  }
  
  ## Now, after re-standardization, set those coefficients that were not selected to exactly zero
  if(m>0)
    fixef[1:m][which.notsel] <- 0
  
  if(m2>0)
    fixef[(m+1):(m + m2)][which.notsel.cat] <- 0
  
  coef <- fixef

  ## Need to adjust intercept if we have performed centering
  if(control$center){
    base.haz[1] <- base.haz[1] - sum(fixef * scale.obj$mu.x)
  }
  
    
  ## list of returnings 
  ret.obj <- list()
  ret.obj$InvFisher <- InvFisher 
  ret.obj$Delta <- Delta 
  ret.obj$baseline <- base.haz
  ret.obj$baseline.cen <- base.haz.cen
  ret.obj$time.vary <- time.vary
  ret.obj$coefficients <- coef
  ret.obj$ranef <- ranef
  ret.obj$nbasis <- nbasis
  ret.obj$spline.degree <- spline.degree
  ret.obj$diff.ord <- diff.ord
  ret.obj$quant.vec <- quant.vec
  ret.obj$penal <- penal
  ret.obj$Q_long <- Q
  ret.obj$Q <- Q1
  ret.obj$Q_long_base <- Q.base
  ret.obj$Q_base <- Q1.base
  ret.obj$Q_long_vary <- Q.vary
  ret.obj$Q_vary <- Q1.vary
  ret.obj$D <- D
  ret.obj$fix.names <- fix.names
  ret.obj$vary.names <- vary.names
  ret.obj$real.factor.names <- real.factor.names
  ret.obj$lev.fac <- lev.fac
  ret.obj$dim.rndeff <- s
  ret.obj$number.effects <- m
  ret.obj$iter <- ll+1
  ret.obj$time.grid <- time.grid
  ret.obj$exact <- exact
  ret.obj$X <- X.both
  ret.obj$Phi.big <- Phi.big
  ret.obj$ada.w.metric <- ada.w.metric
  ret.obj$ada.w.cat <- ada.w.cat
  ret.obj$knots <- knots
  ret.obj$rnd <- rnd
  ret.obj$B.unpen.fact <- B.unpen.fact
  ret.obj$Pen <- Pen
  ret.obj$m <- m
  ret.obj$m2 <- m2
  ret.obj$m.vary <- m.vary
  ret.obj$center <- control$center
  ret.obj$standardize <- control$standardize
  ret.obj$event <- delta
  ret.obj$event.time <- y
  ret.obj$fix <- fix
  ret.obj$vary.coef <- vary.coef
  ret.obj$data <- data
  ret.obj$fix.coef.cat <- fix.coef.cat
  ret.obj$fix.coef <- fix.coef
  ret.obj$cat.present.fix <- cat.present.fix
  ret.obj$loop.ind <- loop.ind
  ret.obj$loglik <- loglik
  ret.obj$df <- df
  ret.obj$xi.max <- xi.max
  ret.obj$adaptive.weights <- adaptive.weights
  ret.obj$AIC <- AIC
  ret.obj$BIC <- BIC
  return(ret.obj)
}

######### Methods

print.coxlasso <- function(x, ...)
{
  cat("Call:\n")
  print(x$call)
  cat("\nFixed Effects:\n")
  if(!is.null(x$coefficients)){ 
    cat("\nCoefficients:\n") 
    print(x$coefficients) }else{ 
      cat("\nNo time-constant effects included!\n")
    }
  
  cat("\nSmooth (time-varying) Effects:\n")
  if(!is.null(x$time.vary))
  {
    print(x$fix.names)
  }else{
    cat("\nNo time-varying effects included!\n") 
  }  
  
  if(!is.null(x$rnd))
  {  
    cat("\nRandom Effects:\n")
    if(x$dim.rndeff==1) cat("\nStdDev:\n") else cat("\nCov:\n") 
    print(x$StdDev)
  }else{
    cat("\nNo random effects included!\n")
  }
}

####

summary.coxlasso <- function(object, ...)
{
  coef <- object$coefficients
  baseline.coef <- object$baseline
  time.vary <- object$time.vary
  ranef <- object$ranef
  res <- list(call=object$call,
              coefficients=coef,baseline=baseline.coef,
              time.vary=time.vary,ranef=ranef,StdDev=object$StdDev,rnd=object$rnd,
              vary.names=object$vary.names,dim.rndeff=object$dim.rndeff)
  class(res) <- "summary.coxlasso"
  res
}

####

print.summary.coxlasso <- function(x, ...)
{
  cat("Call:\n")
  print(x$call)
  cat("\nFixed Effects:\n")
  if(!is.null(x$coefficients)){ 
    cat("\nCoefficients:\n") 
    print(x$coefficients) }else{ 
      cat("\nNo time-constant effects included!\n")
    }
  
  cat("\nSmooth (time-varying) Effects:\n")
  if(!is.null(x$time.vary))
  {
    print(x$vary.names)
    cat("\nTime-varying effects coefficients:\n")
    print(x$time.vary)
  }else{
    cat("\nNo time-varying effects included!\n") 
  }  
  
  if(!is.null(x$rnd))
  {  
    cat("\nRandom Effects Variance Components:\n")
    if(x$dim.rndeff==1) cat("\nStdDev:\n") else cat("\nCov:\n") 
    print(x$StdDev)
    cat("\nRandom Effects:\n")
    print(x$ranef)
  }else{
    cat("\nNo random effects included!\n")
  }
}


plot.coxlasso <- function(x, which.comp = NULL, main = NULL, ylim = NULL, ...)
{
  help.list <- list(...)
  
  if(is.null(help.list))
    help.list<-list()
  
  show.dummy <- FALSE
  if(!is.null(help.list$show.cens))
    show.dummy <- TRUE
  
  if(is.null(help.list$n.grid))
    help.list$n.grid <- 1000
  
  time.seq <- seq(0,max(x$time.grid),length.out=help.list$n.grid)

  Design<-bs.design.asy(time.seq, xl = 0, xr = max(time.seq), spline.degree=x$spline.degree, 
                        nbasis = x$nbasis, quant = x$quant.vec)
  Design  <-cbind(Design%*%x$B.unpen.fact,Design%*%x$Pen)
  
  name.vec <- c("baseline",x$vary.names)
  m<-length(name.vec)
  
  smooth.effects <- Design%*%matrix(c(x$baseline.cen,x$time.vary),x$nbasis,m)
  
  ## transform baseline by exp(.)
  ## AG: 16.11.2023: This is absolutely NOT compatible with the confidence bands that are calculated later on!!!
  smooth.effects[,1] <- exp(smooth.effects[,1])
  
  ## get variance
  var.smooth <- matrix(nrow = help.list$n.grid, ncol = m)
  for(r in 1:m){ 
    # same as diag(B(t) Cov(beta) B(t)')
    var.smooth[,r] <- rowSums((Design %*% x$InvFisher[((r-1) *x$nbasis +1): (r*x$nbasis), ((r-1) *x$nbasis +1): (r*x$nbasis)]) * Design) 
  }
  
  
  if(!exists("plot.data"))
    plot.data <- TRUE
  
  if(is.null(main))
    main<-""
  
  if(is.null(which.comp))
    which.comp<-1:m
  
  p<-length(which.comp)
  if(p>9)
    stop("Too many smooth functions! Please specify at maximum nine.")
  
  a<-ceiling(sqrt(p))
  b<-round(sqrt(p))
  if(b==0)
    b<-1
  
  max.event <- max(x$event.time[x$event==1])
  
  par(mfrow=c(a,b),mar=c(5, 5.5, 4, 2) + 0.1)
  
  x.lim <- NULL
  if(!show.dummy)
    x.lim <- c(0,max.event)
  
  for(i in which.comp)
  {
    if(!show.dummy)
    {
      if(is.null(ylim))
        ylim <- c(min(smooth.effects[,i][time.seq<=max.event]),max(smooth.effects[,i][time.seq<=max.event]))
    }  
    plot(x=time.seq,y=smooth.effects[,i],type="l", lwd=2, ylab=name.vec[i],
         main=main,xlab="time",xlim=x.lim,ylim=ylim,...) 
    lines(x = time.seq, y = smooth.effects[,i] +  qnorm(0.975) * sqrt(var.smooth[,i]), lty = "dashed")
    lines(x = time.seq, y = smooth.effects[,i] - qnorm(0.975) * sqrt(var.smooth[,i]), lty = "dashed")
    
    if(plot.data)
      rug(jitter(x$event.time[x$event==1]))
    
  }
}


####

predict.coxlasso <- function(object, newdata=NULL, time.grid=NULL, id.var = NULL, surv.ma = FALSE, ...)
{
  if(is.null(newdata))
    newdata <- object$data 

    if(!is.data.frame(newdata))
    stop("newdata needs to be a data frame!")
  
  data <- newdata  
  
  ## make sure that those variables that have been factors are now also factors with same levels
  if(object$cat.present.fix)
  {
    data[object$real.factor.names] <- lapply(data[object$real.factor.names],as.factor)
    for(j in 1:length(object$real.factor.names))
      attr(data[[object$real.factor.names[j]]],'levels') <- object$lev.fac[[j]]
  }
  
  if(length(object$fix[[2]])==3)
  {
  time.name <- as.character(object$fix[[2]][[2]])
  event.name <- as.character(object$fix[[2]][[3]])
  }else{
  time.name <- as.character(object$fix[[2]][[3]])
  event.name <- as.character(object$fix[[2]][[4]])
  }

  if(is.null(time.grid))
    time.grid <- seq(0,max(data[,time.name]),length.out=1000)
  time.length <- length(time.grid)
  
  data[,time.name][is.na(data[,time.name])] <- max(time.grid)
  data[,event.name][is.na(data[,event.name])] <- 1

  Resp <- model.response(model.frame(object$fix, data))
  
  if(ncol(Resp)==2)
  {
    y0 <- rep(0,nrow(Resp))  
    y <- Resp[,1]
    delta <- Resp[,2]
  }else{
    y0 <- Resp[,1]  
    y <- Resp[,2]
    delta <- Resp[,3]
  }
  N <- length(y)
  
  Design<-bs.design.asy(time.grid, xl = 0, xr = max(time.grid), spline.degree=object$spline.degree, 
                        nbasis = object$nbasis, quant = object$quant.vec)
  Design  <-cbind(Design%*%object$B.unpen.fact,Design%*%object$Pen)
  
  smooth.effects <- Design%*%matrix(c(object$baseline.cen,object$time.vary),object$nbasis,object$m.vary+1)

  ## design matrix metric covariates
  X <- model.matrix(object$fix.coef, data)
  metric.names <- attr(terms(object$fix.coef),"term.labels")
  
  ## design matrix categorical covariates
  X2 <- model.matrix(object$fix.coef.cat, data)
  cat.names <- attr(terms(object$fix.coef.cat),"term.labels")
  
  ## remove intercept from both design parts
  X <- X[,-1,drop=F]
  X2 <- X2[,-1,drop=F]
  X.both <- cbind(X,X2)
  
  if(ncol(X.both)!=(object$m+object$m2))
    stop("Some variables are missing in newdata!")

  ## set up random formula
  if(!is.null(object$rnd))
  {  
    rnd.len<-length(object$rnd)
    
    rndformula <- as.character(object$rnd)
    
    trmsrnd <- terms(object$rnd[[1]])
    newrndfrml <- "~ -1"
    newrndfrml <- paste(newrndfrml,  if(attr(trmsrnd, "intercept")) names(object$rnd)[1] else "", sep=" + ")
    
    if(length(attr(trmsrnd, "variables"))>1)
    {
      newrndfrml <- paste(newrndfrml,  
                          paste(sapply(attr(trmsrnd,"term.labels"), function(lbl){
                            paste(lbl, names(object$rnd)[1], sep=":")
                          }), collapse=" + "), sep="+") 
    }
    
    data[,colnames(data)==(names(object$rnd)[1])] <- as.factor(as.character(data[,colnames(data)==(names(object$rnd)[1])]))
    subject.names<-names(object$rnd)  
    
    if(length(levels(data[,colnames(data)==(names(object$rnd)[1])]))>1)
    {
    W_start <- model.matrix(formula(newrndfrml), data)
    
    rnlabels <- terms(formula(newrndfrml))
    random.labels <- attr(rnlabels,"term.labels")
    id.levels <- levels(data[,colnames(data)==(names(object$rnd)[1])])
    n <- length(id.levels)
    s <- dim(W_start)[2]/n
    
    if(s>1)
    {
      W<-W_start[,seq(from=1,to=1+(s-1)*n,by=n)]
      for (i in 2:n)
        W<-cbind(W,W_start[,seq(from=i,to=i+(s-1)*n,by=n)])
    }else{
      W<-W_start
    }
  
    XW <- cbind(X.both,W)
      
    ## check if subjects match
    ranef <- rep(0,ncol(W)); names(ranef) <- colnames(W)
    ranef.index <- is.element(names(object$ranef),names(ranef))
    ranef[names(object$ranef[ranef.index])] <- object$ranef[ranef.index]
    ## time-constant part of eta
    eta.vec <- XW %*% c(object$coefficients,ranef)
    }else{
      n <- 1;id.levels <- data[1,colnames(data)==(names(object$rnd)[1])]
      eta.vec <- rep(object$ranef[is.element(names(object$ranef),paste0(names(object$rnd)[1],id.levels))],nrow(data))
      
      if(length(eta.vec)==0)
        eta.vec <- rep(0,nrow(data))
      
      if(!is.null(object$coefficients))
        eta.vec <- eta.vec + X.both %*% object$coefficients 
    }
  ## which is the cluster variable?  
  which.id <- colnames(data)==(names(object$rnd)[1])  
  }else{
    if(is.null(id.var))
      id.var <- "ID.articifical"
      #stop("Set id.var argument to specify which is the ID variable (defining which observations belong to the same subject). If not existing, define it!")
    subject.names <- id.var
    n <- nrow(data)
    id.levels <- 1:n
    data$ID.articifical <- id.levels
    which.id <- colnames(data)==id.var

    ## create RE objects such that code still works
    W <- NULL  
    Q1 <- NULL  
    dimr <- 0
    s <- 0
    eta.vec <- X.both %*% object$coefficients 
  }
  
  ## calculate survival and hazard function, if surv.ma = TRUE
  
  if(surv.ma)
  {
    ### time varying coef
    U <- model.matrix(object$vary.coef, data)
    m.vary <- ncol(U)-1
    
    T.mat <- matrix(NA,time.length,N)
    event.ma <- cbind(y0,y)
    T.mat[1:length(time.grid),1:N] <- apply(event.ma,1,function(z,time.grid){(time.grid<z[2] & time.grid>=z[1])},time.grid=time.grid)
    
    
    id.haz.ma <- numeric()
    id.haz.list <- list()
    counter <- 0
    numb.subj.vec <- numeric()
    ## this is a tricky part: it could be "length(id.levels[i]) > 1" due to random effects
    ## and then there could be several rows due to changing covariates for each "individual"
    for(i in 1:n)
    {
      index.help <- (data[,which.id]==id.levels[i])
      T.mat.help <- T.mat[,index.help,drop=FALSE]
      ## find out the number of individuals
      numb.subj <- sum(T.mat.help[1,])
      U.help <- U[index.help,,drop=FALSE]
      eta.vec.help <- eta.vec[index.help]
      
      index.vec.help <- which(T.mat.help[1,]==TRUE)
      id.help <- rep(c(1:numb.subj),times=c(diff(index.vec.help),ncol(T.mat.help)-tail(index.vec.help,1)+1))
  
      for(io in 1:numb.subj)
      {  
        id.haz.list[[counter+io]] <- rep(NA,time.length)
        
        U.help2 <- U.help[(id.help==io),,drop=FALSE]
        eta.vec.help2 <- eta.vec.help[(id.help==io)]
        T.mat.help2 <- T.mat.help[,(id.help==io),drop=FALSE]
        for(j in 1:nrow(U.help2))
          id.haz.list[[counter+io]][T.mat.help2[,j]] <- (eta.vec.help2[j]+apply((t(smooth.effects)*U.help2[j,]),2,sum))[T.mat.help2[,j]]
    
        id.haz.ma <- cbind(id.haz.ma, exp(id.haz.list[[counter+io]]))
      } 
      counter <- counter + numb.subj
      numb.subj.vec[i] <- numb.subj
    }
    
    if(all(numb.subj.vec==1))
    {  
    colnames(id.haz.ma) <- paste0(subject.names,id.levels)
    }else{
      label.numb <- sequence(numb.subj.vec)
      colnames(id.haz.ma) <- paste0(rep(paste0(subject.names,id.levels),times=numb.subj.vec),".",label.numb)
    }
    surv.ma <- exp(-(apply(id.haz.ma,2,cumsum)))
    ## finalize the hazard and survival matrices
    surv.ma <- rbind(1,surv.ma[-nrow(surv.ma),,drop=F])
    row.names(id.haz.ma) <- row.names(surv.ma) <- round(time.grid,digits=3)
  }
  
 ret.obj <- list()
 ret.obj$eta <- eta.vec
 ret.obj$time.grid <- time.grid
 if(surv.ma)
 {
   ret.obj$haz <- id.haz.ma
   ret.obj$survival <- surv.ma
 }   
 return(ret.obj)  
}

