# cross-validation function

cv.coxlasso <- function(fix, rnd = NULL, vary.coef = NULL, n.folds = 10, xi = NULL, 
                        data, adaptive.weights = NULL, print.fold = TRUE, print.xi = FALSE, 
                        len.xi = 100, lgrid = TRUE, ran.seed = 1909, xi.factor = 1.01, min.fold = 4,
                        pass.on.start = TRUE, control = list(print.iter = FALSE))
{
  # set up
  set.seed(ran.seed)
  N<-dim(data)[1]
  ind<-sample(N,N)
  
  # n.folds
  kk <- n.folds
  nk <- floor(N/kk)

  ## make sure that usually covariates are centered & standardized
  if(is.null(control$standardize)) 
    control$standardize <- TRUE

  if(is.null(control$center)) 
    control$center <- FALSE
  
  ## prepare to later extend control list
  control.list <- control
  
  
  ## determine a rough value for exact for coxFL to make it fast
  design.temp <- model.frame(fix, data)
  Resp <- model.response(design.temp)
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
  
  y.max <- max(y)  
 
  ####### try coxFL, otherwise do ridge until estimation possible#####################

  coxridge.obj <- try(coxFL(fix=fix,rnd=rnd,vary.coef=vary.coef,data=data,control=control.list), silent = TRUE)
  if(inherits(coxridge.obj, "try-error")) {
    xi.ridge <- 10^(-6)
     while(inherits(coxridge.obj, "try-error")) {
       xi.ridge <- xi.ridge * 10
       coxridge.obj <- try(coxridge(fix=fix, rnd = rnd, vary.coef = vary.coef, 
                                   xi.ridge = xi.ridge, data=data, control = control.list), silent = TRUE)  
      }
    } 
    
  Delta.start <- as.matrix(t(rep(0, ncol(coxridge.obj$Delta))))
  adaptive.weights <- abs(coxridge.obj$coef.stand)
  names(adaptive.weights) <- NULL
  
  if(!is.null(rnd))
    Q.start <- 0.1

  if(is.null(xi))
  {
    control.list.temp <- control.list
    control.list.temp$start <- NULL
    
    if(!is.null(rnd))
      control.list.temp$q_start <- NULL
    
    try.grid <- TRUE
   
    ##### take xi.max from ridge regression as starting candidate
    xi.max <-  coxridge.obj$xi.max*0.5
    
    k <- 0 
    
    while(try.grid) # until all coef become 0 
    {
      #print(paste("k-try: ",k,sep=""))
      xi.max <- xi.max*2
      coxlasso.obj <- try(coxlasso(fix=fix, rnd = rnd, vary.coef = vary.coef, 
                                   adaptive.weights = adaptive.weights, 
                                   data = data, xi = xi.max, control = control.list.temp), silent=TRUE)
      if (!inherits(coxlasso.obj, "try-error"))
        try.grid <- !all(coef(coxlasso.obj)==0)
      
      k <- k + 1
    } 
    ## increase xi.max once again to be sure on all sets
    xi.max <- xi.max * xi.factor
    
    if (lgrid == TRUE) 
    {
       xi.grid <- exp(seq(log(xi.max),log(1e-4),length.out = len.xi))
    }else{ 
      xi.grid <- seq(xi.max,0,length.out = len.xi)
    }
   }  
  
  ## set up matrix for CV error
  loglik.mat <-  loglik.mat.sd <- matrix(NA,ncol=kk,nrow=length(xi.grid))

  icept.vec <- rep(NA, kk)
  #### loop over the folds  
  for (i in 1:kk)
  {
    if(print.fold)
      print(paste("CV fold ", i,sep=""))
    
    if (i < kk)
    {
      indi <- ind[(i-1)*nk+(1:nk)]
    }else{
      indi <- ind[((i-1)*nk+1):N]
    }
    
    data.train <- data[-indi,]
    data.test <- data[indi,]
    
    y0.train <- y0[-indi]
    y0.test <- y0[indi]
    
    y.train <- y[-indi]
    y.test <- y[indi]
    y.max.test <- max(y.test)
    
    delta.train <- delta[-indi]
    delta.test <- delta[indi]

    icept <- log(sum(delta.train)/sum(y.train-y0.train))

    icept.vec[i] <- icept
    ### in first row of Delta.temp all 0 except intercept of baseline hazard (as done in coxFL)
    ### just specify icept 
    Delta.temp <- Delta.start
    Delta.temp[1,1] <- icept 

    if(!is.null(rnd))
      Q.temp <- Q.start
    
    coxridge.obj <- try(coxFL(fix=fix,rnd=rnd,vary.coef=vary.coef,data=data.train,control=control), silent = TRUE)
    if (inherits(coxridge.obj, "try-error")) {
      xi.ridge.temp <- 10^(-6)
      while (inherits(coxridge.obj, "try-error")) {
        xi.ridge.temp <- xi.ridge.temp * 10
        coxridge.obj <- try(coxridge(fix=fix, rnd = rnd, vary.coef = vary.coef, 
                                     xi.ridge = xi.ridge.temp, data=data.train, control = control), silent = TRUE)  
      }
    } 

    adaptive.weights.temp <- abs(coxridge.obj$coef.stand)
    names(adaptive.weights.temp) <- NULL
    
    ## loop over xi grid
    for(j in 1:length(xi.grid))
    {
      if(print.xi)
        print(paste("xi: ", j , sep=""))
      
      ## extend control list
      if(pass.on.start)
      { 
        control.list$start <- Delta.temp[j,]
        }else{
        control.list$start <- Delta.temp[1,]
     }  

      if(!is.null(rnd))
      {
        if(pass.on.start)
        {  
          control.list$q_start <- Q.temp[j]
        }else{
          control.list$q_start <- Q.temp[1]
        }  
        
      }
      
      coxlasso.obj <- try(coxlasso(fix=fix, rnd = rnd, vary.coef = vary.coef, 
                          data = data.train, xi = xi.grid[j], 
                          adaptive.weights = adaptive.weights.temp,  
                          control = control.list), silent=TRUE)
      
      if(!inherits(coxlasso.obj,"try-error"))
      {  
          coef.temp <- coxlasso.obj$Delta[coxlasso.obj$iter,]
          Delta.temp<-rbind(Delta.temp,coef.temp)

        if(!is.null(rnd))
          Q.temp<-c(Q.temp,coxlasso.obj$Q_long[[coxlasso.obj$iter]])
        
        status.test <- delta.test
        
        ## note that predict.coxlasso function also already includes the RE prediction, if REs are present; 
        lin.pred <- predict(coxlasso.obj,newdata=data.test)$eta
        
        ## matrix if vary.coef is present; otherwise only intercept vec
        U.test <- model.matrix(coxlasso.obj$vary.coef, data.test)
        
        ## the function bs.design.asy comes from helpers.R
        Phi.test <- bs.design.asy(y.test, xl = 0, xr = y.max, 
                                  spline.degree=coxlasso.obj$spline.degree, nbasis = coxlasso.obj$nbasis,
                                  quant = coxlasso.obj$quant.vec)
        Phi.test  <-cbind(Phi.test%*%coxlasso.obj$B.unpen.fact,Phi.test%*%coxlasso.obj$Pen)
        
        Phi.test.big<-matrix(0,nrow(Phi.test),(coxlasso.obj$m.vary+1)*coxlasso.obj$nbasis)
        for(jj in 1:(coxlasso.obj$m.vary+1))
          Phi.test.big[,((jj-1)*coxlasso.obj$nbasis+1):(jj*coxlasso.obj$nbasis)] <- Phi.test * U.test[,jj]
        
        time.grid.test <- coxlasso.obj$time.grid
        
        Phi.grid <- bs.design.asy(time.grid.test, xl = 0, xr = y.max, 
                              spline.degree=coxlasso.obj$spline.degree, nbasis = coxlasso.obj$nbasis,
                              quant = coxlasso.obj$quant.vec)
        Phi.grid   <-cbind(Phi.grid %*%coxlasso.obj$B.unpen.fact,Phi.grid %*%coxlasso.obj$Pen)
        
        Phi.grid.big <- t(apply(Phi.grid,1,rep,coxlasso.obj$m.vary+1))
        
        ## the function int.approx comes from helpers.R
        part.1 <- status.test * (Phi.test.big %*% c(coxlasso.obj$baseline.cen,coxlasso.obj$time.vary) + lin.pred)
        part.2 <- (apply(cbind(y.test,U.test), 1,int.approx,time.grid=time.grid.test,B=Phi.grid.big,nbasis=coxlasso.obj$nbasis,
          alpha=c(coxlasso.obj$baseline.cen,coxlasso.obj$time.vary)) * exp(lin.pred))
        part.3 <- part.1-part.2
        
        loglik.mat.sd[j,i] <- sd(part.3) 
        loglik.mat[j,i] <- mean(part.3)

      }else{
        Delta.temp<-rbind(Delta.temp,Delta.temp[nrow(Delta.temp),])
        
        if(!is.null(rnd))
          Q.temp<-c(Q.temp,Q.temp[length(Q.temp)])
        
      }
    }
  }
 cv.error.mat <- -loglik.mat  
 
 ## only those xi values are taken into account where min.fold folds are not NA
 xi.notok <- apply(cv.error.mat, 1, isnt.na.fct)<min.fold
 
 cv.error.vec <- rowMeans(cv.error.mat, na.rm = TRUE) # for each xi on grid averaged by kk

 loglik.mat.sd <- loglik.mat.sd^2
 
 ## now pool over folds and then devide by overall sample size as we want var(X_bar)
 cv.error.sd <- sqrt(apply(loglik.mat.sd, 1, sd, na.rm = TRUE)/n.folds/N)
 
 cv.error.vec[xi.notok] <- Inf
 cv.error.sd[xi.notok] <- NA
 
 ind.opt <- which.min(cv.error.vec)
 xi.opt <- xi.grid[ind.opt]
   
 ## calculate stronger penalty parameter similar to glmnet
 ind.1se <- which.max(cv.error.vec <= cv.error.vec[ind.opt] + cv.error.sd[ind.opt])
 
 xi.1se <- xi.grid[ind.1se]
 
 ret.obj <- list()
 ret.obj$xi.opt <- xi.opt
 ret.obj$xi.1se <- xi.1se
 ret.obj$ind.opt <- ind.opt
 ret.obj$ind.1se <- ind.1se
 ret.obj$xi <- xi.grid
 ret.obj$cv.error <- cv.error.vec
 ret.obj$cv.error.mat <- cv.error.mat
 ret.obj$cv.error.sd <- cv.error.sd
 ret.obj$adaptive.weights <- adaptive.weights
 ret.obj$icept.vec <- icept.vec
 ret.obj$coxridge.obj <- coxridge.obj
 class(ret.obj) <- "cv.coxlasso"
 return(ret.obj) # for now give optimal xi, later call coxlasso with optimal xi

}


#################

plot.cv.coxlasso <- function(x, include.xi.1se = TRUE, ylim = NULL, ...){
  up <- x$cv.error + x$cv.error.sd
  lo <- x$cv.error - x$cv.error.sd

  if(is.null(ylim))
    ylim <- c(min(lo,na.rm = TRUE)-2, max(up,na.rm = TRUE)+2)

  plot(x$xi, x$cv.error, xlab=expression(xi), 
       ylab="CV error (negative log-lik)", ylim = ylim, ...)
  error.bars(x$xi, up, lo, width=0.01, col="darkgrey", ylim = ylim,lty=1)
  abline(v=x$xi.opt,lty=3)
  xi.opt.round <- round(x$xi.opt,digits=4)
  axis(3,at=x$xi.opt,labels=substitute(paste(xi["opt"]," = ",nn), list(nn=xi.opt.round)))
  if(include.xi.1se){
  abline(v=x$xi.1se,lty=3)
  xi.1se.round <- round(x$xi.1se,digits=4)
  axis(3,at=x$xi.1se,labels=substitute(paste(xi["1se"]," = ",nn), list(nn=xi.1se.round)))
  }
}

  