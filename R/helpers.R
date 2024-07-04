## help functions

bs.design.asy<-function(x, xl, xr, spline.degree, nbasis, comp = NULL, quant = NULL)
{
  
  ## generate a B-Spline-Matrix with equidistant knots (code by Thomas Kneib & Andreas Groll):
  ## x are the positions where spline to be evaluated
  ## xl, xr intervall boundaries where spline functions are relevant
  xmin<-xl-(xr-xl)/100
  xmax<-xr+(xr-xl)/100
  dx<-(xmax-xmin)/(nbasis-3)
  knots<-seq(xmin-spline.degree*dx,xmax+spline.degree*dx,by=dx)
  if(!is.null(quant))
    knots[3+(1:(nbasis-2))] <- quant
  B<-splines::spline.des(knots,x,spline.degree+1,outer.ok = TRUE)$design
  if(is.null(comp))
  {
    return(B)
  }else{
    return(B[,comp])
  }
}

########

bs.design<-function(x, xl, xr, spline.degree, nbasis, comp = NULL)
{
  
  ## generate a B-Spline-Matrix with equidistant knots (code by Thomas Kneib & Andreas Groll):
  ## x are the positions where spline to be evaluated
  ## xl, xr intervall boundaries where spline functions are relevant
  xmin<-xl-(xr-xl)/100
  xmax<-xr+(xr-xl)/100
  dx<-(xmax-xmin)/(nbasis-3)
  knots<-seq(xmin-spline.degree*dx,xmax+spline.degree*dx,by=dx)
  B<-splines::spline.des(knots,x,spline.degree+1,outer.ok = TRUE)$design
  if(is.null(comp))
  {
    return(B)
  }else{
    return(B[,comp])
  }
}

########


mirror <- function(x,low.tri)
{
  x[low.tri] <- t(x)[low.tri]
  return(x)
}

#############

penal.fct <- function(x,K)  sqrt(t(x) %*% (K %*% x))

#############

is.na.fct <- function(x)  sum(is.na(x))

isnt.na.fct <- function(x)  sum(!is.na(x))


############# function for scaling a design matrixs
scaleDesign <- function(X, index.vec, center, standardize)
{
    ## Which are the non-penalized parameters?
    any.notpen    <- any(is.na(index.vec))
    inotpen.which <- which(is.na(index.vec))
    nrnotpen      <- length(inotpen.which)
    
    ## Index vector of the penalized parameter groups
    if(any.notpen){
      ipen <- index.vec[-inotpen.which]
      ipen.which <- split((1:ncol(X))[-inotpen.which], ipen)
    }else{
      ipen <- index.vec
      ipen.which <- split((1:ncol(X)), ipen)
    }
    
    
    if(center){
      mu.x <- apply(as.matrix(X), 2, mean)
      X <- sweep(as.matrix(X), 2, mu.x)
    }else{
      mu.x <- NULL
    }
    
    ## Standardize the design matrix -> blockwise orthonormalization
    if(standardize){
      ##warning("...Using standardized design matrix.\n")
      stand        <- blockstand(X, ipen.which, inotpen.which)
      X            <- stand$x
      scale.pen    <- stand$scale.pen
      scale.notpen <- stand$scale.notpen
    }else{
      scale.pen    <- NULL
      scale.notpen <- NULL
    }
    
    ret.obj <- list()
    ret.obj$X <- X
    ret.obj$mu.x <- mu.x
    ret.obj$inotpen.which <- inotpen.which
    ret.obj$ipen.which <- ipen.which
    ret.obj$scale.pen <- scale.pen
    ret.obj$scale.notpen <- scale.notpen
    ret.obj$any.notpen <- any.notpen
    return(ret.obj)
}

############# Author: Lukas Meier, Date:  4 Aug 2006, 08:50

blockstand <- function(x, ipen.which, inotpen.which)
{
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
    if(length( scale.pen[[j]]) ==1)
    {
      if(scale.pen[[j]]<0)
      {  
        scale.pen[[j]] <- scale.pen[[j]] * (-1)
        x.ort[,ind] <- x.ort[,ind] *(-1)
      }
    }
  }
  list(x = x.ort, scale.pen = scale.pen, scale.notpen = scale.notpen)
}

#############

factor.test <- function(x)  substr(x,1,9)=="as.factor"

#############

remove.fac <- function(x){
  x1 <- gsub("as.factor","",x)
  substr(x1, 2, nchar(x1)-1)
  
} 
#############

penal.fct.inv <- function(x,K,c.app)  (as.numeric(t(x) %*% (K %*% x))+c.app)^(-0.5)

#############

eucl.norm <- function(x)  sqrt(sum(x^2))

#############

nelem <- function(x) nlevels(factor(x))

#############

count.levels <- function(x)  nlevels(x)

#############

center.fct <- function(fit.vec,mean.vec,sd.vec,nbasis,m, standardize = TRUE, center = TRUE)
{
  if(!center)
    mean.vec <- rep(0,length(mean.vec))
  if(!standardize)
    sd.vec <- rep(1,length(sd.vec))
  
    Delta.base.retrans <- c(1,-mean.vec/sd.vec) %*% t(matrix(fit.vec[1:(nbasis*(m+1))],nbasis,m+1))
  return(Delta.base.retrans)
}

############

int.approx <- function(z,time.grid,B,nbasis,alpha)
{
  index <- time.grid<z[1]
  diff(time.grid[1:2]) * sum(exp(B[index,] %*% (alpha * rep(z[2:length(z)],each=nbasis))))
}

############ error bars for cv.coxlasso

error.bars <-function(z, upper, lower, width = 0.02, ylim, ...){
  xlim <- range(z)
  barw <- diff(xlim) * width
  segments(z, upper, z, lower, ylim, ...)
  segments(z - barw, upper, z + barw, upper, ylim, ...)
  segments(z - barw, lower, z + barw, lower, ylim, ...)
  range(upper, lower, na.rm = T)
}

############ Standardize variables wiht n instead of (n-1) as denominator

mysd <- function(y) sqrt(sum((y-mean(y))^2)/length(y))

