\name{coxlasso}
\alias{coxlasso}
\docType{package}
\title{
A LASSO approach for Cox Frailty Models.}
\description{A LASSO approach for Cox Frailty Models based on the Cox full likelihood is provided.
}
\details{
The \code{coxlasso} algorithm is designed to investigate
the effect structure in the Cox frailty model, which is a
widely used model that accounts for heterogeneity in time-to-event data.
Since in survival models one has to account for possible variation of
the effect strength over time, some features can incorporated with time-varying effects. 


The penalty is depending on the LASSO tuning parameter \eqn{\xi}, which has to be determined by a suitable technique, e.g. by K-fold cross validation.




\tabular{ll}{
Package: \tab pencoxfrail\cr
Type: \tab Package\cr
Version: \tab 1.1.2\cr
Date: \tab 2023-08-25\cr
License: \tab GPL-2\cr
LazyLoad: \tab yes\cr
}
for loading a dataset type data(nameofdataset)
}

\usage{
coxlasso(fix=formula, rnd=NULL, vary.coef=NULL, data, xi, 
              adaptive.weights = NULL, control = list())
}     
\arguments{
  \item{fix}{a two-sided linear formula object describing the LASSO-penalized
    fixed (time-constant) effects part of the model, with the response on the left of a
    \code{~} operator and the terms, separated by \code{+} operators, on
    the right. The response must be a survival object as returned by the \code{\link[survival]{Surv}} function.}  
  \item{rnd}{a two-sided linear formula object describing the
    random-effects part of the model, with the grouping factor on the left of a
    \code{~} operator and the random terms, separated by \code{+} operators, on
    the right.Default is NULL, so no random effects are present.}
  \item{vary.coef}{a one-sided linear formula object describing the
    time-varying effects part of the model, with the time-varying terms, separated by \code{+} operators,
    on the right side of a \code{~} operator.Default is NULL, so no time-varying effects are incorporated.}
  \item{data}{the data frame containing the variables named in the three preceding
    \code{formula} arguments.}
  \item{xi}{the LASSO-penalty parameter that controls the strenght of the penalty term.
  The optimal penalty parameter is a tuning parameter of the procedure that has to be determined, 
 e.g. by K-fold cross validation (see \code{\link{cv.coxlasso}} for details or the quick demo for an example).} 
  \item{adaptive.weights}{for the LASSO-penalized fixed effects a vector of adaptive weights can be passed to the procedure. If no adaptive weights are specified, an unpenalized model (i.e. \eqn{\xi=0}) is fitted by the \code{\link{coxFL}} function and the obtained estimates are used as adaptive weights (see value section).} 
  \item{control}{a list of control values for the estimation algorithm to replace the default values returned by the function \code{\link{coxlassoControl}}. Defaults to an empty list.}
}


\value{Generic functions such as \code{print}, \code{predict}, \code{plot} and \code{summary} have methods to show the results of the fit.

The \code{predict} function uses also estimates of random effects for prediction, if possible (i.e. for known subjects of the grouping factor). 
Either the survival stepfunction or the baseline hazard (not cumulative!) can be calculated by specifying one of two possible methods: \code{method=c("hazard","survival")}. By default, for each new subject in \code{new.data} an individual stepfunction is calculated on a pre-specified time grid, also accounting for covariate changes over time. Alternatively, for \code{new.data} a single vector of a specific (time-constant) covariate combination can be specified.

Usage:  \code{
predict(coxlasso.obj,new.data,time.grid,method=c("hazard","survival"))
}     


The \code{plot} function plots all time-varying effects, including the baseline hazard. 

   \item{call}{a list containing an image of the \code{coxlasso} call that produced the object.}  
     \item{baseline}{a vector containing the estimated B-spline coefficients of the baseline hazard.
     If the covariates corresponding to the time-varying effects are centered (and standardized, see \code{\link{coxlassoControl}}), the coefficients are transformed back to the original scale.} 
          \item{time.vary}{a vector containing the estimated B-spline coefficients of all time-varying effects.
     If the covariates corresponding to the time-varying effects are standardized (see \code{\link{coxlassoControl}}) 
     the coefficients are transformed back to the original scale.} 
  \item{coefficients}{a vector containing the estimated fixed effects.}
  \item{ranef}{a vector containing the estimated random effects.}
  \item{Q}{a scalar or matrix containing the estimates of the random effects standard deviation or variance-covariance parameters, respectively.}
  \item{Delta}{a matrix containing the estimates of fixed and random effects (columns) for each iteration (rows) of the main algorithm (i.e. before the final re-estimation step is performed, see details).}
  \item{Q_long}{a list containing the estimates of the random effects variance-covariance parameters for each iteration of the main algorithm.}
  \item{iter}{number of iterations until the main algorithm has converged.}
  \item{adaptive.weights}{if not given as an argument by the user, a two-column matrix of adaptive weights is calculated by the \code{\link{coxFL}} function; the first column contains the weights \eqn{w_{\Delta,k}}{w_k}, the second column the weights \eqn{v_{k}}{v_k} from \eqn{\xi\cdot J(\zeta,\alpha)}{\xi*J(\zeta,\alpha)}.}
  \item{knots}{vector of knots used in the B-spline representation.}
  \item{Phi.big}{large B-spline design matrix corresponding to the baseline hazard and all time-varying effects. For the time-varying effects, the B-spline functions (as a function of time) have already been multiplied with their associated covariates.}
  \item{time.grid}{the time grid used in when approximating the (Riemann) integral involved in the model's full likelihood.}
  \item{m}{number of metric covariates with time-varying effects.}
  \item{m2}{number of categorical covariates with time-varying effects.}
}



\author{
Andreas Groll \email{groll@statistik.tu-dortmund.de} \cr Maike Hohberg \email{mhohber@uni-goettingen.de}
}

\references{
Groll, A., T. Hastie and G. Tutz (2017). 
Selection of Effects in Cox Frailty Models by
Regularization Methods. \emph{Biometrics} 73(3): 846-856.
}


\seealso{
\code{\link{coxlassoControl},\link{cv.coxlasso},\link{coxFL},\link[survival]{Surv},\link[survival]{pbc}}
}
\examples{
\dontrun{
data(lung)

# remove NAs
lung <- lung[!is.na(lung$inst),]

# transform inst into factor variable
lung$inst <- as.factor(lung$inst)

# just for illustration, create factor with only three ph.ecog classes
lung$ph.ecog[is.na(lung$ph.ecog)] <- 2
lung$ph.ecog[lung$ph.ecog==3] <- 2
lung$ph.ecog <- as.factor(lung$ph.ecog)

fix.form <- as.formula("Surv(time, status) ~ 1 + age + ph.ecog + sex")

lasso.obj <- coxlasso(fix=fix.form, data=lung, xi=10,
                control=list(print.iter=TRUE, exact = 1))
coef(lasso.obj)             
                

# now add random institutional effect
lasso.obj2 <- coxlasso(fix=fix.form, rnd = list(inst=~1), 
              data=lung, xi=10,control=list(print.iter=TRUE, exact = 1))
coef(lasso.obj2)             
# print frailty Std.Dev.
print(lasso.obj2$Q)
# print frailties
print(lasso.obj2$ranef)


# now fit a time-varying effect for age
fix.form <- as.formula("Surv(time, status) ~ 1 + ph.ecog + sex")
vary.coef <- as.formula("~ age")

lasso.obj3 <- coxlasso(fix=fix.form,vary.coef=vary.coef,  
              data=lung, xi=10,control=list(print.iter=TRUE))
summary(lasso.obj3)

# show fit
plot(lasso.obj3)

# predict survival curve of new subject, institution 1 and up to time 300
pred.obj <- predict(lasso.obj2, newdata=data.frame(inst=1, time=NA, status=NA, age=26,
              ph.ecog=2,sex=1), time.grid=seq(0,300,by=1))

# plot predicted hazard function
plot(pred.obj$time.grid,pred.obj$haz,type="l",xlab="time",ylab="hazard")

# plot predicted survival function
plot(pred.obj$time.grid,pred.obj$survival,type="l",xlab="time",ylab="survival")

## specify a larger new data set
new.data <- data.frame(inst=c(1,1,6), time=c(20,40,200), 
      status=c(NA,NA,NA), age=c(26,26,54), ph.ecog=c(0,0,2),sex=c(1,1,1))

## as here no frailties have been specified, id.var needs to be given!
pred.obj2 <- predict(lasso.obj3, newdata=new.data,id.var = "inst")

# plot predicted hazard functions (for the available time intervals)
plot(pred.obj2$time.grid[!is.na(pred.obj2$haz[,1])],
      pred.obj2$haz[,1][!is.na(pred.obj2$haz[,1])],
      type="l",xlab="time",ylab="hazard",xlim=c(0,200),
      ylim=c(0,max(pred.obj2$haz,na.rm=T)))
lines(pred.obj2$time.grid[!is.na(pred.obj2$haz[,3])],
      pred.obj2$haz[,3][!is.na(pred.obj2$haz[,3])],
      col="red",lty=2,)

# plot predicted survival functions (for the available time intervals)
plot(pred.obj2$time.grid[!is.na(pred.obj2$survival[,1])],
    pred.obj2$survival[,1][!is.na(pred.obj2$survival[,1])],
    type="l",xlab="time",ylab="hazard",xlim=c(0,200),
    ylim=c(0,max(pred.obj2$survival,na.rm=T)))
lines(pred.obj2$time.grid[!is.na(pred.obj2$survival[,3])],
      pred.obj2$survival[,3][!is.na(pred.obj2$survival[,3])],
      col="red",lty=2,)


# see also demo("coxlasso-lung")
}}

\concept{Lasso}
\concept{Shrinkage}
\concept{Variable selection}
\concept{Cox Frailty Model}


