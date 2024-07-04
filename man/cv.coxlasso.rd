\name{cv.coxlasso}
\alias{cv.coxlasso}
\concept{cv.coxlasso}
\title{Cross-validation for coxlasso}

\description{performs k-fold cross-validation for \code{coxlasso}, produces a plot, and returns a value for the LASSO tuning parameter \eqn{\xi}.  
}
\details{
The function runs \code{coxlasso} over a grid of values \eqn{\xi} for each training data set with one fold omitted. 

For each run, the value for the full likelihood is calculated and the average for each \eqn{\xi} on the grid is computed over the folds. The function choses the \eqn{\xi} that maximizes this likelihood value as the optimal tuning parameter value.  

}



\usage{
cv.coxlasso(fix, rnd = NULL, vary.coef = NULL, n.folds = 10, xi = NULL,
            data, adaptive.weights = NULL, print.fold = TRUE, print.xi = FALSE,
            len.xi = 100, lgrid = TRUE, ran.seed = 1909, xi.factor = 1.01, min.fold = 4,
            pass.on.start = TRUE, control = list(print.iter = FALSE))
}    

\arguments{
  \item{fix}{a two-sided linear formula object describing the fixed (time-constant) effects part of the model, with the response on the left of a
    \code{~} operator and the terms, separated by \code{+} operators, on
    the right. The response must be a survival object as returned by the \code{\link[survival]{Surv}} function.}  
  \item{rnd}{a two-sided linear formula object describing the
    random-effects part of the model, with the grouping factor on the left of a
    \code{~} operator and the random terms, separated by \code{+} operators, on
    the right. Default is NULL, so no random effects are present.}
  \item{vary.coef}{a one-sided linear formula object describing the
    time-varying effects part of the model, with the time-varying terms, separated by \code{+} operators,
    on the right side of a \code{~} operator. Default is NULL, so no time-varying effects are incorporated.}
  \item{n.folds}{number of folds. Default is 10.} 
  \item{xi}{Optional user-supplied xi sequence; default is NULL, and cv.coxlasso chooses its own sequence}
  \item{data}{the data frame containing the variables named in the three preceding
    \code{formula} arguments.}
 \item{adaptive.weights}{for the LASSO-penalized fixed effects a vector of adaptive weights can be passed to the procedure. If no adaptive weights are specified, an unpenalized model (i.e. \eqn{\xi=0}) is fitted by the \code{\link{coxFL}} function and the obtained estimates are used as adaptive weights (see value section).} 
 \item{print.fold}{Should folds of CV be printed? Default is yes.}
 \item{print.xi}{Should current \eqn{\xi} value be printed? Default is no.}
 \item{len.xi}{Length of \eqn{\xi} grid. Default is 100.}
 \item{lgrid}{Logical; shall a logarithmized grid version for the penalty parameter be used? Default is TRUE.}
 \item{ran.seed}{Random seed number to be set. Default is 1909, the year of birth of Borussia Dortmund football club.}
 \item{xi.factor}{A factor which increases xi.max once again to be sure that xi is large enough on all sets. Default is 1.01}
 \item{min.fold}{Only those xi values are taken into account where at least min.fold folds are not NA. Default is 4.}
 \item{pass.on.start}{Shall starting values be passed onthroughout estimation? Default is TRUE}
  \item{control}{a list of control values for the estimation algorithm to replace the default values returned by the function \code{\link{coxlassoControl}}. Default is \code{print.iter = FALSE}.}
}


\value{The function returns a list \code{"cv.coxlasso"} which includes:    

  \item{cv.error}{a vector of mean CV error (i.e., negative likelihood) values for each \eqn{\xi} on the grid averaged over the folds.}
  \item{xi.opt}{a scalar value of \eqn{\xi} associated with the smallest CV error.}
  \item{xi.1se}{largest value of \eqn{\xi} such that error is within 1 standard error of the minimum.}


The \code{plot} function plots the values of \eqn{\xi} against the corresponding CV error (i.e., negative likelihood) values.   
   
}



\author{
Andreas Groll  \email{groll@statistik.tu-dortmund.de} \cr
Maike Hohberg \email{mhohber@uni-goettingen.de}
}

\references{To appear soon.
}


\seealso{
\code{\link{coxlasso}, \link{coxlassoControl}, \link{coxFL}, \link[survival]{Surv}, \link[survival]{pbc}}
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

# find optimal tuning paramater
cv.coxlasso.obj <- cv.coxlasso(fix = fix.form, data = lung, n.folds = 5)

# estimate coxlasso model with optimal xi
lasso.obj <- coxlasso(fix=fix.form, data=lung, xi=cv.coxlasso.obj$xi.opt,
                control=list(print.iter=TRUE))

coef(lasso.obj)             
                

# see also demo("coxlasso-lung")
}}

\concept{Lasso}
\concept{Shrinkage}
\concept{Variable selection}
\concept{Cox Frailty Model}
\concept{Cross validation}



