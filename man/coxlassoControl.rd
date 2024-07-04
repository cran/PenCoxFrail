\name{coxlassoControl}
\alias{coxlassoControl}
\concept{coxlassoControl}
\title{Control Values for \code{coxlasso} fit}
\description{
  The values supplied in the function call replace the defaults and a list with all possible arguments is returned. The returned list is used as the \code{control} argument to the \code{coxlasso} function.
}

\usage{

coxlassoControl(start = NULL, index=NULL, q_start = NULL, conv.eps = 1e-3, 
                          standardize = TRUE, center = FALSE,
                          smooth=list(nbasis = 6, penal = 1e+2), 
                          print.iter = FALSE, max.iter = 100, c.app = 1e-6,  
                          exact = NULL, xr = NULL, eps = 1e-2, quant.knots = TRUE,...)
} 
    
\arguments{
  \item{start}{a vector of suitable length containing starting values for the spline-coefficients of the baseline hazard and the time-varying effects, followed by the fixed and random effects. The correct ordering is important. Default is a vector full of zeros.}
   \item{index}{vector which defines the grouping of the variables. Components sharing the same number build a group and factor variables get a single number (and are automatically treated as a group). Non-penalized coefficients are marked with NA.}
  \item{q_start}{a scalar or matrix of suitable dimension, specifying starting values for the random-effects variance-covariance matrix. Default is a scalar 0.1 or diagonal matrix with 0.1 in the diagonal, depending on the dimension of the random effects.}
  \item{conv.eps}{controls the speed of convergence. Default is 1e-3.}
  \item{center}{logical. If true, the covariates corresponding to the time-varying effects will be
    centered. Default is FALSE (and centering is only recommended if really necessary; it can also have a strong effect on the baseline hazard, in particular, if a strong penalty is selected).}
  \item{standardize}{logical. If true, the covariates corresponding to the fixed effects will be
    scaled to a variance equal to one. Default is TRUE.}
      \item{smooth}{a list specifying the number of basis functions \code{nbasis} (used for the baseline hazard and all time-varying effects) and the smoothness penalty parameter \code{penal}, which is only applied to the baseline hazard. All time-varying effects are penalized by the specific double-penalty \eqn{\xi\cdot J(\zeta,\alpha)}{\xi*J(\zeta,\alpha)} (see \code{\link{coxlasso}}), which is based on the overall penalty parameter \eqn{\xi} (specified in the main function \code{\link{coxlasso}}) and on the weighting between the two penalty parts \eqn{\zeta}. The degree of the B-splines is fixed to be three (i.e. cubic splines).}
    \item{print.iter}{logical. Should the number of iterations be printed? Default is FALSE.}
  \item{max.iter}{the number of iterations for the final Fisher scoring re-estimation procedure. Default is 200.}
    \item{c.app}{The parameter controlling the exactness of the quadratic approximations of the penalties. Default is 1e-6.}
\item{exact}{controls the exactness of the (Riemann) integral approximations. If not set by the user to a specific value, it will
be automatically chosen such that the time interval [0,t_max] will be divided in about 1000 equal sized Riemann bars.}
\item{xr}{maximal time point that is regarded. Default is NULL and the maximal event or censoring time point in the data is used.}
\item{eps}{Small epsilon that controls which fixed effects are set to zero: parameters with an absolute value smaller than epsilon are taken to be zero. Default is 1e-2.}
\item{quant.knots}{Shall the knots be defined based on quantiles? Default is TRUE.}.
\item{...}{Futher arguments to be passed.}
}

\value{
  a list with components for each of the possible arguments.
}

\author{
Andreas Groll \email{groll@math.lmu.de}
}

\seealso{
  \code{\link{coxlasso}}
}

\examples{
# Use different weighting of the two penalty parts
# and lighten the convergence criterion 
coxlassoControl(c.app = 1e-5, conv.eps=1e-3)
}
