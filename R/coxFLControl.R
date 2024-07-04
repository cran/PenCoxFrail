### control
coxFLControl<-function(start = NULL, index=NULL, q_start = NULL, conv.eps = 1e-3, 
                          standardize = FALSE, center = FALSE, 
                          smooth=list(nbasis = 6, penal = 1e+2), 
                          print.iter = FALSE, max.iter = 100,   
                          exact = NULL, xr = NULL, eps = 1e-3, ...)
{                       
  list(start = start, index=index, q_start = q_start, conv.eps = conv.eps, standardize = standardize, center = center,
       smooth = smooth, print.iter = print.iter, 
       max.iter = max.iter, exact = exact, xr = xr, eps = eps, ...)
}
