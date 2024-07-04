###############
library(PenCoxFrail)

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

set.seed(1909)
cv.obj <- cv.coxlasso(fix = fix.form, n.folds = 5, 
                        data = lung, control = list(print.iter = FALSE))
  
## show result of CV
plot(cv.obj)

## now fit xi grid on whole data set

Delta.start <- NULL

for(j in 1:length(cv.obj$xi))
{
  print(paste("xi: ", j , sep=""))
  lasso.obj <- try(coxlasso(fix=fix.form, data=lung, xi=cv.obj$xi[j],
                            control = list(start=Delta.start[j-1,])),silent=TRUE)
  if(class(lasso.obj)!="try-error")
  {  
    Delta.start <- rbind(Delta.start,lasso.obj$Delta[lasso.obj$iter,])
  }else{
    Delta.start <- rbind(Delta.start,Delta.start[nrow(Delta.start),])
  }  
}

## fit final lasso model with xi.opt
lasso.obj.1 <- coxlasso(fix=fix.form, data=lung, 
                        xi=cv.obj$xi.opt, control = list(Delta.start[cv.obj$ind.opt-1,]))

summary(lasso.obj.1)


## fit final lasso model with xi.1se; here special case: no starting value, as cv.obj$ind.1se = 1;
## hence, Delta.start[cv.obj$ind.1se-1,] makes no sense here
lasso.obj.2 <- coxlasso(fix=fix.form, data=lung, xi=cv.obj$xi.1se)

summary(lasso.obj.2)


## show coefficient paths
xi.opt.round <- round(cv.obj$xi.opt,digits=4)
xi.1se.round <- round(cv.obj$xi.1se,digits=4)

quartz()
par(cex=1.7,mar=c(5,5,3.5,3))
plot(cv.obj$xi,Delta.start[,7],ylim=c(-1.5,1.5),type="l",
     lwd=2,ylab=expression(hat(beta[j])),xlab=expression(xi))
abline(h=0,lty=3,lwd=2,col="grey")
for(j in 8:ncol(Delta.start)){
  lines(cv.obj$xi,Delta.start[,j],lwd=2)
}
axis(3,at=cv.obj$xi.opt,labels=substitute(paste(xi["opt"]," = ",nn), list(nn=xi.opt.round)))
abline(v=cv.obj$xi.opt,lty=2,lwd=2,col=2)
axis(3,at=cv.obj$xi.1se,labels=substitute(paste(xi["1se"]," = ",nn), list(nn=xi.1se.round)))
abline(v=cv.obj$xi.1se,lty=2,lwd=2,col=2)



