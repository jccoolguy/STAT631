####---------------------------####
####    CDF Transformation     ####
####---------------------------####
set.seed(9132024)
n = 2048; 
par(mfrow =c(3,2), pty = "s", mgp = c(2,0.8,0), cex.lab = 1.25)
x = rnorm(n);  ## simulate N(0,1)
y = rexp(n);   ## simulate Exp(1)
hist(x, "scott", main ="x~N(0,1)", xlab = "");
hist(y, "scott", main = "y~Exp(1)", xlab = "")
head(cbind(x,y),4)

ux = pnorm(x); uy = pexp(y) ## CDF transform
hist(pnorm(x),"scott", main = "ux: Normal CDF of x");
hist(pexp(y),"scott", main = "uy: Exp CDF of y")

head(cbind(qnorm(ux), qexp(uy)),4);  ## quantile transform back to x and y

hist(qchisq(ux,6), "scott", main = "chisq quantile of ux"); ## chisquare
hist(qt(uy,3),"scott", main = "t quantile of uy");  ## t_3


set.seed(8343143)
u = runif(n) ### simulate uniform random variables
par(mfrow =c(3,2), pty = "s", mgp = c(2,0.8,0), cex.lab = 1.25)
hist(u, "scott", main = "u~Unif(0,1)", xlab = "")

## quantile transform to 5 different distributions
hist(qnorm(u,1,2), "scott", main = "normal quantile of u")
hist(qt(u,3),"scott", main = "t quantile of u")
hist(qchisq(u, 5),"scott", main = "chisquare quantile of u")
hist(qlnorm(u),"scott", main = "lognormal quantile of u")
hist(qweibull(u, 1),"scott", main = "Weibull quantile of u")

####---------------------------####
####     Gaussian Copula       ####
####---------------------------####
S = matrix(c(1,0.7,0.7,1),2); S
library(MASS); 
set.seed(63290)
X = mvrnorm(n,c(0,0), S); x1 = X[,1]; x2 = X[,2]
cat("Bivariate Normal with rho  = 0.7\n"); S

par(mfcol = c(3,2), pty = "s", mgp = c(2,0.8,0), cex.lab = 1.25)
hist(x1, main = "x1~N(0,1)", xlab = "")
hist(x2, main = "x2~N(0,1)", xlab = "")
plot(x1,x2, main = "Bivariate Normal")
mtext(paste("rho.hat =", round(cor(x1,x2),4)), font = 2, cex = 0.8)

u1 = pnorm(x1); u2 = pnorm(x2); 
hist(pnorm(x1),"scott", main = "u1: Normal CDF of x1")
hist(pnorm(x2),"scott", main = "u2: Normal CDF Of x2")
plot(u1,u2)
title(main = "Gaussian Copula")
mtext(paste("rho.hat =", round(cor(u1,u2),4)), font = 2, cex = 0.8)


v1 = qt(u1, df = 2); v2 = qt(u2, df = 4);
hist(qt(u1, df = 2),"scott", main = "v1~t_2")
hist(qt(u2, df = 4),"scott", main = "v2~t_4")
plot(v1,v2)
mtext(paste("rho.hat =", round(cor(v1,v2),4)), font = 2, cex = 0.8)

library(fGarch)
w1 = qsstd(u1, nu = 3, xi = 0.75); w2 = qsstd(u2, nu = 6, xi = 1.5);
hist(qsstd(u1, nu = 3, xi = 0.75),"scott", main = "w1: skewed t, df = 3")
hist(qsstd(u2, nu = 6, xi = 1.5),"scott", main = "w2: skewed t, df = 6")
plot(w1,w2)
mtext(paste("rho.hat =", round(cor(w1,w2),4)), font = 2, cex = 0.8)

####---------------------------####
#####         t-Copula         ####
####---------------------------####
####    Simulate Bivariate t   ####
set.seed(20283)
df = 3
chi = rchisq(n,df)
y1 = x1/sqrt(chi/df); y2 = x2/sqrt(chi/df); Y = cbind(y1,y2)

####---------------------------####
par(mfcol = c(3,2), pty = "s", mgp = c(2,0.8,0), cex.lab = 1.25)
hist(y1,"scott",  main = "y1~t_3", xlab = "")
hist(y2,"scott", main = "y2~t_3", xlab = "")
plot(Y, main = "Bivariate t, DF = 3")
mtext(paste("rho.hat =", round(cor(y1,y2),4)), font = 2, cex = 0.8)

ut1 = pt(y1,df= 3); ut2 = pt(y2, df =3)
hist(pt(y1, df = 3),"scott", main = "t_3 CDF of y1")
hist(pt(y2, df =3), "scott", main = "t_3 CDF Of y2")
plot(ut1,ut2, main = "t-Copula")
mtext(paste("rho.hat =", round(cor(ut1,ut2),4)), font = 2, cex = 0.8)

####---------------------------####
####    Plot both Copulas      ####
####---------------------------####
par(mfrow = c(2,2), pty = "s", mgp = c(2,0.8,0), cex.lab = 1.25)
plot(x1,x2, main = "Bivariate Normal")
mtext(paste("rho.hat =", round(cor(x1,x2),4)), font = 2, cex = 0.8)
plot(Y, main = "Bivariate t, DF = 3")
mtext(paste("rho.hat =", round(cor(y1,y2),4)), font = 2, cex = 0.8)
plot(u1,u2)
title(main = "Gaussian Copula")
mtext(paste("rho.hat =", round(cor(u1,u2),4)), font = 2, cex = 0.8)
plot(ut1,ut2, main = "t-Copula")
mtext(paste("rho.hat =", round(cor(ut1,ut2),4)), font = 2, cex = 0.8)

####---------------------------#####
## t marginals with different df's
par(mfcol = c(3,2), pty = "s", mgp = c(2,0.8,0), cex.lab = 1.25)
vt1 = qt(ut1, df = 2); vt2 = qt(ut2, df = 4);
hist(qt(ut1, 2),"scott", main = "t quantile of ut1, df = 2")
hist(qt(ut2, df = 4),"scott", main = "t quantile of ut2, df = 4")
plot(vt1,vt2, main = paste("cor.hat = ", round(cor(vt1,vt2),4)))

##  skewed t marginals
wt1 = qsstd(ut1, nu = 3, xi = 0.75); wt2 = qsstd(ut2, nu = 6, xi = 1.5);
hist(qsstd(ut1, nu = 3, xi = 0.75),"scott", main = "skewed t quantile of ut1")
hist(qsstd(ut2, nu = 6, xi = 1.5),"scott", main = "skewed t quantile of ut")
plot(wt1,wt2, main = paste("cor.hat = ", round(cor(wt1,wt2),4)))


####------------------------------####
#### Gaussian and t copulas with  ####
####      the same marginals      ####
####------------------------------####
par(mfrow = c(3,2), pty = "s", font.lab = 2, mgp = c(2,0.8,0), cex.lab = 1.25)
plot(u1,u2, main = "Gaussian Copula\n rho = 0.7")
plot(ut1,ut2, main = "t-Copula\n rho = 0.7, df = 3")

xlim = range(c(v1,vt1)); ylim = range(c(v2,vt2))
plot(v1,v2, xlab = "t_2", ylab = "t_4", xlim = xlim, ylim = ylim)
mtext("Dependence: Gaussian Copula", font = 2, cex = 0.8)
plot(vt1,vt2, xlab = "t_2", ylab = "t_4", xlim = xlim, ylim = ylim)
mtext("Dependence: t-Copula", font = 2, cex = 0.8)

xlim = range(c(w1,wt1)); ylim = range(c(w2,wt2))
plot(w1,w2, xlab = "skewed t_3", ylab = "skewed t_6", xlim = xlim, ylim = ylim)
mtext("Dependence: Gaussian Copula", font = 2, cex = 0.8)
plot(wt1,wt2, xlab = "skewed t_3", ylab = "skewed t_6", xlim = xlim, ylim = ylim)
mtext("Dependence: t-Copula", font = 2, cex = 0.8)














