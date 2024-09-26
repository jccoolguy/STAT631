######################################
##			R code for HW 4			##
######################################
##	get data
library(quantmod);
syb = c("PARA","CMCSA");  d = length(syb)
yt = c()
for(i in 1:d){
	getSymbols(syb[i], from = "2011-01-01", to = "2024-09-14")
	yt = cbind(yt,weeklyReturn(Ad(get(syb[i])), type ="log"))
}
colnames(yt) = syb
yt = as.matrix(100*yt)  ## convert to % for numerical stability
n = dim(yt)[1]

######################################
##	fit 6 copulas to ut
##	Store the 6 fits into cops
library(copula)
copNames = c("t", "normal", "frank","clayton", "gumbel", "joe")
cops = vector("list",6);names(cops) = copNames

for(i in 1:2){ ## 2 ellipticals
	cops[[i]] = fitCopula(copula = ellipCopula(copNames[i],dim = 2), data = ut, method = "ml")
}
for(i in 3:6){ ## 4 archimedeans 
	cops[[i]] = fitCopula(copula = archmCopula(copNames[i],dim = 2), data = ut, method = "ml")
}

######################################
##	Plot parametric and nonparametric 
##	pdf countours overlaid	
######################################
## Create nonparametric copula and its pdf contours
Udex = ((1:n)-.05)/(n)
Cn = C.n(u=cbind(rep(Udex,n),rep(Udex,each=n)), X =ut)
empCop = expression(contour(Udex, Udex, matrix(Cn, n, n), col = "red3", add = TRUE))

## plot contours
## Adjust the size if my par() does not produce
## 		clear enough plots
## There will be 48 warnings about normal copula
## 		you can ignore them but please hide them,
## 		at beginning of the chunk set ```{r warning=FALSE}
 					
par(mfrow = c(2,3), pty = "s", lwd = 1.1, mar = c(4,4,2,1))
nu = round(cops$t@estimate[2]) ## pCoupla takes only integer df
for(i in 1:2){  ## plot the 2 ellipticals
	contour(ellipCopula(copNames[i], param = cops[[i]]@estimate[1], dim = 2, df = nu), pCopula, main = copNames[i])
	eval(empCop)
}
for(i in 3:6){ ## plot the 4 archimedeans 
	contour(archmCopula(copNames[i], param = cops[[i]]@estimate, dim = 2), pCopula, main = copNames[i])
	eval(empCop)
}


