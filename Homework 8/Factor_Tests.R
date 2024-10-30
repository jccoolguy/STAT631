#### The R code for factor models 
#### 1. wald.fun() ## 2. lrt.fun() ## 3.cov.diag.test()


#### --------------------------------------------------------------------
#### 			wald test with least squares estimates
#### --------------------------------------------------------------------
#### est: a vector of OLS estimates
#### est.var: the variance of est
#### n: the sample size
#### p: number of regressors in the linear factor model
#### W: optional, the weight matrix
#### return value--a vector of Wald stat and its pvalue and degreesof freedom

wald.fun = function(est, est.var, n, p, W = diag(length(est))){
	df1 = dim(W)[2]
	if( df1 > length(est) | df1 > dim(W)[1]) stop("Error in dimension of W")
	if(is.null(W) == F & det(t(W)%*%W) < 10e-7) stop("W is not Full rank")
	df2 = n - df1 - p
	west = t(W)%*%est
	wSw  = t(W)%*%est.var%*%W
	test = df2/(df1*n)*t(west)%*%solve(wSw)%*%west
	test = as.vector(test)
	return(c(Wald = test, p.value = 1 - pf(test, df1, df2), df1 = df1, df2 = df2))
}


#### -------------------------------------------------
#### 			Likwlihood Ratio Test 
#### ------------------------------------------------
#### est.var: the variance of est
#### sig: residual variance matrix  
#### sig0: residual covariance under H0
#### n: the sample size


lrt.fun = function(sig, sig0, n){
	if(diff(dim(sig)) !=0 | diff(dim(sig0)) !=0 |dim(sig)[1] != dim(sig0)[1]) stop("Incorrect input matrices.")
	df = dim(sig)[1]
	test = (n-df/2-2)*(log(det(sig0))-log(det(sig)))
	return(c(LRT = test, p.value = 1-pchisq(test,df), df = df))
}

#### -------------------------------------------------
#### 		Testing a matrix is block diagonal
#### --------------------------------------------------
#### Sig: the matrix to be tested
#### Ns: a vector of block dimensions
#### n: the sample size
#### p: number of regressors in the linear factor model

cov.diag.test = function(Sig, Ns, n,  p){
	N = dim(Sig)[1]
	if(sum(Ns)!= N) stop("Sum of Ns is not equal to dimension of Sig!")
	if(abs(det(Sig)) < 10e-7) stop("Sig is singular!")
	nu = n - p - 1;  ## DF of Sig
	a2 = N^2 -sum(Ns^2); a3 = N^3-sum(Ns^3) 
	k = length(Ns)
	if(k == N){
		ttl = "*** Testing if the matrix is diagonal ***"
		test = -(nu-(2*N + 5)/6)*log(det(Sig)/prod(diag(Sig)))
	}
	if(k < N){
		ttl = "*** Testing if the matrix is block diagonal ***"
		i.1 = cumsum(Ns)
		i.0 = i.1 - Ns + 1
		dets = sapply(1:k, function(u){is = i.0[u]:i.1[u]; det(as.matrix(Sig[is,is]))})
		c = 1 - 1/(6*a2*nu)*(2*a3 + 3*a2)
		test = -nu*c*log(det(Sig)/prod(dets))
	}
	out = c(LRT = test, p.value = 1 - pchisq(test,a2/2), df = a2/2);
	cat(ttl)
	cat("\nLRT -statistic:", round(out[1],4),"\t p-value:", round(out[2],4), "\tDF:", a2/2); 
	invisible(out)	
}
