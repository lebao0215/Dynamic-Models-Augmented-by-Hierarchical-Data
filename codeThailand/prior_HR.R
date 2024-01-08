########################################################################################
# Define prior(x) which calculates the prior density of x
# Define sample.prior(n) which draws n samples from the prior distribution
########################################################################################
## 7 lines of changes marked by ## new ##, Le Bao, 09/14/2012

t0.low = 1969.5;	t0.high = 1990.5			## new ##
# the mean and standard deviation for the prior distribution of addtional parameters
x.mean = c(20,  0.42, 0.46, 0.17, -0.68, -0.038)		## generalized ##
x.sd =   c(4.5, 0.23, 0.12, 0.07, 0.24, 0.009)			## generalized ##
# x.sd[-1] = x.sd[-1]*5
X = result$resample;	x.mean.hr = X[,1]
for (j in 1:length(x.mean))	x.mean.hr = cbind(x.mean.hr, (X[,j+1]*tau + x.mean[j]) / (tau+1))
x.sd.hr = c(1, x.sd * sqrt((2*tau+1)/(tau^2+2*tau+1)))

if (region=="Urban"){
	x.mean.hr = cbind(x.mean.hr, rep(0.11, nrow(X)));	x.sd.hr = c(x.sd.hr, 0.04)
}
if (region=="Rural"){
	x.mean.hr = cbind(x.mean.hr, rep(0.17, nrow(X)));	x.sd.hr = c(x.sd.hr, 0.05)
}

# calculate the prior density
prior <- function(x){
	if (is.vector(x)){
		value = 0
		for (i in 1:nrow(x.mean.hr))
		value = value + prod(dnorm(x.mean.hr[i,], x, x.sd.hr))
		value = value/nrow(x.mean.hr)
	}
	if (!is.vector(x))
	{
		value = NULL
		for (j in 1:nrow(x)){
			value.j = 0
			for (i in 1:nrow(x.mean.hr))
			value.j = value.j + prod(dnorm(x.mean.hr[i,], x[j,], x.sd.hr))
			value.j = value.j/nrow(x.mean.hr)
			value = c(value, value.j)
		}
	}
	return(value)
}

# draw samples from the prior distribution
sample.prior <- function(n){
	input = NULL
	for (i in 1:nrow(X)){
		input.i = NULL
		for (j in 1:ncol(X))
		input.i = cbind(input.i, rnorm(n, x.mean.hr[i,j], x.sd.hr[j]))
		input = rbind(input, input.i)
	}
	input[,1] = round(input[,1]);	input[,2] = round(input[,2])
	return(input)
}
