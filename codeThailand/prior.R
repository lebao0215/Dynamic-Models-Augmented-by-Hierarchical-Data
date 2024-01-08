########################################################################################
# Define prior(x) which calculates the prior density of x
# Define sample.prior(n) which draws n samples from the prior distribution
########################################################################################
## 7 lines of changes marked by ## new ##, Le Bao, 09/14/2012

t0.low = 1969.5;	t0.high = 1990.5			## new ##
# the mean and standard deviation for the prior distribution of addtional parameters
x.mean = c(20,  0.42, 0.17, 0.46, -0.68, -0.038)		## generalized ##
x.sd =   c(4.5, 0.23, 0.07, 0.12, 0.24, 0.009)			## generalized ##

# calculate the prior density
prior <- function(x){
	value = 1
	if (is.vector(x))
		value = prod(dnorm(x[-1], x.mean, x.sd))			## new ##
	if (!is.vector(x))
		for (i in 1:length(x.mean)) 					## new ##
		value = value * dnorm(x[,i+1], x.mean[i], x.sd[i]) 	## new ##
	return(value)
}

# draw samples from the prior distribution
sample.prior <- function(n){
	input = t0 = round(runif(n, t0.low, t0.high))			# starting time of epidemic
	for (i in 1:length(x.mean))
	input=cbind(input, rnorm(n, x.mean[i], x.sd[i]))
	input[,2] = round(input[,2])						## new ##
	return(input)
}
