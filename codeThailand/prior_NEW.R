########################################################################################
# Define prior(x) which calculates the prior density of x
# Define sample.prior(n) which draws n samples from the prior distribution
########################################################################################

r0.low = 1/11.5;	r0.high = 5			# lower and upper bounds for r0
y0.low = 0.1^13;	y0.high = 0.0025		# lower and upper bounds for y0
x.mean = c(0.16, -0.23, -1.03, -0.03)
x.sd =   c(0.07,  0.10,  0.32,  0.58)
prior <- function(x){
	if (is.vector(x))
		value = prod(dnorm(x[3:6], x.mean, x.sd))

	if (!is.vector(x))
		value = dnorm(x[,3], x.mean[1], x.sd[1]) * dnorm(x[,4], x.mean[2], x.sd[2]) *
			dnorm(x[,5], x.mean[3], x.sd[3]) * dnorm(x[,6], x.mean[4], x.sd[4])
	return(value)
}

sample.prior <- function(n){
	y0 = runif(n, log(y0.low), log(y0.high))		# log initial pulse at t0
	r0 = runif(n, log(r0.low), log(r0.high))		# log r(t) at t0
	input=cbind(y0, r0)
	for (i in 1:4)
	input=cbind(input, rnorm(n, x.mean[i], x.sd[i]))
	return(input)
}

