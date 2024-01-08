setwd('MeanShiftY0/')
dirData <- 'SSA/'; dirFcns <- 'fcns/'; dirCode <- 'code/'

#Load data
source(paste(dirFcns,'GetData.R',sep=''))
#Load key projection parameters
source(paste(dirFcns,'GetParameters.R',sep=''))
#Load population arrays, treatment arrays, storage arrays for model integrator, and indices for calculating equilibrium prior for force of infection paramter (r)
source(paste(dirFcns,'GenObjects.R',sep=''))
#Load disease model function
source(paste(dirFcns,'fnFlexIntegrator2011a.R',sep=''))
#Load likelihood function and prior function
source(paste(dirCode,'likelihood.R',sep=''))
source(paste(dirCode,'prior.R',sep=''))
source(paste(dirCode,'IMIS.opt.R',sep=''))
library(mvtnorm)

source(paste(dirCode,'prior_New.R',sep=''))
resample_X = sample.prior(3000)

prev = inc = rt = NULL
for (i in 1:3000){
	xx = fnRMS2011(resample_X[i,], turnoverOff=TRUE)
	prev = rbind(prev, xx$annualFittingPrevs)
	inc = rbind(inc, xx$annualIncid)
	rt = rbind(rt, xx$annualR)
}

par(mex=0.7, mar = c(4.5, 4.5, 3.1, 1.1))
plot(apply(prev, 2, quantile, 0.5, na.rm=T)*100~seq(beginYear, endYear),
xlab = "Year", ylab = "Prevalence(%)", main="Samples From The Prior Distribution", 
lty=2, cex.axis=1.2, cex.lab=1.2, cex.main=1.5, lwd=4, xlim = c(1970, endYear), ylim = c(0,100))
for (i in 1:3000)
lines(prev[i,]*100~seq(beginYear, endYear),type="l", lty=2, lwd=1, col="gray") 
lines(apply(prev, 2, quantile, 0.5, na.rm=T)*100~seq(beginYear, endYear),type="l", lty=1, lwd=4, col="black") 
lines(apply(prev, 2, quantile, 0.025, na.rm=T)*100~seq(beginYear, endYear),type="l", lty=2, lwd=2, col="black") 
lines(apply(prev, 2, quantile, 0.975, na.rm=T)*100~seq(beginYear, endYear),type="l", lty=2, lwd=2, col="black") 


country.list <- c("Botswana", "Kenya", "Uganda", "Lesotho", "Liberia",
			"Namibia", "Nigeria", "Tanzania", "Ethiopia")
region.list <- c('Urban', 'Rural')
resample.all = NULL
for (country in country.list)	for (region in region.list){{
load(paste(country,region,".RData",sep=""));	resample.all = rbind(resample.all, result$resample);	}}

apply(resample.all,2,mean)
cov(resample.all)


