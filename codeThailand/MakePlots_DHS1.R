source(paste(dirFcns,'GetData.R',sep=''))
source(paste(dirFcns,'GetDHS.R',sep=''))
source(paste(dirFcns,'GetParameters.R',sep=''))
source(paste(dirFcns,'GenObjects.R',sep=''))
source(paste(dirFcns,'fnFlexIntegrator2011.R',sep=''))
source(paste(dirCode,'likelihood.R',sep=''))
source(paste(dirCode,'prior.R',sep=''))
source(paste(dirCode,'IMIS.opt.R',sep=''))
source(file = "code/plot_prevNew.R")

load(paste("DHS/Incidence",senario,"/", country, region,".RData",sep=""))
year_min = quantile(result$resample[,1], 0.15)
prev.X = result$prev; 			inc.X = result$inc[,(year_min-beginYear+1):(endYear-beginYear+1)]
load(paste("DHS/EPP/", country, region,".RData",sep=""))
prev.EPP = result$prev; 	inc.EPP = result$inc[,(year_min-beginYear+1):(endYear-beginYear+1)]

############################################################
####			Simulation
############################################################
beta.tri = 0.025;		Omega.tri = 150;		T.tri =450
sample.trinomial <- function(theta){
	n = theta[1]
	prev = theta[2]
	inc = theta[3]
	p = rep(0,3)	# negative, non-recent, recent	
	p[1] = 1-prev
	p[3] = (1-prev)*inc*(Omega.tri-beta.tri*T.tri)/365 + beta.tri*prev
	p[2] = prev-p[3]
	n.tri = sample(1:3, n, replace = TRUE, prob = p)
	return(n.tri)
}
N.tri = matrix(0, dim(result$inc)[1], 3)
for (i in 1:dim(result$inc)[1]){
	set.seed(1000+i)
	X = sample.trinomial(c(5000, result$prev[i,35], result$inc[i,35]))
	if (senario=="F")		X = sample.trinomial(c(5000, result$prev[i,35], result$inc[i,35]/2))
	if (senario=="G")		X = sample.trinomial(c(5000, result$prev[i,35], result$inc[i,35]*2))
	N.tri[i,1] = length(which(X==1))
	N.tri[i,2] = length(which(X==2))
	N.tri[i,3] = length(which(X==3))
}
n.tri = round(apply(N.tri, 2, mean))

if (senario=="B")		beta.tri = 0.015
if (senario=="C")		beta.tri = 0.035
if (senario=="D")		Omega.tri = 130
if (senario=="E")		Omega.tri = 170
likelihood.tri <- function(theta){
	l = -100000000
	prev.tri = theta[1]
	inc.tri = theta[2]
	p1 = 1-prev.tri
	p3 = (1-prev.tri)*inc.tri*(Omega.tri-beta.tri*T.tri)/365 + beta.tri*prev.tri
	p2 = prev.tri-p3
	if (p1>0 & p2>0 & p3>0 & p1<1 & p2<1 & p3<1)
	l = n.tri[1]*log(p1) + n.tri[2]*log(p2) + n.tri[3]*log(p3)
	return(-l)
}
theta1.hat = 1-n.tri[1]/sum(n.tri)
theta2.hat = (n.tri[3]-beta.tri*(n.tri[2]+n.tri[3]))/n.tri[1]/(Omega.tri-beta.tri*T.tri)*450
theta.ini = c(theta1.hat, theta2.hat)
likelihood.tri(theta.ini)
optimizer = optim(theta.ini, likelihood.tri, method = "BFGS", hessian = TRUE,
control = list(parscale = rep(0.001,2), maxit = 500))
theta.opt = 100*optimizer$par
theta.cov = solve(optimizer$hessian)
theta.cor = theta.cov[1,2]/sqrt(theta.cov[1,1]*theta.cov[2,2])
theta.sd = 100*sqrt(diag(theta.cov))

# pdf(file=paste("Plot/", country, region, "Senario", senario, ".pdf", sep=""), width = 12, height = 4.5)
jpeg(file=paste("JPG/DHS1/", country, region, "Senario", senario, ".jpg", sep=""), width = 900, height = 400, pointsize = 12, quality = 75)
par(mex=0.7, mar = c(4.5, 4.5, 3, 0.5), mfrow=c(1,2))

## Prevalence
max_prev = max(ceiling(max(apply(result$prev, 2, quantile, 0.975, na.rm=T), na.rm=T)*20)*5, ceiling(max(na_sites, na.rm=T)*20)*5)
plot(na_sites[1,]*100~seq(data_start_yr,data_end_yr), xlab = "Year", lty=2, ylab = "Prevalence(%)", main=paste(country,region,"Fixed beta and Omega"), 
cex.axis=1.2, cex.lab=1.2, cex.main=1.5, xlim = c(1970, endYear), ylim = c(0,max_prev), col="red", lwd=2, type = "n")
plot.clinic()
lines(apply(prev.EPP, 2, quantile, 0.5, na.rm=T)*100~seq(beginYear, endYear), lwd=4, col="black")
lines(apply(prev.EPP, 2, quantile, 0.025, na.rm=T)*100~seq(beginYear, endYear),type="l", lty=2, lwd=2, col="black") 
lines(apply(prev.EPP, 2, quantile, 0.975, na.rm=T)*100~seq(beginYear, endYear),type="l", lty=2, lwd=2, col="black") 
lines(apply(prev.X, 2, quantile, 0.5, na.rm=T)*100~seq(beginYear, endYear), lwd=4, col="blue")
lines(apply(prev.X, 2, quantile, 0.025, na.rm=T)*100~seq(beginYear, endYear),type="l", lty=2, lwd=2, col="blue") 
lines(apply(prev.X, 2, quantile, 0.975, na.rm=T)*100~seq(beginYear, endYear),type="l", lty=2, lwd=2, col="blue") 
segments(2005, theta.opt[1]-2*theta.sd[1], 2005, theta.opt[1]+2*theta.sd[1], col = "brown", lty=1, lwd=2)
points(2005, theta.opt[1],  col = "brown", lwd = 2, pch = 16, cex=1.7)

# legend(1970, max_prev, c("EPP w/o TRI", "EPP with TRI", "TRI prior"), cex=1.2, pch = c(NA,NA,16), 
# text.col = c("black", "blue", "brown"), col = c("black", "blue", "brown"), lty = c(1, 1, 0))
legend(1970, max_prev, c("EPP w/o TRI", "EPP with TRI", "TRI prior", "NPBS data"), cex=1.2, pch = c(NA,NA,16,16), 
text.col = c("black", "blue", "brown", "red"), col = c("black", "blue", "brown", "red"), lty = c(1, 1, 0, 0))
points(NPBS$rate*100~NPBS$year, col = "red", lwd = 2, pch = 16, cex=1.7)

## Incidence
max_inc = ceiling(max(apply(result$inc, 2, quantile, 0.975, na.rm=T), na.rm=T)*120)
plot(na_sites[1,]*100~seq(data_start_yr,data_end_yr), xlab = "Year", lty=2, ylab = "Incidence(%)", main=paste(country,region,"Fixed beta and Omega"), 
cex.axis=1.2, cex.lab=1.2, cex.main=1.5, xlim = c(1970, endYear), ylim = c(0,max_inc), col="red", lwd=2, type = "n")
lines(apply(inc.EPP, 2, quantile, 0.5, na.rm=T)*100~seq(year_min, endYear), lwd=4, col="black")
lines(apply(inc.EPP, 2, quantile, 0.025, na.rm=T)*100~seq(year_min, endYear),type="l", lty=2, lwd=2, col="black") 
lines(apply(inc.EPP, 2, quantile, 0.975, na.rm=T)*100~seq(year_min, endYear),type="l", lty=2, lwd=2, col="black") 
lines(apply(inc.X, 2, quantile, 0.5, na.rm=T)*100~seq(year_min, endYear), lwd=4, col="blue")
lines(apply(inc.X, 2, quantile, 0.025, na.rm=T)*100~seq(year_min, endYear),type="l", lty=2, lwd=2, col="blue") 
lines(apply(inc.X, 2, quantile, 0.975, na.rm=T)*100~seq(year_min, endYear),type="l", lty=2, lwd=2, col="blue") 
segments(2005, theta.opt[2]-2*theta.sd[2], 2005, theta.opt[2]+2*theta.sd[2], col = "brown", lty=1, lwd=2)
points(2005, theta.opt[2],  col = "brown", lwd = 2, pch = 16, cex=1.7)
abline(h=0, lty=2, lwd=1, col="red")
legend(1970, max_inc, c("EPP w/o TRI", "EPP with TRI", "TRI prior"), text.col = c("black", "blue", "brown"), cex=1.2,
col = c("black", "blue", "brown"), lty = c(1, 1, 0), pch = c(NA,NA,16)) 
legend(1970, max_inc/2, c(paste(expression(beta), "=", beta.tri), paste(expression(Omega), "=", Omega.tri)), cex=1.2,lty = c(0, 0)) 
dev.off()


