source(paste(dirFcns,'GetData.R',sep=''))
source(paste(dirFcns,'GetParameters.R',sep=''))
source(paste(dirFcns,'GenObjects.R',sep=''))
source(paste(dirFcns,'fnFlexIntegrator2011a.R',sep=''))
source(paste(dirCode,'likelihood.R',sep=''))
source(paste(dirCode,'prior.R',sep=''))
source(paste(dirCode,'IMIS.opt.R',sep=''))
source(file = "code/plot_prevNew.R")

load(paste(model2, country, region,".RData",sep=""));	t0=1975
prev.A = result$prev; 	inc.A = result$inc
load(paste(model1, country, region,".RData",sep=""))

## Prevalence
pdf(file=paste(model2, country, region, ".pdf", sep=""), width = 7.5, height = 6)
par(mex=0.7, mar = c(4.5, 4.5, 3.1, 2.5))
max_prev = max(ceiling(max(apply(result$prev, 2, quantile, 0.975, na.rm=T), na.rm=T)*20)*5, ceiling(max(na_sites, na.rm=T)*20)*5)
plot(na_sites[1,]*100~seq(data_start_yr,data_end_yr), xlab = "Year", lty=2, ylab = "Prevalence(%)", main=paste(country,region), 
cex.axis=1.2, cex.lab=1.2, cex.main=1.5, xlim = c(1970, endYear), ylim = c(0,max_prev), col="red", lwd=2, type = "n")
legend(1970, max_prev, c("R-flex model", "EPP model", "data average"), col = c("blue", "black", "red"), 
lty = c(1, 1, 1), text.col = c("blue", "black", "red"), cex=1.4)
plot.clinic()
lines(apply(result$prev, 2, quantile, 0.5, na.rm=T)*100~seq(beginYear, endYear), lwd=4)
lines(apply(result$prev, 2, quantile, 0.025, na.rm=T)*100~seq(beginYear, endYear),type="l", lty=2, lwd=2, col="black") 
lines(apply(result$prev, 2, quantile, 0.975, na.rm=T)*100~seq(beginYear, endYear),type="l", lty=2, lwd=2, col="black") 
lines(apply(prev.A, 2, quantile, 0.5, na.rm=T)*100~seq(beginYear, endYear), lwd=4, col="blue")
lines(apply(prev.A, 2, quantile, 0.025, na.rm=T)*100~seq(beginYear, endYear),type="l", lty=2, lwd=2, col="blue") 
lines(apply(prev.A, 2, quantile, 0.975, na.rm=T)*100~seq(beginYear, endYear),type="l", lty=2, lwd=2, col="blue") 
index = which(apply(na_sites, 2, median, na.rm=T)>=0)
lines(apply(na_sites, 2, mean, na.rm=T)[index]*100~seq(data_start_yr,data_end_yr)[index], lwd=2, col="red")
dev.off()

## Incidence
pdf(file=paste(model2, country, region, "Inc.pdf", sep=""), width = 7.5, height = 6)
par(mex=0.7, mar = c(4.5, 4.5, 3.1, 2.5))
max_inc = ceiling(max(apply(result$inc, 2, quantile, 0.975, na.rm=T), na.rm=T)*120)
plot(na_sites[1,]*100~seq(data_start_yr,data_end_yr), xlab = "Year", lty=2, ylab = "Incidence(%)", main=paste(country,region), 
cex.axis=1.2, cex.lab=1.2, cex.main=1.5, xlim = c(1970, endYear), ylim = c(0,max_inc), col="red", lwd=2, type = "n")
lines(apply(result$inc, 2, quantile, 0.5, na.rm=T)*100~seq(beginYear, endYear), lwd=2)
lines(apply(result$inc, 2, quantile, 0.025, na.rm=T)*100~seq(beginYear, endYear),type="l", lty=2, lwd=2) 
lines(apply(result$inc, 2, quantile, 0.975, na.rm=T)*100~seq(beginYear, endYear),type="l", lty=2, lwd=2)
lines(apply(inc.A, 2, quantile, 0.5, na.rm=T)*100~seq(beginYear, endYear), lwd=2, col="blue")
lines(apply(inc.A, 2, quantile, 0.025, na.rm=T)*100~seq(beginYear, endYear),type="l", lty=2, lwd=2, col="blue") 
lines(apply(inc.A, 2, quantile, 0.975, na.rm=T)*100~seq(beginYear, endYear),type="l", lty=2, lwd=2, col="blue")
legend(1970, max_inc, c("R-flex model", "EPP model"), col = c("blue", "black"), 
lty = c(1, 1), text.col = c("blue", "black"), cex=1.4)
abline(h=0, lty=2, lwd=1, col="red")
dev.off()


