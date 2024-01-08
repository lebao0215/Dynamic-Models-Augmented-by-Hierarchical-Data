source(paste(dirFcns,'GetData.R',sep=''))
source(paste(dirFcns,'GetParameters.R',sep=''))
source(paste(dirFcns,'GenObjects.R',sep=''))
source(paste(dirFcns,'fnFlexIntegrator2011a.R',sep=''))
source(paste(dirCode,'likelihood.R',sep=''))
source(paste(dirCode,'prior.R',sep=''))
source(paste(dirCode,'IMIS.opt.R',sep=''))
source(file = "code/plot_prevNew.R")
load(paste(model, country,region,".RData",sep=""))

EPP.csv <- read.csv(paste(country,"Urban.csv",sep=""), header=F) 

pdf(file=paste(model, country, region, ".pdf", sep=""), width = 15, height = 4)
par(mex=0.7, mar = c(4.5, 4.5, 3.1, 2.5), mfrow=c(1,3))
## Prevalence
max_prev = max(ceiling(max(apply(result$prev, 2, quantile, 0.975, na.rm=T), na.rm=T)*24)*5, ceiling(max(na_sites, na.rm=T)*24)*5)
plot(na_sites[1,]*100~seq(data_start_yr,data_end_yr), xlab = "Year", lty=2, ylab = "Prevalence(%)", main=paste(country,region), 
cex.axis=1.2, cex.lab=1.3, cex.main=1.5, xlim = c(1970, endYear), ylim = c(0,max_prev), col="red", lwd=2, type = "n")
legend(1970, max_prev, c("fitted prevalence", "data average"), col = c("black", "red"), lty = c(1, 1), text.col = c("black", "red"), cex=1.4)
plot.clinic()
lines(apply(result$prev, 2, quantile, 0.5, na.rm=T)*100~seq(beginYear, endYear), lwd=4)
lines(apply(result$prev, 2, quantile, 0.025, na.rm=T)*100~seq(beginYear, endYear),type="l", lty=2, lwd=2, col="black") 
lines(apply(result$prev, 2, quantile, 0.975, na.rm=T)*100~seq(beginYear, endYear),type="l", lty=2, lwd=2, col="black") 
index = which(apply(na_sites, 2, median, na.rm=T)>=0)
lines(apply(na_sites, 2, mean, na.rm=T)[index]*100~seq(data_start_yr,data_end_yr)[index], lwd=2, col="red")

lines(as.numeric(EPP.csv[22,-(1:8)])~seq(beginYear+7, endYear), lwd=3, col="blue")

## Infection Rate r(t)
max_rt = ceiling(max(apply(result$rt, 2, quantile, 0.975, na.rm=T), na.rm=T))
min_rt = floor(min(apply(result$rt, 2, quantile, 0.025, na.rm=T), na.rm=T)*1000)/1000

plot(na_sites[1,]*100~seq(data_start_yr,data_end_yr), xlab = "Year", lty=2, ylab = "r(t) on log scale", main=paste(country,region), 
log="y", cex.axis=1.2, cex.lab=1.4, cex.main=1.5, xlim = c(1970, endYear), ylim = c(min_rt, max_rt), col="red", lwd=2, type = "n")
lines(apply(result$rt, 2, quantile, 0.5, na.rm=T)~seq(beginYear, endYear), lwd=3)
lines(apply(result$rt, 2, quantile, 0.025, na.rm=T)~seq(beginYear, endYear),type="l", lty=2, lwd=3, col="black") 
lines(apply(result$rt, 2, quantile, 0.975, na.rm=T)~seq(beginYear, endYear),type="l", lty=2, lwd=3, col="black") 
abline(v=data_end_yr+0.5, lty=2, lwd=1)

## AIDS Death Rate u(t)
max_mu = ceiling(max(apply(result$mu, 2, quantile, 0.975, na.rm=T), na.rm=T)*120)+2
plot(na_sites[1,]*100~seq(data_start_yr,data_end_yr), xlab = "Year", lty=2, ylab = NA, main=paste(country,region), 
cex.axis=1.2, cex.lab=1.2, cex.main=1.5, xlim = c(1970, endYear), ylim = c(0, max_mu), col="red", lwd=2, type = "n")
lines(apply(result$mu, 2, quantile, 0.5, na.rm=T)*100~seq(beginYear, endYear), lwd=3)
lines(apply(result$mu, 2, quantile, 0.025, na.rm=T)*100~seq(beginYear, endYear),type="l", lty=2, lwd=3, col="black") 
lines(apply(result$mu, 2, quantile, 0.975, na.rm=T)*100~seq(beginYear, endYear),type="l", lty=2, lwd=3, col="black") 
abline(v=data_end_yr+0.5, lty=2, lwd=1)

## Incidence
max_inc = ceiling(max(apply(result$inc, 2, quantile, 0.975, na.rm=T), na.rm=T)*120)
# plot(na_sites[1,]*100~seq(data_start_yr,data_end_yr), xlab = "Year", lty=2, ylab = "Incidence(%)", main=paste(country,region), 
# cex.axis=1.2, cex.lab=1.2, cex.main=1.5, xlim = c(1970, endYear), ylim = c(0,max_inc), col="red", lwd=2, type = "n")
axis(4, col = "red", col.axis="red", lwd = 2, at=seq(0,max_inc,2)*max_mu/max_inc, labels=seq(0,max_inc,2))
lines(apply(result$inc, 2, quantile, 0.5, na.rm=T)*100*max_mu/max_inc~seq(beginYear, endYear), lwd=3, col="red")
lines(apply(result$inc, 2, quantile, 0.025, na.rm=T)*100*max_mu/max_inc~seq(beginYear, endYear),type="l", lty=2, lwd=3, col="red") 
lines(apply(result$inc, 2, quantile, 0.975, na.rm=T)*100*max_mu/max_inc~seq(beginYear, endYear),type="l", lty=2, lwd=3, col="red")
legend(1970, max_mu, c("AIDS Death Rate (%)", "Incidence (%)"), col = c("black", "red"), lty = c(1, 1), text.col = c("black", "red"), cex=1.4)

lines(as.numeric(EPP.csv[24,-(1:8)]/(EPP.csv[25,-(1:8)]-EPP.csv[23,-(1:8)]))*100~seq(beginYear+7, endYear), lwd=3, col="blue")

dev.off()

