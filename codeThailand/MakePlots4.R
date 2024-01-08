source(paste(dirFcns,'GetData.R',sep=''))
source(paste(dirFcns,'GetParameters.R',sep=''))
source(paste(dirFcns,'GenObjects.R',sep=''))
source(paste(dirFcns,'fnFlexIntegrator2011a.R',sep=''))
source(paste(dirCode,'likelihood.R',sep=''))
source(paste(dirCode,'prior.R',sep=''))
source(paste(dirCode,'IMIS.opt.R',sep=''))
source(file = "code/plot_prevNew.R")
load(paste(model, country,region,".RData",sep=""))

# dRt = result$resample[,3]*(result$resample[,4]-result$rt) - result$resample[,5]*result$prev + result$resample[,6]*result$inc + result$resample[,7]*result$mu
# dRt = result$resample[,3]*(result$resample[,4]-result$rt) - result$resample[,5]*result$prev + result$resample[,6]*result$inc*(1-result$prev) + result$resample[,7]*result$mu*result$prev
# dRt = result$resample[,3]*(result$resample[,4]-result$rt) - result$resample[,5]*result$prev + result$resample[,6]*result$inc
dRt = result$resample[,3]*(result$resample[,4]-result$rt) - result$resample[,5]*result$prev + result$resample[,6]*result$mu



pdf(file=paste(model, country, region, "Diag.pdf", sep=""), width = 12, height = 16)
par(mex=0.7, mar = c(4.5, 4.5, 3.1, 1.1), mfrow=c(4,3))

## Infection Rate r(t)
max_rt = ceiling(max(apply(result$rt, 2, quantile, 0.975, na.rm=T), na.rm=T))
plot(na_sites[1,]*100~seq(data_start_yr,data_end_yr), xlab = "Year", lty=2, ylab = "Infection Rate r(t)", 
main=paste("beta = ", round(mean(result$resample[,3]),3), "   sd = ", round(sd(result$resample[,3]),3)), 
cex.axis=1.2, cex.lab=1.2, cex.main=1.5, xlim = c(1970, endYear), ylim = c(0, max_rt), col="red", lwd=2, type = "n")
lines(apply(result$rt, 2, quantile, 0.5, na.rm=T)~seq(beginYear, endYear), lwd=3)
lines(apply(result$rt, 2, quantile, 0.025, na.rm=T)~seq(beginYear, endYear),type="l", lty=2, lwd=3, col="black") 
lines(apply(result$rt, 2, quantile, 0.975, na.rm=T)~seq(beginYear, endYear),type="l", lty=2, lwd=3, col="black") 
abline(v=data_end_yr+0.5, lty=2, lwd=1)
abline(h=0, lty=2, lwd=1, col="red")

## Reproductive Rate R0(t)
max_R0 = ceiling(max(apply(result$rt/result$mu, 2, quantile, 0.975, na.rm=T), na.rm=T))
plot(na_sites[1,]*100~seq(data_start_yr,data_end_yr), xlab = "Year", lty=2, ylab = "R0(t)", main=paste(country,region), 
cex.axis=1.2, cex.lab=1.2, cex.main=1.5, xlim = c(1970, endYear), ylim = c(0, min(max_R0,10)), col="red", lwd=2, type = "n")
lines(apply(result$rt/result$mu, 2, quantile, 0.5, na.rm=T)~seq(beginYear, endYear), lwd=3)
lines(apply(result$rt/result$mu, 2, quantile, 0.025, na.rm=T)~seq(beginYear, endYear),type="l", lty=2, lwd=3, col="black") 
lines(apply(result$rt/result$mu, 2, quantile, 0.975, na.rm=T)~seq(beginYear, endYear),type="l", lty=2, lwd=3, col="black") 
abline(v=data_end_yr+0.5, lty=2, lwd=1)
abline(h=1, lty=2, lwd=1, col="red")

# plot(apply(result$rt/result$mu, 2, quantile, 0.5, na.rm=T), apply(dRt, 2, quantile, 0.5, na.rm=T), type="l", lwd=3,
# xlab = "Reproductive Rate R0(t)", ylab = "Infection Rate r(t)")

plot(apply(result$rt, 2, quantile, 0.5, na.rm=T), apply(dRt, 2, quantile, 0.5, na.rm=T), type="l", lwd=3,
xlab = "Infection Rate r(t)", ylab = "Change of Infection Rate r(t)",
main=paste("beta = ", round(mean(result$resample[,4]),3), "   sd = ", round(sd(result$resample[,4]),3)), ylim=c(-0.05, 0.02))
points(tail(apply(result$rt, 2, quantile, 0.5, na.rm=T),5), tail(apply(dRt, 2, quantile, 0.5, na.rm=T),5), col="red", lwd=3)
abline(h=0, lty=2, lwd=1, col="red")

## Prevalence
max_prev = max(ceiling(max(apply(result$prev, 2, quantile, 0.975, na.rm=T), na.rm=T)*20)*5, ceiling(max(na_sites, na.rm=T)*20)*5)
plot(na_sites[1,]*100~seq(data_start_yr,data_end_yr), xlab = "Year", lty=2, ylab = "Prevalence(%)", main=paste(country,region), 
cex.axis=1.2, cex.lab=1.2, cex.main=1.5, xlim = c(1970, endYear), ylim = c(0,max_prev), col="red", lwd=2, type = "n")
legend(1970, max_prev, c("fitted prevalence", "data average"), col = c("black", "red"), lty = c(1, 1), text.col = c("black", "red"))
plot.clinic()
lines(apply(result$prev, 2, quantile, 0.5, na.rm=T)*100~seq(beginYear, endYear), lwd=4)
lines(apply(result$prev, 2, quantile, 0.025, na.rm=T)*100~seq(beginYear, endYear),type="l", lty=2, lwd=2, col="black") 
lines(apply(result$prev, 2, quantile, 0.975, na.rm=T)*100~seq(beginYear, endYear),type="l", lty=2, lwd=2, col="black") 
index = which(apply(na_sites, 2, median, na.rm=T)>=0)
lines(apply(na_sites, 2, mean, na.rm=T)[index]*100~seq(data_start_yr,data_end_yr)[index], lwd=2, col="red")

plot(apply(result$prev, 2, quantile, 0.5, na.rm=T)*100, apply(result$rt, 2, quantile, 0.5, na.rm=T), type="l", lwd=3,
xlab = "Prevalence(%)", ylab = "Infection Rate r(t)", 
main=paste("beta = ", round(mean(result$resample[,5]),3), "   sd = ", round(sd(result$resample[,5]),3)))

plot(apply(result$prev, 2, quantile, 0.5, na.rm=T)*100, apply(dRt, 2, quantile, 0.5, na.rm=T), type="l", lwd=3,
xlab = "Prevalence(%)", ylab = "Change of Infection Rate r(t)", 
main=paste("beta = ", round(mean(result$resample[,5]),3), "   sd = ", round(sd(result$resample[,5]),3)), ylim=c(-0.05, 0.02))
points(tail(apply(result$prev, 2, quantile, 0.5, na.rm=T)*100,5), tail(apply(dRt, 2, quantile, 0.5, na.rm=T),5), col="red", lwd=3)
abline(h=0, lty=2, lwd=1, col="red")

## Incidence
max_inc = ceiling(max(apply(result$inc, 2, quantile, 0.975, na.rm=T), na.rm=T)*120)
plot(na_sites[1,]*100~seq(data_start_yr,data_end_yr), xlab = "Year", lty=2, ylab = "Incidence(%)", main=paste(country,region), 
cex.axis=1.2, cex.lab=1.2, cex.main=1.5, xlim = c(1970, endYear), ylim = c(0,max_inc), col="red", lwd=2, type = "n")
lines(apply(result$inc, 2, quantile, 0.5, na.rm=T)*100~seq(beginYear, endYear), lwd=2)
lines(apply(result$inc, 2, quantile, 0.025, na.rm=T)*100~seq(beginYear, endYear),type="l", lty=2, lwd=2) 
lines(apply(result$inc, 2, quantile, 0.975, na.rm=T)*100~seq(beginYear, endYear),type="l", lty=2, lwd=2)
abline(v=data_end_yr+0.5, lty=2, lwd=1)
abline(h=0, lty=2, lwd=1, col="red")

plot(apply(result$inc, 2, quantile, 0.5, na.rm=T)*100, apply(result$rt, 2, quantile, 0.5, na.rm=T), type="l", lwd=3,
xlab = "Incidence(%)", ylab = "Infection Rate r(t)", 
main=paste("beta = ", round(mean(result$resample[,6]),3), "   sd = ", round(sd(result$resample[,6]),3)))

plot(apply(result$inc, 2, quantile, 0.5, na.rm=T)*100, apply(dRt, 2, quantile, 0.5, na.rm=T), type="l", lwd=3,
xlab = "Incidence(%)", ylab = "Change of Infection Rate r(t)", 
main=paste("beta = ", round(mean(result$resample[,6]),3), "   sd = ", round(sd(result$resample[,6]),3)), ylim=c(-0.05, 0.02))
points(tail(apply(result$inc, 2, quantile, 0.5, na.rm=T)*100,5), tail(apply(dRt, 2, quantile, 0.5, na.rm=T),5), col="red", lwd=3)
abline(h=0, lty=2, lwd=1, col="red")

## AIDS Death Rate u(t)
max_mu = ceiling(max(apply(result$mu, 2, quantile, 0.975, na.rm=T), na.rm=T)*100)
plot(na_sites[1,]*100~seq(data_start_yr,data_end_yr), xlab = "Year", lty=2, ylab = "AIDS Death Rate u(t) (%)", main=paste(country,region), 
cex.axis=1.2, cex.lab=1.2, cex.main=1.5, xlim = c(1970, endYear), ylim = c(0, max_mu), col="red", lwd=2, type = "n")
lines(apply(result$mu, 2, quantile, 0.5, na.rm=T)*100~seq(beginYear, endYear), lwd=3)
lines(apply(result$mu, 2, quantile, 0.025, na.rm=T)*100~seq(beginYear, endYear),type="l", lty=2, lwd=3, col="black") 
lines(apply(result$mu, 2, quantile, 0.975, na.rm=T)*100~seq(beginYear, endYear),type="l", lty=2, lwd=3, col="black") 
abline(v=data_end_yr+0.5, lty=2, lwd=1)
abline(h=0, lty=2, lwd=1, col="red")

plot(apply(result$mu, 2, quantile, 0.5, na.rm=T)*100, apply(result$rt, 2, quantile, 0.5, na.rm=T), type="l", lwd=3,
xlab = "AIDS Death Rate u(t)", ylab = "Infection Rate r(t)")

plot(apply(result$mu, 2, quantile, 0.5, na.rm=T)*100, apply(dRt, 2, quantile, 0.5, na.rm=T), type="l", lwd=3,
xlab = "AIDS Death Rate u(t) (%)", ylab = "Change of Infection Rate r(t)", 
ylim=c(-0.05, 0.02))#, main=paste("beta = ", round(mean(result$resample[,7]),3), "   sd = ", round(sd(result$resample[,7]),3)))
points(tail(apply(result$mu, 2, quantile, 0.5, na.rm=T)*100,5), tail(apply(dRt, 2, quantile, 0.5, na.rm=T),5), col="red", lwd=3)
abline(h=0, lty=2, lwd=1, col="red")

dev.off()

jpeg(file=paste(model, country, region, "Pairs.jpg", sep=""), width = 800, height = 800)
# X_input = as.data.frame(result$resample[,3:7])
# names(X_input) = c("beta0", "beta1", "beta2", "beta3", "beta4")
X_input = as.data.frame(result$resample[,3:6])
names(X_input) = c("beta0", "beta1", "beta2", "beta3")
pairs(X_input, main = paste(country, region, "Pairs Plot"), lower.panel=panel.smooth, upper.panel=panel.cor, 
cex = 1, bg="light blue", diag.panel=panel.hist, cex.labels = 2, font.labels=2)
dev.off()

