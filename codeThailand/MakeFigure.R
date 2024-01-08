country = "Nigeria"
dirData <- 'Data/'; dirFcns <- 'fcnsNigeria/'; dirCode <- 'codeNigeria/'
option.truncate = 0
source(paste(dirFcns,'GetParameters.R',sep=''))    # Load key projection parameters
source(paste(dirFcns,'GenObjects.R',sep=''))		# Load population arrays, treatment arrays, storage arrays for model integrator
data.id = "Test"
source(paste(dirFcns,'GetData.R',sep=''))			# Load data
source(paste(dirFcns,'GetDHS.R',sep=''))

load(file = paste('Result/',country, site,"Train.RData",sep=""))
prev.old = result$prev

par(mex=0.7, mar = c(4.5, 4.5, 5.5, 2.5))
source(file = paste(dirCode, "plot_prevNew.R", sep=""))
max_prev = max(data$P_table*100, na.rm=T)
max_prev = max(c(max_prev, apply(prev.old,2,quantile,0.975,na.rm=T)*100), na.rm=T)
na_sites = data$P_table;  isna_sites = is.na(na_sites)
nr_sites = dim(na_sites)[1];    nr_years = dim(na_sites)[2]


plot(prev.old[1,]*100~seq(1970,2015), 
     xlab = "Year", lty=2, ylab = "Prevalence(%)", cex.axis=1.2, cex.lab=1.2, cex.main=1.5, 
     main=paste(country, site), xlim = c(1970, 2015), ylim = c(0,max_prev), col="red", lwd=2, type = "n")
for (i in 1:nrow(prev.old))
  lines(prev.old[i,]*100~seq(1970,2015), lwd=1, lty=2, col="gray")
plot.clinic()

lines(apply(prev.old,2,mean,na.rm=T)*100~seq(1970,2015), lwd=4, col="black")
lines(apply(prev.old,2,quantile,0.025,na.rm=T)*100~seq(1970,2015), lwd=2, lty=2, col="black")
lines(apply(prev.old,2,quantile,0.975,na.rm=T)*100~seq(1970,2015), lwd=2, lty=2, col="black")

# lines(apply(prev.new,2,mean,na.rm=T)*100~seq(1970,2015), lwd=4, col="blue")
# lines(apply(prev.new,2,quantile,0.025,na.rm=T)*100~seq(1970,2015), lwd=2, lty=2, col="blue")
# lines(apply(prev.new,2,quantile,0.975,na.rm=T)*100~seq(1970,2015), lwd=2, lty=2, col="blue")

legend(1970, max_prev, c("Original EPP", "Imputed Model", "Data Average"), text.col = c("black", "blue", "red"), cex=1.2,
       col = c("black", "blue", "red"), lty = c(1, 1, 1)) 
points(NPBS$rate*100~NPBS$year, col = "red", lwd = 2, pch = 16, cex=1.7)

  