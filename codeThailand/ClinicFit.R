country = "BELARUS"
region.list = c("IDU", "SW", "MSM");	Turnover = c(8, 5, 10000)
dirData <- paste('Le Bao/', country, "/", sep="")		# the data folder
dirFcns <- 'fcns/'	# the deterministic functions folder
dirCode <- 'code/'	# the Bayesian melding folder
source(paste(dirFcns,'GetParameters.R',sep=''))		# Load key projection parameters
source(paste(dirFcns,'GenObjects.R',sep=''))		# Load population arrays, treatment arrays, storage arrays for model integrator
DATA = list()
for (data.k in 1:length(region.list)){
	region = region.list[data.k]
	source(paste(dirFcns,'GetData.R',sep=''))			# Load data
	DATA[[data.k]] = dataOut
}

Nst = dataOut$data$N_table
sample_size_pop = spectrumObjects$p15to49[(data_start_yr-min(popData[,1])+1):(data_end_yr-min(popData[,1])+1)]
xst = (dataOut$data$P_table*dataOut$data$N_table+0.5) / (dataOut$data$N_table+1)
Wst = qnorm(xst)
var_eps = 2*pi*xst*(1-xst)/dataOut$data$N_table*exp(Wst^2)
alpha = 0.58; beta = 93

row.names = c(row.names, paste(country, region))

for (option.hr in 1:3){

ptm.opt = proc.time()
if (option.epp==1){
	load(paste(model1, country, region,".RData",sep=""))
	prev.model = result$prev
	tstart = 1970
}
if (option.epp==2){
	load(paste(model2, country, region,".RData",sep=""))
	prev.model = result$prev
	par.current = round(apply(result$resample[,3:5], 2, mean),2)
	tstart = 1970
}
if (option.epp==3){
	load(paste(model3, country, region,".RData",sep=""))
	prev.model = prev.A
	tstart = 1970
}

set.seed(1000)
xx = array(0, c(nr_sites,dim(prev.model)[2],dim(prev.model)[1]))
for (i in 1:dim(prev.model)[1]){
	prev = prev.model[i,]
	compare_prev = prev[(data_start_yr-tstart+1):(data_end_yr-tstart+1)]
	rho = compare_prev[1:nr_years]
	tr_diff = matrix(NA,nr_sites, nr_years)	# difference between transformed prev and transformed obs
	var_su = prev.random = matrix(NA, nr_sites, length(prev))
	b.mean = b.var = b.random = rep(NA, nr_sites)

	for (site in 1:nr_sites){
		# sites is matrix with observed prevalence (row = site)
		site_index = which(isna_sites[site,]==F)
		xst_site = xst[site,site_index];	Wst_site = Wst[site,site_index];	var_site = var_eps[site,site_index]
		tr_diff[site,site_index] = Wst_site - qnorm( (rho[site_index]*sample_size_pop[site_index] +0.5)/(sample_size_pop[site_index]+1))

		b.var[site] = 1/(sum(1/var_eps[site,site_index]))
		b.mean[site] = sum(tr_diff[site,site_index]/var_eps[site,site_index])*b.var[site]
		b.temp = rnorm(100, b.mean[site], sqrt(b.var[site]))
		b.density = (b.temp^2/2+1/beta)^(-alpha-0.5)
		b.density = b.density/sum(b.density) 
		b.index = sample(100, 1, replace = TRUE, prob = b.density)			# Draw resamples
		b.random[site] = b.temp[b.index]

		Wst_site = b.random[site] + qnorm(prev)
		xst_site = pnorm(Wst_site)
		Nst_site = rep(NA, endYear-beginYear+1)
		Nst_site[(data_start_yr-tstart+1):(data_end_yr-tstart+1)] = Nst[site,]
		Nst_site[is.na(Nst_site)] = 300
		var_su[site,] = 2*pi*xst_site*(1-xst_site)/Nst_site*exp(Wst_site^2)
		site_index = which(prev>0)
		prev.random[site,site_index] = pnorm(qnorm(prev[site_index]) + b.random[site] + rnorm(length(site_index), 0, sqrt(var_su[site,site_index])) )		
	}
	xx[,,i] = prev.random
}

clinic.upper = clinic.med = clinic.lower = NULL
for (site in 1:nr_sites){
	clinic.lower = rbind(clinic.lower, apply(xx[site,,], 1, quantile, 0.025, na.rm=T))
	clinic.med = rbind(clinic.med, apply(xx[site,,], 1, quantile, 0.5, na.rm=T))
	clinic.upper = rbind(clinic.upper, apply(xx[site,,], 1, quantile, 0.975, na.rm=T))
}
clinic.width = clinic.upper[,(data_start_yr-tstart+1):(data_end_yr-tstart+1)]-clinic.lower[,(data_start_yr-tstart+1):(data_end_yr-tstart+1)]

sites_new = xst
A = sites_new > clinic.lower[,(data_start_yr-tstart+1):(data_end_yr-tstart+1)] & sites_new < clinic.upper[,(data_start_yr-tstart+1):(data_end_yr-tstart+1)]
B = !is.na(sites_new)
ss2 = mean((sites_new[which(sites_new>=0)]-clinic.med[,(data_start_yr-tstart+1):(data_end_yr-tstart+1)][which(sites_new>=0)])^2)
mae = mean(abs(sites_new[which(sites_new>=0)]-clinic.med[,(data_start_yr-tstart+1):(data_end_yr-tstart+1)][which(sites_new>=0)]))

if (option.epp==1){
	width.epp = c(width.epp, mean(clinic.width))
	cover.epp = c(cover.epp, sum(A,na.rm=T)/sum(B,na.rm=T))
	mae.epp = c(mae.epp, mae)
	time.epp = c(time.epp, time.used[3])
}
if (option.epp==2){
	size.flex = c(size.flex, sum(B,na.rm=T))
	width.rflex = c(width.rflex, mean(clinic.width))
	cover.rflex = c(cover.rflex, sum(A,na.rm=T)/sum(B,na.rm=T))
	mae.rflex = c(mae.rflex, mae)
	time.rflex = c(time.rflex, time.used[3])
	par.rflex = rbind(par.rflex, apply(result$resample, 2, mean))
	cor.rflex = rbind(cor.rflex, as.vector(cor(result$resample)))
}
if (option.epp==3){
	width.rstoch = c(width.rstoch, mean(clinic.width))
	cover.rstoch = c(cover.rstoch, sum(A,na.rm=T)/sum(B,na.rm=T))
	mae.rstoch = c(mae.rstoch, mae)
	time.rstoch = c(time.rstoch, time.used[3])
}
	mae.current = c(mae.current, mae)
	ptm.use = (proc.time() - ptm.opt)[3]
	print(paste("time used=", round(ptm.use,1), "seconds"))
} # option.epp


if (0){
pdf(file=paste(country_data, "Clinic", option.proj, ".pdf", sep=""), width = 6, height = ceiling(dim(sites_old)[1]/3)*2)
par(mar = c(3, 3, 2, 1), mfrow=c(ceiling(dim(sites_old)[1]/3),3))
prev_max = max(c(prev_max, max(sites_old,na.rm=T)))*110
for (site in 1:nr_sites){
	plot(sites_old[site,]*100~as.vector(years), lty=2, xlab=NA, ylab=NA, main=paste("Clinic", site),
	cex.axis=1.5, cex.lab=2, cex.main=1.5, col="red", lwd=2, 
	xlim = c(min(which(sites_old[site,]>=0))+first_year_obs-3, last_year_obs+1), ylim=c(0,prev_max))
	lines(clinic.med1[site,]*100~seq(1970, END_YEAR), lwd=3)
	lines(clinic.lower1[site,]*100~seq(1970, END_YEAR), lty=2, lwd=2)
	lines(clinic.upper1[site,]*100~seq(1970, END_YEAR), lty=2, lwd=2)
	if (option.proj==1)	abline(v=last_year_obs-4.5, lty=2, lwd=2)
	lines(clinic.med2[site,]*100~seq(1975, END_YEAR), lwd=3, col="blue")
	lines(clinic.lower2[site,]*100~seq(1975, END_YEAR), lty=2, lwd=2, col="blue")
	lines(clinic.upper2[site,]*100~seq(1975, END_YEAR), lty=2, lwd=2, col="blue")
	# legend(last_year_obs-4, max(sites_old[site,],na.rm=T)*150, col=c("black", "blue", "red"), 
	# lwd=3, lty=c(1,1,3), text.col=c("black", "blue", "red"), cex=1.5, bty="o",
	# legend=c("EPP Model", "Flexible Model", "Clinic Data"))
}
dev.off()
}












