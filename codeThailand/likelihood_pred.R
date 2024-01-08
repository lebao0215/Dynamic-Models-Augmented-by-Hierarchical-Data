########################################################################################
# Define the probit_unequalvar likelihood
########################################################################################
## 4 lines of changes marked by ## new ##, Le Bao, 09/14/2012

# Define the parameters in the prior for sigma02
alpha = 0.58;	beta = 93
na_sites = dataOut$data$P_table;	isna_sites = is.na(na_sites)

if (option.truncate==1){
  data.year = which(apply(!is.na(na_sites), 2, sum)>0)
  length.truncate = ceiling(length(data.year)/3)
  index.truncate = (length.truncate+1):(length(data.year)-length.truncate)
  na_sites[,index.truncate] = NA;  isna_sites = is.na(na_sites)
  region.data = apply(!is.na(na_sites), 1, sum)
  na_sites = na_sites[region.data>0,];  isna_sites = is.na(na_sites)
  dataOut$data$N_table = dataOut$data$N_table[region.data>0,]
  dataOut$data$I_table = dataOut$data$I_table[region.data>0,]
  dataOut$data$P_table = dataOut$data$P_table[region.data>0,]
}

nr_sites = dim(na_sites)[1];		nr_years = dim(na_sites)[2]

r0.low = 1/11.5;	r0.high = 5			# lower and upper bounds for r0
y0.low = 0.1^13;	y0.high = 0.0025		# lower and upper bounds for y0
t0.low = 1969.5;	t0.high = 1990.5		## new ##

likelihood.log <- function(x){
Like = -10^10

# Calculate the variance of estimated prevalance for each ANC data
Nst = dataOut$data$N_table		# The sample size for each clinic over time is a matrix, each row is one clinic
sample_size_pop = spectrumObjects$p15to49[(data_start_yr-min(popData[,1])+1):(data_end_yr-min(popData[,1])+1)]
xst = (dataOut$data$P_table*dataOut$data$N_table+0.5) / (dataOut$data$N_table+1)
Wst = qnorm(xst)
var_eps = 2*pi*xst*(1-xst)/dataOut$data$N_table*exp(Wst^2)
# Calculate the empirical variance of clinic effects, var_bi, which is used as the upper bound of integration
if (nr_sites>1){
	Wst.mean = apply(Wst,2,mean,na.rm=T)
	bi = NULL
	for (i in 1:nr_sites)	bi = c(bi, mean(Wst[i,]-Wst.mean, na.rm=T)^2)
	var_bi = mean(bi, na.rm=T)
	site_index = which(!is.na(bi))
}

if (x[1]>=t0.low & x[1]<=t0.high & x[2]>=1 & x[7]<=0){	## new ##

	result = fnRMS2011(x,turnoverOff=TRUE)		## new ##
	result$annualFittingPrevs[is.na(result$annualFittingPrevs)] = 0		## new ##
	# Compare_prev has the model output for the years of which we have observed prevalence
	rho = result$annualFittingPrevs[(data_start_yr-beginYear+1):(data_end_yr-beginYear+1)]
	# tr_diff: difference between transformed prev and transformed obs
	tr_diff = NULL
	for (i in 1:nr_sites){
		if (length(NPBS)>1)	tr_diff = rbind(tr_diff, Wst[i,] + x[8] - qnorm( (rho*sample_size_pop+0.5)/(sample_size_pop+1) ))
		if (length(NPBS)==1)	tr_diff = rbind(tr_diff, Wst[i,] - qnorm( (rho*sample_size_pop+0.5)/(sample_size_pop+1) ))
	}


	# the change of log r(t) before the data period
	r.pre.data = log(result$annualR[6:(data_start_yr-beginYear+3)])
	r.change = r.pre.data[-1] - r.pre.data[-length(r.pre.data)]

	check.output = max(result$annualFittingPrevs[1:(1980-beginYear+1)], na.rm=T)<0.1 &		# constraint on prevalence before 1980
				result$annualFittingPrevs[data_end_yr-beginYear+1]>1e-4 & 			# to speed up calculation
				result$annualFittingPrevs[data_end_yr-beginYear+1]<0.99 &			# to speed up calculation
				max(r.change, na.rm=T)<=0								# r(t) can not increase in pre-data period

	integrand = function(sigma02){
		pdf_site = rep(NA,nr_sites); result = rep(NA,length(sigma02))
		for (i in 1:length(sigma02)){
			for (site in site_index){
				nr_years_site = sum(!isna_sites[site,])
				Sigma = sigma02[i]*matrix(1,nr_years_site, nr_years_site) + 
					diag(x = var_eps[site,!isna_sites[site,]], ncol = nr_years_site, nrow = nr_years_site)
				pdf_site[site] = -0.5*log(det(Sigma)) - 0.5*nr_years_site*log(2*pi) - 
				0.5*t(tr_diff[site,!isna_sites[site,]])%*%solve(Sigma)%*%tr_diff[site,!isna_sites[site,]] 
			}
			result[i] = exp( sum(pdf_site[site_index]) + dgamma(1/sigma02[i], alpha, 1/beta, log=T) )
		}
		return(result)
	}

	# calculate the log likelihood when there is multiple clinics
	if (nr_sites>=2 & check.output)
		Like = log(exp(-600)+integrate(integrand, lower = 0.000000000000001, 
			upper = var_bi, subdivisions = 1000, stop.on.error = FALSE)$value)

	# calculate the log likelihood when there is a single clinic
	if (nr_sites==1 & check.output)
		Like = sum(dnorm(tr_diff,0,sqrt(var_eps),log=T),na.rm=T)

	# to solve numeric issues when the likelihood is too small, inflate the log likelihood by the number of observations
	# Like = Like + sum(!isna_sites)*2

	# Now include NPBS information
	if (length(NPBS)>1){
		if (option.truncate == 1 & length(NPBS$rate)>1)		NPBS$rate[-1] = NA
		NPBS_diff = qnorm(NPBS$rate)-qnorm(result$annualFittingPrevs[NPBS$year-beginYear+1])
		Like = Like + sum(dnorm(NPBS_diff,0,sqrt(var_NPBS),log=T),na.rm=T)
	}
}
return(Like)
}


