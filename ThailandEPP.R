# run full model independently
setwd("/storage/work/lub14/Thailand")
country = "Thailand"
dirData <- 'Data/'; dirFcns <- 'fcnsThailand/'; dirCode <- 'codeThailand/'
data.name = HRG = "FSWb"
Area.name = c("Central", "North", "Northeast", "South")

site = 1    # focus on one particular state
option.truncate = 0
source(paste(dirFcns,'GetParameters.R',sep=''))  	# Load key projection parameters
source(paste(dirFcns,'GenObjects.R',sep=''))		# Load population arrays, treatment arrays, storage arrays for model integrator
# data.id = data.name
data.id = paste(data.name, "Aug" ,sep="")
source(paste(dirFcns,'GetData.R',sep=''))			# Load data

source(paste(dirFcns,'fnFlexIntegrator2011.R',sep=''))	# remove bias parameter and add turnover
source(paste(dirCode,'prior.R',sep=''))                 # remove bias parameter
source(paste(dirCode,'likelihood.R',sep=''))            # remove NPBS
xxx = c(1983.00, 22.00, 0.40, 0.43, 0.25, -0.64, -0.04)
print(likelihood.log(xxx, dataOut, data.name, sample_size_pop))
source(paste(dirCode,'IMIS.opt.R',sep=''))
library(mvtnorm)
set.seed(1000)
ptm <- proc.time()
result = IMIS.opt(1000, 3000, 120)
time.used = round((proc.time() - ptm)/60)
print(time.used)	# 16 min
print(round(apply(result$resample, 2, median),2))

save(result, time.used, file = paste('Result/', country, data.name, site, "_Aug", ".RData",sep=""))
