library(brms)
library(bayesplot)
library(splines)

country = "Thailand"
dirData <- 'Data/'; dirFcns <- 'fcnsThailand/'; dirCode <- 'codeThailand/'
data.name = HRG = "FSWb"
name.main = "Indirect Sex Workers"
label.name = c("(a)", "(b)", "(c)", "(d)", "(e)", "(f)", "(g)", "(h)")
Area.name = c("Central", "North", "Northeast", "South")

# load data
Data <- read.csv(paste(dirData, country, data.name, ".csv" ,sep=""), header=TRUE, fill=TRUE, colClasses="character")
# make variables numeric
colnames(data)
data <- as.data.frame((apply(Data,2,as.numeric)))
# check which Years have data
data.filled.Year = which(apply(data, 2, sum, na.rm=T)>0)
data.filled.Year = 1983 + data.filled.Year[-1] 
data.filled.Year
# data years
data.Year <- 1985:2020

# create data frame to use
createDF <- function(inData, year.Range) {
  HIV.p = HIV.n = Year = t1 = Site = Area = NULL
  for (i in 1:(nrow(inData)/2)){
    # postive cases at site i
    HIV.i = round(inData[2*i-1,-1]/100 * inData[2*i,-1])
    HIV.p = c(HIV.p, as.numeric(HIV.i))
    HIV.n = c(HIV.n, as.numeric(inData[2*i,-1]))
    t1 = c(t1, year.Range-1990)
    Year = c(Year, year.Range) 
    Site = c(Site, rep(i, length(year.Range)))
    Area = c(Area, rep(inData[2*i-1,1], length(year.Range)))
  }
  Site = as.factor(Site);    Area = as.factor(Area); 
  outData=data.frame(HIV.p,HIV.n,t1,Site,Area,Year)  
  return(outData)
}
iData = createDF(data, data.Year)
iData.complete = iData[!is.na(iData$HIV.p),]
Num.Year = length(unique(iData.complete$t1)) # 9 unique years
t1.min = min(iData.complete$t1)
t1.max = max(iData.complete$t1)
t1.n = t1.max-t1.min+1
  
prior.ns <- prior(normal(0, 0.1), class = "sd", group = "Area")
time_taken <- system.time(
  fit.full <- brm(HIV.p | trials(HIV.n) ~  s(t1, bs="ps",k=4) + (ns(t1, df=3)|Area) + (1|Site), data = iData.complete,
                  family = binomial("probit"), prior = prior.ns,
                  control = list(adapt_delta = 0.98, max_treedepth=20))
    )
Time.full = time_taken["elapsed"]
save(fit.full, Time.full, file = paste("Result/", country, data.name, "Full.RData" ,sep=""))

for (Area.i in unique(iData$Area)){
  source("KL.R")
  save(K.all, KL.all, Time.all, Posterior_mean.all,
       Posterior_mean, Posterior_cov, Posterior_samples,
       file = paste('Result/Area', country, data.name, Area.i, ".RData",sep=""))
}

# load data
Data <- read.csv(paste(dirData, country, data.name, ".csv" ,sep=""), header=TRUE, fill=TRUE, colClasses="character")
# make variables numeric
data <- as.data.frame((apply(Data,2,as.numeric)))
colnames(data) = NULL

Area.KL.all = KL.min = Time.CV = NULL
K.all = c(25, 50, 100, 200, 400, 800, 1600, 3200)
for (Area.i in 1:4){
  load(paste('Result/Area', country, data.name, Area.i, ".RData",sep=""))
  Area.KL.all = rbind(Area.KL.all, KL.all)
  KL.min = c(KL.min, which(KL.all==min(KL.all)))
  print(Time.all)
  Time.CV = rbind(Time.CV, as.numeric(Time.all))
  
  data.i = matrix(NA, 2, ncol(data))
  data.i[1,1] = Area.i
  n.impute = K.all[which(KL.all==min(KL.all))]
  #cbind(seq(t1.min,t1.max,by=1)+1990, y.impute, n.impute)
  data.i[1, (t1.min+1990-1984+1):(t1.max+1990-1984+1)] = Posterior_mean.all[which(KL.all==min(KL.all)),]*100
  data.i[2, (t1.min+1990-1984+1):(t1.max+1990-1984+1)] = n.impute
  data = rbind(data, data.i[1,], data.i[2,])
}
print(cbind(1:4, KL.min))
round(Area.KL.all,2)

round(sum(Time.CV)/60/4,2)

data[is.na(data)] = ""
write.csv(rbind(c("State", 1985:2020), data), file=paste(dirData, country, data.name, "Aug.csv" ,sep=""), row.names=FALSE)
