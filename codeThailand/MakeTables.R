load(paste(model,country,region,".RData",sep=""))

result$resample = exp(result$resample)
resample.50 = rbind(resample.50, apply(result$resample, 2, median))
resample.025 = rbind(resample.025, apply(result$resample, 2, quantile, 0.025))
resample.975 = rbind(resample.975, apply(result$resample, 2, quantile, 0.975))
# if (country!="Liberia" | region!="Rural") 
resample.all = rbind(resample.all, result$resample)

xx = NULL
for (i in 1:dim(result$resample)[2])	xx = cbind(xx, result$resample[,i]- mean(result$resample[,i]))
if (country!="Liberia" | region!="Rural") resample.centered = rbind(resample.centered, xx)
print(paste(country, region))
print(round(cor(xx), 2))

