xx = result$resample
xx[,5] = result$resample[,5]/result$resample[,6]
xx[,4] = result$resample[,4]/result$resample[,6]
xx[,3] = result$resample[,3] - result$resample[,4]*result$resample[,5]/result$resample[,6]
index = which(abs(xx[,3])<quantile(abs(xx[,3]),0.95) & 
abs(xx[,4])<quantile(abs(xx[,4]),0.95) & abs(xx[,5])<quantile(abs(xx[,5]),0.95))
pairs(xx[index,])

xx = result$resample
xx[,4] = result$resample[,4]+result$resample[,6]
xx[,3] = result$resample[,3]+result$resample[,5]
pairs(xx)


range(result$resample[,3])
range(result$resample[,4])
range(result$resample[,5])
range(result$resample[,6])
pairs(result$resample)
