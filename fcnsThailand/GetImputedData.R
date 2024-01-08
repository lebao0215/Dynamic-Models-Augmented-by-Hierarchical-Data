#Create data objects
data.imputed <- read.csv(paste(dirData,country,"Impute", data.name, seed.k, ".csv" ,sep=""), header=TRUE, fill=TRUE, colClasses="character", na=NA)
data.site = data.imputed[which(data.imputed$State==site),-1]
data.site[1,] = as.numeric(data.site[1,]); data.site[2,] = as.numeric(data.site[2,])
new_start_year = min(which(data.site[1,]>0))+1985-1
new_end_year = max(which(data.site[1,]>0))+1985-1

index.year = (data_start_yr-new_start_year+1):(data_end_yr-new_start_year+1)
n.site = nrow(data$N_table)
data.site[2,] = as.numeric(data.site[2,]) / sum(as.numeric(data.site[2,]),na.rm=T)*n.impute
data.site[2,] = as.numeric(data.site[2,])
data.new = data
data.new$N_table = data.new$I_table = data.new$P_table = matrix(NA, n.site+1, sum(!is.na(data.site[1,])))
data.new$N_table[1:n.site,index.year] = data$N_table
data.new$N_table[1+n.site,] = as.numeric(data.site[2,(new_start_year-1985+1):(new_end_year-1985+1)])

data.new$P_table[1:n.site,index.year] = data$P_table
data.new$P_table[1+n.site,] = as.numeric(data.site[1,(new_start_year-1985+1):(new_end_year-1985+1)])

data.new$I_table = data.new$P_table*data.new$N_table

index.NA = which(data.new$N_table<=1)
if (length(index.NA)>0)   data.new$N_table[index.NA] = data.new$P_table[index.NA] = data.new$I_table[index.NA] = NA

if (length(index.NA)<ncol(data.new$N_table)){
  data = data.new
  data_start_yr <- new_start_year
  data_end_yr <- new_end_year
}
