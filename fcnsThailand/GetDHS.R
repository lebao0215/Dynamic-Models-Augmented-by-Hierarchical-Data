#Function to load and reorganize ANC surveillance data from csv file
fnGetNPBS <- function(country,site){ 

 #Reorganize data
 data <- read.csv(paste(dirData,country,"NPBS.csv",sep=""), header=TRUE, fill=TRUE, colClasses="character")
 year = as.numeric(substr(colnames(data)[-1],start=2,stop=5))[1:2]
 
 rate = as.numeric(data[site,2:3])
 size = as.numeric(data[site,4:5])

 out <- list(year, rate/100, size)
 names(out) <- c("year", "rate", "size")
 return(out)
}

NPBS <- fnGetNPBS(country,site)
var_NPBS = 2*pi*NPBS$rate*(1-NPBS$rate)/NPBS$size*exp(qnorm(NPBS$rate)^2)
