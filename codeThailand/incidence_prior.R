  load(paste("EPP/", country, region,".RData",sep=""))
	# prev.EPP = result$prev; 	inc.EPP = result$inc[,(year_min-beginYear+1):(endYear-beginYear+1)]

	beta.tri = 0.025;		Omega.tri = 150;		T.tri =450
	sample.trinomial <- function(theta){
		n = theta[1]
		prev = theta[2]
		inc = theta[3]
		p = rep(0,3)	# negative, non-recent, recent	
		p[1] = 1-prev
		p[3] = (1-prev)*inc*(Omega.tri-beta.tri*T.tri)/365 + beta.tri*prev
		p[2] = prev-p[3]
		n.tri = sample(1:3, n, replace = TRUE, prob = p)
		return(n.tri)
	}
	N.tri = matrix(0, dim(result$inc)[1], 3)
	for (i in 1:dim(result$inc)[1]){
		set.seed(1000+i)
		X = sample.trinomial(c(sample_size_TRI, result$prev[i,46], result$inc[i,46]))
		if (senario=="F")		X = sample.trinomial(c(sample_size_TRI, result$prev[i,46], result$inc[i,46]/2))
		if (senario=="G")		X = sample.trinomial(c(sample_size_TRI, result$prev[i,46], result$inc[i,46]*2))
		N.tri[i,1] = length(which(X==1))
		N.tri[i,2] = length(which(X==2))
		N.tri[i,3] = length(which(X==3))
	}
	n.tri = round(apply(N.tri, 2, mean))
	if (senario=="B")		beta.tri = 0.015
	if (senario=="C")		beta.tri = 0.035	
	if (senario=="D")		Omega.tri = 130
	if (senario=="E")		Omega.tri = 170

	likelihood.tri <- function(theta){
		l = -10000
		prev.tri = theta[1]
		inc.tri = theta[2]
		p1 = 1-prev.tri
		p3 = (1-prev.tri)*inc.tri*(Omega.tri-beta.tri*T.tri)/365 + beta.tri*prev.tri
		p2 = prev.tri-p3
		if (p1>0 & p2>0 & p3>0 & p1<1 & p2<1 & p3<1)
		l = n.tri[1]*log(p1) + n.tri[2]*log(p2) + n.tri[3]*log(p3)
		return(-l)
	}
	theta1.hat = 1-n.tri[1]/sum(n.tri)
	theta2.hat = (n.tri[3]-beta.tri*(n.tri[2]+n.tri[3]))/n.tri[1]/(Omega.tri-beta.tri*T.tri)*450
	theta.ini = c(theta1.hat, theta2.hat)
	likelihood.tri(theta.ini)
	optimizer = optim(theta.ini, likelihood.tri, method = "BFGS", hessian = TRUE,
	control = list(parscale = rep(0.001,2), maxit = 500))
	theta.opt = 100*optimizer$par
	theta.cov = solve(optimizer$hessian)
	theta.cor = theta.cov[1,2]/sqrt(theta.cov[1,1]*theta.cov[2,2])
	theta.sd = 100*sqrt(diag(theta.cov))
