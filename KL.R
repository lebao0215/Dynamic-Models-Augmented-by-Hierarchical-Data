iData.Area = iData.complete[which(iData.complete$Area==Area.i),]
new_data = data.frame(t1=seq(t1.min, t1.max, by=1), 
                      Area=rep(as.numeric(Area.i),t1.n), 
                      Site=rep(as.numeric(iData.Area$Site[1]),t1.n),
                      HIV.n = rep(300, t1.n))
Posterior_samples = posterior_predict(fit.full, newdata = new_data)
Posterior_mean = apply(Posterior_samples/300, 2, mean)
Posterior_cov = cov(Posterior_samples/300)
Posterior_mean.all = Posterior_mean

K.all = c(25, 50, 100, 200, 400, 800, 1600, 3200)
KL.all = Time.all.G = Time.all = NULL

for (K in K.all){
  iData.ghost = cbind(round(Posterior_mean*K), rep(K, length(Posterior_mean)), t1.min:t1.max, 
                      rep(max(as.numeric(iData.Area$Site))+1,length(Posterior_mean)), 
                      rep(Area.i,length(Posterior_mean)), 1990+(t1.min:t1.max))
  colnames(iData.ghost) = colnames(iData.Area)
  iData.Area.new = rbind(iData.Area, iData.ghost)
  iData.Area.new$HIV.n = as.numeric(iData.Area.new$HIV.n)
  iData.Area.new$HIV.p = as.numeric(iData.Area.new$HIV.p)
  iData.Area.new$t1 = as.numeric(iData.Area.new$t1)
  iData.Area.new$Site = as.factor(iData.Area.new$Site)
  time_taken <- system.time(
    fit.ghost <- brm(HIV.p | trials(HIV.n) ~  s(t1, bs="ps",k=4) + (1|Site), data = iData.Area.new,
                family = binomial("probit"), control = list(adapt_delta = 0.9, max_treedepth=10))
  )
  Time.all = c(Time.all, time_taken["elapsed"])
  new_data = data.frame(t1=seq(t1.min, t1.max, by=1), 
                        Area=rep(Area.i,t1.n), Site=rep(as.numeric(iData.Area.new$Site[1]),t1.n),
                        HIV.n = rep(300, t1.n))
  
  Posterior_ghost = posterior_predict(fit.ghost, newdata = new_data)
  Posterior_mean.K = apply(Posterior_ghost/300, 2, mean)
  Posterior_cov.K = cov(Posterior_ghost/300)
  Posterior_mean.all = rbind(Posterior_mean.all, Posterior_mean.K)
  
  KL = sum(diag(solve(Posterior_cov.K)%*%Posterior_cov)) - length(Posterior_mean) + 
    t(Posterior_mean.K-Posterior_mean)%*%solve(Posterior_cov.K)%*%(Posterior_mean.K-Posterior_mean) +
    sum(log(eigen(Posterior_cov.K)$values)) - sum(log(eigen(Posterior_cov)$values))
  KL.all = c(KL.all, KL)
  if (K==25){
    fit.K = fit.ghost
    KL.min = KL
  }
  if (K>25 & KL<KL.min){
    fit.K = fit.ghost
    KL.min = KL
  }
}

print(KL.all)
print(Time.all)
time_taken <- system.time(
  fit.ind <- brm(HIV.p | trials(HIV.n) ~  s(t1, bs="ps",k=4) + (1|Site), data = iData.Area,
                   family = binomial("probit"), control = list(adapt_delta = 0.9, max_treedepth=10))
)
Time.ind = time_taken["elapsed"]

