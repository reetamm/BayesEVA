rm(list = ls())
library(extRemes)

load('LA.RData')

tmp = aggregate(tmax,list(years),FUN=max)
y=tmp[,2]

fit1 = fevd(x=y,type = 'GEV')
plot(fit1)
fit2 = fevd(x=y,type = 'GEV',method = 'Bayesian',iter = 11000)
plot(fit2,burn.in=1000)
fit3 = fevd(x=x,data=tmp,threshold = 95,type = 'GP')
fit4 = fevd(x=x,data=tmp,threshold = 95,type = 'GP',method = 'Bayesian')

View(fevd)
