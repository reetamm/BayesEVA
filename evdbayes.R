rm(list = ls())
library(evd)
library(evdbayes)
load('HCDN/HCDN_annual_max.RData')

nmiss = apply(Y,1,function(x)sum(is.na(x)))
which(nmiss[HUC02=='10L']<22)
y = Y[HUC02=='10L',][19,]
# hist(y)
# load('precip.Rda')
# y <- as.matrix(maxs[ , 5:ncol(maxs)])[1,]
# summary(y)
mle.fit = fgev(y)
mle.fit
mat <- diag(c(10000, 100, 0.25))
pr.strm <- prior.norm(mean = c(0,0,0), cov = mat, trendsd = 0)
n <- 15000 ; t0 <- mle.fit$estimate ; s <- mle.fit$std.err ; b <- 5000
# s[2] = log(s[2]); t0[2] = log(t0[2])
rn.post <- posterior(n, t0, pr.strm, "gev", data = y, burn = b,psd=s)
summary(rn.post)

par(mfrow=c(1,3))
plot(1:10001,rn.post[,1],'l')
plot(1:10001,rn.post[,2],'l')
plot(1:10001,rn.post[,3],'l')
par(mfrow=c(1,1))
