rm(list = ls())
library(evd)
library(evmix)
library(evdbayes)
load('LA.RData')

tmp = aggregate(tmax,list(years),FUN=max)
y=tmp[,2]
# hist(y)
# load('precip.Rda')
# y <- as.matrix(maxs[ , 5:ncol(maxs)])[1,]
# summary(y)
mle.fit = fgev(y)
mle.fit
mat <- diag(c(10000, 100, 0.25))
pr.strm <- prior.norm(mean = c(0,0,0), cov = mat, trendsd = 0)
n <- 11000 ; t0 <- mle.fit$estimate ; s <- mle.fit$std.err ; b <- 1000
# s[2] = log(s[2]); t0[2] = log(t0[2])
rn.post <- posterior(n, t0, pr.strm, "gev", data = y, burn = b,psd=s)
summary(rn.post)

par(mfrow=c(1,3))
plot(1:10001,rn.post[,1],'l')
plot(1:10001,rn.post[,2],'l')
plot(1:10001,rn.post[,3],'l')
par(mfrow=c(1,1))


thresh = aggregate(tmax,list(years),FUN=quantile,0.95)
summary(thresh[,2])
y = tmax[tmax>95]
length(y)
mle.fit = fgpd(y,u=9)
mle.fit$mle
mle.fit$se

t0 <- c(95,mle.fit$mle) ; s <- c(.33,mle.fit$se)
rn.post <- posterior(n, t0, pr.strm, "gpd", data = y, burn = b,psd=s)

summary(rn.post)

par(mfrow=c(1,3))
plot(1:10001,rn.post[,1],'l')
plot(1:10001,rn.post[,2],'l')
plot(1:10001,rn.post[,3],'l')
par(mfrow=c(1,1))
