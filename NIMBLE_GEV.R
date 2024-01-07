rm(list = ls())
library(methods)  # otherwise new() not being found - weird
library(nimble)
library(extRemes)
library(coda)
library(postpack)

dgev <- nimbleFunction(
    run = function(x = double(0), mu = double(0), sigma = double(0), xi = double(0), 
                   log = integer(0, default = 0)) {
        std <- (x - mu) / sigma
        if(1 + xi*std <= 0)
            logProb <- -Inf
        if(1 + xi*std > 0)
            logProb <- -log(sigma) - (1+1/xi) * log(1+xi*std) - (1+xi*std)^(-1/xi)
        returnType(double(0))
        if(log)  return(logProb)
        if(!log) return(exp(logProb))
    })

rgev <- nimbleFunction(
    run = function(n = integer(0), mu = double(0), sigma = double(0), xi = double(0)) {
        returnType(double(0))
        if(n != 1) print("rgev only allows n = 1; using n = 1.")
        u <- runif(1)
        if(xi == 0)
            x = -log(-log(u))
        if(xi != 0)
            x = ((-log(u))^(-xi)-1)/xi
        return (x*sigma + mu)
    })

data(Fort)
FortMax <- aggregate(Prec ~ year, data = Fort, max)
gev_model <- fevd(FortMax$Prec, type = 'GEV')
# gev_model$results$par[2] = log(gev_model$results$par[2])

code <- nimbleCode({
    # likelihood
    for(i in 1:n) {
        y[i] ~ dgev(mu, sig, xi)
    }
    # priors
    mu ~ dnorm(0, 10)
    log(sig) ~ dnorm(-1, 1)
    xi ~ dnorm(0, .5)
})

inits <- as.list(gev_model$results$par)
names(inits) <- c('mu', 'log_sig', 'xi')

output <- nimbleMCMC(code = code, data = list(y = FortMax$Prec),
                     constants = list(n = length(FortMax$Prec)),
                     inits = inits, niter = 11000, nburnin = 1000,
                     nchains = 2, summary = T, setSeed = 303*c(1:2))
output$summary
output$samples$chain1[,1] = exp(output$samples$chain1[,1])
output$samples$chain2[,1] = exp(output$samples$chain2[,1])
all.chains = do.call(rbind,output$samples)
dim(all.chains)
# all.chains[,1] = exp(all.chains[,1])
apply(all.chains,2,mean)
apply(all.chains,2,sd)
apply(all.chains,2,quantile,c(0.025,0.975))
# par(mfrow=c(1,3))
# plot(1:10000,output$samples[,1],'l')
# plot(1:10000,output$samples[,2],'l')
# plot(1:10000,output$samples[,3],'l')
# par(mfrow=c(1,1))


coda_samples = post_convert(output$samples)
plot(coda_samples)
geweke.diag(coda_samples)
gelman.diag(coda_samples)
effectiveSize(coda_samples)

par(mfrow=c(1,3))
plot(1:10000,all.chains[1:10000,2],'l',col='blue',xlab = 'Iteration',ylab = expression(mu))
lines(1:10000,all.chains[10001:20000,2],col='red')
plot(1:10000,all.chains[1:10000,1],'l',col='blue',xlab = 'Iteration',ylab = expression(sigma))
lines(1:10000,all.chains[10001:20000,1],col='red')
plot(1:10000,all.chains[1:10000,3],'l',col='blue',xlab = 'Iteration',ylab = expression(xi))
lines(1:10000,all.chains[10001:20000,3],col='red')
par(mfrow=c(1,1))
