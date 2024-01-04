rm(list = ls())
library(methods)  # otherwise new() not being found - weird
library(nimble)
library(extRemes)

dgpd <- nimbleFunction(
run = function(x = double(0), mu = double(0), sigma = double(0), xi = double(0), 
               log = integer(0, default = 0)) {
    std <- (x - mu)/sigma
    r  <- log(1/sigma)
    id <- xi == 0
    good <- (std > 0 & ((1 + xi * std) > 0)) | is.na(std)
    if(id & good) 
        logProb <- r - std
    if(!id & good) 
        logProb <- r - (1/xi + 1) * log(1 + xi * std)
    if(!good) 
        logProb <- -Inf
    returnType(double(0))
    if(log)  return(logProb)
    if(!log) return(exp(logProb))
    })

rgpd <- nimbleFunction(
    run = function(n = integer(0), mu = double(0), sigma = double(0), xi = double(0)) {
        returnType(double(0))
        if(n != 1) print("rgpd only allows n = 1; using n = 1.")
        u <- runif(1)
        if (xi == 0) {
            x = mu - sigma * log(u)
        } else {
            x =  mu + sigma * (u^(-xi) - 1)/xi
        }
        return (x)
    })

assign('dgpd', dgpd, .GlobalEnv)
data(Fort)
FortMax <- aggregate(Prec ~ year, data = Fort, max)

gev_model <- fevd(FortMax$Prec, type = 'GEV')
gpd_model = fevd(x=Prec,data=Fort,threshold = 0.23,type = 'GP')
FortGPD = Fort$Prec[Fort$Prec>0.23]
code <- nimbleCode({
    # likelihood
    for(i in 1:n) {
        y[i] ~ dgpd(mu, sigma, xi)
    }
    
    # priors
    # mu ~ dflat()
    sigma ~ dunif(0, 100)
    xi ~ dunif(-1, 1)
})

inits <- as.list(gpd_model$results$par)
names(inits) <- c('sigma', 'xi')
constants = list(n = length(FortGPD), mu=0.23)

model <- nimbleModel(code, data = list(y = FortGPD), inits = inits,
                     constants = constants)
output <- nimbleMCMC(code = code, data = list(y = FortGPD),
                     constants = constants,
                     inits = inits, niter = 11000, nburnin = 1000,
                     nchains = 2, summary = T, setSeed = 303*c(1:2))
output$summary

coda_samples = post_convert(output$samples)
plot(coda_samples)
geweke.diag(coda_samples)
gelman.diag(coda_samples)
effectiveSize(coda_samples)
