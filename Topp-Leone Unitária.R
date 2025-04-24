#log-likelihood of the unit gamma
UTP<-expression(
  log(2) + log(mu)+log(1-y) + mu*log(y) - log(y) + mu * log(2-y) - log(2-y)
)
m1UTP<-D(UTP,"mu")
print(m1UTP)

UTP<-function (mu.link = "logit")
{
  mstats <- checklink("mu.link", "UTP", substitute(mu.link),
                      c("logit", "probit", "cloglog", "cauchit", "log", "own"))

  structure(list(family = c("UTP", "Unit-Topp-Leone"),
                 parameters = list(mu = TRUE),
                 nopar = 1,
                 type = "Continuous",
                 mu.link = as.character(substitute(mu.link)),
                 mu.linkfun = mstats$linkfun,
                 mu.linkinv = mstats$linkinv,
                 mu.dr = mstats$mu.eta,
                 dldm = function(y, mu) {
                   dldm <- eval(m1UTP)
                   dldm
                 },
                 d2ldm2 = function(y,mu) {
                   dldm <- eval(m1UTP)
                   d2ldm2 <- -dldm * dldm
                   d2ldm2 <- ifelse(d2ldm2 < -1e-15, d2ldm2,-1e-15)
                   d2ldm2
                 },
                 G.dev.incr = function(y, mu, w, ...) -2 * log(dUTP(y=y, mu=mu)),
                 rqres = expression(
                   rqres(pfun = "pUTP", type = "Continuous", y = y, mu = mu)
                 ),
                 mu.initial = expression(mu <- rep(mean(y),length(y))),
                 mu.valid = function(mu) all(mu > 0 & mu < 1),
                 y.valid = function(y) all(y > 0 & y < 1)
  ),
  class = c("gamlss.family", "family"))
}




#------------------------------------------------------------------------------------------
# density function
dUTP <- function(y, mu = 0.7) {
  if (any(mu <= 0) | any(mu > 1)) stop(paste("mu must be between 0 and 1", "\n", ""))
  if (any( y < 0) | any(y > 1)) stop(paste("x must be between 0 and 1", "\n", ""))
  
  fy1 <- 2 * mu * (1 - y) * (y * (2 - y))^(mu - 1)
  fy1
}

# integrate(dUTP,0,1) # checking the pdf
#------------------------------------------------------------------------------------------
# cumulative distribution function
pUTP<-function(q, mu = 0.7){
  if (any(mu <= 0) | any(mu >= 1)) stop(paste("mu must be between 0 and 1", "\n", ""))
  if (any(q <= 0) | any(q >= 1)) stop(paste("x must be between 0 and 1", "\n", ""))
  cdf1 <- (q * (2 - q))^mu
  cdf1
}
# pUTP(.5)
# integrate(dUTP,0,.5) # checking the cdf with the pdf
#------------------------------------------------------------------------------------------
# quantile function
qUTP<-function(u,mu)
{
  q<-(1 - sqrt(1 - u^(1 / mu)))
  q
}
# u=pUTP(.5)
# qUTP(u,mu=.7) # checking the qf with the cdf
#------------------------------------------------------------------------------------------
# inversion method for random generation
rUTP<-function(n,mu)
{
  u<- runif(n)
  y<- qUTP(u,mu=mu)
  y
}

# Checking the results
library(gamlss)




set.seed(10)
n<-1000
# Case 1: without regressors
mu_true<-.523431

mu_result<-c()
for (i in 1:1000) {
  y<-rUTP(n,mu_true)
  fit1<-gamlss(y~1, family="UTP", trace = F)
  logit_link<-make.link("logit")
  mu_result[i]<-logit_link$linkinv(fit1$mu.coefficients)

  }
result1<- matrix(c(mu_true, mean(mu_result)), 2,1)
colnames(result1)<-c("mu")
rownames(result1)<-c("true value","mean")
print(round(result1,3))


# Case 2: with regressors

n <- 1000
beta0 <- 0.8
beta1 <- 0.3

beta0_result <- beta1_result <- numeric(5000)

for (i in 1:10000){
  x <- runif(n) 
  eta <- beta0 + beta1 * x
  mu <- exp(eta) / (1 + exp(eta))
  
  y <- rUTP(n, mu)
  
  fit <- gamlss(y ~ x, family = "UTP", trace = FALSE)
  beta0_result[i] <- fit$mu.coefficients[1]
  beta1_result[i] <- fit$mu.coefficients[2]
}

result2 <- matrix(
  c(beta0, mean(beta0_result),
    beta1, mean(beta1_result)),
  2, 2, byrow = TRUE
)
colnames(result2) <- c("true", "mean_estimate")
rownames(result2) <- c("beta0", "beta1")

print(round(result2, 3))
