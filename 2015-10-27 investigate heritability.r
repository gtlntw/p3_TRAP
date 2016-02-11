pe_var <- 1.447215
p_dis=0.1

#simulate the mean
mean(rbinom(10000, 1, plogis(rnorm(10000, qlogis(p_dis), pe_var))))
var(rbinom(10000, 1, plogis(rnorm(10000, qlogis(p_dis), pe_var))))


#analytical mean
integrate(function(x) plogis(x)*dnorm(x, qlogis(p_dis), pe_var), lower = -Inf, upper = Inf)
p <- integrate(function(x) plogis(x)*dnorm(x, qlogis(p_dis), pe_var), lower = -Inf, upper = Inf)$value
p*(1-p)

integrate(function(x) plogis(x)*(1-plogis(x))*dnorm(x, qlogis(p_dis), pe_var), lower = -Inf, upper = Inf)$value +
integrate(function(x) plogis(x)^2*dnorm(x, qlogis(p_dis), pe_var), lower = -Inf, upper = Inf)$value -
integrate(function(x) plogis(x)*dnorm(x, qlogis(p_dis), pe_var), lower = -Inf, upper = Inf)$value^2

######################################################
pe_var <- 0.5  #1
e_var <- 0.5 #1
beta0 <- -2.564219 #-3.417
p_dis=0.01

#simulate the mean and variance
mean(rbinom(10000, 1, plogis(rnorm(10000, beta0, pe_var + e_var))))
var(rbinom(10000, 1, plogis(rnorm(10000, beta0, pe_var + e_var))))


#analytical mean
optimize(function(y) (p_dis-integrate(function(x) plogis(x)*dnorm(x, y, pe_var + e_var), lower = -Inf, upper = Inf)$value)^2, c(-10,10))

integrate(function(x) plogis(x)*dnorm(x, beta0, pe_var + e_var), lower = -Inf, upper = Inf)$value

integrate(function(x) plogis(x)*dnorm(x, beta0, pe_var + e_var), lower = -Inf, upper = Inf)
p <- integrate(function(x) plogis(x)*dnorm(x, beta0, pe_var + e_var), lower = -Inf, upper = Inf)$value
p*(1-p)

integrate(function(x) plogis(x)*(1-plogis(x))*dnorm(x, beta0, pe_var + e_var), lower = -Inf, upper = Inf)$value +
integrate(function(x) plogis(x)^2*dnorm(x, beta0, pe_var + e_var), lower = -Inf, upper = Inf)$value -
integrate(function(x) plogis(x)*dnorm(x, beta0, pe_var + e_var), lower = -Inf, upper = Inf)$value^2

