#power comparison between analytical and simulation
##load project related functions
source('function.r')
result <- read.csv("2014-08-07 family compare to sibpair - check for power inconsistency.csv", header=T)
result_0.01 <- subset(result, f==0.01)
result_0.2 <- subset(result, f==0.2)

cex=1.2
lwd=3
affect <- c(0,1,1,1) #with two generations
p_dis=0.01 #baselinbe is 0.01, the other option is 0.1
alpha=0.05
n1 <- 1000

curve(sapply(x, function(x) lyapunov(affect=affect, n_family=n1/2, f=0.01, mu=x, p_dis=0.01, delta=1, alpha=alpha, stat=F)[3]),1,2, lwd=lwd, ylab="power", main="f=0.01")
lines(seq(1,2,by=0.1), tapply(result_0.01$p_value<alpha, result_0.01$r, mean), col=2, lty=2)

#check the distribution of test statistics under the null
hist(subset(result_0.01$T, result_0.01$r==1), freq=F, xlab="Z", main="Dist. of test statistics underl the null")
curve(dnorm, add=T)
#check the distribution of test statistics under r=1.5
hist(subset(result_0.01$T, result_0.01$r==1.5), freq=F, xlab="Z", main="Dist. of test statistics underl r=1.5")
curve(dnorm(x, mean=2.123511), add=T)

#debugging and compare with the previous resuls
result_old <- read.csv("2014-08-05 family compare to sibpair - Lyapunov CLT power comparison.csv", header=T)
result_old_0.01 <- subset(result_old, f==0.01)
tapply(result_old_0.01$sd, result_old_0.01$r, mean)

result <- read.csv("2014-08-07 family compare to sibpair - check for power inconsistency.csv", header=T)
result_0.01 <- subset(result, f==0.01)
tapply(result_0.01$sd, result_0.01$r, mean)
