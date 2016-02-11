#power comparison between analytical and simulation
##read in f=0.01, n_family=500 with two generations
source('function.r')
result <- read.csv("2014-08-07 family compare to sibpair - check for power inconsistency.csv", header=T)
result_0.01 <- subset(result, f==0.01)
##read in f=0.01, n_family=400 with three generations
result_3rd <- read.csv("2014-08-19 power of family with three generation.csv", header=T)
result_3rd_0.01 <- subset(result_3rd, f==0.01)


cex=1.2
lwd=3
affect <- c(0,1,1,1) #with two generations
p_dis=0.01 #baselinbe is 0.01, the other option is 0.1
alpha=0.05
n1 <- 1000 #no. of sibpairs

plot(unique(result_0.01$r), tapply(result_0.01$p_value<alpha, result_0.01$r, mean), col=1, lty=1, lwd=lwd, type="l", xlab="mean relative risk", ylab="power")
lines(unique(result_3rd_0.01$r), tapply(result_3rd_0.01$p_value<alpha, result_3rd_0.01$r, mean), col=2, lty=1,lwd=lwd)
legend("bottomright", c("two generations", "three generations"), col=c(1,2), lty=1)
