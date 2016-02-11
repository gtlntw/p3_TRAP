#power comparison between analytical and simulation
##load project related functions
source('function.r')

result <- read.csv("2014-11-03 answer for the discrepency btw old and new simulations.csv", header=T)

cex=1.2
lwd=3
affect <- c(0,1,1,1) #with two generations
p_dis=0.01 #baselinbe is 0.01, the other option is 0.1
alpha=0.05
n1 <- 1000

curve(sapply(x, function(x) lyapunov.conditional(affect=affect, n_family=n1/2, f=0.0102, mu=x, p_dis=0.05, delta=1, alpha=alpha, stat=F)[2]),1,2, lwd=lwd, ylab="power", main="f=0.01, n_family=500", ylim=c(0,1))
lines(seq(1, 2, length=11), tapply(result$p.value<alpha, result$r, mean), col=2, lty=2, lwd=lwd)
legend("topleft",c("Analytical", "Simulated"), lwd=lwd, lty=1:2, col=1:2)
