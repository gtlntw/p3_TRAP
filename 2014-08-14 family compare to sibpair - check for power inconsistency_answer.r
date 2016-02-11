#power comparison between analytical and simulation
##load project related functions
source('function.r')
#n=500
result <- read.csv("2014-08-07 family compare to sibpair - check for power inconsistency.csv", header=T)
result_0.01 <- subset(result, f==0.01)
result_0.2 <- subset(result, f==0.2)

cex=1.2
lwd=3
affect <- c(0,1,1,1) #with two generations
p_dis=0.01 #baselinbe is 0.01, the other option is 0.1
alpha=0.05
n1 <- 1000

curve(sapply(x, function(x) lyapunov(affect=affect, n_family=n1/2, f=0.01, mu=x, p_dis=0.01, delta=1, alpha=alpha, stat=F)[3]),1,2, lwd=lwd, ylab="power", main="f=0.01 n_fam=500")
lines(seq(1,2,by=0.1), tapply(result_0.01$p_value<alpha, result_0.01$r, mean), col=2, lty=2, lwd=lwd)

#n=2000
result <- read.csv("2014-08-12 family compare to sibpair - check for power inconsistency_n2000.csv", header=T)
result_0.01 <- subset(result, f==0.01)
result_0.2 <- subset(result, f==0.2)

cex=1.2
lwd=3
affect <- c(0,1,1,1) #with two generations
p_dis=0.01 #baselinbe is 0.01, the other option is 0.1
alpha=0.05
n1 <- 4000

curve(sapply(x, function(x) lyapunov(affect=affect, n_family=n1/2, f=0.01, mu=x, p_dis=0.01, delta=1, alpha=alpha, stat=F)[3]),1,2, lwd=lwd, ylab="power", main="f=0.01, n_family=2000")
lines(seq(1, 1.6, length=11), tapply(result_0.01$p_value<alpha, result_0.01$r, mean), col=2, lty=2, lwd=lwd)


#correct power comparison
##load project related functions
source('function.r')


affect <- c(0,1,1,1) #with two generations
affect_3 <- c(0,1,1,1,1) #with third generation
p_dis=0.01 #baselinbe is 0.01, the other option is 0.1
alpha=10^-6
fn <- lyapunov #can handle the family with third generation


#1. power comparison across family, sibpair and case-control
cex=1.2
lwd=3
dev.new(width=10, height=5)
par(mfrow=c(1,3), cex=cex, mar=c(3.1, 0.1, 1.6, 0.1), oma=c(0,3,0,1), lwd=2)
n1 <- 1000
sigma2 <- 0
curve(sapply(x, function(x) power.sas(mu=x, sigma2=sigma2, f=0.01, n_sb=n1, alpha=alpha)), 1, 4.5, col="black", ylim=c(0,1),xlab="",ylab="", lwd=lwd)
mtext("Power", 2, line=2, cex=cex, adj=0.5)
mtext("(a)", 3, cex=cex, line=0, adj=-0.3)
mtext("f=0.01", 3, cex=cex, line=0.2, adj=0.5)
curve(sapply(x, function(x) fn(affect=affect, n_family=n1/2, f=0.01, mu=x, p_dis=0.01, delta=1, alpha=alpha, stat=F)[3]), col="seagreen3", add=TRUE, lwd=lwd)
#curve(sapply(x, function(x) fn(affect=affect_3, n_family=400, f=0.01, mu=x, p_dis=0.01, delta=1, alpha=alpha, stat=F)[3]), col="blue", add=TRUE, lwd=lwd)
curve(sapply(x, function(x) power_comp.sas(mu=x, sigma2=sigma2, f=0.01, n_pair=n1, out="A", alpha=alpha)), col="darkorange", add=TRUE, lwd=lwd)
grid(NA, NULL)
par(yaxt="n")
curve(sapply(x, function(x) power.sas(mu=x, sigma2=sigma2, f=0.05, n_sb=n1, alpha=alpha)), 1, 4.5, col="black", ylim=c(0,1), lwd=lwd)
mtext("f=0.05", 3, cex=cex, line=0.2, adj=0.5)
curve(sapply(x, function(x) fn(affect=affect, n_family=n1/2, f=0.05, mu=x, p_dis=0.05, delta=1, alpha=alpha, stat=F)[3]), col="seagreen3", add=TRUE, lwd=lwd)
#curve(sapply(x, function(x) fn(affect=affect_3, n_family=400, f=0.05, mu=x, p_dis=0.01, delta=1, alpha=alpha, stat=F)[3]), col="blue", add=TRUE, lwd=lwd)
curve(sapply(x, function(x) power_comp.sas(mu=x, sigma2=sigma2, f=0.05, n_pair=n1, out="A", alpha=alpha)), col="darkorange", add=TRUE, lwd=lwd)
grid(NA, NULL)
mtext("Mean Relative Risk", 1, outer=FALSE, line=2, cex=cex)
curve(sapply(x, function(x) power.sas(mu=x, sigma2=sigma2, f=0.2, n_sb=n1, alpha=alpha)), 1, 4.5, col="black", ylim=c(0,1), lwd=lwd)
mtext("f=0.2", 3, cex=cex, line=0.2, adj=0.5)
curve(sapply(x, function(x) fn(affect=affect, n_family=n1/2, f=0.2, mu=x, p_dis=0.01, delta=1, alpha=alpha, stat=F)[3]), col="seagreen3", add=TRUE, lwd=lwd)
#curve(sapply(x, function(x) fn(affect=affect_3, n_family=400, f=0.2, mu=x, p_dis=0.01, delta=1, alpha=alpha, stat=F)[3]), col="blue", add=TRUE, lwd=lwd)
curve(sapply(x, function(x) power_comp.sas(mu=x, sigma2=sigma2, f=0.2, n_pair=n1, out="A", alpha=alpha)), col="darkorange", add=TRUE, lwd=lwd)
grid(NA, NULL)
legend("right", c("Family", "Sibpair", "Conventional CC"), lty=1, col=c("seagreen3", "black", "darkorange"), cex=1, lwd=lwd, bty="n")
#legend("right", c("Family", "Family w. 3rd", "Sibpair", "Conventional CC"), lty=1, col=c("seagreen3", "blue", "black", "darkorange"), cex=1, lwd=lwd, bty="n")
