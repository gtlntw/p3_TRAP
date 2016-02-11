source("function.r")
result <- read.csv("2014-05-12 family compare to sibpair - simulation result.csv", header=T)

##calculate the simulated false positive rate
alpha=0.05
sim_result <- with(result, 
{
  tapply(true<alpha, list(f, r), mean)
}
     )
x <- seq(1,4, by=0.2)

##draw the analytical result first
affect <- c(0,1,1,1) #with two generations
p_dis=0.01 #baselinbe is 0.01, the other option is 0.1

fn <- power.family.5 #can handle the family with third generation

cex=1.2
lwd=3
dev.new(width=10, height=5)
par(mfrow=c(1,3), cex=cex, mar=c(3.1, 0.1, 1.6, 0.1), oma=c(0,3,0,1), lwd=2)
n1 <- 1000
sigma2 <- 0
curve(sapply(x, function(x) power.sas(mu=x, sigma2=sigma2, f=0.01, n_sb=n1, alpha=alpha)), 1, 3, col="black", ylim=c(0,1),xlab="",ylab="", lwd=lwd)
mtext("Power", 2, line=2, cex=cex, adj=0.5)
mtext("(a)", 3, cex=cex, line=0, adj=-0.3)
mtext("f=0.01", 3, cex=cex, line=0.2, adj=0.5)
curve(sapply(x, function(x) fn(affect=affect, p_dis=p_dis, mu=x, f=0.01, n_family=n1/2, alpha=alpha)), col="seagreen3", add=TRUE, lwd=lwd)
lines(x, sim_result[1,], col="seagreen3", lwd=lwd, lty=2)
grid(NA, NULL)
par(yaxt="n")
curve(sapply(x, function(x) power.sas(mu=x, sigma2=sigma2, f=0.05, n_sb=n1, alpha=alpha)), 1, 3, col="black", lwd=lwd)
mtext("f=0.05", 3, cex=cex, line=0.2, adj=0.5)
curve(sapply(x, function(x) fn(affect=affect, p_dis=p_dis, mu=x, f=0.05, n_family=n1/2, alpha=alpha)), col="seagreen3", add=TRUE, lwd=lwd)
lines(x, sim_result[2,], col="seagreen3", lwd=lwd, lty=2)
grid(NA, NULL)
mtext("Mean Relative Risk", 1, outer=FALSE, line=2, cex=cex)
curve(sapply(x, function(x) power.sas(mu=x, sigma2=sigma2, f=0.2, n_sb=n1, alpha=alpha)), 1, 3, col="black", lwd=lwd)
mtext("f=0.2", 3, cex=cex, line=0.2, adj=0.5)
curve(sapply(x, function(x) fn(affect=affect, p_dis=p_dis, mu=x, f=0.2, n_family=n1/2, alpha=alpha)), col="seagreen3", add=TRUE, lwd=lwd)
lines(x, sim_result[3,], col="seagreen3", lwd=lwd, lty=2)
grid(NA, NULL)
legend("right", c("Family", "Sibpair", "Family Simulated"), lty=c(1,1,2), col=c("seagreen3", "black", "seagreen3"), cex=1, lwd=lwd, bty="n")
