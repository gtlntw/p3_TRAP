##load project related functions
source('function.r')

##Figure 3.a
#power curve var=0
cex=1.2
lwd=3
dev.new(width=9, height=5)
par(mfrow=c(1,3), cex=cex, mar=c(3.1, 0.1, 1.6, 0.1), oma=c(0,3,0,1), lwd=2)
n1 <- 1000
sigma2 <- 0
curve(sapply(x, function(x) power.sas(mu=x, sigma2=sigma2, f=0.01, n_sb=n1)), 1, 4.5, col="black", ylim=c(0,1),xlab="",ylab="", lwd=lwd)
mtext("Power", 2, line=2, cex=cex, adj=0.5)
mtext("(a)", 3, cex=cex, line=0, adj=-0.3)
mtext("f=0.01", 3, cex=cex, line=0.2, adj=0.5)
curve(sapply(x, function(x) power.family(mu=x, f=0.01, n_family=n1/2)), col="seagreen3", add=TRUE, lwd=lwd)
curve(sapply(x, function(x) power_comp.sas(mu=x, sigma2=sigma2, f=0.01, n_pair=n1, out="A")), col="darkorange", add=TRUE, lwd=lwd)
grid(NA, NULL)
par(yaxt="n")
curve(sapply(x, function(x) power.sas(mu=x, sigma2=sigma2, f=0.05, n_sb=n1)), 1, 4.5, col="black", lwd=lwd)
mtext("f=0.05", 3, cex=cex, line=0.2, adj=0.5)
curve(sapply(x, function(x) power.family(mu=x, f=0.05, n_family=n1/2)), col="seagreen3", add=TRUE, lwd=lwd)
curve(sapply(x, function(x) power_comp.sas(mu=x, sigma2=sigma2, f=0.05, n_pair=n1, out="A")), col="darkorange", add=TRUE, lwd=lwd)
grid(NA, NULL)
mtext("Mean Relative Risk", 1, outer=FALSE, line=2, cex=cex)
curve(sapply(x, function(x) power.sas(mu=x, sigma2=sigma2, f=0.2, n_sb=n1)), 1, 4.5, col="black", lwd=lwd)
mtext("f=0.2", 3, cex=cex, line=0.2, adj=0.5)
curve(sapply(x, function(x) power.family(mu=x, f=0.2, n_family=n1/2)), col="seagreen3", add=TRUE, lwd=lwd)
curve(sapply(x, function(x) power_comp.sas(mu=x, sigma2=sigma2, f=0.2, n_pair=n1, out="A")), col="darkorange", add=TRUE, lwd=lwd)
grid(NA, NULL)
legend("right", c("Family", "Sibpair", "Conventional CC"), lty=1, col=c("seagreen3", "black", "darkorange"), cex=1, lwd=lwd, bty="n")
