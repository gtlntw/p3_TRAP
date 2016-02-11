#correct power comparison
##load project related functions
source('../function.r')


affect <- c(0,1,1,1) #with two generations
affect_3 <- c(0,1,1,1,1) #with third generation
p_dis=0.01 #baselinbe is 0.01, the other option is 0.1
alpha=2.5*10^-6
fn <- lyapunov #can handle the family with third generation


#1. power comparison across family, sibpair and case-control -- horizontal view view
cex=1.2
lwd=4
dev.new(width=10, height=5)
par(mfrow=c(1,3), cex=cex, mar=c(3.1, 0.1, 1.6, 0.1), oma=c(0,3,0,1), lwd=2)
n1 <- 1000
sigma2 <- 0
curve(sapply(x, function(x) fn(affect=affect, n_family=n1/2, f=0.01, mu=x, p_dis=0.01, delta=1, alpha=alpha, stat=F)[3]), 1, 4.5, ylim=c(0,1),xlab="",ylab="", lwd=lwd)
curve(sapply(x, function(x) power_comp.sas(mu=x, sigma2=sigma2, f=0.01, n_pair=n1, out="A", alpha=alpha)), col="darkorange", add=TRUE, lwd=lwd)
mtext("Power", 2, line=2, cex=cex, adj=0.5)
mtext("f=0.01", 3, cex=cex, line=0.2, adj=0.5)
par(yaxt="n")
curve(sapply(x, function(x) fn(affect=affect, n_family=n1/2, f=0.05, mu=x, p_dis=0.05, delta=1, alpha=alpha, stat=F)[3]), 1, 4.5, ylim=c(0,1),xlab="",ylab="", lwd=lwd)
curve(sapply(x, function(x) power_comp.sas(mu=x, sigma2=sigma2, f=0.05, n_pair=n1, out="A", alpha=alpha)), col="darkorange", add=TRUE, lwd=lwd)
mtext("f=0.05", 3, cex=cex, line=0.2, adj=0.5)
mtext("Mean Relative Risk", 1, outer=FALSE, line=2, cex=cex)
curve(sapply(x, function(x) fn(affect=affect, n_family=n1/2, f=0.2, mu=x, p_dis=0.01, delta=1, alpha=alpha, stat=F)[3]), 1, 4.5, ylim=c(0,1),xlab="",ylab="", lwd=lwd)
curve(sapply(x, function(x) power_comp.sas(mu=x, sigma2=sigma2, f=0.2, n_pair=n1, out="A", alpha=alpha)), col="darkorange", add=TRUE, lwd=lwd)
mtext("f=0.2", 3, cex=cex, line=0.2, adj=0.5)
legend("bottomright", c("Family", "Case-Control"), lty=1, col=c("black", "darkorange"), cex=1.1, lwd=lwd)


#2. mu varying -- power comparison across family, sibpair and case-control -- vertical view
cex=1.2
lwd=4
dev.new(width=5, height=10)
par(mfrow=c(3,1), cex=cex, mar=c(0.8, 0.1, 0.8, 0.1), oma=c(2.5,3,1.5,1))
n1 <- 1000
sigma2 <- 0
par(xaxt="n")
curve(sapply(x, function(x) fn(affect=affect, n_family=n1/2, f=0.01, mu=x, p_dis=0.01, delta=1, alpha=alpha, stat=F)[3]), 1, 4.5, ylim=c(0,1),xlab="",ylab="", lwd=lwd)
mtext("f=0.01", 3, cex=1.5, line=0.2, adj=0.5)
curve(sapply(x, function(x) power_comp.sas(mu=x, sigma2=sigma2, f=0.01, n_pair=n1, out="A", alpha=alpha)), col="darkorange", add=TRUE, lwd=lwd)
par(xaxt="n")
curve(sapply(x, function(x) fn(affect=affect, n_family=n1/2, f=0.05, mu=x, p_dis=0.05, delta=1, alpha=alpha, stat=F)[3]),1, 4.5, ylim=c(0,1),xlab="",ylab="", lwd=lwd)
mtext("f=0.05", 3, cex=1.5, line=0.2, adj=0.5)
mtext("Power", 2, line=2, cex=1.5, adj=0.5)
curve(sapply(x, function(x) power_comp.sas(mu=x, sigma2=sigma2, f=0.05, n_pair=n1, out="A", alpha=alpha)), col="darkorange", add=TRUE, lwd=lwd)
par(xaxt="s")
curve(sapply(x, function(x) fn(affect=affect, n_family=n1/2, f=0.2, mu=x, p_dis=0.01, delta=1, alpha=alpha, stat=F)[3]), 1, 4.5, ylim=c(0,1),xlab="",ylab="", lwd=lwd)
mtext("f=0.2", 3, cex=1.5, line=0.2, adj=0.5)
curve(sapply(x, function(x) power_comp.sas(mu=x, sigma2=sigma2, f=0.2, n_pair=n1, out="A", alpha=alpha)), col="darkorange", add=TRUE, lwd=lwd)
mtext("Mean Relative Risk", 1, outer=FALSE, line=2, cex=cex)
legend("bottomright", c("Family", "Case-Control"), lty=1, col=c("black","darkorange"), cex=1, lwd=lwd)

#3. sigma^2 varying -- power comparison across family, sibpair and case-control -- vertical view
cex=1.2
lwd=4
dev.new(width=5, height=10)
par(mfrow=c(3,1), cex=cex, mar=c(0.8, 0.1, 0.8, 0.1), oma=c(2.5,3,1.5,1))
n1 <- 1000
mu=1.2
par(xaxt="n")
curve(sapply(x, function(x) fn(affect=affect, n_family=n1/2, f=0.01, mu=mu, sigma2=x, p_dis=0.01, delta=1, alpha=alpha, stat=F)[3]), 1, 4.5, ylim=c(0,1),xlab="",ylab="", lwd=lwd)
mtext("f=0.01", 3, cex=1.5, line=0.2, adj=0.5)
curve(sapply(x, function(x) power_comp.sas(mu=mu, sigma2=x, f=0.01, n_pair=n1, out="A", alpha=alpha)), col="darkorange", add=TRUE, lwd=lwd)
par(xaxt="n")
curve(sapply(x, function(x) fn(affect=affect, n_family=n1/2, f=0.05, mu=mu, sigma2=x, p_dis=0.05, delta=1, alpha=alpha, stat=F)[3]),1, 4.5, ylim=c(0,1),xlab="",ylab="", lwd=lwd)
mtext("f=0.05", 3, cex=1.5, line=0.2, adj=0.5)
mtext("Power", 2, line=2, cex=1.5, adj=0.5)
curve(sapply(x, function(x) power_comp.sas(mu=mu, sigma2=x, f=0.05, n_pair=n1, out="A", alpha=alpha)), col="darkorange", add=TRUE, lwd=lwd)
par(xaxt="s")
curve(sapply(x, function(x) fn(affect=affect, n_family=n1/2, f=0.2, mu=mu, sigma2=x, p_dis=0.01, delta=1, alpha=alpha, stat=F)[3]), 1, 4.5, ylim=c(0,1),xlab="",ylab="", lwd=lwd)
mtext("f=0.2", 3, cex=1.5, line=0.2, adj=0.5)
curve(sapply(x, function(x) power_comp.sas(mu=mu, sigma2=x, f=0.2, n_pair=n1, out="A", alpha=alpha)), col="darkorange", add=TRUE, lwd=lwd)
mtext("Mean Relative Risk", 1, outer=FALSE, line=2, cex=cex)
legend("right", c("Family", "Case-Control"), lty=1, col=c("black","darkorange"), cex=1, lwd=lwd)

#4. sigma^2 varying -- power comparison across family, sibpair and case-control -- vertical view
cex=1.2
lwd=4
dev.new(width=4.5, height=4.5)
par(mfrow=c(1,1), cex=cex, mar=c(0.8, 0.1, 0.8, 0.1), oma=c(2.5,3,1.5,1))
n1 <- 1000
mu=1.2
curve(sapply(x, function(x) fn(affect=affect, n_family=n1/2, f=0.01, mu=mu, sigma2=x, p_dis=0.01, delta=1, alpha=alpha, stat=F)[3]), 1, 4.5, ylim=c(0,0.5),xlab="",ylab="", lwd=lwd)
mtext("f=0.01, mean relative risk=1.2", 3, cex=1.5, line=0.2, adj=0.5)
mtext("Variance of Relative Risk", 1, outer=FALSE, line=2, cex=cex)
mtext("Power", 2, line=2, cex=1.5, adj=0.5)
curve(sapply(x, function(x) power_comp.sas(mu=mu, sigma2=x, f=0.01, n_pair=n1, out="A", alpha=alpha)), col="darkorange", add=TRUE, lwd=lwd)
legend("topleft", c("Family", "Case-Control"), lty=1, col=c("black","darkorange"), cex=1.3, lwd=lwd)
