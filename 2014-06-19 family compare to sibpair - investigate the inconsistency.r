source("function.r")
result <- read.csv("allResults.csv", header=T)

##calculate the simulated false positive rate
alpha=0.05
sim_result <- with(result, 
{
  tapply(true<alpha, list(f, r), mean)
}
     )
x <- seq(1,4, by=0.2)

f.1 <- subset(result, f==0.01)
f.5 <- subset(result, f==0.05)
f.20 <- subset(result, f==0.20)

mean(f.1$avg)
sd(f.1$avg)
hist(f.1$avg)
mean(f.1$e_avg)
mean(f.1$sd^2)
mean(f.1$e_sd^2)
power.family.5(affect=c(0,1,1,1), p_dis=0.01, mu=1, f=0.01, n_family=500, alpha=0.05)
mean(f.1$p_value <0.05)
hist(f.1$p_value)
mean(f.1$n_test.data)

mean(f.5$avg)
hist(f.5$avg)
mean(f.5$e_avg)
mean(f.5$sd^2)
sd(f.5$sd^2)
mean(f.5$e_sd^2)
power.family.5(affect=c(0,1,1,1), p_dis=0.01, mu=1, f=0.05, n_family=500, alpha=0.05)
mean(f.5$p_value <0.05)
hist(f.5$p_value)
mean(f.5$n_test.data)

mean(f.20$avg)
sd(f.20$avg)
hist(f.20$avg)
mean(f.20$e_avg)
mean(f.20$sd^2)
sd(f.20$sd^2)
mean(f.20$e_sd^2)
power.family.5(affect=c(0,1,1,1), p_dis=0.01, mu=1, f=0.20, n_family=500, alpha=0.05)
mean(f.20$p_value <0.05)
hist(f.20$p_value)
mean(f.20$n_test.data)

#n increases from 500 to 2500
result <- read.csv("2014-06-19 family compare to sibpair - investigate the inconsistency_n2500.csv", header=T)
f.1 <- subset(result, f==0.01)
f.5 <- subset(result, f==0.05)
mean(f.1$true<0.05)
mean(f.5$true<0.05)

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
