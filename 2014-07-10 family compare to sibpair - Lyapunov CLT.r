source("function.r")

lyapunov(fam.str=c(2,2), affect=c(c(0,1),c(1,1)), n_family=500, f=0.01, mu=1, p_dis=0.01, delta=1, stat=TRUE)

n_family <- 50000
fn <- Vectorize(function(x) lyapunov(fam.str=c(2,2), affect=c(c(0,1),c(1,1)), n_family=x, f=0.01, mu=1, p_dis=0.01, delta=1, stat=FALSE)[2])
curve(fn, 500, n_family, xlab="Sample Size", ylab="Lyapunov Condition", ylim=c(0, 0.35))
fn <- Vectorize(function(x) lyapunov(fam.str=c(2,2), affect=c(c(0,1),c(1,1)), n_family=x, f=0.05, mu=1, p_dis=0.01, delta=1, stat=FALSE)[2])
curve(fn, 500, n_family, add=TRUE, col=2)
fn <- Vectorize(function(x) lyapunov(fam.str=c(2,2), affect=c(c(0,1),c(1,1)), n_family=x, f=0.20, mu=1, p_dis=0.01, delta=1, stat=FALSE)[2])
curve(fn, 500, n_family, add=TRUE, col=4, lwd=3)
legend("topright", c("f=0.01", "f=0.05", "f=0.20"), col=c(1,2,4), lty=1)
