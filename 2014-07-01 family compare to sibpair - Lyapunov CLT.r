source("function.r")

#n=500
result <- read.csv("2014-07-01 family compare to sibpair - Lyapunov CLT_n500.csv", header=T)
f.1 <- subset(result, f==0.01)
f.5 <- subset(result, f==0.05)
f.20 <- subset(result, f==0.20)

mean(f.1$T)
hist(f.1$T)
mean(f.1$p_value <0.05)
hist(f.1$p_value)
mean(f.1$Lyapunov)

mean(f.5$T)
hist(f.5$T)
mean(f.5$p_value <0.05)
hist(f.5$p_value)
mean(f.5$Lyapunov)

mean(f.20$T)
hist(f.20$T)
mean(f.20$p_value <0.05)
hist(f.20$p_value)
mean(f.20$Lyapunov)


#n=1000
result <- read.csv("2014-07-01 family compare to sibpair - Lyapunov CLT_n1000.csv", header=T)
f.1 <- subset(result, f==0.01)
f.5 <- subset(result, f==0.05)
f.20 <- subset(result, f==0.20)

mean(f.1$T)
hist(f.1$T)
mean(f.1$p_value <0.05)
hist(f.1$p_value)
mean(f.1$Lyapunov)

mean(f.5$T)
hist(f.5$T)
mean(f.5$p_value <0.05)
hist(f.5$p_value)
mean(f.5$Lyapunov)

mean(f.20$T)
hist(f.20$T)
mean(f.20$p_value <0.05)
hist(f.20$p_value)
mean(f.20$Lyapunov)


#n=2000
result <- read.csv("2014-07-01 family compare to sibpair - Lyapunov CLT_n2000.csv", header=T)
f.1 <- subset(result, f==0.01)
f.5 <- subset(result, f==0.05)
f.20 <- subset(result, f==0.20)

mean(f.1$T)
hist(f.1$T)
mean(f.1$p_value <0.05)
hist(f.1$p_value)
mean(f.1$Lyapunov)

mean(f.5$T)
hist(f.5$T)
mean(f.5$p_value <0.05)
hist(f.5$p_value)
mean(f.5$Lyapunov)

mean(f.20$T)
hist(f.20$T)
mean(f.20$p_value <0.05)
hist(f.20$p_value)
mean(f.20$Lyapunov)
