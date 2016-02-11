##load project related functions
source('function.r')

#generate data
mu <- 3
gene.data.family(m=mu, var=0, f=0.01, SRR=5000)

#simulation
result <- replicate(1000, {
  n_sample <- 500
  alpha <- 0.05
  data.temp=NULL
  temp <- data.family[sample(nrow(data.family), n_sample),]
  attach(temp)
  data.temp <- temp[which(!(H1_f==0 & H2_f==0 & H1_m==0 & H2_m==0) & !(H1_f==1 & H2_f==1 & H1_m==1 & H2_m==1)),]
  detach(temp)
  T <- data.temp$T_stat
  E.T <- mean(T)
  VAR.T <- var(T)
  #abs((T-E.T)/sqrt(VAR.T))
  
  c(E.T, VAR.T, length(T))
  }
)

pdf("myplot.pdf")
hist(result[1,])
dev.off()

var(result[1,])
mean(result[2,])/mean(result[3,])
mean(result[2,]/(result[3,]))
stat[5]/(stat[6]*n_family)
  
#setting
alpha=10^-6
n_family <- 500
#expectation for the family setting
stat <- power.family.4(fam.str=c(2,2), affect=c(c(0,1),c(1,1)), n_family=500, f=0.01, mu=mu, p_dis=0.01, stat=T)
E_N <- stat[1]
sigma2_N <- stat[2]
p_N_n0 <- stat[3]
crit.pt.L <- qnorm(alpha/2, mean=E_N, sd=sqrt(sigma2_N/(p_N_n0*n_family)))
crit.pt.H <- qnorm(1-alpha/2, mean=E_N, sd=sqrt(sigma2_N/(p_N_n0*n_family)))
#observed mean
E_A <- result[1,]
sigma2_A <- result[2,]
n_sample <- result[3,]
#calculate the test statistics
round(power.family.4(fam.str=c(2,2), affect=c(c(0,1),c(1,1)), n_family=500, f=0.01, mu=mu, p_dis=0.01), 6)
power <- pnorm(crit.pt.L, mean=E_A, sd=sqrt(sigma2_A/(n_sample))) + pnorm(crit.pt.H, mean=E_A, sd=sqrt(sigma2_A/(n_sample)), lower=F)
mean(power)

mean(E_A)
stat[4]
mean(sigma2_A)
stat[5]
mean(n_sample)
stat[6]*n_family