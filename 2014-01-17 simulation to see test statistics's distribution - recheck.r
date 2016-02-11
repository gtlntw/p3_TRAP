##load project related functions
source('function.r')

#generate data
gene.data.family(m=2, var=0, f=0.01, SRR=1000)
#simulation
result <- replicate(1000, {
  n_sample <- 500
  alpha <- 0.05
  data.temp=NULL
  temp <- data.family[sample(nrow(data.family), n_sample),]
  attach(temp)
  data.temp <- temp[which(!(H1_f==0 & H2_f==0 & H1_m==0 & H2_m==0) & !(H1_f==1 & H2_f==1 & H1_m==1 & H2_m==1)),]
  detach(temp)
  T <- data.temp$T
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