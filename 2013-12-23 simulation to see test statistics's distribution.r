result <- replicate(1000, {
  n_sample <- 500
  alpha <- 0.05
  data.temp=NULL
  temp <- data.family[sample(nrow(data.family), n_sample),]
  with(temp,
  data.temp <- temp[which(H1_f==0 & H2_f==0 & H1_m==0 & H2_m==0) & !(H1_f==1 & H2_f==1 & H1_m==1 & H2_m==1),]
  )
  T <- data.temp$T
  E.T <- mean(T)
  VAR.T <- var(T)
  abs((T-E.T)/sqrt(VAR.T))
  }
)
