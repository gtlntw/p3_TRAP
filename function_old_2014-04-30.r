##Theoratical allele frequency: E(HS|AAr, S=1,2) and E(HN|AAr, S=0,1)on a chromosome
EH <- function(mu=1, sigma2=0, f=0.01) {
  mu <- mu
  sigma2 <- sigma2
  p <- f
  P_AA <- 0.25*(1+p*(mu-1))^4 + 0.5*(1+p*(mu-1))^2*(1+p*(mu^2+sigma2-1)) + 0.25*(1+p*(mu^2+sigma2-1))^2
  
  E_S_C <- (0.5*(1+p*(mu-1))^2*(1+p*(mu^2+sigma2-1))) / P_AA + 2*(0.25*(1+p*(mu^2+sigma2-1))^2) / P_AA
  E_NS_C <- 4*(0.25*(1+p*(mu-1))^4)/P_AA + 2*0.5*(1+p*(mu-1))^2*(1+p*(mu^2+sigma2-1))/P_AA
  
  P_HS1 <- ((mu^2+sigma2)*((1-p)+mu*p)^2*p*0.5+(mu^2+sigma2)*2*p*(1-p)*0.25)/P_AA
  P_HS2 <- (mu^2+sigma2)^2*p^2*0.25/P_AA
  E_HS <- P_HS1 + 2*P_HS2
  
  P_HN1 <- (mu*4*p*(1-p)^3*0.25 + (mu*(1-p)+mu*(mu^2+sigma2)*p)*2*p*(1-p)*0.5)/P_AA
  P_HN2 <- (mu^2*6*p^2*(1-p)^2*0.25 + (mu^2*(1-p)+mu^2*(mu^2+sigma2)*p)*p^2*0.5)/P_AA
  P_HN3 <- (mu^3*4*p^3*(1-p)*0.25)/P_AA
  P_HN4 <- (mu^4*p^4*0.25)/P_AA
  E_HN <- P_HN1 + 2*P_HN2 + 3*P_HN3 + 4*P_HN4
  
  c(E_HS/E_S_C, E_HN/E_NS_C)
}
EH(mu=5, sigma2=0, f=0.02)

## E(N_C|AAr) E(N_S_C|AAr) E(N_NS_C|AAr) no. of chromosomes
E_C <- function(mu=1, sigma2=0, f=0.01) {
  mu <- mu
  sigma2 <- sigma2
  p <- f
  P_AA <- 0.25*(1+p*(mu-1))^4 + 0.5*(1+p*(mu-1))^2*(1+p*(mu^2+sigma2-1)) + 0.25*(1+p*(mu^2+sigma2-1))^2
  
  E_S_C <- (1*0.5*(1+p*(mu-1))^2*(1+p*(mu^2+sigma2-1)) + 2*0.25*(1+p*(mu^2+sigma2-1))^2)/P_AA
  E_NS_C <- (4*(0.25*(1+p*(mu-1))^4) + 2*0.5*(1+p*(mu-1))^2*(1+p*(mu^2+sigma2-1)))/P_AA
  
  E_C <- (4*(0.25*(1+p*(mu-1))^4) + 3*0.5*(1+p*(mu-1))^2*(1+p*(mu^2+sigma2-1)) + 2*0.25*(1+p*(mu^2+sigma2-1))^2)/P_AA
  
  c(E_S_C, E_NS_C, E_C)
}

power.sas <- function(mu=1.2, sigma2=3, f=0.01, n_sb=50, alpha=10^-6) {
  e_c <- E_C(mu=mu, sigma2=sigma2, f=f) ## expected no. of chromosomes per sibpair
  N <- n_sb*e_c[3] #no. of  total chromosomes
  n1 <- n_sb*e_c[1]  ##no. of shared chromosomes
  n2 <- n_sb*e_c[2]  #no. of non-shared chromosomes
  
  P_A <- EH(mu=mu, sigma2=sigma2, f=f) ##proportion under HA
  p1 <- P_A[1]
  p2 <- P_A[2]
  p0 <- 0
  w1 <- n1/N
  w2 <- n2/N
  
  pnorm( ((p2-p1-p0)*(N*w1*w2)^.5 -  qnorm(1-alpha/2)*((w1*p1+w2*p2)*(1-w1*p1-w2*p2))^.5)/(w2*p1*(1-p1)+w1*p2*(1-p2))^.5 ) +
    pnorm( (-(p2-p1-p0)*(N*w1*w2)^.5 -  qnorm(1-alpha/2)*((w1*p1+w2*p2)*(1-w1*p1-w2*p2))^.5)/(w2*p1*(1-p1)+w1*p2*(1-p2))^.5) 
}
power.sas(mu=3.0, sigma2=0, f=0.01, n_sb=1000, alpha=10^-6)

##use chisquare distribution with noncentrality parameter to calculate power
power.sas.ncp <- function(mu=1.2, sigma2=3, f=0.01, n_sb=50, alpha=10^-6) {
  e_c <- E_C(mu=mu, sigma2=sigma2, f=f) ## expected no. of chromosomes per sibpair
  N <- n_sb*e_c[3] #no. of  total chromosomes
  n1 <- n_sb*e_c[1]  ##no. of shared chromosomes
  n2 <- n_sb*e_c[2]  #no. of non-shared chromosomes
  
  P_A <- EH(mu=mu, sigma2=sigma2, f=f) ##proportion under HA
  p1 <- P_A[1]
  p2 <- P_A[2]
  p0 <- 0
  w1 <- n1/N
  w2 <- n2/N
  
  pchisq((-qnorm(1-alpha/2)*((w1*p1+w2*p2)*(1-w1*p1-w2*p2))^.5/(w2*p1*(1-p1)+w1*p2*(1-p2))^.5)^2, df=1, ncp=(-(p2-p1-p0)*(N*w1*w2)^.5/(w2*p1*(1-p1)+w1*p2*(1-p2))^.5)^2, lower=F)
}
power.sas.ncp(mu=3.0, sigma2=0, f=0.01, n_sb=1000, alpha=10^-6)

##sample size calculation -- doesn't really use
samplesize.sas <- function(mu=1.2, sigma2=3, f=0.01, power=0.8, alpha=10^-6) {
  e_c <- E_C(mu=mu, sigma2=sigma2, f=f) ## expected no. of chromosomes per sibpair
  P_A <- EH(mu=mu, sigma2=sigma2, f=f) ##proportion under HA
  p1 <- P_A[1]
  p2 <- P_A[2]
  p0 <- 0
  
  obj_fn <- function(n_sb) {
    N <- n_sb*e_c[3] #no. of  total chromosomes
    n1 <- n_sb*e_c[1]  ##no. of shared chromosomes
    n2 <- n_sb*e_c[2]  #no. of non-shared chromosomes
    w1 <- n1/N
    w2 <- n2/N
    
    (pnorm( ((p2-p1-p0)*(N*w1*w2)^.5 -  qnorm(1-alpha/2)*((w1*p1+w2*p2)*(1-w1*p1-w2*p2))^.5)/(w2*p1*(1-p1)+w1*p2*(1-p2))^.5 ) +
       pnorm( (-(p2-p1-p0)*(N*w1*w2)^.5 -  qnorm(1-alpha/2)*((w1*p1+w2*p2)*(1-w1*p1-w2*p2))^.5)/(w2*p1*(1-p1)+w1*p2*(1-p2))^.5) 
     - power)^2
  }
  
  init_n <- (qnorm(1-alpha)*((0.33*p1+0.67*p2)*(1-0.33*p1-0.67*p2))^0.5 + qnorm(power)*(0.67*p1*(1-p1)+0.33*p2*(1-p2))^0.5)^2 / (0.33*0.67*(p2-p1-p0)^2)
  ceiling(optimize(obj_fn, c(1,init_n))$minimum)
  #ceiling(optim(1000, obj_fn)$par)
}
samplesize.sas(mu=5, sigma2=0, f=0.01, power=0.8, alpha=10^-6)


##Theoratical allele frequency: Unaffected and Affected, E(H|A) and E(H|AAr) on a chromosome for comparison
EH_comp <- function(mu=1, sigma2=0, f=0.01, out="A") {
  mu <- mu
  sigma2 <- sigma2
  p <- f
  P_A <- (1+p*(mu-1))^2
  P_AA <- 0.25*(1+p*(mu-1))^4 + 0.5*(1+p*(mu-1))^2*(1+p*(mu^2+sigma2-1)) + 0.25*(1+p*(mu^2+sigma2-1))^2
  if(out=="A") {
    return(c(EH_A <- (mu*2*p*(1-p)+2*mu^2*p^2)/P_A/2, p))
  }
  else {
    return(c(EH_AA <- (0.25*2*p*(1-p)*(mu*(1+p*(mu-1))^2) + 0.5*2*p*(1-p)*(0.5*((1-p)*mu+p*mu^2)+0.5*((1-p)*(mu^2+sigma2)+p*mu*(mu^2+sigma2))) 
                       + 0.25*2*p*(1-p)*(mu^2+sigma2) + 2*(0.25*p^2*mu^2*(1+p*(mu-1))^2+0.5*p^2*((1-p)*mu*(mu^2+sigma2)+p*mu^2*(mu^2+sigma2))+0.25*p^2*(mu^2+sigma2)^2) )/P_AA/2,p))
  }  		   
}

EH_comp(mu=5, sigma2=0, f=0.02, out="A")


power_comp.sas <- function(mu=1.2, sigma2=3, f=0.01, n_pair=50, alpha=10^-6, out="A") {
  N <- 4*n_pair
  P_A <- EH_comp(mu=mu, sigma2=sigma2, f=f, out=out) ##proportion under HA
  p1 <- P_A[1]
  p2 <- P_A[2]
  p0 <- 0
  w1 <- 0.5
  w2 <- 0.5
  
  pnorm( ((p2-p1-p0)*(N*w1*w2)^.5 -  qnorm(1-alpha/2)*((w1*p1+w2*p2)*(1-w1*p1-w2*p2))^.5)/(w2*p1*(1-p1)+w1*p2*(1-p2))^.5 ) +
    pnorm( (-(p2-p1-p0)*(N*w1*w2)^.5 -  qnorm(1-alpha/2)*((w1*p1+w2*p2)*(1-w1*p1-w2*p2))^.5)/(w2*p1*(1-p1)+w1*p2*(1-p2))^.5) 
}
power_comp.sas(mu=3.0, sigma2=0, f=0.01, n_pair=1000, alpha=10^-6, out="A")

samplesize_comp.sas <- function(mu=1.2, sigma2=3, f=0.01, power=0.8, alpha=10^-6, out="A") {
  P_A <- EH_comp(mu=mu, sigma2=sigma2, f=f, out=out) ##proportion under HA
  p1 <- P_A[1]
  p2 <- P_A[2]
  p0 <- 0
  w1 <- 0.5
  w2 <- 0.5
  
  obj_fn <- function(n_pair) {
    
    N <- 4*n_pair    
    (pnorm( ((p2-p1-p0)*(N*w1*w2)^.5 -  qnorm(1-alpha/2)*((w1*p1+w2*p2)*(1-w1*p1-w2*p2))^.5)/(w2*p1*(1-p1)+w1*p2*(1-p2))^.5 ) +
       pnorm( (-(p2-p1-p0)*(N*w1*w2)^.5 -  qnorm(1-alpha/2)*((w1*p1+w2*p2)*(1-w1*p1-w2*p2))^.5)/(w2*p1*(1-p1)+w1*p2*(1-p2))^.5) 
     - power)^2
  }
  
  init_n <- (qnorm(1-alpha)*((0.5*p1+0.5*p2)*(1-0.5*p1-0.5*p2))^0.5 + qnorm(power)*(0.5*p1*(1-p1)+0.5*p2*(1-p2))^0.5)^2 / (0.5*0.5*(p2-p1-p0)^2)
  ceiling(optimize(obj_fn, c(1,init_n))$minimum)
  #ceiling(optim(1000, obj_fn)$par)
  #c(init_n, optimize(obj_fn, c(1,init_n)))
}
samplesize_comp.sas(mu=5, sigma2=0, f=0.01, power=0.8, alpha=10^-6)



### Using Gamma distribution to generate the efffect of haplotype
##simulate a case/control
gene.data <- function(m=1, var=0, f=0.01, SRR=5, p_dis=0.01, pop=FALSE) {
  adj <- ifelse(var==0 & m==1, 1, ifelse(var==0, m, 1))
  c = (m - 1) #shift to accommondate gamma distribution
  beta <- ifelse(var==0 & m==1, 0, var/c) #prevent not a number error when var==0 & m==1
  alpha = c/beta
  KL <- (1+f*(m-1))^2 #contribution from locus
  KLKLR <- 0.25*(1+f*(m-1))^4 + 0.5*(1+f*(m-1))^2*(1+f*(m^2+var-1)) + 0.25*(1+f*(m^2+var-1))^2
  KG <- p_dis/KL #contribution from other genome
  SRR <- SRR #sibling relaive risk
  KGKGR <- SRR*p_dis*p_dis/KLKLR #implement the heriatbility from other locus given SRR
  
  ##generate case and control individuals from population
  if(pop){ #skip individual part
    n_pop <- 5000000
    H1 <- rbinom(n_pop,1,f) #if the first haplotype carries risk raviant
    H2 <- rbinom(n_pop,1,f) #if the second haplotype carries risk raviant
    H1_rr <- ifelse(H1 == 1, rgamma(length(H1), alpha, scale=beta) + adj, 1) #RR of the first haplotype
    H2_rr <- ifelse(H2 == 1, rgamma(length(H2), alpha, scale=beta) + adj, 1) #RR of the second haplotype
    penetrance <- H1_rr*H2_rr*KG #penetrance of disease given haplotypes
    penetrance <- ifelse(penetrance>1, 1, penetrance)
    dis <- rbinom(length(penetrance),1,penetrance) #disease status
    data.pop <<- data.frame(H1=H1, H2=H2, H1_rr=H1_rr, H2_rr=H2_rr, penetrance=penetrance, dis=dis)
  }
  
  ##generate sibpairs
  n_family <- 50000000
  H1_f <- rbinom(n_family,1,f) #if the first haplotype carries risk raviant
  H2_f <- rbinom(n_family,1,f) #if the second haplotype carries risk raviant
  H1_rr_f <- ifelse(H1_f == 1, rgamma(length(H1_f), alpha, scale=beta) + adj, 1) #RR of the first haplotype
  H2_rr_f <- ifelse(H2_f == 1, rgamma(length(H2_f), alpha, scale=beta) + adj, 1) #RR of the second haplotype
  H1_m <- rbinom(n_family,1,f) #if the first haplotype carries risk raviant
  H2_m <- rbinom(n_family,1,f) #if the second haplotype carries risk raviant
  H1_rr_m <- ifelse(H1_m == 1, rgamma(length(H1_m), alpha, scale=beta) + adj, 1) #RR of the first haplotype
  H2_rr_m <- ifelse(H2_m == 1, rgamma(length(H2_m), alpha, scale=beta) + adj, 1) #RR of the second haplotype
  
  H1_inh_s <- sample(c(1,2),n_family, replace=TRUE) 
  H1_s <- ifelse(H1_inh_s==1, H1_f, H2_f)
  H1_rr_s <- ifelse(H1_inh_s==1, H1_rr_f, H2_rr_f)
  
  H2_inh_s <- sample(c(1,2),n_family, replace=TRUE) 
  H2_s <- ifelse(H2_inh_s==1, H1_m, H2_m)
  H2_rr_s <- ifelse(H2_inh_s==1, H1_rr_m, H2_rr_m)
  
  #penetrance_s <- H1_rr_s*H2_rr_s #penetrance of disease given haplotypes
  #penetrance_s <- ifelse(penetrance_s>1, 1, penetrance_s)
  #dis_s <- rbinom(length(penetrance_s),1,penetrance_s) #disease status
  
  H1_inh <- sample(c(1,2), n_family, replace=TRUE) 
  H1 <- ifelse(H1_inh==1, H1_f, H2_f)
  H1_rr <- ifelse(H1_inh==1, H1_rr_f, H2_rr_f)
  
  H2_inh <- sample(c(1,2), n_family, replace=TRUE) 
  H2 <- ifelse(H2_inh==1, H1_m, H2_m)
  H2_rr <- ifelse(H2_inh==1, H1_rr_m, H2_rr_m)
  
  rm(H1_f, H2_f, H1_rr_f, H2_rr_f, H1_m, H2_m, H1_rr_m, H2_rr_m) #save memory
  
  S0 <- (H1_inh_s != H1_inh) & (H2_inh_s != H2_inh)
  S2 <- (H1_inh_s == H1_inh) & (H2_inh_s == H2_inh)
  S1 <- ((H1_inh_s == H1_inh) & (H2_inh_s != H2_inh)) | ((H1_inh_s != H1_inh) & (H2_inh_s == H2_inh))
  S11 <- ((H1_inh_s == H1_inh) & (H2_inh_s != H2_inh))
  S12 <- ((H1_inh_s != H1_inh) & (H2_inh_s == H2_inh))
  
  S <- (S0==TRUE)*0 + (S1==TRUE)*1 + (S2==TRUE)*2
  HS.mis <- HS <- (S1==TRUE&S11)*(H1) + (S1==TRUE&S12)*(H2) + (S2==TRUE)*(H1+H2)
  HN.mis <- HN <- (S0==TRUE)*(H1+H2+H1_s+H2_s) + (S1==TRUE&S11)*(H2+H2_s) + (S1==TRUE&S12)*(H1+H1_s)
  HS.mis <- ifelse(S1==TRUE & HS==0 & HN==2, 1, HS.mis) 
  HN.mis <- ifelse(S1==TRUE & HS==0 & HN==2, 0, HN.mis) 
  
  penetrance <- H1_rr*H2_rr*H1_rr_s*H2_rr_s*KGKGR #penetrance of disease given haplotypes of both siblings
  penetrance <- ifelse(penetrance>1, 1, penetrance)
  dis <- rbinom(length(penetrance),1,penetrance) #disease status of both affected
  
  #data.family <- data.frame(H1_f, H2_f, H1_rr_f, H2_rr_f,   #for varification
  #                          H1_m, H2_m, H1_rr_m, H2_rr_m,
  #                          H1_inh, H1, H1_rr, H2_inh, H2, H2_rr, penetrance, dis,
  #                          H1_inh_s, H1_s, H1_rr_s, H2_inh_s, H2_s, H2_rr_s, penetrance_s, dis_s)
  
  idx <- which(dis==1)
  
  data.family <<- data.frame(H1=H1[idx], H1_rr=H1_rr[idx], H2=H2[idx], H2_rr=H2_rr[idx], dis=dis[idx],
                             H1_s=H1_s[idx], H1_rr_s=H1_rr_s[idx], H2_s=H2_s[idx], H2_rr_s=H2_rr_s[idx],
                             S0=S0[idx], S1=S1[idx], S2=S2[idx], S11=S11[idx], S12=S12[idx],
                             HS=HS[idx],HN=HN[idx],S=S[idx], HS.mis=HS.mis[idx],HN.mis=HN.mis[idx])
}

##EM algorithm for imputation
EM <- function(data) {
  #initialization
  pn.init <- sum(subset(data, S==0)$HN.mis)/(4*sum(data$S==0)) #probability of rare variant on shared chromosome
  pn.cur <- ifelse(pn.init==0, runif(1), pn.init)
  ps.init <- sum(subset(data, S==2)$HS.mis)/(2*sum(data$S==2)) #probability of rare variant on non-shared chromosome
  ps.cur <- ifelse(ps.init==0, runif(1), ps.init)
  kn <- sum(subset(data, S==0 | (S==1 & !(HS.mis==1 &HN.mis==0)))$HN.mis) #known non-shared variants (On ibd 0 or single variants on ibd 1)
  ks <- sum(subset(data, (S==1 & !(HS.mis==1 &HN.mis==0)) | S==2)$HS.mis) #known shared variants  (On ibd 2 or more than two variants on ibd 1)
  cn <- 4*sum(data$S==0) + 2*sum(data$S==1) # total number of shared chromosomes
  cs <- sum(data$S==1) + 2*sum(data$S==2) #total number of non-shared chromosomes
  u <- sum(data$S1==TRUE & data$HS.mis==1 & data$HN.mis==0) #number of unknown configuration (Double hets in IBD 1)
  delta <- Inf
  iter <- 1
  
  while(delta > 10^-6) {
    #E step
    #us <- u*ps.cur/(pn.cur+ps.cur)
    #un <- u*pn.cur/(pn.cur+ps.cur)
    us <- u* ps.cur*(1-pn.cur)^2 / (ps.cur*(1-pn.cur)^2 + (1-ps.cur)*pn.cur^2)
    un <- u* (1-ps.cur)*pn.cur^2 / (ps.cur*(1-pn.cur)^2 + (1-ps.cur)*pn.cur^2)  
    #M-step
    pn.new <- (kn + 2*un)/cn
    ps.new <- (ks+us)/cs
    #print(c(mu.new, sigma2.new, f.new, cor.factor, optim.result$value))
    
    #check convergence
    delta <- max(abs(pn.cur - pn.new), abs(ps.cur - ps.new))
    pn.cur <- pn.new
    ps.cur <- ps.new
    
    #print(c(pn.cur, ps.cur, iter))
    #iter <- iter + 1
  }
  #c(pn.init, ps.init)
  c(ps.cur*(1-pn.cur)^2 / (ps.cur*(1-pn.cur)^2 + (1-ps.cur)*pn.cur^2))
}

EM2 <- function(data) {
  #initialization
  p.init <- (sum(subset(data, S==0)$HN.mis)+sum(subset(data, S==2)$HS.mis))/((2*sum(data$S==2))+(4*sum(data$S==0))) #probability of rare variant on shared chromosome
  p.cur <- ifelse(p.init==0, runif(1), p.init)
  k <- sum(subset(data, S==0 | (S==1 & !(HS.mis==1 &HN.mis==0)))$HN.mis) + sum(subset(data, (S==1 & !(HS.mis==1 &HN.mis==0)) | S==2)$HS.mis)
  c <- 4*sum(data$S==0) + 2*sum(data$S==1) + sum(data$S==1) + 2*sum(data$S==2) 
  u <- sum(data$S1==TRUE & data$HS.mis==1 & data$HN.mis==0) #number of unknown configuration (Double hets in IBD 1)
  delta <- Inf
  iter <- 1
  
  while(delta > 10^-6) {
    us <- 1-p.cur
    un <- p.cur
    #M-step
    p.new <- (k + us + 2*un)/c
    #print(c(mu.new, sigma2.new, f.new, cor.factor, optim.result$value))
    
    #check convergence
    delta <- max(abs(p.cur - p.new))
    p.cur <- p.new
    
    #print(c(pn.cur, ps.cur, iter))
    #iter <- iter + 1
  }
  #c(pn.init, ps.init)
  c(1-p.cur)
}

MI <- function(corr_fac, count.HS1_HN0_S1, case.count.mis, control.count.mis, n_chr_s, n_chr_ns) {
  diff <- NULL
  var <- NULL
  u <- count.HS1_HN0_S1
  N <- 10
  
  for(i in 1:N) {
    us <- rbinom(1, u, corr_fac)
    uns <- u - us
    
    xs <- case.count.mis - u + us
    xns <- control.count.mis + 2*uns
    
    p1 <- xs/n_chr_s
    p2 <- xns/n_chr_ns
    p <- (xs+xns)/(n_chr_s+n_chr_ns)
    
    diff <- cbind(diff, p1-p2)
    var <- cbind(var, p*(1-p)*(1/n_chr_s+1/n_chr_ns))
  }
  
  TD <- mean(diff)
  VARD <- mean(var) + (1+1/N)*sum((diff-TD)^2)/(N-1)
  pchisq(TD^2/VARD, df=1, lower=F)
}

#Power curve simulation
#count allele on shared and non-shared chromosomes
#sibpair
sim <- function(rep=100000, var=0, f=0.01, sibpair=1000, SRR=5, r=seq(1, 3.5, by=0.1)) {
  sim.result <- NULL
  r <<- r
  for(i in r) {
    print(paste("Generating Data...m=", i, "var=", var, "f=", f))
    gene.data(m=i, var=var, f=f, SRR=SRR)
    mark <<- 1 #reset the current iteration
    print(mark)  #start iteration
    sim <- replicate(rep, {
      ##sibpair internal control
      data.sample <- data.family[sample(1:nrow(data.family), sibpair), ]
      attach(data.sample)
      n_chr_s <- sum(S1==1) + 2*sum(S2==1)
      n_chr_ns <- 4*sum(S0==1) + 2*sum(S1==1)
      
      case.count.mis <- sum(HS.mis)
      control.count.mis <- sum(HN.mis)
      
      #chi-square test on misclassification data
      temp <- matrix(c(case.count.mis, n_chr_s - case.count.mis, control.count.mis, n_chr_ns - control.count.mis), ncol = 2)
      test.mis <- chisq.test(temp, correct=F) ##prevent anti-conservative
      
      
      #chi-square test on corrected data by EM
      corr_fac <- EM(data.sample)
      
      count.HS1_HN0_S1 <- sum(S1==TRUE & HS.mis==1 & HN.mis==0)
      
      case.count.mis.cor <- corr_fac*count.HS1_HN0_S1 + (case.count.mis - count.HS1_HN0_S1)
      control.count.mis.cor <- control.count.mis + 2*(1-corr_fac)*count.HS1_HN0_S1
      
      temp <- matrix(c(case.count.mis.cor, n_chr_s - case.count.mis.cor, control.count.mis.cor, n_chr_ns - control.count.mis.cor), ncol = 2)
      test.mis.cor <- chisq.test(temp, correct=F) ##prevent anti-conservative
      
      #chi-square test on corrected data by Multiple Imputation
      test.mis.cor.MI <- MI(corr_fac, count.HS1_HN0_S1, case.count.mis, control.count.mis, n_chr_s, n_chr_ns)
      
      #power under true model
      case.count.true <- sum(HS)
      control.count.true <- sum(HN)
      
      temp <- matrix(c(case.count.true, n_chr_s - case.count.true, control.count.true, n_chr_ns - control.count.true), ncol = 2)
      test.true <- chisq.test(temp, correct=F) ##prevent anti-conservative
      
      detach(data.sample)
      if(mark%%1000 == 0) print(mark) #print every 1000 iterations
      mark <<- mark + 1
      c(test.mis$p.value, test.mis.cor.MI, test.true$p.value)
    }
    )
    sim.result <- cbind(sim.result, rowMeans(sim<10^-6))
    #sim.result
  }
  sim.result
}

sim_var <- function(rep=100000, mu=1, f=0.01, sibpair=1000, SRR=5, r=seq(1, 3.5, by=0.1)) {
  sim.result <- NULL
  r <<- r
  for(i in r) {
    print(paste("Generating Data...var=", i, "mu=", mu, "f=", f))
    gene.data(m=mu, var=i, f=f, SRR=SRR)
    mark <<- 1 #reset the current iteration
    print(mark)  #start iteration
    sim <- replicate(rep, {
      ##sibpair internal control
      data.sample <- data.family[sample(1:nrow(data.family), sibpair), ]
      attach(data.sample)
      n_chr_s <- sum(S1==1) + 2*sum(S2==1)
      n_chr_ns <- 4*sum(S0==1) + 2*sum(S1==1)
      
      case.count.mis <- sum(HS.mis)
      control.count.mis <- sum(HN.mis)
      
      #chi-square test on misclassification data
      temp <- matrix(c(case.count.mis, n_chr_s - case.count.mis, control.count.mis, n_chr_ns - control.count.mis), ncol = 2)
      test.mis <- chisq.test(temp, correct=F) ##prevent anti-conservative
      
      
      #chi-square test on corrected data by EM
      corr_fac <- EM(data.sample)
      
      count.HS1_HN0_S1 <- sum(S1==TRUE & HS.mis==1 & HN.mis==0)
      
      case.count.mis.cor <- corr_fac*count.HS1_HN0_S1 + (case.count.mis - count.HS1_HN0_S1)
      control.count.mis.cor <- control.count.mis + 2*(1-corr_fac)*count.HS1_HN0_S1
      
      temp <- matrix(c(case.count.mis.cor, n_chr_s - case.count.mis.cor, control.count.mis.cor, n_chr_ns - control.count.mis.cor), ncol = 2)
      test.mis.cor <- chisq.test(temp, correct=F) ##prevent anti-conservative
      
      #chi-square test on corrected data by Multiple Imputation
      test.mis.cor.MI <- MI(corr_fac, count.HS1_HN0_S1, case.count.mis, control.count.mis, n_chr_s, n_chr_ns)
      
      #power under true model
      case.count.true <- sum(HS)
      control.count.true <- sum(HN)
      
      temp <- matrix(c(case.count.true, n_chr_s - case.count.true, control.count.true, n_chr_ns - control.count.true), ncol = 2)
      test.true <- chisq.test(temp, correct=F) ##prevent anti-conservative
      
      detach(data.sample)
      if(mark%%1000 == 0) print(mark) #print every 1000 iterations
      mark <<- mark + 1
      c(test.mis$p.value, test.mis.cor.MI, test.true$p.value)
    }
    )
    sim.result <- cbind(sim.result, rowMeans(sim<10^-6))
    #sim.result
  }
  sim.result
}




###############Interaction model########################################
inter_model <- function(p=c(0.01, 0.05), trr_alpha=2, trr_beta=6.2, gamma=1, bl=0.0000001) {
  p <- p #Allele frequencies
  
  phi <- function(x) ifelse(x==1, p[1], 1-p[1])
  phj <- function(x) ifelse(x==1, p[1], 1-p[1])
  pgs <- function(x) ifelse(x==1, p[2], 1-p[2])
  pgt <- function(x) ifelse(x==1, p[2], 1-p[2])
  
  trr_alpha <- trr_alpha #target relative risk
  trr_beta <- trr_beta
  gamma <- gamma
  bl <- bl
  PofS <- c(0.25, 0.5, 0.25)
  
  ##new obj function
  obj <- function(x, trr_alpha, trr_beta, gamma) {
    beta_L <- x[1] 
    beta_G <- x[2]
    
    i=0;j=0;k=0;l=0;
    for(hj in 0:1) {
      for(gs in 0:1) {
        for(gt in 0:1){
          i=i+(beta_G*gamma)^(gs+gt)*gamma^(hj*(gs+gt))*pgs(gs)*pgt(gt)
          j=j+(beta_G)^(gs+gt)*gamma^(hj*(gs+gt))*pgs(gs)*pgt(gt)
        }
      }
      k=k+beta_L^(hj)*phj(hj)*i
      l=l+beta_L^(hj)*phj(hj)*j
    }
    rr_alpha=(beta_L*k)/l
    
    i=0;j=0;k=0;l=0;
    for(gt in 0:1) {
      for(hi in 0:1) {
        for(hj in 0:1){
          i=i+(beta_L*gamma)^(hi+hj)*gamma^(gt*(hi+hj))*phi(hi)*phj(hj)
          j=j+(beta_L)^(hi+hj)*gamma^(gt*(hi+hj))*phi(hi)*phj(hj)
        }
      }
      k=k+beta_G^(gt)*pgt(gt)*i
      l=l+beta_G^(gt)*pgt(gt)*j
    }
    rr_beta=(beta_G*k)/l
    
    return(abs(rr_alpha-trr_alpha)+abs(rr_beta - trr_beta))
  }
  #obj(c(1.671843,1.93), 2, 2, 2)
  
  result <- optim(c(1,1), obj, trr_alpha=trr_alpha, trr_beta=trr_beta, gamma=gamma)
  result$par
  
  #solution of beta_L and beta_G
  alpha <- result$par[1]
  beta <- result$par[2]
  
  #penetrance of one individual P(D|A,B)
  Pgt <- array(NA, c(3,3))
  for(a in 0:2) {
    for (b in 0:2) {
      Pgt[a+1,b+1] <- bl*alpha^a*beta^b*gamma^(a*b)
      if(Pgt[a+1,b+1]>1) Pgt[a+1,b+1] <- 1
      #print(Pgt[a+1,b+1])
    }
  }
  
  # P(i, j|A,B,S,loc)  given the first sibling and sharing status what is the genotype of the second sibling
  Ptr <- array(NA, c(3,3,3,2)) 
  for (loc in 0:1) {
    for (a in 0:2) {
      for (i in 0:2) {
        if(a==i) {
          Ptr[a+1, i+1, 2+1, loc+1]=1;
        }
        else
          Ptr[a+1, i+1, 2+1, loc+1]=0;
        Ptr[a+1, i+1, 0+1, loc+1]=p[loc+1]^i*(1-p[loc+1])^(2-i)
      }
      Ptr[a+1, 1+1, 0+1, loc+1]=2*Ptr[a+1, 1+1, 0+1, loc+1]
    }
    Ptr[2+1, 0+1, 1+1, loc+1]=0
    Ptr[0+1, 2+1, 1+1, loc+1]=0
    Ptr[2+1, 1+1, 1+1, loc+1]=1-p[loc+1]
    Ptr[2+1, 2+1, 1+1, loc+1]=p[loc+1]
    Ptr[1+1, 2+1, 1+1, loc+1]=0.5*p[loc+1]
    Ptr[1+1, 1+1, 1+1, loc+1]=0.5
    Ptr[1+1, 0+1, 1+1, loc+1]=0.5*(1-p[loc+1])
    Ptr[0+1, 1+1, 1+1, loc+1]=p[loc+1]
    Ptr[0+1, 0+1, 1+1, loc+1]=1-p[loc+1]
  }
  
  #penetrance of P(DD|A,B)
  Pab <- array(0, c(3,3)) 
  for(a in 0:2){
    for(b in 0:2){
      for(i in 0:2){
        for(j in 0:2){
          Pab[a+1, b+1]=Pab[a+1, b+1]+Pgt[i+1, j+1]*(PofS[0+1]*Ptr[a+1, i+1, 0+1, 0+1]+PofS[1+1]*Ptr[a+1, i+1, 1+1, 0+1]+PofS[2+1]*Ptr[a+1, i+1, 2+1, 0+1])*(PofS[0+1]*Ptr[b+1, j+1, 0+1, 1+1]+PofS[1+1]*Ptr[b+1, j+1, 1+1, 1+1]+PofS[2+1]*Ptr[b+1, j+1, 2+1, 1+1])
        }
      }
      Pab[a+1, b+1]=Pab[a+1, b+1]*Pgt[a+1, b+1]
    }
  }
  
  #probability of observed genotype P(A,B)
  Pf <- array(0, c(3,3))
  norm <- 0
  for(a in 0:2){
    for(b in 0:2) {
      Pf[a+1, b+1]=p[0+1]^a*(1-p[0+1])^(2-a)*p[1+1]^b*(1-p[1+1])^(2-b)
      if(a==1) Pf[a+1, b+1]=2*Pf[a+1, b+1]
      if(b==1) Pf[a+1, b+1]=2*Pf[a+1, b+1]
      norm = norm + Pf[a+1, b+1]*Pab[a+1, b+1]
    }
  }
  
  #prob. of P(A|DD)
  prob <- array(0, 3)
  for(a in 0:2){
    for(b in 0:2){
      prob[a+1] = prob[a+1]+Pf[a+1, b+1]*Pab[a+1, b+1]/norm
    }
  }
  
  #prevalence
  prev <- 0
  for(a in 0:2){
    for(b in 0:2){
      prev = prev+Pf[a+1, b+1]*Pgt[a+1, b+1]
    }
  }
  
  #sibling relative risk
  srr <- 0
  for(a in 0:2){
    for(b in 0:2){
      srr = srr+Pf[a+1, b+1]*Pab[a+1, b+1]
    }
  }
  prev_dd <- srr ##P(DD)
  srr=srr/prev^2
  
  #print(c(alpha, beta, srr))
  
  #baseline
  prob_b<-array(0, 3)
  norm_b <- 0
  for(a in 0:2){
    for(b in 0:2){
      norm_b = norm_b + Pf[a+1, b+1]*Pgt[a+1, b+1] ##P(D)
      prob_b[a+1]=prob_b[a+1] + Pf[a+1, b+1]*Pgt[a+1, b+1]; ##P(D,A)
    }
  }
  
  for(a in 0:2){
    prob_b[a+1]=prob_b[a+1]/norm_b;  ##P(A|D)
    #print "$a $prob[$a] $prob_b[$a]\n";
  }
  
  Ex_sib=prob[1+1]+2*prob[2+1]; Ex_unrelated=prob_b[1+1]+2*prob_b[2+1];
  #print "Expectations sib-pair $Ex unrelated $Ex_b\t";
  
  #Calculate P(AA|vns1, vns2,vs) summing over P(b,b')
  PAAs <- array(0, c(3,3,3));
  for(vns1 in 0:2){
    for(vns2 in 0:2){
      for(vs in 0:2){
        if(vns1+vs<3 & vns2+vs<3)
        {
          for(b in 0:2){
            for(i in 0:2){
              
              ah=vns1+vs;
              ah2=vns2+vs;
              p_sib=0;
              for(S in 0:2){
                p_sib=p_sib+PofS[S+1]*Ptr[b+1, i+1, S+1, 1+1]*dbinom(b,2,p[1+1]); ##contribution from loci B
                #print "tt $S $PofS[$S] $Ptr[$b, $i, $S, 1]\n";
              }
              #print "$b $i $p_sib\n";
              PAAs[vns1+1, vns2+1, vs+1]=PAAs[vns1+1, vns2+1, vs+1]+Pgt[ah+1, b+1]*Pgt[ah2+1, i+1]*p_sib;
            }}
          #print "test $vns1 $vns2 $vs $PAAs[vns1+1, vns2+1, vs+1]\n";
        }
        else {PAAs[vns1+1, vns2+1, vs+1]=-1;}
      }}}
  #Calculate the prob of shared/nonshared variants based on S
  Psns <- array(0, c(3,3,3))
  norm <- 0;
  for(vns1 in 0:2){
    for(vns2 in 0:2){
      for(vs in 0:2){
        if(vns1+vs<3 & vns2+vs<3)
        {
          for(S in 0:2){
            Psns[vns1+1, vns2+1, vs+1]=Psns[vns1+1, vns2+1, vs+1]+dbinom(vs,S,p[0+1])*dbinom(vns1,2-S,p[0+1])*dbinom(vns2,2-S,p[0+1])*PofS[S+1];
          }
          Psns[vns1+1, vns2+1, vs+1]=Psns[vns1+1, vns2+1, vs+1]*PAAs[vns1+1, vns2+1, vs+1]; ##P(DD, ns1, ns2, ns)
          norm=norm+Psns[vns1+1, vns2+1, vs+1]  ##P(DD)
        }else  {Psns[vns1+1, vns2+1, vs+1]=-1;}
      }}}
  #Expectation of shared and non-shared variants
  vncount <- 0;vscount <-0;
  for(vns1 in 0:2){
    for(vns2 in 0:2){
      for(vs in 0:2){
        if(vns1+vs<3 & vns2+vs<3)
        {
          Psns[vns1+1, vns2+1, vs+1]=Psns[vns1+1, vns2+1, vs+1]/norm;  ##P(vns1, vns2, ns|DD)
          vncount=vncount+(vns1+vns2)*Psns[vns1+1, vns2+1, vs+1]; ##E(vns1+vns2|DD)
          vscount=vscount+vs*Psns[vns1+1, vns2+1, vs+1]; ##E(vs|DD)
          
          #print "$vns1 $vns2 $vs $Psns[vns1+1, vns2+1, vs+1]\n"
        }
      }}}
  
  
  ##calculate expectation of S given DD
  #penetrance of P(DD|A,B,S)
  Pab_S <- array(0, c(3,3,3))
  for(S in 0:2) {
    for(a in 0:2){
      for(b in 0:2){
        for(i in 0:2){
          for(j in 0:2){
            Pab_S[a+1, b+1, S+1]=Pab_S[a+1, b+1, S+1]+Pgt[i+1, j+1]*Ptr[a+1, i+1, S+1, 0+1]*(PofS[0+1]*Ptr[b+1, j+1, 0+1, 1+1]+PofS[1+1]*Ptr[b+1, j+1, 1+1, 1+1]+PofS[2+1]*Ptr[b+1, j+1, 2+1, 1+1])
          }
        }
        Pab_S[a+1, b+1, S+1]=Pab_S[a+1, b+1, S+1]*Pgt[a+1, b+1]
      }
    }
  }
  
  
  #Calculate P(AA|vns1, vns2,vs, S) summing over P(b,b')
  PAAs_S <- array(0, c(3,3,3,3));
  for(S in 0:2) {
    for(vns1 in 0:2){
      for(vns2 in 0:2){
        for(vs in 0:2){
          if(vns1+vs<3 & vns2+vs<3)
          {
            for(b in 0:2){
              for(i in 0:2){
                
                ah=vns1+vs;
                ah2=vns2+vs;
                p_sib=0;
                for(S2 in 0:2){
                  p_sib=p_sib+PofS[S2+1]*Ptr[b+1, i+1, S2+1, 1+1]*dbinom(b,2,p[1+1]); ##contribution from loci B
                  #print "tt $S $PofS[$S] $Ptr[$b, $i, $S, 1]\n";
                }
                #print "$b $i $p_sib\n";
                PAAs_S[vns1+1, vns2+1, vs+1, S+1]=PAAs_S[vns1+1, vns2+1, vs+1, S+1]+Pgt[ah+1, b+1]*Pgt[ah2+1, i+1]*p_sib;
              }}
            #print "test $vns1 $vns2 $vs $PAAs_S[vns1+1, vns2+1, vs+1]\n";
          }
          else {PAAs_S[vns1+1, vns2+1, vs+1, S+1]=-1;}
        }}}
  }
  
  #P(S|AA)
  Ps <- array(0, 3)
  norm <- 0;
  for(S in 0:2){
    for(vns1 in 0:2){
      for(vns2 in 0:2){
        for(vs in 0:2){
          if(vns1+vs<3 & vns2+vs<3)
          {
            Ps[S+1]=Ps[S+1]+PAAs_S[vns1+1, vns2+1, vs+1, S+1]*dbinom(vs,S,p[0+1])*dbinom(vns1,2-S,p[0+1])*dbinom(vns2,2-S,p[0+1])*PofS[S+1];##P(DD, ns1, ns2, ns)
          }
        }}}
    norm=norm+Ps[S+1]  ##P(DD)
  }
  #Expectation of shared and non-shared chromosome
  ES <- 0;ENS <-0;
  for(S in 0:2){
    Ps[S+1]=Ps[S+1]/norm;  ##P(S|DD)
    ES=ES+S*Ps[S+1]; ##E(vns1+vns2|DD)
    ENS=ENS+(2-S)*2*Ps[S+1]; ##E(vs|DD)
  }
  return(list(EH=c(vscount/ES, vncount/ENS), E_c=c(ES, ENS, ES+ENS), alpha=alpha, beta=beta, prevalence=prev, srr=srr, Ex_sib=Ex_sib, Ex_unrelated=Ex_unrelated))
}

power_inter.sas <- function(p=c(0.01, 0.05), trr_alpha=2, trr_beta=6.2, gamma=1, n_sb=50, alpha=10^-6) {
  result <- inter_model(p, trr_alpha, trr_beta, gamma)
  e_c <- result$E_c  ## expected no. of chromosomes per sibpair
  N <- n_sb*e_c[3] #no. of  total chromosomes
  n1 <- n_sb*e_c[1]  ##no. of shared chromosomes
  n2 <- n_sb*e_c[2]  #no. of non-shared chromosomes
  
  P_A <- result$EH ##proportion under HA
  p1 <- P_A[1]
  p2 <- P_A[2]
  p0 <- 0
  w1 <- n1/N
  w2 <- n2/N
  
  pnorm( ((p2-p1-p0)*(N*w1*w2)^.5 -  qnorm(1-alpha/2)*((w1*p1+w2*p2)*(1-w1*p1-w2*p2))^.5)/(w2*p1*(1-p1)+w1*p2*(1-p2))^.5 ) +
    pnorm( (-(p2-p1-p0)*(N*w1*w2)^.5 -  qnorm(1-alpha/2)*((w1*p1+w2*p2)*(1-w1*p1-w2*p2))^.5)/(w2*p1*(1-p1)+w1*p2*(1-p2))^.5) 
}
power_inter.sas(p=c(0.01, 0.05), trr_alpha=2, trr_beta=6.2, gamma=1, n_sb=1000, alpha=10^-6)

power_inter_comp.sas <- function(p=c(0.01, 0.05), trr_alpha=2, trr_beta=6.2, gamma=1, n_sb=50, alpha=10^-6, out="A") {
  result <- inter_model(p, trr_alpha, trr_beta, gamma)
  N <- n_sb*4 #no. of  total chromosomes
  n1 <- n_sb*2  ##no. of case chromosomes
  n2 <- n_sb*2  #no. of control chromosomes
  
  P_A <- ifelse(out=="A", result$Ex_unrelated / 2,  result$Ex_sib / 2)##proportion under HA
  p1 <- P_A
  p2 <- p[1]
  p0 <- 0
  w1 <- n1/N
  w2 <- n2/N
  
  pnorm( ((p2-p1-p0)*(N*w1*w2)^.5 -  qnorm(1-alpha/2)*((w1*p1+w2*p2)*(1-w1*p1-w2*p2))^.5)/(w2*p1*(1-p1)+w1*p2*(1-p2))^.5 ) +
    pnorm( (-(p2-p1-p0)*(N*w1*w2)^.5 -  qnorm(1-alpha/2)*((w1*p1+w2*p2)*(1-w1*p1-w2*p2))^.5)/(w2*p1*(1-p1)+w1*p2*(1-p2))^.5) 
}
power_inter_comp.sas(p=c(0.01, 0.05), trr_alpha=2, trr_beta=6.2, gamma=1, n_sb=1000, alpha=10^-6)


###########
#only consider 1 variant in founders
power.family <- function(fam.str=c(2,2), affect=c(c(0,1),c(1,1)), IBD=c(c(1,3), c(3, 2)), n_family=500, f=0.01, mu=1, alpha=10^-6) {
  P_A_F <- function(F, mu) {
    ifelse(F==1, mu, 1)
  }
  P_AA_H <- function(HS1, HS2, mu) {
    mu^HS1*(mu^2)^HS2
  }
  
  #specify the family structure
  n_I <-  fam.str[1]
  n_II <- fam.str[2]
  
  #affected status
  affect <- affect  
  
  #specify the IBD status
  IBD <- IBD
  
  
  #test statistics and calculate the power
  sce <- matrix(c(1,0,0,0,
                  0,1,0,0,
                  0,0,1,0,
                  0,0,0,1), nrow=4)
  p_A_sum=0
  p_A <- rep(NA, 4)
  p_N_sum=0
  p_N <- rep(NA, 4)
  S <- rep(NA, 4)
  for(i in 1:4) {
    HS1 <- HS2 <- 0
    F1 <- sum(sce[i,1:2])
    F2 <- sum(sce[i,3:4])
    HS1 <- ifelse(length(which(IBD==which(sce[i,]==1)))==1, 1, 0)
    HS2 <- ifelse(length(which(IBD==which(sce[i,]==1)))==2, 1, 0)
    S[i] <- F1*affect[1] + F2*affect[2] + HS1 + HS2*2
    p_A[i] <- P_A_F(F2, mu)*P_AA_H(HS1,HS2, mu) #(1/2)^HS1
    p_N[i] <- P_A_F(F2, 1)*P_AA_H(HS1,HS2, 1)
    p_A_sum <- p_A_sum+p_A[i]
    p_N_sum <- p_N_sum+p_N[i]
  } 
  E_A <- sum(S*p_A/p_A_sum)
  E_N <- sum(S*p_N/p_N_sum)
  sigma2_N <- sum((S-E_N)^2*p_N/p_N_sum)
  sigma2_A <- sum((S-E_A)^2*p_A/p_A_sum)
  
  crit.pt.L <- qnorm(alpha/2, mean=E_N, sd=sqrt(sigma2_N/n_family))
  crit.pt.H <- qnorm(1-alpha/2, mean=E_N, sd=sqrt(sigma2_N/n_family))
  power <- pnorm(crit.pt.L, mean=E_A, sd=sqrt(sigma2_A/n_family)) + pnorm(crit.pt.H, mean=E_A, sd=sqrt(sigma2_A/n_family), lower=F)
  #c(S, p_A/p_A_sum, E_A, E_N, sigma2_A, sigma2_N, power)  
  power
}
power.family(mu=2)

##consider 0,1,2,3,4 variants in founders
power.family.all <- function(fam.str=c(2,2), affect=c(c(0,1),c(1,1)), n_family=500, f=0.01, mu=1, alpha=10^-6) {
  P_A_F <- function(F, mu) { #prob of being affected for founder
    ifelse(F==1, mu, 1)
  }
  P_AA_H <- function(HS1, HS2, mu, Sharing) { #prob of being double affected siblings
    if(Sharing==1)  result <- mu^HS1*(mu^2)^HS2
    if(Sharing==2)  result <- (mu^2)^HS2
    result
  }
  
  #specify the family structure
  n_I <-  fam.str[1] #no. of the first generation
  n_II <- fam.str[2] #no. of the second generation
  
  #affected status
  affect <- affect  
    
  #initialization
  p_A_sum=0 #sum of prob. under the alternative
  p_N_sum=0 #sum of prob. under the null
  E_A=0 #expectation under the alternative
  E_N=0 #expectation under the null
  p_S_A=rep(0,7) #prob. of S under the alternative
  p_S_N=rep(0,7) #prob. of S under the null S=0 is places at [1]
  p_IBD = c(0.25, 0.5, 0.25) #PROB.of IBD
  # update <- function() { #update the sum
  #  S <- h1*(1 %in% IBD_Str[i,])*affect[1] + h2*(2 %in% IBD_Str[i,])*affect[1] + h3*(3 %in% IBD_Str[i,])*affect[2] + h4*(4 %in% IBD_Str[i,])*affect[2] + HS1 + HS2*2
    #p_A <- P_A_F(F2, mu)*P_AA_H(HS1,HS2, mu) #(1/2)^HS1
    #p_N <- P_A_F(F2, 1)*P_AA_H(HS1,HS2, 1)
  #  p_S_A[S+1] <<- p_S_A[S+1] + p_A 
  #  p_S_N[S+1] <<- p_S_N[S+1] + p_N
  #  p_A_sum <<- p_A_sum+p_A
  #  p_N_sum <<- p_N_sum+p_N
    
  #  E_A <<- E_A + S*p_A
  #  E_N <<- E_N + S*p_N
  #}
  for(IBD in 0:2) {
    for(h1 in 0:1) { #father's first haplotype
      for(h2 in 0:1) { #father's second haplotype
        for(h3 in 0:1) { #mother's first haplotype
          for(h4 in 0:1) { #mother's second haplotype
            #initialization
            sce <- c(h1,h2,h3,h4) #founder's haplotype
            p_A <- 0
            p_N <- 0
            S <- NA
            
            HS1 <- HS2 <- 0 #shared variant
            F1 <- sum(sce[1:2]) #founder: father
            F2 <- sum(sce[3:4]) #founder: mothrt
            FB <- sum(sce) #founder's haplotype
            p_f <- f^(FB)*(1-f)^(4-FB) #allele frequency of founders
            
            if(IBD==0) {
              p_A <- P_A_F(F2, mu)*mu^(FB)*p_f*p_IBD[IBD+1] #(1/2)^HS1
              p_N <- P_A_F(F2, 1)*1*p_f*p_IBD[IBD+1]
              
              S <- h1*affect[1] + h2*affect[1] + h3*affect[2] + h4*affect[2] + FB
              #p_A <- P_A_F(F2, mu)*P_AA_H(HS1,HS2, mu) #(1/2)^HS1
              #p_N <- P_A_F(F2, 1)*P_AA_H(HS1,HS2, 1)
              p_S_A[S+1] <- p_S_A[S+1] + p_A 
              p_S_N[S+1] <- p_S_N[S+1] + p_N
              p_A_sum <- p_A_sum+p_A
              p_N_sum <- p_N_sum+p_N
              
              E_A <- E_A + S*p_A
              E_N <- E_N + S*p_N
            }
            #print(IBD)
            if(IBD==1) {
              IBD_Str=matrix(c(1,3,1,4,
                               2,3,2,4,
                               3,1,3,2,
                               4,1,4,2), byrow=T, nrow=4)
              for(i in 1:4) {
                HS1 <- sum(table(match(IBD_Str[i,], which(sce==1)))==1) #how many variant shared once
                HS2 <- sum(table(match(IBD_Str[i,], which(sce==1)))==2) #how many variant shared twice
                p_A <- P_A_F(F2, mu)*P_AA_H(HS1,HS2, mu,IBD)*p_f*p_IBD[IBD+1]*.25 #(1/2)^HS1
                p_N <- P_A_F(F2, 1)*1*p_f*p_IBD[IBD+1]*.25
                
                S <- h1*(1 %in% IBD_Str[i,])*affect[1] + h2*(2 %in% IBD_Str[i,])*affect[1] + h3*(3 %in% IBD_Str[i,])*affect[2] + h4*(4 %in% IBD_Str[i,])*affect[2] + HS1 + HS2*2
                #p_A <- P_A_F(F2, mu)*P_AA_H(HS1,HS2, mu) #(1/2)^HS1
                #p_N <- P_A_F(F2, 1)*P_AA_H(HS1,HS2, 1)
                p_S_A[S+1] <- p_S_A[S+1] + p_A 
                p_S_N[S+1] <- p_S_N[S+1] + p_N
                p_A_sum <- p_A_sum+p_A
                p_N_sum <- p_N_sum+p_N
                
                E_A <- E_A + S*p_A
                E_N <- E_N + S*p_N
              }
            }
            if(IBD==2) {
              IBD_Str=matrix(c(1,3,1,3,
                               1,4,1,4,
                               2,3,2,3,
                               2,4,2,4), byrow=T, nrow=4)
              for(i in 1:4) {
                HS1 <- sum(table(match(IBD_Str[i,], which(sce==1)))==1)
                HS2 <- sum(table(match(IBD_Str[i,], which(sce==1)))==2)
                p_A <- P_A_F(F2, mu)*P_AA_H(HS1,HS2, mu,IBD)*p_f*p_IBD[IBD+1]*.25 #(1/2)^HS1
                p_N <- P_A_F(F2, 1)*1*p_f*p_IBD[IBD+1]*.25
                
                S <- h1*(1 %in% IBD_Str[i,])*affect[1] + h2*(2 %in% IBD_Str[i,])*affect[1] + h3*(3 %in% IBD_Str[i,])*affect[2] + h4*(4 %in% IBD_Str[i,])*affect[2] + HS1 + HS2*2
                #p_A <- P_A_F(F2, mu)*P_AA_H(HS1,HS2, mu) #(1/2)^HS1
                #p_N <- P_A_F(F2, 1)*P_AA_H(HS1,HS2, 1)
                p_S_A[S+1] <- p_S_A[S+1] + p_A 
                p_S_N[S+1] <- p_S_N[S+1] + p_N
                p_A_sum <- p_A_sum+p_A
                p_N_sum <- p_N_sum+p_N
                
                E_A <- E_A + S*p_A
                E_N <- E_N + S*p_N
              }
            }
          }
        }
      }
    }
  }
  
  E_A <- E_A/p_A_sum #normalized by sum
  E_N <- E_N/p_N_sum #normalized by sum
  p_S_A <- p_S_A/p_A_sum #normalized by sum
  p_S_N <- p_S_N/p_N_sum #normalized by sum
  
  S=0:6
  sigma2_N <- sum((S-E_N)^2*p_S_N[S+1]) #variance under the null
  sigma2_A <- sum((S-E_A)^2*p_S_A[S+1]) #variance under the alternative
  crit.pt.L <- qnorm(alpha/2, mean=E_N, sd=sqrt(sigma2_N/n_family))
  crit.pt.H <- qnorm(1-alpha/2, mean=E_N, sd=sqrt(sigma2_N/n_family))
  power <- pnorm(crit.pt.L, mean=E_A, sd=sqrt(sigma2_A/n_family)) + pnorm(crit.pt.H, mean=E_A, sd=sqrt(sigma2_A/n_family), lower=F)
  #c(S, p_A/p_A_sum, E_A, E_N, sigma2_A, sigma2_N, power)  
  power
}
power.family.all(mu=1)

##consider 0,1,2,3,4 variants in founders and dump 0 variant families
power.family.all.2 <- function(fam.str=c(2,2), affect=c(c(0,1),c(1,1)), n_family=500, f=0.01, mu=1, alpha=10^-6) {
  P_A_F <- function(F, mu) { #prob of being affected for founder
    ifelse(F==1, mu, 1)
  }
  P_AA_H <- function(HS1, HS2, mu, Sharing) { #prob of being double affected siblings
    if(Sharing==1)  result <- mu^HS1*(mu^2)^HS2
    if(Sharing==2)  result <- (mu^2)^HS2
    result
  }
  
  #specify the family structure
  n_I <-  fam.str[1] #no. of the first generation
  n_II <- fam.str[2] #no. of the second generation
  
  #affected status
  affect <- affect  
  
  #initialization
  p_A_sum=0 #sum of prob. under the alternative
  p_N_sum=0 #sum of prob. under the null
  E_A=0 #expectation under the alternative
  E_N=0 #expectation under the null
  p_S_A=rep(0,7) #prob. of S under the alternative
  p_S_N=rep(0,7) #prob. of S under the null S=0 is places at [1]
  p_IBD = c(0.25, 0.5, 0.25) #PROB.of IBD
  # update <- function() { #update the sum
  #  S <- h1*(1 %in% IBD_Str[i,])*affect[1] + h2*(2 %in% IBD_Str[i,])*affect[1] + h3*(3 %in% IBD_Str[i,])*affect[2] + h4*(4 %in% IBD_Str[i,])*affect[2] + HS1 + HS2*2
  #p_A <- P_A_F(F2, mu)*P_AA_H(HS1,HS2, mu) #(1/2)^HS1
  #p_N <- P_A_F(F2, 1)*P_AA_H(HS1,HS2, 1)
  #  p_S_A[S+1] <<- p_S_A[S+1] + p_A 
  #  p_S_N[S+1] <<- p_S_N[S+1] + p_N
  #  p_A_sum <<- p_A_sum+p_A
  #  p_N_sum <<- p_N_sum+p_N
  
  #  E_A <<- E_A + S*p_A
  #  E_N <<- E_N + S*p_N
  #}
  for(IBD in 0:2) {
    for(h1 in 0:1) { #father's first haplotype
      for(h2 in 0:1) { #father's second haplotype
        for(h3 in 0:1) { #mother's first haplotype
          for(h4 in 0:1) { #mother's second haplotype
            #initialization
            sce <- c(h1,h2,h3,h4) #founder's haplotype
            p_A <- 0
            p_N <- 0
            S <- NA
            
            HS1 <- HS2 <- 0 #shared variant
            F1 <- sum(sce[1:2]) #founder: father
            F2 <- sum(sce[3:4]) #founder: mothrt
            FB <- sum(sce) #founder's haplotype
            p_f <- f^(FB)*(1-f)^(4-FB) #allele frequency of founders
            
            if(IBD==0) {
              p_A <- P_A_F(F2, mu)*mu^(FB)*p_f*p_IBD[IBD+1] #(1/2)^HS1
              p_N <- P_A_F(F2, 1)*1*p_f*p_IBD[IBD+1]
              
              S <- h1*affect[1] + h2*affect[1] + h3*affect[2] + h4*affect[2] + FB
              #p_A <- P_A_F(F2, mu)*P_AA_H(HS1,HS2, mu) #(1/2)^HS1
              #p_N <- P_A_F(F2, 1)*P_AA_H(HS1,HS2, 1)
              p_S_A[S+1] <- p_S_A[S+1] + p_A 
              p_S_N[S+1] <- p_S_N[S+1] + p_N
              p_A_sum <- p_A_sum+p_A
              p_N_sum <- p_N_sum+p_N
              
              E_A <- E_A + S*p_A
              E_N <- E_N + S*p_N
            }
            #print(IBD)
            if(IBD==1) {
              IBD_Str=matrix(c(1,3,1,4,
                               2,3,2,4,
                               3,1,3,2,
                               4,1,4,2), byrow=T, nrow=4)
              for(i in 1:4) {
                HS1 <- sum(table(match(IBD_Str[i,], which(sce==1)))==1) #how many variant shared once
                HS2 <- sum(table(match(IBD_Str[i,], which(sce==1)))==2) #how many variant shared twice
                p_A <- P_A_F(F2, mu)*P_AA_H(HS1,HS2, mu,IBD)*p_f*p_IBD[IBD+1]*.25 #(1/2)^HS1
                p_N <- P_A_F(F2, 1)*1*p_f*p_IBD[IBD+1]*.25
                
                S <- h1*(1 %in% IBD_Str[i,])*affect[1] + h2*(2 %in% IBD_Str[i,])*affect[1] + h3*(3 %in% IBD_Str[i,])*affect[2] + h4*(4 %in% IBD_Str[i,])*affect[2] + HS1 + HS2*2
                #p_A <- P_A_F(F2, mu)*P_AA_H(HS1,HS2, mu) #(1/2)^HS1
                #p_N <- P_A_F(F2, 1)*P_AA_H(HS1,HS2, 1)
                p_S_A[S+1] <- p_S_A[S+1] + p_A 
                p_S_N[S+1] <- p_S_N[S+1] + p_N
                p_A_sum <- p_A_sum+p_A
                p_N_sum <- p_N_sum+p_N
                
                E_A <- E_A + S*p_A
                E_N <- E_N + S*p_N
              }
            }
            if(IBD==2) {
              IBD_Str=matrix(c(1,3,1,3,
                               1,4,1,4,
                               2,3,2,3,
                               2,4,2,4), byrow=T, nrow=4)
              for(i in 1:4) {
                HS1 <- sum(table(match(IBD_Str[i,], which(sce==1)))==1)
                HS2 <- sum(table(match(IBD_Str[i,], which(sce==1)))==2)
                p_A <- P_A_F(F2, mu)*P_AA_H(HS1,HS2, mu,IBD)*p_f*p_IBD[IBD+1]*.25 #(1/2)^HS1
                p_N <- P_A_F(F2, 1)*1*p_f*p_IBD[IBD+1]*.25
                
                S <- h1*(1 %in% IBD_Str[i,])*affect[1] + h2*(2 %in% IBD_Str[i,])*affect[1] + h3*(3 %in% IBD_Str[i,])*affect[2] + h4*(4 %in% IBD_Str[i,])*affect[2] + HS1 + HS2*2
                #p_A <- P_A_F(F2, mu)*P_AA_H(HS1,HS2, mu) #(1/2)^HS1
                #p_N <- P_A_F(F2, 1)*P_AA_H(HS1,HS2, 1)
                p_S_A[S+1] <- p_S_A[S+1] + p_A 
                p_S_N[S+1] <- p_S_N[S+1] + p_N
                p_A_sum <- p_A_sum+p_A
                p_N_sum <- p_N_sum+p_N
                
                E_A <- E_A + S*p_A
                E_N <- E_N + S*p_N
              }
            }
          }
        }
      }
    }
  }
  
  E_A <- E_A/p_A_sum #normalized by sum
  E_N <- E_N/p_N_sum #normalized by sum
  p_S_A <- p_S_A/p_A_sum #normalized by sum
  p_S_N <- p_S_N/p_N_sum #normalized by sum
  
  S=1:6
  sigma2_N <- sum((S-sum(S*p_S_N[-1]/sum(p_S_N[-1])))^2*p_S_N[-1]/sum(p_S_N[-1])) #variance under the null
  sigma2_A <- sum((S-sum(S*p_S_A[-1]/sum(p_S_N[-1])))^2*p_S_A[-1]/sum(p_S_A[-1])) #variance under the alternative
  c(E_A, E_N, sum(S*p_S_A[-1]/sum(p_S_A[-1])), sum(S*p_S_N[-1]/sum(p_S_N[-1])), sigma2_N, sigma2_A)
  crit.pt.L <- qnorm(alpha/2, mean=E_N, sd=sqrt(sigma2_N/(sum(p_S_N[-1])/sum(p_S_N)*n_family)))
  crit.pt.H <- qnorm(1-alpha/2, mean=E_N, sd=sqrt(sigma2_N/(sum(p_S_N[-1])/sum(p_S_N)*n_family)))
  power <- pnorm(crit.pt.L, mean=E_A, sd=sqrt(sigma2_A/(sum(p_S_A[-1])/sum(p_S_A)*n_family))) + pnorm(crit.pt.H, mean=E_A, sd=sqrt(sigma2_A/(sum(p_S_A[-1])/sum(p_S_A)*n_family)), lower=F)
  #c(S, p_A/p_A_sum, E_A, E_N, sigma2_A, sigma2_N, power)  
  print(power)
}
power.family.all.2(mu=1)

#calculate only 0 and 1 variant in founders and dump 0
power.family.2 <- function(fam.str=c(2,2), affect=c(c(0,1),c(1,1)), n_family=500, f=0.01, mu=1, p_dis=0.01, alpha=10^-6) {
  #calculate KG
  KL <- (1+f*(mu-1))^2
  KG=p_dis/KL

  P_A_F <- function(F, A, mu) { #prob of being affected for founder
    if(A==1 & F==1) return(mu)
    if(A==0 & F==1) return(1-mu*KG)
    return(1)
  }
  
  P_AA_H <- function(HS1, HS2, AA, mu, sharing_IBD) { #prob of being double affected siblings
    if(sum(AA)==2) {
      if(sharing_IBD==0)  result <- mu^HS1
      if(sharing_IBD==1)  result <- mu^HS1*(mu^2)^HS2
      if(sharing_IBD==2)  result <- (mu^2)^HS2
    }
    if(sum(AA)==0) {
      if(sharing_IBD==0)  result <- (1-mu*KG)^HS1
      if(sharing_IBD==1)  result <- (1-mu*KG)^HS1*((1-mu*KG)^2)^HS2
      if(sharing_IBD==2)  result <- ((1-mu*KG)^2)^HS2
    }
    result
  }
  
  #specify the family structure
  n_I <-  fam.str[1] #no. of the first generation
  n_II <- fam.str[2] #no. of the second generation
  
  #affected status
  affect <- affect  
  
  #initialization
  p_A_sum=0 #sum of prob. under the alternative
  p_N_sum=0 #sum of prob. under the null
  p_A_n0_sum=0 #sum of prob. under the alternative non 0 in founders
  p_N_n0_sum=0 #sum of prob. under the null non 0 in founders
  E_A=0 #expectation under the alternative
  E_N=0 #expectation under the null
  p_S_A=rep(0,7) #prob. of S under the alternative
  p_S_N=rep(0,7) #prob. of S under the null S=0 is places at [1]
  p_IBD = c(0.25, 0.5, 0.25) #PROB.of IBD
  # update <- function() { #update the sum
  #  S <- h1*(1 %in% IBD_Str[i,])*affect[1] + h2*(2 %in% IBD_Str[i,])*affect[1] + h3*(3 %in% IBD_Str[i,])*affect[2] + h4*(4 %in% IBD_Str[i,])*affect[2] + HS1 + HS2*2
  #p_A <- P_A_F(F2, mu)*P_AA_H(HS1,HS2, mu) #(1/2)^HS1
  #p_N <- P_A_F(F2, 1)*P_AA_H(HS1,HS2, 1)
  #  p_S_A[S+1] <<- p_S_A[S+1] + p_A 
  #  p_S_N[S+1] <<- p_S_N[S+1] + p_N
  #  p_A_sum <<- p_A_sum+p_A
  #  p_N_sum <<- p_N_sum+p_N
  
  #  E_A <<- E_A + S*p_A
  #  E_N <<- E_N + S*p_N
  #}
  for(IBD in 0:2) {
    founder <- matrix(c(0,0,0,0,
                    1,0,0,0,
                    0,1,0,0,
                    0,0,1,0,
                    0,0,0,1), nrow=5, byrow=T)
    for(h in 1:nrow(founder)) {
      #initialization
      sce <- as.vector(founder[h,]) #founder's haplotype
      h1 <- sce[1]
      h2 <- sce[2]
      h3 <- sce[3]
      h4 <- sce[4]
      p_A <- 0
      p_N <- 0
      S <- NA
      
      HS1 <- HS2 <- 0 #shared variant
      F1 <- sum(sce[1:2]) #founder: father
      F2 <- sum(sce[3:4]) #founder: mothrt
      FB <- sum(sce) #founder's haplotype
      p_f <- f^(FB)*(1-f)^(4-FB) #allele frequency of founders
      
      if(IBD==0) {
        HS1 <- FB
        HS2 <- 0
        p_A <- P_A_F(F1, affect[1], mu)*P_A_F(F2, affect[2], mu)*P_AA_H(HS1,HS2,AA=affect[3:4], mu,IBD)*p_f*p_IBD[IBD+1] #(1/2)^HS1
        p_N <- P_A_F(F1, affect[1], 1)*P_A_F(F2, affect[2], 1)*P_AA_H(HS1,HS2,AA=affect[3:4], 1,IBD)*p_f*p_IBD[IBD+1]
        
        S <- h1*affect[1] + h2*affect[1] + h3*affect[2] + h4*affect[2] + FB
        #p_A <- P_A_F(F2, mu)*P_AA_H(HS1,HS2, mu) #(1/2)^HS1
        #p_N <- P_A_F(F2, 1)*P_AA_H(HS1,HS2, 1)
        p_A_sum <- p_A_sum+p_A
        p_N_sum <- p_N_sum+p_N
        if(h!=1) {
          p_S_A[S+1] <- p_S_A[S+1] + p_A 
          p_S_N[S+1] <- p_S_N[S+1] + p_N
          
          p_A_n0_sum <- p_A_n0_sum + p_A
          E_A <- E_A + S*p_A
          
          p_N_n0_sum <- p_N_n0_sum + p_N
          E_N <- E_N + S*p_N
        }
      }
      #print(IBD)
      if(IBD==1) {
        IBD_Str=matrix(c(1,3,1,4,
                         2,3,2,4,
                         3,1,3,2,
                         4,1,4,2), byrow=T, nrow=4)
        for(i in 1:4) {
          HS1 <- sum(table(match(IBD_Str[i,], which(sce==1)))==1) #how many variant shared once
          HS2 <- sum(table(match(IBD_Str[i,], which(sce==1)))==2) #how many variant shared twice
          p_A <- P_A_F(F1, affect[1], mu)*P_A_F(F2, affect[2], mu)*P_AA_H(HS1,HS2,AA=affect[3:4], mu,IBD)*p_f*p_IBD[IBD+1]*.25 #(1/2)^HS1
          p_N <- P_A_F(F1, affect[1], 1)*P_A_F(F2, affect[2], 1)*P_AA_H(HS1,HS2,AA=affect[3:4], 1,IBD)*p_f*p_IBD[IBD+1]*.25
          
          S <- h1*affect[1] + h2*affect[1] + h3*affect[2] + h4*affect[2] + HS1 + HS2*2
          #p_A <- P_A_F(F2, mu)*P_AA_H(HS1,HS2, mu) #(1/2)^HS1
          #p_N <- P_A_F(F2, 1)*P_AA_H(HS1,HS2, 1)
          p_A_sum <- p_A_sum+p_A
          p_N_sum <- p_N_sum+p_N
          if(h!=1) {
            p_S_A[S+1] <- p_S_A[S+1] + p_A 
            p_S_N[S+1] <- p_S_N[S+1] + p_N
            
            p_A_n0_sum <- p_A_n0_sum + p_A
            E_A <- E_A + S*p_A
            
            p_N_n0_sum <- p_N_n0_sum + p_N
            E_N <- E_N + S*p_N
          }
        }
      }
      if(IBD==2) {
        IBD_Str=matrix(c(1,3,1,3,
                         1,4,1,4,
                         2,3,2,3,
                         2,4,2,4), byrow=T, nrow=4)
        for(i in 1:4) {
          HS1 <- sum(table(match(IBD_Str[i,], which(sce==1)))==1) #how many variant shared once
          HS2 <- sum(table(match(IBD_Str[i,], which(sce==1)))==2) #how many variant shared twice
          p_A <- P_A_F(F1, affect[1], mu)*P_A_F(F2, affect[2], mu)*P_AA_H(HS1,HS2,AA=affect[3:4], mu,IBD)*p_f*p_IBD[IBD+1]*.25 #(1/2)^HS1
          p_N <- P_A_F(F1, affect[1], 1)*P_A_F(F2, affect[2], 1)*P_AA_H(HS1,HS2,AA=affect[3:4], 1,IBD)*p_f*p_IBD[IBD+1]*.25
          
          S <- h1*affect[1] + h2*affect[1] + h3*affect[2] + h4*affect[2] + HS1 + HS2*2
          #p_A <- P_A_F(F2, mu)*P_AA_H(HS1,HS2, mu) #(1/2)^HS1
          #p_N <- P_A_F(F2, 1)*P_AA_H(HS1,HS2, 1)
          p_A_sum <- p_A_sum+p_A
          p_N_sum <- p_N_sum+p_N
          if(h!=1) {
            p_S_A[S+1] <- p_S_A[S+1] + p_A 
            p_S_N[S+1] <- p_S_N[S+1] + p_N
            
            p_A_n0_sum <- p_A_n0_sum + p_A
            E_A <- E_A + S*p_A
            
            p_N_n0_sum <- p_N_n0_sum + p_N
            E_N <- E_N + S*p_N
          }
        }
      }
    }
  }
  
  #
  p_S_A <- p_S_A/p_A_n0_sum #normalized by sum
  p_S_N <- p_S_N/p_N_n0_sum #normalized by sum
  E_A <- E_A/p_A_n0_sum #normalized by sum
  E_N <- E_N/p_N_n0_sum #normalized by sum

  p_A_n0 <- p_A_n0_sum/p_A_sum
  p_N_n0 <- p_N_n0_sum/p_N_sum
  
  S=0:6
  sigma2_N <- sum((S-E_N)^2*p_S_N) #variance under the null
  sigma2_A <- sum((S-E_A)^2*p_S_A) #variance under the alternative
  print(c(E_A, E_N, p_A_n0, p_N_n0, sigma2_A, sigma2_N))
  print(p_S_A)
  print(p_S_N)
  crit.pt.L <- qnorm(alpha/2, mean=E_N, sd=sqrt(sigma2_N/(p_N_n0*n_family)))
  crit.pt.H <- qnorm(1-alpha/2, mean=E_N, sd=sqrt(sigma2_N/(p_N_n0*n_family)))
  power <- pnorm(crit.pt.L, mean=E_A, sd=sqrt(sigma2_A/(p_A_n0*n_family))) + pnorm(crit.pt.H, mean=E_A, sd=sqrt(sigma2_A/(p_A_n0*n_family)), lower=F)
  #c(S, p_A/p_A_sum, E_A, E_N, sigma2_A, sigma2_N, power)  
  print(power)
}
power.family.2(mu=1)


#calculate only 0,1,2,3,4 variant in founders and dump 0 and 4
power.family.3 <- function(fam.str=c(2,2), affect=c(c(0,1),c(1,1)), n_family=500, f=0.01, mu=1, p_dis=0.01, alpha=10^-6) {
  #calculate KG
  KL <- (1+f*(mu-1))^2
  KG=p_dis/KL
  
  P_A_F <- function(F, A, mu) { #prob of being affected for founder
    if(A==1 & F==1) return(mu)
    if(A==0 & F==1) return(1-mu*KG)
    return(1)
  }
  
  P_AA_H <- function(HS1, HS2, AA, mu, sharing_IBD) { #prob of being double affected siblings
    if(sum(AA)==2) {
      if(sharing_IBD==0)  result <- mu^HS1
      if(sharing_IBD==1)  result <- mu^HS1*(mu^2)^HS2
      if(sharing_IBD==2)  result <- (mu^2)^HS2
    }
    if(sum(AA)==1) {
      if(sharing_IBD==0)  {
        if(HS1==0) result <- 1
        if(HS1==1) result <- 0.5*mu+0.5*(1-mu*KG)
        if(HS1==2) result <- 1/6*mu^2 + 4/6*mu*(1-mu*KG) + 1/6*(1-mu*KG)^2
        if(HS1==3) result <- 2/4*mu^2*(1-mu*KG) + 2/4*mu*(1-mu*KG)^2
        if(HS1==4) result <- mu^2*(1-mu*KG)^2
      }
      if(sharing_IBD==1)  {
        if(HS1==0) result <- (mu*(1-mu*KG))^HS2
        if(HS1==1) result <- 0.5*mu*(mu*(1-mu*KG))^HS2 + 0.5*(1-mu*KG)*(mu*(1-mu*KG))^HS2
        if(HS1==2) result <- mu*(1-mu*KG)*(mu*(1-mu*KG))^HS2
      }
      if(sharing_IBD==2)  {
        result <- (mu*(1-mu*KG))^HS2
      }
    }
    if(sum(AA)==0) {
      if(sharing_IBD==0)  result <- (1-mu*KG)^HS1
      if(sharing_IBD==1)  result <- (1-mu*KG)^HS1*((1-mu*KG)^2)^HS2
      if(sharing_IBD==2)  result <- ((1-mu*KG)^2)^HS2
    }
    result
  }
  
  #specify the family structure
  n_I <-  fam.str[1] #no. of the first generation
  n_II <- fam.str[2] #no. of the second generation
  
  #affected status
  affect <- affect  
  
  #initialization
  p_A_sum=0 #sum of prob. under the alternative
  p_N_sum=0 #sum of prob. under the null
  p_A_n0_sum=0 #sum of prob. under the alternative non 0 in founders
  p_N_n0_sum=0 #sum of prob. under the null non 0 in founders
  E_A=0 #expectation under the alternative
  E_N=0 #expectation under the null
  p_S_A=rep(0,10) #prob. of S under the alternative
  p_S_N=rep(0,10) #prob. of S under the null S=0 is places at [1]
  p_IBD = c(0.25, 0.5, 0.25) #PROB.of IBD
  # update <- function() { #update the sum
  #  S <- h1*(1 %in% IBD_Str[i,])*affect[1] + h2*(2 %in% IBD_Str[i,])*affect[1] + h3*(3 %in% IBD_Str[i,])*affect[2] + h4*(4 %in% IBD_Str[i,])*affect[2] + HS1 + HS2*2
  #p_A <- P_A_F(F2, mu)*P_AA_H(HS1,HS2, mu) #(1/2)^HS1
  #p_N <- P_A_F(F2, 1)*P_AA_H(HS1,HS2, 1)
  #  p_S_A[S+1] <<- p_S_A[S+1] + p_A 
  #  p_S_N[S+1] <<- p_S_N[S+1] + p_N
  #  p_A_sum <<- p_A_sum+p_A
  #  p_N_sum <<- p_N_sum+p_N
  
  #  E_A <<- E_A + S*p_A
  #  E_N <<- E_N + S*p_N
  #}
  for(IBD in 0:2) {
    founder <- matrix(c(0,0,0,0,
                        1,0,0,0,
                        0,1,0,0,
                        0,0,1,0,
                        0,0,0,1,
                        1,1,0,0,
                        1,0,1,0,
                        1,0,0,1,
                        0,1,1,0,
                        0,1,0,1,
                        0,0,1,1,
                        1,1,1,0,
                        1,1,0,1,
                        1,0,1,1,
                        0,1,1,1,
                        1,1,1,1), nrow=16, byrow=T)
    for(h in 1:nrow(founder)) {
      #initialization
      sce <- as.vector(founder[h,]) #founder's haplotype
      h1 <- sce[1]
      h2 <- sce[2]
      h3 <- sce[3]
      h4 <- sce[4]
      p_A <- 0
      p_N <- 0
      S <- NA
      
      HS1 <- HS2 <- 0 #shared variant
      F1 <- sum(sce[1:2]) #founder: father
      F2 <- sum(sce[3:4]) #founder: mothrt
      FB <- sum(sce) #founder's haplotype
      p_f <- f^(FB)*(1-f)^(4-FB) #allele frequency of founders
      
      if(IBD==0) {
        HS1 <- FB
        HS2 <- 0
        p_A <- P_A_F(F1, affect[1], mu)*P_A_F(F2, affect[2], mu)*P_AA_H(HS1,HS2,AA=affect[3:4], mu,IBD)*p_f*p_IBD[IBD+1] #(1/2)^HS1
        p_N <- P_A_F(F1, affect[1], 1)*P_A_F(F2, affect[2], 1)*P_AA_H(HS1,HS2,AA=affect[3:4], 1,IBD)*p_f*p_IBD[IBD+1]
        
        S <- h1*affect[1] + h2*affect[1] + h3*affect[2] + h4*affect[2] + FB
        #p_A <- P_A_F(F2, mu)*P_AA_H(HS1,HS2, mu) #(1/2)^HS1
        #p_N <- P_A_F(F2, 1)*P_AA_H(HS1,HS2, 1)
        p_A_sum <- p_A_sum+p_A
        p_N_sum <- p_N_sum+p_N
        if(h!=1 & h!=16) {
          p_S_A[S+1] <- p_S_A[S+1] + p_A 
          p_S_N[S+1] <- p_S_N[S+1] + p_N
          
          p_A_n0_sum <- p_A_n0_sum + p_A
          E_A <- E_A + S*p_A
          
          p_N_n0_sum <- p_N_n0_sum + p_N
          E_N <- E_N + S*p_N
        }
      }
      #print(IBD)
      if(IBD==1) {
        IBD_Str=matrix(c(1,3,1,4,
                         2,3,2,4,
                         3,1,3,2,
                         4,1,4,2), byrow=T, nrow=4)
        for(i in 1:4) {
          HS1 <- sum(table(match(IBD_Str[i,], which(sce==1)))==1) #how many variant shared once
          HS2 <- sum(table(match(IBD_Str[i,], which(sce==1)))==2) #how many variant shared twice
          p_A <- P_A_F(F1, affect[1], mu)*P_A_F(F2, affect[2], mu)*P_AA_H(HS1,HS2,AA=affect[3:4], mu,IBD)*p_f*p_IBD[IBD+1]*.25 #(1/2)^HS1
          p_N <- P_A_F(F1, affect[1], 1)*P_A_F(F2, affect[2], 1)*P_AA_H(HS1,HS2,AA=affect[3:4], 1,IBD)*p_f*p_IBD[IBD+1]*.25
          
          S <- h1*affect[1] + h2*affect[1] + h3*affect[2] + h4*affect[2] + HS1 + HS2*2
          #p_A <- P_A_F(F2, mu)*P_AA_H(HS1,HS2, mu) #(1/2)^HS1
          #p_N <- P_A_F(F2, 1)*P_AA_H(HS1,HS2, 1)
          p_A_sum <- p_A_sum+p_A
          p_N_sum <- p_N_sum+p_N
          if(h!=1 & h!=16) {
            p_S_A[S+1] <- p_S_A[S+1] + p_A 
            p_S_N[S+1] <- p_S_N[S+1] + p_N
            
            p_A_n0_sum <- p_A_n0_sum + p_A
            E_A <- E_A + S*p_A
            
            p_N_n0_sum <- p_N_n0_sum + p_N
            E_N <- E_N + S*p_N
          }
        }
      }
      if(IBD==2) {
        IBD_Str=matrix(c(1,3,1,3,
                         1,4,1,4,
                         2,3,2,3,
                         2,4,2,4), byrow=T, nrow=4)
        for(i in 1:4) {
          HS1 <- sum(table(match(IBD_Str[i,], which(sce==1)))==1) #how many variant shared once
          HS2 <- sum(table(match(IBD_Str[i,], which(sce==1)))==2) #how many variant shared twice
          p_A <- P_A_F(F1, affect[1], mu)*P_A_F(F2, affect[2], mu)*P_AA_H(HS1,HS2,AA=affect[3:4], mu,IBD)*p_f*p_IBD[IBD+1]*.25 #(1/2)^HS1
          p_N <- P_A_F(F1, affect[1], 1)*P_A_F(F2, affect[2], 1)*P_AA_H(HS1,HS2,AA=affect[3:4], 1,IBD)*p_f*p_IBD[IBD+1]*.25
          
          S <- h1*affect[1] + h2*affect[1] + h3*affect[2] + h4*affect[2] + HS2*2
          #p_A <- P_A_F(F2, mu)*P_AA_H(HS1,HS2, mu) #(1/2)^HS1
          #p_N <- P_A_F(F2, 1)*P_AA_H(HS1,HS2, 1)
          p_A_sum <- p_A_sum+p_A
          p_N_sum <- p_N_sum+p_N
          if(h!=1 & h!=16) {
            p_S_A[S+1] <- p_S_A[S+1] + p_A 
            p_S_N[S+1] <- p_S_N[S+1] + p_N
            
            p_A_n0_sum <- p_A_n0_sum + p_A
            E_A <- E_A + S*p_A
            
            p_N_n0_sum <- p_N_n0_sum + p_N
            E_N <- E_N + S*p_N
          }
        }
      }
    }
  }
  
  
  
  #
  print(c(p_A_n0_sum, p_A_sum))
  print(c(p_N_n0_sum, p_N_sum))
  p_S_A <- p_S_A/p_A_n0_sum #normalized by sum
  p_S_N <- p_S_N/p_N_n0_sum #normalized by sum
  E_A <- E_A/p_A_n0_sum #normalized by sum
  E_N <- E_N/p_N_n0_sum #normalized by sum
  
  p_A_n0 <- p_A_n0_sum/p_A_sum
  p_N_n0 <- p_N_n0_sum/p_N_sum
  
  S=0:(length(p_S_A)-1)
  sigma2_N <- sum((S-E_N)^2*p_S_N) #variance under the null
  sigma2_A <- sum((S-E_A)^2*p_S_A) #variance under the alternative
  print(c(E_A, E_N, p_A_n0, p_N_n0, sigma2_A, sigma2_N))
  print(p_S_A)
  print(p_S_N)
  crit.pt.L <- qnorm(alpha/2, mean=E_N, sd=sqrt(sigma2_N/(p_N_n0*n_family)))
  crit.pt.H <- qnorm(1-alpha/2, mean=E_N, sd=sqrt(sigma2_N/(p_N_n0*n_family)))
  power <- pnorm(crit.pt.L, mean=E_A, sd=sqrt(sigma2_A/(p_A_n0*n_family))) + pnorm(crit.pt.H, mean=E_A, sd=sqrt(sigma2_A/(p_A_n0*n_family)), lower=F)
  #c(S, p_A/p_A_sum, E_A, E_N, sigma2_A, sigma2_N, power)  
  print(power)
}
power.family.3(mu=1)

#calculate only 0,1,2,3,4 variant in founders and dump 0 and 4, fix the sibling affected status
power.family.4 <- function(fam.str=c(2,2), affect=c(c(0,1),c(1,1)), n_family=500, f=0.01, mu=1, p_dis=0.01, alpha=10^-6, stat=F) {
  #I go over the possible parents' genotypes and the combination of children's genotypes
  #but not the affected status, this effect will be included in the combination of children genotyps
  
  #calculate KG
  KL <- (1+f*(mu-1))^2
  KG=p_dis/KL
  
  P_A_F <- function(F, A, mu) { #prob of being affected for founder
    if(A==1 & F==1) return(mu)
    if(A==0 & F==1) return(1-mu*KG)
    return(1)
  }
  
  P_AA_H <- function(HS1, HS2, AA, mu, sharing_IBD, H_sib) { #prob of being double affected siblings
    if(sum(AA)==2) {
      if(sharing_IBD==0)  result <- mu^HS1
      if(sharing_IBD==1)  result <- mu^HS1*(mu^2)^HS2
      if(sharing_IBD==2)  result <- (mu^2)^HS2
    }
    if(sum(AA)==1) {
      if(sharing_IBD==0)  {
        #         if(HS1==0) result <- 1
        #         if(HS1==1) result <- 0.5*mu+0.5*(1-mu*KG)
        #         if(HS1==2) result <- 1/6*mu^2 + 4/6*mu*(1-mu*KG) + 1/6*(1-mu*KG)^2
        #         if(HS1==3) result <- 2/4*mu^2*(1-mu*KG) + 2/4*mu*(1-mu*KG)^2
        #         if(HS1==4) result <- mu^2*(1-mu*KG)^2
        result <- prod(ifelse(rep(AA, each=2)==1, mu, (1-mu*KG))^H_sib) #multiple each haplotype's effect which determined by the affected status
      }
      if(sharing_IBD==1)  {
        if(HS1==0) result <- (mu*(1-mu*KG))^HS2
        #         if(HS1==1) result <- 0.5*mu*(mu*(1-mu*KG))^HS2 + 0.5*(1-mu*KG)*(mu*(1-mu*KG))^HS2
        if(HS1==1) result <- ifelse(sum(H_sib[1:2]) > sum(H_sib[3:4]), ifelse(AA[1]==1,mu,1-mu*KG), ifelse(AA[2]==1, mu, 1-mu*KG))*(mu*(1-mu*KG))^HS2 #see if the affected carries the non-shared variant
        if(HS1==2) result <- mu*(1-mu*KG)*(mu*(1-mu*KG))^HS2
      }
      if(sharing_IBD==2)  {
        result <- (mu*(1-mu*KG))^HS2
      }
    }
    if(sum(AA)==0) {
      if(sharing_IBD==0)  result <- (1-mu*KG)^HS1
      if(sharing_IBD==1)  result <- (1-mu*KG)^HS1*((1-mu*KG)^2)^HS2
      if(sharing_IBD==2)  result <- ((1-mu*KG)^2)^HS2
    }
    result
  }
  
  #specify the family structure
  n_I <-  fam.str[1] #no. of the first generation
  n_II <- fam.str[2] #no. of the second generation
  
  #affected status
  affect <- affect  
  
  #initialization
  p_A_sum=0 #sum of prob. under the alternative
  p_N_sum=0 #sum of prob. under the null
  p_A_n0_sum=0 #sum of prob. under the alternative non 0 in founders
  p_N_n0_sum=0 #sum of prob. under the null non 0 in founders
  E_A=0 #expectation under the alternative
  E_N=0 #expectation under the null
  p_S_A=rep(0,10) #prob. of S under the alternative
  p_S_N=rep(0,10) #prob. of S under the null S=0 is places at [1]
  p_IBD = c(0.25, 0.5, 0.25) #PROB.of IBD
  # update <- function() { #update the sum
  #  S <- h1*(1 %in% IBD_Str[i,])*affect[1] + h2*(2 %in% IBD_Str[i,])*affect[1] + h3*(3 %in% IBD_Str[i,])*affect[2] + h4*(4 %in% IBD_Str[i,])*affect[2] + HS1 + HS2*2
  #p_A <- P_A_F(F2, mu)*P_AA_H(HS1,HS2, mu) #(1/2)^HS1
  #p_N <- P_A_F(F2, 1)*P_AA_H(HS1,HS2, 1)
  #  p_S_A[S+1] <<- p_S_A[S+1] + p_A 
  #  p_S_N[S+1] <<- p_S_N[S+1] + p_N
  #  p_A_sum <<- p_A_sum+p_A
  #  p_N_sum <<- p_N_sum+p_N
  
  #  E_A <<- E_A + S*p_A
  #  E_N <<- E_N + S*p_N
  #}
  for(IBD in 0:2) {
    founder <- matrix(c(0,0,0,0,
                        1,0,0,0,
                        0,1,0,0,
                        0,0,1,0,
                        0,0,0,1,
                        1,1,0,0,
                        1,0,1,0,
                        1,0,0,1,
                        0,1,1,0,
                        0,1,0,1,
                        0,0,1,1,
                        1,1,1,0,
                        1,1,0,1,
                        1,0,1,1,
                        0,1,1,1,
                        1,1,1,1), nrow=16, byrow=T)
    for(h in 1:nrow(founder)) {
      #initialization
      sce <- as.vector(founder[h,]) #founder's haplotype
      h1 <- sce[1]
      h2 <- sce[2]
      h3 <- sce[3]
      h4 <- sce[4]
      p_A <- 0
      p_N <- 0
      S <- NA
      
      HS1 <- HS2 <- 0 #shared variant
      F1 <- sum(sce[1:2]) #founder: father
      F2 <- sum(sce[3:4]) #founder: mothrt
      FB <- sum(sce) #founder's haplotype
      p_f <- f^(FB)*(1-f)^(4-FB) #allele frequency of founders
      
      if(IBD==0) {
        IBD_Str=matrix(c(1,2,3,4,
                         1,2,4,3,
                         1,3,2,4,
                         1,3,4,2,
                         1,4,2,3,
                         1,4,3,2,
                         2,1,3,4,
                         2,1,4,3,
                         2,3,1,4,
                         2,3,4,1,
                         2,4,1,3,
                         2,4,3,1,
                         3,1,2,4,
                         3,1,4,2,
                         3,2,1,4,
                         3,2,4,1,
                         3,4,1,2,
                         3,4,2,1,
                         4,1,2,3,
                         4,1,3,2,
                         4,2,1,3,
                         4,2,3,1,
                         4,3,1,2,
                         4,3,2,1), byrow=T, nrow=24)
        for(i in 1:nrow(IBD_Str)) {
          HS1 <- FB
          HS2 <- 0
          H_sib <- sce[IBD_Str[i,]] #variant or not on each haplotype
          
          p_A <- P_A_F(F1, affect[1], mu)*P_A_F(F2, affect[2], mu)*P_AA_H(HS1,HS2,AA=affect[3:4], mu,IBD,H_sib)*p_f*p_IBD[IBD+1]*1/24 #(1/2)^HS1
          p_N <- P_A_F(F1, affect[1], 1)*P_A_F(F2, affect[2], 1)*P_AA_H(HS1,HS2,AA=affect[3:4], 1,IBD,H_sib)*p_f*p_IBD[IBD+1]*1/24
          
          S <- h1*affect[1] + h2*affect[1] + h3*affect[2] + h4*affect[2] + sum(H_sib*rep(affect[3:4], each=2))
          #p_A <- P_A_F(F2, mu)*P_AA_H(HS1,HS2, mu) #(1/2)^HS1
          #p_N <- P_A_F(F2, 1)*P_AA_H(HS1,HS2, 1)
          p_A_sum <- p_A_sum+p_A
          p_N_sum <- p_N_sum+p_N
          if(h!=1 & h!=16) {
            p_S_A[S+1] <- p_S_A[S+1] + p_A 
            p_S_N[S+1] <- p_S_N[S+1] + p_N
            
            p_A_n0_sum <- p_A_n0_sum + p_A
            E_A <- E_A + S*p_A
            
            p_N_n0_sum <- p_N_n0_sum + p_N
            E_N <- E_N + S*p_N
          }
        }
      }
      #print(IBD)
      if(IBD==1) {
        IBD_Str=matrix(c(1,3,1,4,
                         1,4,1,3,
                         2,3,2,4,
                         2,4,2,3,
                         3,1,3,2,
                         3,2,3,1,
                         4,1,4,2,
                         4,2,4,1), byrow=T, nrow=8)
        for(i in 1:nrow(IBD_Str)) {
          H_sib <- sce[IBD_Str[i,]] #variant or not on each haplotype
          
          HS1 <- sum(table(match(IBD_Str[i,], which(sce==1)))==1) #how many variant shared once
          HS2 <- sum(table(match(IBD_Str[i,], which(sce==1)))==2) #how many variant shared twice
          p_A <- P_A_F(F1, affect[1], mu)*P_A_F(F2, affect[2], mu)*P_AA_H(HS1,HS2,AA=affect[3:4], mu,IBD,H_sib)*p_f*p_IBD[IBD+1]*1/8 #(1/2)^HS1
          p_N <- P_A_F(F1, affect[1], 1)*P_A_F(F2, affect[2], 1)*P_AA_H(HS1,HS2,AA=affect[3:4], 1,IBD,H_sib)*p_f*p_IBD[IBD+1]*1/8
          
          S <- h1*affect[1] + h2*affect[1] + h3*affect[2] + h4*affect[2] + sum(H_sib*rep(affect[3:4], each=2))
          #p_A <- P_A_F(F2, mu)*P_AA_H(HS1,HS2, mu) #(1/2)^HS1
          #p_N <- P_A_F(F2, 1)*P_AA_H(HS1,HS2, 1)
          p_A_sum <- p_A_sum+p_A
          p_N_sum <- p_N_sum+p_N
          if(h!=1 & h!=16) {
            p_S_A[S+1] <- p_S_A[S+1] + p_A 
            p_S_N[S+1] <- p_S_N[S+1] + p_N
            
            p_A_n0_sum <- p_A_n0_sum + p_A
            E_A <- E_A + S*p_A
            
            p_N_n0_sum <- p_N_n0_sum + p_N
            E_N <- E_N + S*p_N
          }
        }
      }
      if(IBD==2) {
        IBD_Str=matrix(c(1,3,1,3,
                         1,4,1,4,
                         2,3,2,3,
                         2,4,2,4), byrow=T, nrow=4)
        for(i in 1:4) {
          H_sib <- sce[IBD_Str[i,]] #variant or not on each haplotype
          HS1 <- sum(table(match(IBD_Str[i,], which(sce==1)))==1) #how many variant shared once
          HS2 <- sum(table(match(IBD_Str[i,], which(sce==1)))==2) #how many variant shared twice
          p_A <- P_A_F(F1, affect[1], mu)*P_A_F(F2, affect[2], mu)*P_AA_H(HS1,HS2,AA=affect[3:4], mu,IBD)*p_f*p_IBD[IBD+1]*.25 #(1/2)^HS1
          p_N <- P_A_F(F1, affect[1], 1)*P_A_F(F2, affect[2], 1)*P_AA_H(HS1,HS2,AA=affect[3:4], 1,IBD)*p_f*p_IBD[IBD+1]*.25
          
          S <- h1*affect[1] + h2*affect[1] + h3*affect[2] + h4*affect[2] + sum(H_sib*rep(affect[3:4], each=2))
          #p_A <- P_A_F(F2, mu)*P_AA_H(HS1,HS2, mu) #(1/2)^HS1
          #p_N <- P_A_F(F2, 1)*P_AA_H(HS1,HS2, 1)
          p_A_sum <- p_A_sum+p_A
          p_N_sum <- p_N_sum+p_N
          if(h!=1 & h!=16) {
            p_S_A[S+1] <- p_S_A[S+1] + p_A 
            p_S_N[S+1] <- p_S_N[S+1] + p_N
            
            p_A_n0_sum <- p_A_n0_sum + p_A
            E_A <- E_A + S*p_A
            
            p_N_n0_sum <- p_N_n0_sum + p_N
            E_N <- E_N + S*p_N
          }
        }
      }
    }
  }
  
  
  
  #
  print(c(p_A_n0_sum, p_A_sum))
  print(c(p_N_n0_sum, p_N_sum))
  p_S_A <- p_S_A/p_A_n0_sum #normalized by sum
  p_S_N <- p_S_N/p_N_n0_sum #normalized by sum
  E_A <- E_A/p_A_n0_sum #normalized by sum
  E_N <- E_N/p_N_n0_sum #normalized by sum
  
  p_A_n0 <- p_A_n0_sum/p_A_sum
  p_N_n0 <- p_N_n0_sum/p_N_sum
  
  S=0:(length(p_S_A)-1)
  sigma2_N <- sum((S-E_N)^2*p_S_N) #variance under the null
  sigma2_A <- sum((S-E_A)^2*p_S_A) #variance under the alternative
  print(c(E_A, E_N, p_A_n0, p_N_n0, sigma2_A, sigma2_N))
  print(p_S_A)
  print(p_S_N)
  crit.pt.L <- qnorm(alpha/2, mean=E_N, sd=sqrt(sigma2_N/(p_N_n0*n_family)))
  crit.pt.H <- qnorm(1-alpha/2, mean=E_N, sd=sqrt(sigma2_N/(p_N_n0*n_family)))
  power <- pnorm(crit.pt.L, mean=E_A, sd=sqrt(sigma2_A/(p_A_n0*n_family))) + pnorm(crit.pt.H, mean=E_A, sd=sqrt(sigma2_A/(p_A_n0*n_family)), lower=F)
  #c(S, p_A/p_A_sum, E_A, E_N, sigma2_A, sigma2_N, power) 
  if(stat==T) return(c(E_N, sigma2_N, p_N_n0, E_A, sigma2_A, p_A_n0))
  print(power)
}
power.family.4(fam.str=c(2,2), affect=c(c(0,1),c(1,1)), n_family=500, f=0.01, mu=1, p_dis=0.01)

#calculate only 0,1,2,3,4 variant in founders and dump 0 and 4, fix the sibling affected status
power.family.5 <- function(fam.str=c(2,2,1), affect=c(c(0,1),c(1,1),1), n_family=500, f=0.01, mu=1, p_dis=0.01, alpha=10^-6, stat=F) {
  #I go over the possible parents' genotypes and the combination of children's genotypes
  #but not the affected status, this effect will be included in the combination of children genotyps
  
  #calculate KG
  KL <- (1+f*(mu-1))^2
  KG=p_dis/KL
  
  P_A_F <- function(F, A, mu) { #prob of being affected for founder
    if(is.na(A)) return(1)
    if(A==1 & F>=1) return(mu^F)
    if(A==0 & F>=1) return((1-mu*KG)^F)
    return(1)
  }
  
  P_AA_H <- function(HS1, HS2, AA, mu, sharing_IBD, H_sib) { #prob of being double affected siblings
    if(sum(AA)==2) {
      if(sharing_IBD==0)  result <- mu^HS1
      if(sharing_IBD==1)  result <- mu^HS1*(mu^2)^HS2
      if(sharing_IBD==2)  result <- (mu^2)^HS2
    }
    if(sum(AA)==1) {
      if(sharing_IBD==0)  {
        #         if(HS1==0) result <- 1
        #         if(HS1==1) result <- 0.5*mu+0.5*(1-mu*KG)
        #         if(HS1==2) result <- 1/6*mu^2 + 4/6*mu*(1-mu*KG) + 1/6*(1-mu*KG)^2
        #         if(HS1==3) result <- 2/4*mu^2*(1-mu*KG) + 2/4*mu*(1-mu*KG)^2
        #         if(HS1==4) result <- mu^2*(1-mu*KG)^2
        result <- prod(ifelse(rep(AA, each=2)==1, mu, (1-mu*KG))^H_sib) #multiple each haplotype's effect which determined by the affected status
      }
      if(sharing_IBD==1)  {
        if(HS1==0) result <- (mu*(1-mu*KG))^HS2
        #         if(HS1==1) result <- 0.5*mu*(mu*(1-mu*KG))^HS2 + 0.5*(1-mu*KG)*(mu*(1-mu*KG))^HS2
        if(HS1==1) result <- ifelse(sum(H_sib[1:2]) > sum(H_sib[3:4]), ifelse(AA[1]==1,mu,1-mu*KG), ifelse(AA[2]==1, mu, 1-mu*KG))*(mu*(1-mu*KG))^HS2 #see if the affected carries the non-shared variant
        if(HS1==2) result <- mu*(1-mu*KG)*(mu*(1-mu*KG))^HS2
      }
      if(sharing_IBD==2)  {
        result <- (mu*(1-mu*KG))^HS2
      }
    }
    if(sum(AA)==0) {
      if(sharing_IBD==0)  result <- (1-mu*KG)^HS1
      if(sharing_IBD==1)  result <- (1-mu*KG)^HS1*((1-mu*KG)^2)^HS2
      if(sharing_IBD==2)  result <- ((1-mu*KG)^2)^HS2
    }
    result
  }
  
  #specify the family structure
  n_I <-  fam.str[1] #no. of the first generation
  n_II <- fam.str[2] #no. of the second generation
  
  #affected status
  affect <- affect  
  
  #initialization
  p_A_sum=0 #sum of prob. under the alternative
  p_N_sum=0 #sum of prob. under the null
  p_A_n0_sum=0 #sum of prob. under the alternative non 0 in founders
  p_N_n0_sum=0 #sum of prob. under the null non 0 in founders
  E_A=0 #expectation under the alternative
  E_N=0 #expectation under the null
  p_S_A=rep(0,10) #prob. of S under the alternative
  p_S_N=rep(0,10) #prob. of S under the null S=0 is places at [1]
  p_IBD = c(0.25, 0.5, 0.25) #PROB.of IBD
  # update <- function() { #update the sum
  #  S <- h1*(1 %in% IBD_Str[i,])*affect[1] + h2*(2 %in% IBD_Str[i,])*affect[1] + h3*(3 %in% IBD_Str[i,])*affect[2] + h4*(4 %in% IBD_Str[i,])*affect[2] + HS1 + HS2*2
  #p_A <- P_A_F(F2, mu)*P_AA_H(HS1,HS2, mu) #(1/2)^HS1
  #p_N <- P_A_F(F2, 1)*P_AA_H(HS1,HS2, 1)
  #  p_S_A[S+1] <<- p_S_A[S+1] + p_A 
  #  p_S_N[S+1] <<- p_S_N[S+1] + p_N
  #  p_A_sum <<- p_A_sum+p_A
  #  p_N_sum <<- p_N_sum+p_N
  
  #  E_A <<- E_A + S*p_A
  #  E_N <<- E_N + S*p_N
  #}
  for(IBD in 0:2) {
    
    #create founders' haplotypes
    n_founder=4
    result <- sapply(0:n_founder, function(y) {
      combn(n_founder, y, function(x) {
        temp <- rep(0,n_founder)
        sapply(x, function(x) {
          temp[x] <<- 1
        } )
        temp
      })
    })  
    founder <- matrix(unlist(result), ncol=n_founder, byrow=T)
    for(h in 1:nrow(founder)) {
      #initialization
      sce <- as.vector(founder[h,1:4]) #founder's haplotype
      h1 <- sce[1]
      h2 <- sce[2]
      h3 <- sce[3]
      h4 <- sce[4]
      p_A <- 0
      p_N <- 0
      S <- NA
      
      HS1 <- HS2 <- 0 #shared variant
      F1 <- sum(sce[1:2]) #founder: father
      F2 <- sum(sce[3:4]) #founder: mothrt
      FB <- sum(sce[1:4]) #founder's haplotype
      p_f <- f^(FB)*(1-f)^(4-FB) #allele frequency of founders
      
      if(IBD==0) {
        IBD_Str=matrix(c(1,2,3,4,
                         1,2,4,3,
                         1,3,2,4,
                         1,3,4,2,
                         1,4,2,3,
                         1,4,3,2,
                         2,1,3,4,
                         2,1,4,3,
                         2,3,1,4,
                         2,3,4,1,
                         2,4,1,3,
                         2,4,3,1,
                         3,1,2,4,
                         3,1,4,2,
                         3,2,1,4,
                         3,2,4,1,
                         3,4,1,2,
                         3,4,2,1,
                         4,1,2,3,
                         4,1,3,2,
                         4,2,1,3,
                         4,2,3,1,
                         4,3,1,2,
                         4,3,2,1), byrow=T, nrow=24)
        for(i in 1:nrow(IBD_Str)) {
          HS1 <- FB
          HS2 <- 0
          H_sib <- sce[IBD_Str[i,]] #variant or not on each haplotype
          
          p_A <- P_A_F(F1, affect[1], mu)*P_A_F(F2, affect[2], mu)*P_AA_H(HS1,HS2,AA=affect[3:4], mu,IBD,H_sib)*p_f*p_IBD[IBD+1]*1/nrow(IBD_Str)*P_A_F(H_sib[4], affect[5], mu) #(1/2)^HS1
          p_N <- P_A_F(F1, affect[1], 1)*P_A_F(F2, affect[2], 1)*P_AA_H(HS1,HS2,AA=affect[3:4], 1,IBD,H_sib)*p_f*p_IBD[IBD+1]*1/nrow(IBD_Str)*P_A_F(H_sib[4], affect[5], 1)
          
          S <- h1*affect[1] + h2*affect[1] + h3*affect[2] + h4*affect[2] + sum(H_sib*rep(affect[3:4], each=2)) + ifelse(is.na(affect[5]), 0, H_sib[4]*affect[5])
          #p_A <- P_A_F(F2, mu)*P_AA_H(HS1,HS2, mu) #(1/2)^HS1
          #p_N <- P_A_F(F2, 1)*P_AA_H(HS1,HS2, 1)
          p_A_sum <- p_A_sum+p_A
          p_N_sum <- p_N_sum+p_N
          if(h!=1 & h!=16) {
            p_S_A[S+1] <- p_S_A[S+1] + p_A 
            p_S_N[S+1] <- p_S_N[S+1] + p_N
            
            p_A_n0_sum <- p_A_n0_sum + p_A
            E_A <- E_A + S*p_A
            
            p_N_n0_sum <- p_N_n0_sum + p_N
            E_N <- E_N + S*p_N
          }
        }
      }
      #print(IBD)
      if(IBD==1) {
        IBD_Str=matrix(c(1,3,1,4,
                         1,4,1,3,
                         3,1,4,1,
                         4,1,3,1,
                         2,3,2,4,
                         2,4,2,3,
                         3,2,4,2,
                         4,2,3,2,
                         3,1,3,2,
                         3,2,3,1,
                         1,3,2,3,
                         2,3,1,3,
                         4,1,4,2,
                         4,2,4,1,
                         1,4,2,4,
                         2,4,1,4), byrow=T, nrow=16)
        for(i in 1:nrow(IBD_Str)) {
          H_sib <- sce[IBD_Str[i,]] #variant or not on each haplotype
          
          HS1 <- sum(table(match(IBD_Str[i,], which(sce==1)))==1) #how many variant shared once
          HS2 <- sum(table(match(IBD_Str[i,], which(sce==1)))==2) #how many variant shared twice
          p_A <- P_A_F(F1, affect[1], mu)*P_A_F(F2, affect[2], mu)*P_AA_H(HS1,HS2,AA=affect[3:4], mu,IBD,H_sib)*p_f*p_IBD[IBD+1]*1/nrow(IBD_Str)*P_A_F(H_sib[4], affect[5], mu) #(1/2)^HS1
          p_N <- P_A_F(F1, affect[1], 1)*P_A_F(F2, affect[2], 1)*P_AA_H(HS1,HS2,AA=affect[3:4], 1,IBD,H_sib)*p_f*p_IBD[IBD+1]*1/nrow(IBD_Str)*P_A_F(H_sib[4], affect[5], 1)
          
          S <- h1*affect[1] + h2*affect[1] + h3*affect[2] + h4*affect[2] + sum(H_sib*rep(affect[3:4], each=2)) + ifelse(is.na(affect[5]), 0, H_sib[4]*affect[5])
          #p_A <- P_A_F(F2, mu)*P_AA_H(HS1,HS2, mu) #(1/2)^HS1
          #p_N <- P_A_F(F2, 1)*P_AA_H(HS1,HS2, 1)
          p_A_sum <- p_A_sum+p_A
          p_N_sum <- p_N_sum+p_N
          if(h!=1 & h!=16) {
            p_S_A[S+1] <- p_S_A[S+1] + p_A 
            p_S_N[S+1] <- p_S_N[S+1] + p_N
            
            p_A_n0_sum <- p_A_n0_sum + p_A
            E_A <- E_A + S*p_A
            
            p_N_n0_sum <- p_N_n0_sum + p_N
            E_N <- E_N + S*p_N
          }
        }
      }
      if(IBD==2) {
        IBD_Str=matrix(c(1,3,1,3,
                         1,4,1,4,
                         3,1,3,1,
                         4,1,4,1,
                         2,3,2,3,
                         2,4,2,4,
                         3,2,3,2,
                         4,2,4,2), byrow=T, nrow=8)
        for(i in 1:nrow(IBD_Str)) {
          H_sib <- sce[IBD_Str[i,]] #variant or not on each haplotype
          HS1 <- sum(table(match(IBD_Str[i,], which(sce==1)))==1) #how many variant shared once
          HS2 <- sum(table(match(IBD_Str[i,], which(sce==1)))==2) #how many variant shared twice
          p_A <- P_A_F(F1, affect[1], mu)*P_A_F(F2, affect[2], mu)*P_AA_H(HS1,HS2,AA=affect[3:4], mu,IBD)*p_f*p_IBD[IBD+1]*1/nrow(IBD_Str)*P_A_F(H_sib[4], affect[5], mu) #(1/2)^HS1
          p_N <- P_A_F(F1, affect[1], 1)*P_A_F(F2, affect[2], 1)*P_AA_H(HS1,HS2,AA=affect[3:4], 1,IBD)*p_f*p_IBD[IBD+1]*1/nrow(IBD_Str)*P_A_F(H_sib[4], affect[5], 1)
          
          S <- h1*affect[1] + h2*affect[1] + h3*affect[2] + h4*affect[2] + sum(H_sib*rep(affect[3:4], each=2)) + ifelse(is.na(affect[5]), 0, H_sib[4]*affect[5])
          #p_A <- P_A_F(F2, mu)*P_AA_H(HS1,HS2, mu) #(1/2)^HS1
          #p_N <- P_A_F(F2, 1)*P_AA_H(HS1,HS2, 1)
          p_A_sum <- p_A_sum+p_A
          p_N_sum <- p_N_sum+p_N
          if(h!=1 & h!=16) {
            p_S_A[S+1] <- p_S_A[S+1] + p_A 
            p_S_N[S+1] <- p_S_N[S+1] + p_N
            
            p_A_n0_sum <- p_A_n0_sum + p_A
            E_A <- E_A + S*p_A
            
            p_N_n0_sum <- p_N_n0_sum + p_N
            E_N <- E_N + S*p_N
          }
        }
      }
    }
  }
  
  
  
  #
  print(c(p_A_n0_sum, p_A_sum))
  print(c(p_N_n0_sum, p_N_sum))
  p_S_A <- p_S_A/p_A_n0_sum #normalized by sum
  p_S_N <- p_S_N/p_N_n0_sum #normalized by sum
  E_A <- E_A/p_A_n0_sum #normalized by sum
  E_N <- E_N/p_N_n0_sum #normalized by sum
  
  p_A_n0 <- p_A_n0_sum/p_A_sum
  p_N_n0 <- p_N_n0_sum/p_N_sum
  
  S=0:(length(p_S_A)-1)
  sigma2_N <- sum((S-E_N)^2*p_S_N) #variance under the null
  sigma2_A <- sum((S-E_A)^2*p_S_A) #variance under the alternative
  print(c(E_A, E_N, p_A_n0, p_N_n0, sigma2_A, sigma2_N))
  print(p_S_A)
  print(p_S_N)
  crit.pt.L <- qnorm(alpha/2, mean=E_N, sd=sqrt(sigma2_N/(p_N_n0*n_family)))
  crit.pt.H <- qnorm(1-alpha/2, mean=E_N, sd=sqrt(sigma2_N/(p_N_n0*n_family)))
  power <- pnorm(crit.pt.L, mean=E_A, sd=sqrt(sigma2_A/(p_A_n0*n_family))) + pnorm(crit.pt.H, mean=E_A, sd=sqrt(sigma2_A/(p_A_n0*n_family)), lower=F)
  #c(S, p_A/p_A_sum, E_A, E_N, sigma2_A, sigma2_N, power) 
  if(stat==T) return(c(E_N, sigma2_N, p_N_n0, E_A, sigma2_A, p_A_n0))
  print(power)
}
power.family.5(fam.str=c(2,2), affect=c(c(0,1),c(1,1),1), n_family=500, f=0.01, mu=1, p_dis=0.01)


### Using Gamma distribution to generate the efffect of haplotype
##simulate family
gene.data.family <- function(m=1, var=0, f=0.01, SRR=5, p_dis=0.01) {
  adj <- ifelse(var==0 & m==1, 1, ifelse(var==0, m, 1))
  c = (m - 1) #shift to accommondate gamma distribution
  beta <- ifelse(var==0 & m==1, 0, var/c) #prevent not a number error when var==0 & m==1
  alpha = c/beta
  KL <- (1+f*(m-1))^2 #contribution from locus
  KLKLR <- 0.25*(1+f*(m-1))^4 + 0.5*(1+f*(m-1))^2*(1+f*(m^2+var-1)) + 0.25*(1+f*(m^2+var-1))^2
  KG <- p_dis/KL #contribution from other genome
  SRR <- SRR #sibling relaive risk
  KGKGR <- SRR*p_dis*p_dis/KLKLR #implement the heriatbility from other locus given SRR
  
  ##generate sibpairs
  n_family <- 1000000
  H1_f <- rbinom(n_family,1,f) #if the first haplotype carries risk raviant
  H2_f <- rbinom(n_family,1,f) #if the second haplotype carries risk raviant
  H1_rr_f <- ifelse(H1_f == 1, rgamma(length(H1_f), alpha, scale=beta) + adj, 1) #RR of the first haplotype
  H2_rr_f <- ifelse(H2_f == 1, rgamma(length(H2_f), alpha, scale=beta) + adj, 1) #RR of the second haplotype
  H1_m <- rbinom(n_family,1,f) #if the first haplotype carries risk raviant
  H2_m <- rbinom(n_family,1,f) #if the second haplotype carries risk raviant
  H1_rr_m <- ifelse(H1_m == 1, rgamma(length(H1_m), alpha, scale=beta) + adj, 1) #RR of the first haplotype
  H2_rr_m <- ifelse(H2_m == 1, rgamma(length(H2_m), alpha, scale=beta) + adj, 1) #RR of the second haplotype
  
  H1_inh_s <- sample(c(1,2),n_family, replace=TRUE) 
  H1_s <- ifelse(H1_inh_s==1, H1_f, H2_f)
  H1_rr_s <- ifelse(H1_inh_s==1, H1_rr_f, H2_rr_f)
  
  H2_inh_s <- sample(c(1,2),n_family, replace=TRUE) 
  H2_s <- ifelse(H2_inh_s==1, H1_m, H2_m)
  H2_rr_s <- ifelse(H2_inh_s==1, H1_rr_m, H2_rr_m)
  
  #penetrance_s <- H1_rr_s*H2_rr_s #penetrance of disease given haplotypes
  #penetrance_s <- ifelse(penetrance_s>1, 1, penetrance_s)
  #dis_s <- rbinom(length(penetrance_s),1,penetrance_s) #disease status
  
  H1_inh <- sample(c(1,2), n_family, replace=TRUE) 
  H1 <- ifelse(H1_inh==1, H1_f, H2_f)
  H1_rr <- ifelse(H1_inh==1, H1_rr_f, H2_rr_f)
  
  H2_inh <- sample(c(1,2), n_family, replace=TRUE) 
  H2 <- ifelse(H2_inh==1, H1_m, H2_m)
  H2_rr <- ifelse(H2_inh==1, H1_rr_m, H2_rr_m)
  
  #rm(H1_f, H2_f, H1_rr_f, H2_rr_f, H1_m, H2_m, H1_rr_m, H2_rr_m) #save memory
  
  S0 <- (H1_inh_s != H1_inh) & (H2_inh_s != H2_inh)
  S2 <- (H1_inh_s == H1_inh) & (H2_inh_s == H2_inh)
  S1 <- ((H1_inh_s == H1_inh) & (H2_inh_s != H2_inh)) | ((H1_inh_s != H1_inh) & (H2_inh_s == H2_inh))
  S11 <- ((H1_inh_s == H1_inh) & (H2_inh_s != H2_inh))
  S12 <- ((H1_inh_s != H1_inh) & (H2_inh_s == H2_inh))
  
  S <- (S0==TRUE)*0 + (S1==TRUE)*1 + (S2==TRUE)*2
  HS.mis <- HS <- (S1==TRUE&S11)*(H1) + (S1==TRUE&S12)*(H2) + (S2==TRUE)*(H1+H2)
  HN.mis <- HN <- (S0==TRUE)*(H1+H2+H1_s+H2_s) + (S1==TRUE&S11)*(H2+H2_s) + (S1==TRUE&S12)*(H1+H1_s)
  HS.mis <- ifelse(S1==TRUE & HS==0 & HN==2, 1, HS.mis) 
  HN.mis <- ifelse(S1==TRUE & HS==0 & HN==2, 0, HN.mis) 
  
  #modified to unaffected father, affected mother and both affected children
  penetrance <- (1-H1_rr_f*H2_rr_f*KG)*(H1_rr_m*H2_rr_m*KG)*H1_rr*H2_rr*H1_rr_s*H2_rr_s*KGKGR #penetrance of disease given haplotypes of both siblings
  penetrance <- ifelse(penetrance>1, 1, penetrance)
  dis <- rbinom(length(penetrance),1,penetrance) #disease status of both affected
  
  #data.family <- data.frame(H1_f, H2_f, H1_rr_f, H2_rr_f,   #for varification
  #                          H1_m, H2_m, H1_rr_m, H2_rr_m,
  #                          H1_inh, H1, H1_rr, H2_inh, H2, H2_rr, penetrance, dis,
  #                          H1_inh_s, H1_s, H1_rr_s, H2_inh_s, H2_s, H2_rr_s, penetrance_s, dis_s)
  
  idx <- which(dis==1)
  
  #count the T statistics
  T_stat = (H1_f==1)*((H1_inh_s==1)+(H1_inh==1)) + (H2_f==1)*((H1_inh_s==2)+(H1_inh==2)) + (H1_m==1)*((H2_inh_s==1)+(H2_inh==1)) + (H2_m==1)*((H2_inh_s==2)+(H2_inh==2)) + H1_m + H2_m
  
  data.family <<- data.frame(H1_f=H1_f[idx], H2_f=H2_f[idx], H1_m=H1_m[idx], H2_m=H2_m[idx], 
                             H1=H1[idx], H1_rr=H1_rr[idx], H2=H2[idx], H2_rr=H2_rr[idx], dis=dis[idx],
                             H1_s=H1_s[idx], H1_rr_s=H1_rr_s[idx], H2_s=H2_s[idx], H2_rr_s=H2_rr_s[idx],
                             S0=S0[idx], S1=S1[idx], S2=S2[idx], S11=S11[idx], S12=S12[idx],
                             HS=HS[idx],HN=HN[idx],S=S[idx], HS.mis=HS.mis[idx],HN.mis=HN.mis[idx],T_stat=T_stat[idx])
}
#gene.data.family(m=1, var=0, f=0.01, SRR=1000)
