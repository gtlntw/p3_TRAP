## ... your simualtion codek
##the founder's haplotypes
#import the haplotypes generated by cosi
haplotype <- read.table("out_100k_10k_1kb.hap-1", header=F)
colnames(haplotype) <- c("HAP", "CHROM", paste("SNP", 1:(ncol(haplotype)-2), sep=""))
snp <-read.table("out_100k_10k_1kb.pos-1", header=T)
#make allele 1 is the minor allele, 2 is the common allele
temp.idx <- snp$FREQ1 > snp$FREQ2
temp.freq <- snp$FREQ2
snp$FREQ2[temp.idx] <- snp$FREQ1[temp.idx]
snp$FREQ1[temp.idx] <- temp.freq[temp.idx]
#also change the genotype file
haplotype[,which(temp.idx==T)+2] <- 3 - haplotype[,which(temp.idx==T)+2]

#allele frequency
nrow(snp) #total number of snp
sum(snp$FREQ1 < 0.05) # number of snp with f < 0.05
sum(snp$FREQ1 < 0.01) # number of snp with f < 0.01
sum(snp$FREQ1 == 0.0001) # number of singletons

##assign risk variants and the corresponding effect size (proportional to allele frequency)
# null <- FALSE
n_haplo <- 10000
n_snp <- ncol(haplotype)-2
prevalence <- p_dis
b0_sqrt <- sqrt(prevalence)   #baseline

#set up causual SNPs
#generate risk haplotypes
# risk.variant.id <- c(3, 8,19,21,23,27,44,47,49,50)
risk.variant.id <- c(3,39)
risk.haplo.id <- which(apply(2-haplotype[, risk.variant.id+2], 1, sum)>0)
(risk.haplo.f <- mean(apply(2-haplotype[, risk.variant.id+2], 1, sum)>0)) #carrier haplotype frequency
haplotype.risk <- rep(1, length=nrow(haplotype))
#assign mean relative risk calculate the haplotype variants p(A|h)
haplotype.risk[risk.haplo.id] <- r
mean(haplotype.risk[risk.haplo.id]) #mean relative risk
haplotype.risk <<- haplotype.risk*b0_sqrt

##gene drop simulation for two generations
family_strct.2g3c <- data.frame(family=c(1,1,1,1,1), person=c(1,2,3,4,5), father=c(0,0,1,1,1), 
                                mother=c(0,0,2,2,2), sex=c(1,2,1,1,1), affect=c(1,2,2,2,2)) #1=male, 2=female, 1=unaffected, 2=affected


rep.idx <<- 1
sim_result <- replicate(n_rep, {
  print(rep.idx)
  rep.idx <<- rep.idx + 1
  

  family_generated_2g3c <<- gene_family(family_strct=family_strct.2g3c, n=n_family, haplotype.risk=haplotype.risk) 
  
  #remove the founder from the data
  family_generated_3c <- family_generated_2g3c
  temp.idx <- which(family_generated_2g3c$data_family$person %in% 1:2) #take out founders
  family_generated_3c$data_family <- family_generated_2g3c$data_family[-temp.idx, ]
  family_generated_3c$tran_vec <- family_generated_2g3c$tran_vec[-temp.idx, ]

  #test based on new imputation framework
  result.trap.impute.3c <- family.test.nofounder.impute(data=family_generated_3c, f=39, method="trap")
  result.trap.impute.noaggregate.3c <- family.test.nofounder.impute.noaggregate(data=family_generated_3c, f=39, method="trap")
  result.trafic.impute.3c <- family.test.nofounder.impute(data=family_generated_3c, f=39, method="trafic")
  result.trafic.impute.noaggregate.3c <- family.test.nofounder.impute.noaggregate(data=family_generated_3c, f=39, method="trafic")
  
  
  ##below is the baseline power to compare with to see how much power lost
  result.trap.nofounder.3c <- family.test.nofounder(data=family_generated_3c, f=39)
  result.trafic.ext.3c <- family.test.trafic.ext(data=family_generated_3c, f=39)
   
  
  #only report p.value
  c(result.trap.impute.3c=result.trap.impute.3c, 
    result.trap.impute.noaggregate.3c=result.trap.impute.noaggregate.3c,
    result.trafic.impute.3c=result.trafic.impute.3c,
    result.trafic.impute.noaggregate.3c=result.trafic.impute.noaggregate.3c,
    result.trap.nofounder.3c=result.trap.nofounder.3c[7],
    result.trafic.ext.3c=result.trafic.ext.3c)
})
## Write out your results to a csv file
result.df <- as.data.frame(t(sim_result))
colnames(result.df) <- c("result.trap.impute.3c","result.trap.impute.noaggregate.3c",
                         "result.trafic.impute.3c", "result.trafic.impute.noaggregate.3c",
                         "result.trap.nofounder.3c", "result.trafic.ext.3c")
result.df <- cbind(seed,r,p_dis,risk.haplo.f,n_family,result.df)
write.csv(result.df, paste("res_",r,"_",n_family,"_",seed,".csv",sep=""), row.names=FALSE)
## R.miSniam
## End(Not run)
## Not run:

#Load results with 3g and check the allele frequency
result <- read.csv("2015-05-19 imputation framework for TRAP and TRAFIC_ext fix baseline and test on rarer variant.csv", header=T)
result <- result %>% melt(id.vars=c("seed", "r", "p_dis","rish.haplo.f","n_family"), variable.name="method", value.name="p.value")
result.plot <- result %>% group_by(r, method) %>% 
  summarise(n=n(), power=mean(p.value < 0.05))

#TRAP
filter(result.plot, grepl("trap", method)) %>% ggplot(aes(x=r, y=power, group=method, col=method)) +
  #   geom_point(size=3, alpha=1) +
  geom_line(size=1.2, alpha=0.7) +
  ggtitle("f=0.004, consider effect size of risk haplotypes, TRAP") +
  labs(x="relative risk r") +
  theme_gray(base_size = 20)

#TRAFIC
filter(result.plot, grepl("trafic", method)) %>% ggplot(aes(x=r, y=power, group=method, col=method)) +
  #   geom_point(size=3, alpha=1) +
  geom_line(size=1.2, alpha=0.7) +
  ggtitle("f=0.004, consider effect size of risk haplotypes, TRAFIC") +
  labs(x="relative risk r") +
  theme_gray(base_size = 20)

#TRAP and TRAFIC
result.plot %>% ggplot(aes(x=r, y=power, group=method, col=method)) +
  #   geom_point(size=3, alpha=1) +
  geom_line(size=1.2, alpha=0.7) +
  ggtitle("f=0.004, consider effect size of risk haplotypes, TRAP and TRAFIC") +
  labs(x="relative risk r") +
  theme_gray(base_size = 20)