source("function.r")
#need to take in parameter
## passed in from the command prompt.
parseCommandArgs()
## Will disply the default values on the first run,
print(seed)
print(n_rep)
print(r)
print(n_family) #number of family
print(p_dis) #prevelance

## ... your simualtion code
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
risk.variant.id <- c(3, 8,19,21,23,27,44,47,49,50)
risk.haplo.id <- which(apply(2-haplotype[, risk.variant.id+2], 1, sum)>0)
(rish.haplo.f <- mean(apply(2-haplotype[, risk.variant.id+2], 1, sum)>0)) #carrier haplotype frequency
haplotype.risk <- rep(1, length=nrow(haplotype))
#assign mean relative risk calculate the haplotype variants p(A|h)
haplotype.risk[risk.haplo.id] <- r
mean(haplotype.risk[risk.haplo.id]) #mean relative risk
haplotype.risk <<- haplotype.risk*b0_sqrt

##gene drop simulation for three generations
family_strct.2g2c <- data.frame(family=c(1,1,1,1), person=c(1,2,3,4), father=c(0,0,1,1), 
                                mother=c(0,0,2,2), sex=c(1,2,1,1), affect=c(1,2,2,2)) #1=male, 2=female, 1=unaffected, 2=affected

rep.idx <<- 1
sim_result <- replicate(n_rep, {
  print(rep.idx)
  rep.idx <<- rep.idx + 1
  
#   family_generated_2g2c <- gene_family(family_strct=family_strct.2g2c, n=n_family) 
#   result.complete <- family.test(data=family_generated_2g2c, f=risk.variant.id)
  
  family_generated_2c <- gene_family(family_strct=family_strct.2g2c, n=2*n_family) #double the family size to compensate for the lose of founders
  temp.idx <- which(family_generated_2c$data_family$person==1:2) #take out founders
  family_generated_2c$data_family <- family_generated_2c$data_family[-temp.idx, ]
  family_generated_2c$tran_vec <- family_generated_2c$tran_vec[-temp.idx, ]
  result.nofounder <- family.test.nofounder(data=family_generated_2c, f=risk.variant.id)
  result.trafic <- trafic.test(family_generated_2c, risk.variant.id)

  c(result.nofounder, result.trafic)
})
## Write out your results to a csv file
write.csv(data.frame(seed=seed, r=r, finail.test.stat=sim_result[1,], sum_e=sim_result[2,], se=sim_result[3,], 
                     mean_observed=sim_result[4,], mean_mean=sim_result[5,], mean_var=sim_result[6,],
                     p.value=sim_result[7,], n_info_family=sim_result[8,],
                     trafic.true=sim_result[9,], trafic.haplotype=sim_result[10,], trafic.genotype=sim_result[11,]),
          paste("res_",r,"_",n_family,"_",seed,".csv",sep=""), row.names=FALSE)

#plot for the results under the null
#result only 2c using TRAP and TRAFIC
result.2c <- read.csv("2015-02-27 figure 7. compare 2g2c and 2conly.csv", stringsAsFactors=FALSE)
result.2g2c <- read.csv("2015-02-27 figure 7_2g2c.csv", stringsAsFactors=FALSE)
#
dev.new(width=6, height=5)
lwd=3
par(mfrow=c(1,1), mar=c(4, 4, 2, 1), oma=c(0,0,0,0))
alpha=0.05
plot(sort(unique(result.2g2c$r)), tapply(result.2g2c$p.value<alpha, result.2g2c$r, mean), col="black", lty=1, lwd=lwd, type="l", xlab="mean relative risk", ylab="power", ylim=c(0,1))
lines(sort(unique(result.2c$r)), tapply(result.2c$p.value<alpha, result.2c$r, mean), col="seagreen2", lty=1,lwd=lwd)
lines(sort(unique(result.2c$r)), tapply(result.2c$trafic.haplotype<alpha, result.2c$r, mean), col="royalblue", lty=1,lwd=lwd)
#lines(unique(result_three2nd_0.1_0.01$r), tapply(result_three2nd_0.1_0.01$p_value<alpha, result_three2nd_0.1_0.01$r, mean), col=5, lty=1,lwd=lwd)
legend("left", c("TRAP 2g2c", "TRAP 2c", "TRAFIC 2c"), col=c("black","seagreen2","royalblue"), lty=1, lwd=lwd)

