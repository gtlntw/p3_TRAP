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
print(null)
print(family_strct)

## ... your simualtion code
##the founder's haplotypes
#import the haplotypes generated by cosi
haplotype <- read.table("../out_100k_10k_1kb.hap-1", header=F)
colnames(haplotype) <- c("HAP", "CHROM", paste("SNP", 1:(ncol(haplotype)-2), sep=""))
snp <-read.table("../out_100k_10k_1kb.pos-1", header=T)
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


##gene drop simulation for three generations
family_strct.2g2c <- data.frame(family=c(1,1,1,1), person=c(1,2,3,4), father=c(0,0,1,1), 
                                mother=c(0,0,2,2), sex=c(1,2,1,1), affect=c(1,2,2,2)) #1=male, 2=female, 1=unaffected, 2=affected
family_strct.3g <- data.frame(family=c(1,1,1,1,1,1), person=c(1,2,3,4,5,6), father=c(0,0,1,1,0,4), 
                              mother=c(0,0,2,2,0,5), sex=c(1,2,1,1,2,1), affect=c(1,2,2,2,1,2)) #1=male, 2=female, 1=unaffected, 2=affected
family_strct.2g3c <- data.frame(family=c(1,1,1,1,1), person=c(1,2,3,4,5), father=c(0,0,1,1,1), 
                                mother=c(0,0,2,2,2), sex=c(1,2,1,1,1), affect=c(1,2,2,2,2)) #1=male, 2=female, 1=unaffected, 2=affected

null==FALSE

rep.idx <<- 1
sim_result <- replicate(n_rep, {
  ##generate risk haplotypes
  #   risk.variant.id <- sort(c(sample(which((snp$FREQ1 < 0.05)==T & snp$FREQ1 != 1/n_haplo ), 5), sample(which(snp$FREQ1 == 1/n_haplo), 5))) #randomly pick 5 variants with f < 0.05 and 5 singletons as casual 
  risk.variant.id <- c(3, 8,19,21,23,27,44,47,49,50)
  risk.haplo.id <- which(apply(2-haplotype[, risk.variant.id+2], 1, sum)>0)
  (rish.haplo.f <- mean(apply(2-haplotype[, risk.variant.id+2], 1, sum)>0)) #carrier haplotype frequency
  
  haplotype.risk <- rep(1, length=nrow(haplotype))
  #assign mean relative risk calculate the haplotype variants p(A|h)
  if(null==FALSE) haplotype.risk[risk.haplo.id] <- r
  mean(haplotype.risk) #mean relative risk
  haplotype.risk <<- haplotype.risk*b0_sqrt
  
  print(rep.idx)
  rep.idx <<- rep.idx + 1
  family_generated <- gene_family(family_strct=eval(parse(text=family_strct)), n=n_family) 
  if(family_strct=="family_strct.3g") { #delete the 6th person
    temp.idx <- which(family_generated$data_family$person==5) #take out 5th person
    family_generated$data_family <- family_generated$data_family[-temp.idx, ]
    family_generated$tran_vec <- family_generated$tran_vec[-temp.idx, ]
  }
  family.test(data=family_generated, f=risk.variant.id)
})
## Write out your results to a csv file
write.csv(data.frame(seed=seed, r=r, finail.test.stat=sim_result[1,], sum_e=sim_result[2,], se=sim_result[3,], 
                     mean_observed=sim_result[4,], mean_mean=sim_result[5,], mean_var=sim_result[6,],
                     p.value=sim_result[7,], n_info_family=sim_result[8,]),
          paste("res_",r,"_",n_family,"_",seed,".csv",sep=""), row.names=FALSE)


#plot for the results under the null
alpha=0.05
lwd=3
result.2g2c <- read.csv("2014-10-19 figure 7_2g2c.csv", stringsAsFactors=FALSE)
result.2g3c <- read.csv("2014-10-19 figure 7_2g3c.csv", stringsAsFactors=FALSE)
result.3g_5ppl <- read.csv("2014-10-19 figure 7_3g_5ppl.csv", stringsAsFactors=FALSE)

#
plot(sort(unique(result.2g2c$r)), tapply(result.2g2c$p.value<alpha, result.2g2c$r, mean), col=1, lty=1, lwd=lwd, type="l", xlab="mean relative risk", ylab="power", ylim=c(0,1))
lines(sort(unique(result.2g3c$r)), tapply(result.2g3c$p.value<alpha, result.2g3c$r, mean), col=2, lty=1,lwd=lwd)
lines(sort(unique(result.3g_5ppl$r)), tapply(result.3g_5ppl$p.value<alpha, result.3g_5ppl$r, mean), col=4, lty=1,lwd=lwd)
#lines(unique(result_three2nd_0.1_0.01$r), tapply(result_three2nd_0.1_0.01$p_value<alpha, result_three2nd_0.1_0.01$r, mean), col=5, lty=1,lwd=lwd)
legend("topleft", c("two generations 2:2:0", "three generations 2:2:1", "three 2nd generations 2:3:0"), col=c(1,2, 4), lty=1)

