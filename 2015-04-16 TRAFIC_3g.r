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

##gene drop simulation for two generations
family_strct.2g2c <- data.frame(family=c(1,1,1,1), person=c(1,2,3,4), father=c(0,0,1,1), 
                                mother=c(0,0,2,2), sex=c(1,2,1,1), affect=c(1,2,2,2)) #1=male, 2=female, 1=unaffected, 2=affected
family_strct.3g <- data.frame(family=c(1,1,1,1,1,1,1), person=c(1,2,3,4,5,6,7), father=c(0,0,1,1,0,4,4), 
                                mother=c(0,0,2,2,0,5,5), sex=c(1,2,1,1,2,1,1), affect=c(1,2,1,2,1,2,1)) #1=male, 2=female, 1=unaffected, 2=affected

rep.idx <<- 1
sim_result <- replicate(n_rep, {
  print(rep.idx)
  rep.idx <<- rep.idx + 1
  
  #simulation for two-generation families
  family_generated_2g2c <- gene_family(family_strct=family_strct.2g2c, n=n_family, haplotype.risk=haplotype.risk) 
  result.trap.2g2c <- family.test(data=family_generated_2g2c, f=risk.variant.id)
  result.trafic.ext.2g2c <- family.test.trafic.ext(data=family_generated_2g2c, f=risk.variant.id)
    
  #simulation for three-generation families
  family_generated_3g <- gene_family(family_strct=family_strct.3g, n=n_family, haplotype.risk=haplotype.risk) 
  result.trap.3g <- family.test(data=family_generated_3g, f=risk.variant.id)
  result.trafic.ext.3g <- family.test.trafic.ext(data=family_generated_3g, f=risk.variant.id)
  
  #only report p.value
  c(result.trap.2g2c[7], result.trafic.ext.2g2c, result.trap.3g[7], result.trafic.ext.3g)
  
})
## Write out your results to a csv file
write.csv(data.frame(seed=seed, r=r, result.trap.2g2c=sim_result[1,], result.trafic.ext.2g2c=sim_result[2,],
                     result.trap.3g=sim_result[3,], result.trafic.ext.3g=sim_result[4,]),
          paste("res_",r,"_",n_family,"_",seed,".csv",sep=""), row.names=FALSE)

#Load results with 3g and calculate the power
result.3g <- read.csv("2015-04-16 TRAFIC_3g.csv", header=T)
result.3g <- melt(result.3g, id.vars=c("seed", "r"), variable.name="method", value.name="p.value")
result.3g.power <- result.3g %>% group_by(r, method) %>% 
  summarise(n=n(), power=mean(p.value<0.05))

#Load results with 2g2c and calculate the power
result.2g2c <- read.csv("2015-03-27 TRAFIC_ext.csv", header=T)
result.2g2c <- melt(result.2g2c, id.vars=c("seed", "r"), variable.name="method", value.name="p.value")
result.2g2c.power <- result.2g2c %>% filter(method %in% c("result.2g2c", "result.trafic.ext")) %>% group_by(r, method) %>% 
  summarise(n=n(), power=mean(p.value<0.05))

#combine 3g and 2g2c results
result.combined <- rbind(result.2g2c.power, result.3g.power)

ggplot(result.combined, aes(x=r, y=power, col=method)) +
  #   geom_point(size=3, alpha=1) +
  geom_line(size=1.2, alpha=0.7) +
  ggtitle("f=0.0102, consider effect size of risk haplotypes") +
  labs(x="relative risk r") +
  theme_gray(base_size = 20) +
  scale_color_discrete(name="Method", 
                       breaks=c("result.2g2c", "omnibus", "result.trafic.ext"),
                       labels=c("TRAP 2g2c", "Omnibus", "TRAFIC 2g2c"))


