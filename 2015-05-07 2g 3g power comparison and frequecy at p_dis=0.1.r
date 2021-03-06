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
family_strct.2g.3a.1u <- data.frame(family=c(1,1,1,1), person=c(1,2,3,4), father=c(0,0,1,1), 
                                mother=c(0,0,2,2), sex=c(1,2,1,1), affect=c(1,2,2,2)) #1=male, 2=female, 1=unaffected, 2=affected
family_strct.3g.3a.4u <- data.frame(family=c(1,1,1,1,1,1,1), person=c(1,2,3,4,5,6,7), father=c(0,0,1,1,0,4,4), 
                                 mother=c(0,0,2,2,0,5,5), sex=c(1,2,1,1,2,1,1), affect=c(1,2,1,2,1,2,1)) #1=male, 2=female, 1=unaffected, 2=affected
family_strct.3g.2a.5u <- data.frame(family=c(1,1,1,1,1,1,1), person=c(1,2,3,4,5,6,7), father=c(0,0,1,1,0,4,4), 
                                 mother=c(0,0,2,2,0,5,5), sex=c(1,2,1,1,2,1,1), affect=c(1,2,1,1,1,2,1)) #1=male, 2=female, 1=unaffected, 2=affected
family_strct.3g.3a.5u <- data.frame(family=c(1,1,1,1,1,1,1,1), person=c(1,2,3,4,5,6,7,8), father=c(0,0,1,1,0,4,4,4), 
                                mother=c(0,0,2,2,0,5,5,5), sex=c(1,2,1,1,2,1,1,1), affect=c(1,2,1,2,1,2,1,1)) #1=male, 2=female, 1=unaffected, 2=affected

rep.idx <<- 1
sim_result <- list()
replicate(n_rep, {
  print(rep.idx)
  
  #simulation for two-generation families
  family_generated_2g.3a.1u <- gene_family(family_strct=family_strct.2g.3a.1u, n=n_family, haplotype.risk=haplotype.risk) 
  result.trap.2g.3a.1u <- family.test(data=family_generated_2g.3a.1u, f=risk.variant.id)
  result.trafic.ext.2g.3a.1u <- family.test.trafic.ext(data=family_generated_2g.3a.1u, f=risk.variant.id)
  result.maf.2g.3a.1u <- family.maf(data=family_generated_2g.3a.1u, f=risk.variant.id)
  
  #simulation for three-generation families
  #three affected and four unaffected
  family_generated_3g.3a.4u <- gene_family(family_strct=family_strct.3g.3a.4u, n=n_family, haplotype.risk=haplotype.risk) 
  result.trap.3g.3a.4u <- family.test(data=family_generated_3g.3a.4u, f=risk.variant.id)
  result.trafic.ext.3g.3a.4u <- family.test.trafic.ext(data=family_generated_3g.3a.4u, f=risk.variant.id)
  result.maf.3g.3a.4u <- family.maf(data=family_generated_3g.3a.4u, f=risk.variant.id)
  #two affected and five unaffected
  family_generated_3g.2a.5u <- gene_family(family_strct=family_strct.3g.2a.5u, n=n_family, haplotype.risk=haplotype.risk) 
  result.trap.3g.2a.5u <- family.test(data=family_generated_3g.2a.5u, f=risk.variant.id)
  result.trafic.ext.3g.2a.5u <- family.test.trafic.ext(data=family_generated_3g.2a.5u, f=risk.variant.id)
  result.maf.3g.2a.5u <- family.maf(data=family_generated_3g.2a.5u, f=risk.variant.id)
  #three affected and five unaffected
  family_generated_3g.3a.5u <- gene_family(family_strct=family_strct.3g.3a.5u, n=n_family, haplotype.risk=haplotype.risk)
  result.trap.3g.3a.5u <- family.test(data=family_generated_3g.3a.5u, f=risk.variant.id)
  result.trafic.ext.3g.3a.5u <- family.test.trafic.ext(data=family_generated_3g.3a.5u, f=risk.variant.id)
  result.maf.3g.3a.5u <- family.maf(data=family_generated_3g.3a.5u, f=risk.variant.id)
  #report maf
  sim_result[[rep.idx]] <<- rbind(
    data.frame(family.strct="2g.3a.1u", result.maf.2g.3a.1u, trap=result.trap.2g.3a.1u[7], trafic=result.trafic.ext.2g.3a.1u),
    data.frame(family.strct="3g.3a.4u", result.maf.3g.3a.4u, trap=result.trap.3g.3a.4u[7], trafic=result.trafic.ext.3g.3a.4u),
    data.frame(family.strct="3g.2a.5u", result.maf.3g.2a.5u, trap=result.trap.3g.2a.5u[7], trafic=result.trafic.ext.3g.2a.5u),
    data.frame(family.strct="3g.3a.5u", result.maf.3g.3a.5u, trap=result.trap.3g.3a.5u[7], trafic=result.trafic.ext.3g.3a.5u)
  )
  rep.idx <<- rep.idx + 1
})

sim_result <- do.call(rbind, sim_result)

## Write out your results to a csv file
#becareful that power for each iteration has three rows; needs to process to keep only one
result.df <- as.data.frame(sim_result)
colnames(result.df) <- c("family.strct","haplotype_on_affect",
                         "noncarrier","carrier",
                         "total","carrier.prop", "trap", "trafic")
result.df <- cbind(seed,r,p_dis,rish.haplo.f,n_family,result.df)
write.csv(result.df, paste("res_",r,"_",n_family,"_",seed,".csv",sep=""), row.names=FALSE)

## R.miSniam
## End(Not run)
## Not run:

##allele ferquency on shared chromosomes
#Load results with 3g and check the allele frequency
result <- read.csv("2015-05-07 2g 3g power comparison and frequecy at p_dis=0.1.csv", header=T)
result$family.strct <- factor(result$family.strct, levels=c("2g.3a.1u", "3g.3a.4u",
                                                            "3g.2a.5u", "3g.3a.5u"))
result$haplotype_on_affect <- factor(result$haplotype_on_affect)

# result <- melt(id.vars=c("seed", "r", "p_dis","rish.haplo.f","n_family", ""), variable.name="", value.name="p.value")
result.plot <- result %>% group_by(r, family.strct, haplotype_on_affect) %>% 
  summarise(n=n(), non.carrier.mean=mean(noncarrier), carrier.mean=mean(carrier),
            maf=mean(carrier.prop))

#maf on shared 1,2,3 chromosome
ggplot(result.plot, aes(x=r, y=maf, group=haplotype_on_affect, col=haplotype_on_affect)) +
  #   geom_point(size=3, alpha=1) +
  geom_line(size=1.2, alpha=0.7) +
  ggtitle("f=0.0102, consider effect size of risk haplotypes") +
  labs(x="relative risk r") +
  theme_gray(base_size = 20) + 
  facet_grid(.~family.strct)


##power comparison between 
result <- read.csv("2015-05-07 2g 3g power comparison and frequecy at p_dis=0.1.csv", header=T)
result$family.strct <- factor(result$family.strct, levels=c("2g.3a.1u", "3g.3a.4u",
                                                            "3g.2a.5u", "3g.3a.5u"))
result$haplotype_on_affect <- factor(result$haplotype_on_affect)
result <- result %>% select(-noncarrier, -carrier, -total, -carrier.prop) %>% 
  filter(haplotype_on_affect == 1) %>% melt(id.vars=c("seed", "r", "rish.haplo.f", "n_family", "family.strct", "haplotype_on_affect", "p_dis"),
                                            variable.name="method", value.name="p.value")

# result <- melt(id.vars=c("seed", "r", "p_dis","rish.haplo.f","n_family", ""), variable.name="", value.name="p.value")
result.plot <- result %>% group_by(r, family.strct, method) %>% 
  summarise(n=n(), power=mean(p.value < 0.05))

#power for trap
ggplot(filter(result.plot, method=="trap"), aes(x=r, y=power, group=family.strct, col=family.strct)) +
  #   geom_point(size=3, alpha=1) +
  geom_line(size=1.2, alpha=0.7) +
  ggtitle("f=0.0102, consider effect size of risk haplotypes") +
  labs(x="relative risk r") +
  theme_gray(base_size = 20)

#power for trafic
ggplot(filter(result.plot, method=="trafic"), aes(x=r, y=power, group=family.strct, col=family.strct)) +
  #   geom_point(size=3, alpha=1) +
  geom_line(size=1.2, alpha=0.7) +
  ggtitle("f=0.0102, consider effect size of risk haplotypes") +
  labs(x="relative risk r") +
  theme_gray(base_size = 20)

#power for trap and trafic
ggplot(result.plot, aes(x=r, y=power, col=family.strct, lty=method)) +
  #   geom_point(size=3, alpha=1) +
  geom_line(size=1.2, alpha=0.7) +
  ggtitle("f=0.0102, consider effect size of risk haplotypes") +
  labs(x="relative risk r") +
  theme_gray(base_size = 20)