source("function.r")
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
null <- TRUE
n_haplo <- 10000
n_snp <- ncol(haplotype)-2
prevalence <- 0.05
b0_sqrt <- sqrt(prevalence)   #baseline
temp.idx <- sort(c(sample(which((snp$FREQ1 < 0.05)==T & snp$FREQ1 != 1/n_haplo ), 5), sample(which(snp$FREQ1 == 1/n_haplo), 5))) #randomly pick 5 variants with f < 0.05 and 5 singletons as casual 
b <- rep(1, length=(ncol(haplotype)-2))
if(null==FALSE) b[temp.idx] <- abs(log10(snp$FREQ1[temp.idx])) #effect size is log10 of its allele frequency
#calculate the haplotype variants p(A|h)
haplotype.risk <- apply(2-haplotype[, -c(1:2)], 1, function(x) prod(b^x)*b0_sqrt)
#calculate the probability of drawing given the affected status of founder p(h|A)
# haplotype.draw_prob <- (haplotype.risk %o% haplotype.risk)/n_haplo^2/prevalence
mean(b) #mean relative risk
mean(apply(2-haplotype[, temp.idx+2], 1, sum)>0) #carrier haplotype frequency

##gene drop simulation for three generations
family_strct.3g <- data.frame(family=c(1,1,1,1,1,1), person=c(1,2,3,4,5,6), father=c(0,0,1,1,0,4), 
                                   mother=c(0,0,2,2,0,5), sex=c(1,2,1,1,2,1), affect=c(1,2,2,2,1,2)) #1=male, 2=female, 1=unaffected, 2=affected
family_strct.2g3c <- data.frame(family=c(1,1,1,1,1), person=c(1,2,3,4,5), father=c(0,0,1,1,1), 
                              mother=c(0,0,2,2,2), sex=c(1,2,1,1,1), affect=c(1,2,2,2,2)) #1=male, 2=female, 1=unaffected, 2=affected


##simulation with 1000 replications
seed=1000
f=0.01
n_pair=1000
rep=2
rep.idx <<- 1
sim_result <- replicate(rep, {
  print(rep.idx)
  rep.idx <<- rep.idx + 1
  case_control_generated <- gene_case_control(n=n_pair)
  case_control.test(data=case_control_generated, f=f)
})

write.csv(data.frame(f=f, prop.1=sim_result[1,], prop.2=sim_result[2,], p.value=sim_result[3,]),
          paste("res_",f,"_",n_pair,"_",seed,"_case_control",".csv",sep=""), row.names=FALSE)

rep.idx <<- 1
sim_result <- replicate(n_rep, {
  print(rep.idx)
  rep.idx <<- rep.idx + 1
  family_generated <- gene_family(family_strct=family_strct.2g3c, n=400) 
  family.test(data=family_generated, f=f)
})

##3generations but consider only 5 ppl
rep.idx <<- 1
sim_result <- replicate(2, {
  print(rep.idx)
  rep.idx <<- rep.idx + 1
  family_generated <- gene_family(family_strct=family_strct.3g, n=400)
  temp.idx <- which(family_generated$data_family$person==5) #take out 5th person
  family_generated$data_family <- family_generated$data_family[-temp.idx, ]
  family_generated$tran_vec <- family_generated$tran_vec[-temp.idx, ]
#   family_generated$data_family[which(family_generated$data_family$person==6), "person"] <- 5
  family.test(data=family_generated, f=f)
})



#plot for the results under the null
alpha=0.05
result.cc_null <- read.csv("2014-09-30 simulation for family data 1kb region_null_case_control.csv", stringsAsFactors=FALSE)
result.family_null <- read.csv("2014-09-30 simulation for family data 1kb region_null_family.csv", stringsAsFactors=FALSE)
result.2g3c_null <- read.csv("2014-10-03 simulation for family data 1kb region more generations and children_2g3c_null.csv", stringsAsFactors=FALSE)
# result.3g_null <- read.csv("2014-10-03 simulation for family data 1kb region more generations and children_3g_null.csv", stringsAsFactors=FALSE)
result.3g_5ppl_null <- read.csv("2014-10-03 simulation for family data 1kb region more generations and children_3g_5ppl_null.csv", stringsAsFactors=FALSE)
(freq <- c(mean((result.cc_null$p.value < alpha)), mean((result.family_null$p.value < alpha)), mean((result.2g3c_null$p.value < alpha)),  mean((result.3g_5ppl_null$p.value < alpha))))
xx <- barplot(freq, names=c("case-control", "family 2:2:0", "family 2:3:0", "family 2:2:1"), ylim=c(0,1), ylab="power")
text(x = xx, y = freq, label=freq, pos = 3, cex = 1, col = "red")

#plot for the results under the alternative
alpha=0.05
result.cc <- read.csv("2014-09-30 simulation for family data 1kb region_case_control.csv", stringsAsFactors=FALSE)
result.family <- read.csv("2014-09-30 simulation for family data 1kb region_family.csv", stringsAsFactors=FALSE)
result.2g3c <- read.csv("2014-10-03 simulation for family data 1kb region more generations and children_2g3c.csv", stringsAsFactors=FALSE)
# result.3g <- read.csv("2014-10-03 simulation for family data 1kb region more generations and children_3g.csv", stringsAsFactors=FALSE)
result.3g_5ppl <- read.csv("2014-10-03 simulation for family data 1kb region more generations and children_3g_5ppl.csv", stringsAsFactors=FALSE)
(freq <- c(mean((result.cc$p.value < alpha)), mean((result.family$p.value < alpha)), mean((result.2g3c$p.value < alpha)), mean((result.3g_5ppl$p.value < alpha))))
xx <- barplot(freq, names=c("case-control", "family 2:2:0", "family 2:3:0", "family 2:2:1"), ylim=c(0,1), ylab="power")
text(x = xx, y = freq, label=freq, pos = 3, cex = 1, col = "red")

#effect size
plot(seq(0.0001, 0.05, by=0.001), -log10(seq(0.0001, 0.05, by=0.001)), type="l", xlab="MAF", ylab="Reltative Risk", ylim=c(1,4.5), lwd=2)

