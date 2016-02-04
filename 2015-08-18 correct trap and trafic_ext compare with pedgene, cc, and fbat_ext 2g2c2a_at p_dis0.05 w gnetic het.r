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
# risk.variant.id <- sort(sample(which(snp$FREQ1>=0.001), 3))
risk.variant.id <- risk.variant #2 for common 7 for rare 39 for super rare
risk.haplo.id <- which(apply(2-as.matrix(haplotype[, risk.variant.id+2]), 1, sum)>0)
(risk.haplo.f <- mean(apply(2-as.matrix(haplotype[, risk.variant.id+2]), 1, sum)>0)) #carrier haplotype frequency
b <- rep(1, length=(ncol(haplotype)-2)) #initialize every haplotype is neutral
b[risk.variant.id] <- r #effect size is equal for every risk variant
#calculate the haplotype variants p(A|h)
haplotype.risk <<- apply(2-haplotype[, -c(1:2)], 1, function(x) prod(b^x))
#calculate the probability of drawing given the affected status of founder p(h|A)
# haplotype.draw_prob <- (haplotype.risk %o% haplotype.risk)/n_haplo^2/prevalence
(mean(haplotype.risk[risk.haplo.id])) #mean relative risk
(var(haplotype.risk[risk.haplo.id])) #variance relative risk
haplotype.risk <<- haplotype.risk*b0_sqrt

##gene drop simulation for two generations
family_strct.2g2c2a <- data.frame(family=c(1,1,1,1), person=c(1,2,3,4), father=c(0,0,1,1), 
                                mother=c(0,0,2,2), sex=c(1,2,1,1), affect=c(1,1,2,2)) #1=male, 2=female, 1=unaffected, 2=affected

##load library for pedgene and fbat_ext
##load library for pedgene and fbat_ext
.libPaths(c("/net/frodo/home/khlin/R/x86_64-pc-linux-gnu-library/3.2/", .libPaths()))
# library(Matrix, lib.loc = "../R/x86_64-pc-linux-gnu-library/3.2/")
# library(CompQuadForm, lib.loc = "../R/x86_64-pc-linux-gnu-library/3.2/")
# library(survey,lib.loc = "../R/x86_64-pc-linux-gnu-library/3.2/")
# library(quadprog,lib.loc = "../R/x86_64-pc-linux-gnu-library/3.2/")
# library(kinship2,lib.loc = "../R/x86_64-pc-linux-gnu-library/3.2/")
library(pedgene)
library(methods)


rep.idx <<- 1
sim_result <- replicate(n_rep, {
  print(rep.idx)
  rep.idx <<- rep.idx + 1

  family_generated_2g2c2a <<- gene_family(family_strct=family_strct.2g2c2a, n=n_family, haplotype.risk=haplotype.risk) 
  
  #add variant 3 since FB-SKAT can only work with 2+ in a gene
  family_generated_2g2c2a_diploid <- hap2dip(data=family_generated_2g2c2a, risk.variant.id=risk.variant, save.file=T)
  
  ##test based on new imputation framework
  result.trap.2g2c2a <- family.test(data=family_generated_2g2c2a, f=risk.variant)$p.value
  result.trafic_ext.2g2c2a <- family.test.trafic.ext(data=family_generated_2g2c2a, f=risk.variant)$p.value
  result.pedgene.2g2c2a <- pedgene(ped=family_generated_2g2c2a_diploid$ped, geno=family_generated_2g2c2a_diploid$geno)
  result.pedgene.vc.2g2c2a <- result.pedgene.2g2c2a$pgdf$pval.kernel
  result.pedgene.burden.2g2c2a <- result.pedgene.2g2c2a$pgdf$pval.burden
  
  system(paste("./FB-SKAT/FB-SKAT_v2.1 data_",seed,".ped variant_pass.txt genes.txt weight.txt results_", seed, "_ mendelian_errors_",seed,"_ 10000 1 0.01 0 0", sep=""))
  system(paste("./FB-SKAT/FB-SKAT_v2.1 data_",seed,".ped variant_pass.txt genes.txt weight.txt results_", seed, "_ mendelian_errors_",seed,"_ 10000 1 0.01 1 0", sep=""), ignore.stdout = TRUE)
  x <- read.table(paste("results_", seed,"_1.000000_0.010000_0.txt", sep="")) #variance component
  result.fbskat.vc.2g2c2a <- 1-pchisq(x[,4],x[,5])
  x <- read.table(paste("results_", seed,"_1.000000_0.010000_1.txt", sep="")) #burden
  result.fbskat.burden.2g2c2a=1-pchisq(x[,4],x[,5])
  
  #case-control
  data.cc <- gene_case_control(n_case_control_pair=2*n_family)
  result.cc <- case_control.test(data=data.cc, f=risk.variant)[3]
  
  #only report p.value
  c(result.trap.2g2c2a, result.trafic_ext.2g2c2a, 
    result.pedgene.vc.2g2c2a, result.pedgene.burden.2g2c2a,
    result.fbskat.vc.2g2c2a, result.fbskat.burden.2g2c2a, 
    result.cc)
})
## Write out your results to a csv file
result.df <- as.data.frame(t(sim_result))
colnames(result.df) <- c("result.trap.2g2c2a", "result.trafic_ext.2g2c2a", 
                         "result.pedgene.vc.2g2c2a", "result.pedgene.burden.2g2c2a",
                         "result.fbskat.vc.2g2c2a", "result.fbskat.burden.2g2c2a",
                         "result.cc")
result.df <- cbind(seed,r,p_dis,risk.variant=paste(risk.variant,collapse="-"),risk.haplo.f,n_family,result.df)
write.csv(result.df, paste("res_",r,"_",n_family,"_",seed,".csv",sep=""), row.names=FALSE)
## R.miSniam
## End(Not run)
## Not run:

###########################################
##passed in from the command prompt.
parseCommandArgs()
## Will disply the default values on the first run,
##initialize
seed=1000
#number of replications
n_rep=250
#number of family
n_family=1000
#prevalence
p_dis=0.05
# print(null) #null=FALSE
#family structure
# print(family_strct) #family_strct='family_strct.2g2c2a'

print(getwd())

##command line
parallel <- function(...) {
  names.args <- names(list(...))
  value.args <- as.vector(list(...))
  ##commands to run
  for(i in 1:6) {
    rfile="mainSim.R"
    cmd <- paste("R --vanilla --args seed ", seed, " ", paste(names.args, value.args, collapse=" "),
                 " < ", rfile, " > ", "mainSim_", paste(value.args, collapse="_"), ".Rout", seed, " 2>&1", sep="")
    print(cmd)
    cat(cmd, file=fileConn, sep="\n") #write jobs to text
    #     writeLines(cmd, fileConn) #write jobs to text
    
    #add seed number
    seed<<-seed+1
  }
}
##clean-up & initialization
system("rm *.csv")
system("rm *.Rout*")
system("rm *.out*")
system("rm *.sh")
system("rm result*.txt")
system("rm mendelian*.txt")
system("rm data*.ped")
fileConn<-file("jobs.txt", "w")

##create your jobs here
for(i in seq(1,2, length.out = 11)) {
  parallel(n_rep=n_rep, r=i, n_family=n_family, p_dis=p_dis, risk.variant="\"c(4,16,42)\"") #common 0.2135
}
for(i in seq(1,3, length.out = 11)) {
  parallel(n_rep=n_rep, r=i, n_family=n_family, p_dis=p_dis, risk.variant="\"c(4,15,32)\"") #rare 0.0116
}
for(i in seq(1,4, length.out = 11)) {
  parallel(n_rep=n_rep, r=i, n_family=n_family, p_dis=p_dis, risk.variant="\"c(18,25,47)\"") #super rare 0.0062
}

##write jobs.txt
close(fileConn)
##submit the jobs.txt using runslurm.pl
system("runslurm.pl -replace -logfile jobs.log -copts \"--time=72:00:00 --mem=4096\" jobs.txt")
##submit jobs.txt via runbatch.pl i.e. mosic mode
# dir=paste("/net/frodo", substr(getwd(), 0, 500), sep="")
# cmd <- paste("runbatch.pl -replace -logfile jobs.log -concurrent 40 -copts \"-E", dir," -j2,4-8\" jobs.txt &", sep="")
# print(cmd)
# system(cmd)
###################################################
# The color-blind friednly palette
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
# The palette with black:
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
#check power using true, imputed founder carrier, minimum offspring carrier
#three versions -- command and rare and super rare #2 for common, 7 for rare, 39 for super rare
result <- read.csv("2015-08-18 correct trap and trafic_ext compare with pedgene, cc, and fbat_ext 2g2c2a_at p_dis0.05 w gnetic het.csv", header=T)
names(result)[4] <- "risk.variant"
result <- result %>% melt(id.vars=c("seed", "r", "p_dis","risk.variant","risk.haplo.f","n_family"), variable.name="method", value.name="p.value")
result.plot <- result %>% group_by(risk.variant, risk.haplo.f, r, method) %>% 
  summarise(n=n(), power=mean(p.value<2.5*10^-6, na.rm=T))
result.plot$risk.haplo.f <- factor(result.plot$risk.haplo.f, labels=c("f=0.0062","f=0.0116","f=0.2135"))

#only burden test
pd <- position_dodge(0.0)
filter(result.plot, !grepl("vc", method)) %>% ggplot(aes(x=r, y=power, ymax=max(power), group=method, col=method)) +
  #   geom_point(size=3, alpha=1) +
  facet_wrap(~risk.haplo.f, ncol=3, scale="free_x") +
  geom_line(size=1.2, alpha=0.7, position=pd) +
  geom_point(size=1.2, position=pd) +
#   ggtitle("f=0.202, consider effect size of risk haplotypes, TRAP") +
#   ggtitle("f=0.0178, consider effect size of risk haplotypes, TRAP") +
#   ggtitle("f=0.0039, consider effect size of risk haplotypes, TRAP") +
  labs(x="relative risk r") +
  scale_y_continuous(limits=c(0,1)) +
  theme_gray(base_size = 20) +
  theme(legend.position="bottom", panel.background = element_rect(fill = 'grey85')) +
  scale_color_manual(values=cbbPalette)

#both burden and variance component test
pd <- position_dodge(0.0)
filter(result.plot) %>% ggplot(aes(x=r, y=power, ymax=max(power), group=method, col=method)) +
  #   geom_point(size=3, alpha=1) +
  facet_wrap(~risk.haplo.f, ncol=3, scale="free_x") +
  geom_line(size=1.2, alpha=0.7, position=pd) +
  geom_point(size=1.2, position=pd) +
  #   ggtitle("f=0.202, consider effect size of risk haplotypes, TRAP") +
  #   ggtitle("f=0.0178, consider effect size of risk haplotypes, TRAP") +
  #   ggtitle("f=0.0039, consider effect size of risk haplotypes, TRAP") +
  labs(x="relative risk r") +
  scale_y_continuous(limits=c(0,1)) +
  theme_gray(base_size = 20) +
  theme(legend.position="right", panel.background = element_rect(fill = 'grey85')) +
  scale_color_manual(values=cbbPalette)
