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
risk.variant.id <- c(risk.variant) #2 for common 7 for rare 39 for super rare
risk.haplo.id <- which(apply(2-as.matrix(haplotype[, risk.variant.id+2]), 1, sum)>0)
(risk.haplo.f <- mean(apply(2-as.matrix(haplotype[, risk.variant.id+2]), 1, sum)>0)) #carrier haplotype frequency
haplotype.risk <- rep(1, length=nrow(haplotype))
#assign mean relative risk calculate the haplotype variants p(A|h)
haplotype.risk[risk.haplo.id] <- r
mean(haplotype.risk[risk.haplo.id]) #mean relative risk
haplotype.risk <<- haplotype.risk*b0_sqrt

##gene drop simulation for two generations
family_strct.2g2c <- data.frame(family=c(1,1,1,1), person=c(1,2,3,4), father=c(0,0,1,1), 
                                mother=c(0,0,2,2), sex=c(1,2,1,1), affect=c(1,2,2,2)) #1=male, 2=female, 1=unaffected, 2=affected


rep.idx <<- 1
sim_result <- replicate(n_rep, {
  print(rep.idx)
  rep.idx <<- rep.idx + 1

  family_generated_2g2c <<- gene_family(family_strct=family_strct.2g2c, n=n_family, haplotype.risk=haplotype.risk) 
  
  #remove the founder from the data
  family_generated_2c <- family_generated_2g2c
  family_generated_founder <- list()
  temp.idx <- which(family_generated_2g2c$data_family$person %in% 1:2) #take out founders
  family_generated_founder$data_family <- family_generated_2g2c$data_family[temp.idx, ]
  family_generated_founder$tran_vec <- family_generated_2g2c$tran_vec[temp.idx, ]
  family_generated_2c$data_family <- family_generated_2g2c$data_family[-temp.idx, ]
  family_generated_2c$tran_vec <- family_generated_2g2c$tran_vec[-temp.idx, ]

  family_generated_2g2c_subsample <- family_generated_2g2c
  subsample.idx <- 1:((n_family*4)/2)
  family_generated_2g2c_subsample$data_family <- family_generated_2g2c$data_family[subsample.idx, ]
  family_generated_2g2c_subsample$tran_vec <- family_generated_2g2c$tran_vec[subsample.idx, ]
  
#   family_generated_2c_subsample <- family_generated_2c
#   subsample.idx <- 1:((n_family*2)/2)
#   family_generated_2c_subsample$data_family <- family_generated_2c$data_family[subsample.idx, ]
#   family_generated_2c_subsample$tran_vec <- family_generated_2c$tran_vec[subsample.idx, ]
    
  ##test to compare
  result.trafic.2c <- trafic.test(data=family_generated_2c, f=risk.variant.id, true.only = T)$p.value
#   result.trafic.2c.sub <- trafic.test(data=family_generated_2c_subsample, f=risk.variant.id, true.only = T)$p.value
  result.trafic.ext.2g2c <- family.test.trafic.ext(data=family_generated_2g2c_subsample, f=risk.variant.id, nofounderpheno=T)
  result.trap.2g2c <- family.test.nofounderphenotype(data=family_generated_2g2c_subsample, f=risk.variant.id)[7]
  
  
  ##generate case-control samples
  data.cc <- gene_case_control(n_case_control_pair=n_family)
  result.cc <- case_control.test(data=data.cc, f=risk.variant)[3]
  
  #only report p.value
  c(result.trafic.2c, result.trafic.ext.2g2c, result.trap.2g2c, result.cc)
})
## Write out your results to a csv file
result.df <- as.data.frame(t(sim_result))
colnames(result.df) <- c("result.trafic.2c", "result.trafic.ext.2g2c", "result.trap.2g2c", "result.cc")
result.df <- cbind(seed,r,p_dis,risk.variant,risk.haplo.f,n_family,result.df)
write.csv(result.df, paste("res_",r,"_",n_family,"_",seed,".csv",sep=""), row.names=FALSE)
## R.miSniam
## End(Not run)
## Not run:

####################commands to run jobs in parallel
##passed in from the command prompt.
parseCommandArgs()
## Will disply the default values on the first run,
##initialize
seed=1000
#number of replications
n_rep=500
#number of family
n_family=1000
#prevalence
p_dis=0.15
# print(null) #null=FALSE
#family structure
# print(family_strct) #family_strct='family_strct.2g2c'

print(getwd())

##command line
parallel <- function(...) {
  names.args <- names(list(...))
  value.args <- as.vector(list(...))
  ##commands to run
  for(i in 1:8) {
    rfile="mainSim.R"
    cmd <- paste("R --vanilla --args seed ", seed, " ", paste(names.args, value.args, collapse=" "),
                 " < ", rfile, " > ", "mainSim_", paste(value.args, collapse="_"), ".Rout", seed, " 2>&1", sep="")
    print(cmd)
    writeLines(cmd, fileConn) #write jobs to text
    
    #add seed number
    seed<<-seed+1
  }
}
##clean-up & initialization
system("rm *.csv")
system("rm *.Rout*")
system("rm *.out*")
fileConn<-file("jobs.txt", "w")

##create your jobs here
for(i in seq(1,2, length.out = 11)) {
  parallel(n_rep=n_rep, r=i, n_family=n_family, p_dis=p_dis, risk.variant=2) #common
}
for(i in seq(1,2.5, length.out = 11)) {
  parallel(n_rep=n_rep, r=i, n_family=n_family, p_dis=p_dis, risk.variant=7) #rare
}
for(i in seq(1,3.5, length.out = 11)) {
  parallel(n_rep=n_rep, r=i, n_family=n_family, p_dis=p_dis, risk.variant=39) #super rare
}

##write jobs.txt
close(fileConn)
##submit the jobs.txt using runslurm.pl
system("runslurm.pl -replace -logfile jobs.log -copts \"--time=16:00:00 --mem=512\" jobs.txt")
##submit jobs.txt via runbatch.pl i.e. mosic mode
# dir=paste("/net/frodo", substr(getwd(), 0, 500), sep="")
# cmd <- paste("runbatch.pl -replace -logfile jobs.log -concurrent 40 -copts \"-E", dir," -j2,4-8\" jobs.txt &", sep="")
# print(cmd)
# system(cmd)
###################################################
#check power using true, imputed founder carrier, minimum offspring carrier
#three versions -- command and rare and super rare #2 for common, 7 for rare, 39 for super rare
result <- read.csv("p1 2015-06-20 trafic, trafic_ext, trap, cc power comparison.csv", header=T)
result <- result %>% melt(id.vars=c("seed", "r", "p_dis","risk.variant","risk.haplo.f","n_family"), variable.name="method", value.name="p.value")
result.plot <- result %>% group_by(risk.variant, risk.haplo.f, r, method) %>% 
  summarise(n=n(), power=mean(p.value<2.5*10^-4, na.rm=T))
result.plot$risk.haplo.f <- factor(result.plot$risk.haplo.f, labels=c("f=0.0039","f=0.0178","f=0.202"))
result.plot$method <- factor(result.plot$method, label=c("TRAFIC", "TRAFIC_EXT", "TRAP", "CC"))

#CSG presentation
pd <- position_dodge(0.0)
filter(result.plot, risk.variant!=7) %>% ggplot(aes(x=r, y=power, ymax=max(power), group=method, col=method)) +
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
  theme(legend.position="bottom", panel.background = element_rect(fill = 'grey85'))

#temp keep
pd <- position_dodge(0.0)
filter(result.plot, risk.variant %in% c(2,39))%>% 
  ggplot(aes(x=r, y=power, ymax=max(power), group=method, col=method)) +
  facet_grid(.~risk.haplo.f,labeller = label_bquote(expr=f == .(x)), scale="free_x") + 
  #   geom_point(size=3, alpha=1) +
  geom_line(size=1.2, alpha=0.7, position=pd) +
  geom_point(size=1.2, position=pd) +
  #   ggtitle("f=0.202, consider effect size of risk haplotypes, TRAP") +
  #   ggtitle("f=0.0178, consider effect size of risk haplotypes, TRAP") +
  #   ggtitle("f=0.0039, consider effect size of risk haplotypes, TRAP") +
  labs(x="relative risk r") +
  coord_cartesian(ylim=c(0, 1))+
  theme_gray(base_size = 20)