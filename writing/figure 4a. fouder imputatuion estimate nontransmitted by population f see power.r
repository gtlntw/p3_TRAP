## ... your simualtion codek
source("/net/frodo/home/khlin/p3/function.r")

#need to take in parameter
## passed in from the command prompt.
parseCommandArgs()
## Will disply the default values on the first run,
print(seed)
print(n_rep)
print(r)
print(n_family) #number of family
print(p_dis) #prevelance
# print(risk.variant)
print(family_strct)
print(f) #snps with f to consider

##1. read in cosi haplotype
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
n_snp <- nrow(snp)

##2. gives structure and 
family_strct_ped <- select_family_strct(family_strct)

##3.load library for pedgene and fbat_ext
.libPaths(c("/net/frodo/home/khlin/R/x86_64-pc-linux-gnu-library/3.2/", .libPaths()))
library(pedgene)
library(methods)
library(SKAT)
n_cc <- length(family_strct_ped$person)/2*n_family

##4. effect size
# risk.variant.id <- c(5,15,16,19,20,27,31,35,43,49)
(n_snp_f <- sum(snp$FREQ1<f))
(n_snp_f_causal <- round(n_snp_f*0.1)) #10 percent of snp with maf < f are causal
# risk.variant.id <<- sort(sample(which(0.0001<snp$FREQ1 & snp$FREQ1<f), round(n_snp_f_causal)))
# risk.haplo.id <- which(apply(2-as.matrix(haplotype[, risk.variant.id+2]), 1, sum)>0)
risk.variant.id <- switch(as.character(f),
                          "0.01"=c(4,37,38),
                          "0.05"=c(12,20,37,44),
                          "0.2"=c(1,4,22,26)
)
risk.haplo.f <<- mean(apply(2-as.matrix(haplotype[, risk.variant.id+2]), 1, sum)>0) #carrier haplotype frequency
cat("risk.haplo.f = ", risk.haplo.f, "\n")
Beta <- rep(0, length.out=50)
Beta[risk.variant.id] <- log(r) #effect size of OR=2

##5 .create genes.txt weight.txt variant_pass.txt for fb-skat
write(paste("geneA 1", n_snp_f_causal), paste0("genes_",seed,".txt"))
write(rep("1", n_snp_f_causal), paste0("weight_",seed,".txt"), sep="\n")
write(rep("1", n_snp_f_causal), paste0("variant_pass_",seed,".txt"), sep="\n")

sim_result <- list()
for(i in 1:n_rep) {
  print(i)
  sim.fail <- F
  
  family_generated <<- gene_family_pe(family_strct=family_strct_ped, n_family=n_family, p_dis=p_dis, Beta=Beta)
  
  #remove the founder from the data
  family_generated_c <- family_generated
  temp.idx <- which(family_generated$data_family$person %in% 1:2) #take out founders
  family_generated_c$data_family <- family_generated$data_family[-temp.idx, ]
  family_generated_c$tran_vec <- family_generated$tran_vec[-temp.idx, ]
  
  ##test based on new imputation framework
  result.trap.3c.noimpute <- family.test(data=family_generated, f=risk.variant.id, nofounderphenotype=T)$p.value
  result.trap.3c.impute.founder.pop.f <- family.test.nofounder.impute(data=family_generated_c, f=risk.variant.id, sample.f=F, pop.f.off=0)
  result.trap.3c.impute.founder.pop.f.off.10 <- family.test.nofounder.impute(data=family_generated_c, f=risk.variant.id, sample.f=F, pop.f.off=10)
  result.trap.3c.impute.founder.pop.f.off.20 <- family.test.nofounder.impute(data=family_generated_c, f=risk.variant.id, sample.f=F, pop.f.off=20)
  result.trap.3c.impute.founder.pop.f.off.50 <- family.test.nofounder.impute(data=family_generated_c, f=risk.variant.id, sample.f=F, pop.f.off=50)
  result.trap.3c.impute.founder.pop.f.off.n10 <- family.test.nofounder.impute(data=family_generated_c, f=risk.variant.id, sample.f=F, pop.f.off=-10)
  result.trap.3c.impute.founder.pop.f.off.n20 <- family.test.nofounder.impute(data=family_generated_c, f=risk.variant.id, sample.f=F, pop.f.off=-20)
  result.trap.3c.impute.founder.pop.f.off.n50 <- family.test.nofounder.impute(data=family_generated_c, f=risk.variant.id, sample.f=F, pop.f.off=-50)
  result.trap.3c.impute.founder.sample.f <- family.test.nofounder.impute(data=family_generated_c, f=risk.variant.id, sample.f=T)
  
  #skip return when errors occured
  if(sim.fail==F) {
    #only report p.value
    sim_result[[i]] <- data.frame(result.trap.3c.noimpute, result.trap.3c.impute.founder.pop.f,
                                  result.trap.3c.impute.founder.pop.f.off.10,
                                  result.trap.3c.impute.founder.pop.f.off.20, result.trap.3c.impute.founder.pop.f.off.50,
                                  result.trap.3c.impute.founder.pop.f.off.n10,
                                  result.trap.3c.impute.founder.pop.f.off.n20, result.trap.3c.impute.founder.pop.f.off.n50,
                                  result.trap.3c.impute.founder.sample.f) 
  } else {
    warning("error occured!!")
    next
  }
}
##remove junk files produced by fbskat
system(paste("rm results_",seed,"*.txt", sep=""))
system(paste("rm mendelian_errors_",seed,"*.txt", sep=""))
system(paste("rm data_",seed,"*.ped", sep=""))
system(paste("rm genes_",seed,"*.txt", sep=""))
system(paste("rm weight_",seed,"*.txt", sep=""))
system(paste("rm variant_pass_",seed,"*.txt", sep=""))
## Write out your results to a csv file
result.df <- do.call(rbind, sim_result)
result.df <- cbind(seed,f,r,p_dis,risk.variant.id=paste(c(risk.variant.id), collapse = "_"),risk.haplo.f,n_family,family_strct,result.df)
write.csv(result.df, paste("res_",seed,"_",r,"_",f,"_",n_family,".csv",sep=""), row.names=FALSE)
## R.miSniam
## End(Not run)
## Not run:

###############################################################################
####################commands to run jobs in parallel
##initialize
opts['seed']=1000 #starting seed number
opts['n_rep_total']=1000 #total number of replications
opts['n_rep']=50 #number of replications in each parallele job
opts['n_ite']=opts['n_rep_total']/opts['n_rep'] #number of parallele jobs
opts['n_family']=1000 #number of family
opts['p_dis']=0.1 #prevalence

######################
#1.0. log the start time
######################
tgt = '{outputDir}/start.runmake.{jobName}.OK'.format(**opts)
dep = ''
cmd = ['[ ! -f {outputDir}/runmake_{jobName}_time.log ] && echo > {outputDir}/runmake_{jobName}_time.log; date | awk \'{{print "Simulations pipeline\\n\\nstart: "$$0}}\' >> {outputDir}/runmake_{jobName}_time.log'.format(**opts)]
makeJob('local', tgt, dep, cmd)

opts["exclude"] = "--exclude=../exclude_node.txt"
opts["param"] = "--time=1-12:0 {exclude}".format(**opts) #indicate this is a quick job
######################
#1.1. run simulations by calling mainSim.R
######################
inputFilesOK = []
opts['f'] = 0.01 #super rare
opts['family_strct'] = '\"2g.3a.2u\"' #family structure
for i in numpy.arange(1,2.4,0.1):
	for j in range(opts['n_ite']):
		opts['r'] = i
		tgt = 'callmainSim_{seed}.OK'.format(**opts)
		inputFilesOK.append(tgt)
		dep = ''
		cmd = ['R --vanilla --args seed {seed} n_rep {n_rep} r {r} f {f} n_family {n_family} p_dis {p_dis} family_strct {family_strct} < mainSim.R > mainSim{seed}_{n_rep}_{r}_{n_family}_{p_dis}_{f}_{family_strct}.Rout 2>&1'.format(**opts)]
		makeJob(opts['launchMethod'], tgt, dep, cmd)
		opts['seed'] += 1	

opts['f'] = 0.05 #less rare
for i in numpy.arange(1,1.8,0.1):
	for j in range(opts['n_ite']):
		opts['r'] = i
		tgt = 'callmainSim_{seed}.OK'.format(**opts)
		inputFilesOK.append(tgt)
		dep = ''
		cmd = ['R --vanilla --args seed {seed} n_rep {n_rep} r {r} f {f} n_family {n_family} p_dis {p_dis} family_strct {family_strct} < mainSim.R > mainSim{seed}_{n_rep}_{r}_{n_family}_{p_dis}_{f}_{family_strct}.Rout 2>&1'.format(**opts)]
		makeJob(opts['launchMethod'], tgt, dep, cmd)
		opts['seed'] += 1	

opts['f'] = 0.20 #common
for i in numpy.arange(1,1.6,0.1):
	for j in range(opts['n_ite']):
		opts['r'] = i
		tgt = 'callmainSim_{seed}.OK'.format(**opts)
		inputFilesOK.append(tgt)
		dep = ''
		cmd = ['R --vanilla --args seed {seed} n_rep {n_rep} r {r} f {f} n_family {n_family} p_dis {p_dis} family_strct {family_strct} < mainSim.R > mainSim{seed}_{n_rep}_{r}_{n_family}_{p_dis}_{f}_{family_strct}.Rout 2>&1'.format(**opts)]
		makeJob(opts['launchMethod'], tgt, dep, cmd)
		opts['seed'] += 1	

######################
#1.2. combine the result
######################
tgt = 'pasteResults.OK'
dep = ' '.join(inputFilesOK)
cmd = ['python paste_mainSim_results.py']
makeJob('local', tgt, dep, cmd)

######################
#1.3. copy the results to a folder with name $jobName
######################
tgt = 'backupResult.OK'
dep = 'pasteResults.OK'
cmd = ['cp * {jobName}/'.format(**opts)]
makeJob('local', tgt, dep, cmd)

######################
#1.4. log the end time
######################
tgt = '{outputDir}/end.runmake.{jobName}.OK'.format(**opts)
dep = 'pasteResults.OK'
cmd = ['date | awk \'{{print "\\n\\nend: "$$0}}\' >> {outputDir}/runmake_{jobName}_time.log'.format(**opts)]
makeJob('local', tgt, dep, cmd)

###################################################
# The color-blind friednly palette
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
# The palette with black:
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
#check power using true, imputed founder carrier, minimum offspring carrier
#three versions -- command and rare and super rare #2 for common, 7 for rare, 39 for super rare
result <- read.csv("figure 4a. fouder imputatuion estimate nontransmitted by population f see power.csv", header=T)
result <- result %>% gather(key="method", value="p.value", -c(1:8))
result.plot <- result %>% group_by(f, risk.haplo.f, r, method) %>% 
  summarise(n=n(), power=mean(p.value<0.05, na.rm=T))
result.plot$method <- factor(result.plot$method, levels=c("result.trap.3c.noimpute","result.trap.3c.impute.founder.pop.f","result.trap.3c.impute.founder.pop.f.off.10",
																													"result.trap.3c.impute.founder.pop.f.off.20","result.trap.3c.impute.founder.pop.f.off.50",
																													"result.trap.3c.impute.founder.pop.f.off.n10","result.trap.3c.impute.founder.pop.f.off.n20",
																													"result.trap.3c.impute.founder.pop.f.off.n50","result.trap.3c.impute.founder.sample.f"),
														 labels=c("w.o. imputation", "pop.f", "pop.f off 10%",
                                                "pop.f off 20%", "pop.f off 50%", "pop.f off -10%"
                                                , "pop.f off -20%", "pop.f off -50%", "sample.f"))

#only include 20%off
pd <- position_dodge(0.0)
filter(result.plot, !grepl("conditional|10%|5%|50%", method)) %>% ggplot(aes(x=r, y=power, ymax=max(power), group=method, col=method)) +
  #   geom_point(size=3, alpha=1) +
  facet_wrap(~f, ncol=3, scale="free_x") +
  geom_line(size=1.2, alpha=0.7, position=pd) +
  geom_point(size=1.2, position=pd) +
  #   ggtitle("f=0.202, consider effect size of risk haplotypes, TRAP") +
  #   ggtitle("f=0.0178, consider effect size of risk haplotypes, TRAP") +
  #   ggtitle("f=0.0039, consider effect size of risk haplotypes, TRAP") +
  labs(x="odds ratio r") +
  scale_y_continuous(limits=c(0,1)) +
  theme_gray(base_size = 20) +
  theme(legend.position="right", panel.background = element_rect(fill = 'grey85')) +
  scale_color_manual(values=cbbPalette)

#everything
pd <- position_dodge(0.0)
filter(result.plot, !grepl("conditional", method)) %>% ggplot(aes(x=r, y=power, ymax=max(power), group=method, col=method)) +
  #   geom_point(size=3, alpha=1) +
  facet_wrap(~f, ncol=3, scale="free_x") +
  geom_line(size=1.2, alpha=0.7, position=pd) +
  geom_point(size=1.2, position=pd) +
  #   ggtitle("f=0.202, consider effect size of risk haplotypes, TRAP") +
  #   ggtitle("f=0.0178, consider effect size of risk haplotypes, TRAP") +
  #   ggtitle("f=0.0039, consider effect size of risk haplotypes, TRAP") +
  labs(x="odds ratio r") +
  scale_y_continuous(limits=c(0,1)) +
  theme_gray(base_size = 20) +
  theme(legend.position="right", panel.background = element_rect(fill = 'grey85'))
