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

rep.idx <<- 1
sim_result <- replicate(n_rep, {
  print(rep.idx)
  rep.idx <<- rep.idx + 1
  
  family_generated <<- gene_family_pe(family_strct=family_strct_ped, n_family=n_family, p_dis=p_dis, Beta=Beta)
  
  #add variant 3 since FB-SKAT can only work with 2+ in a gene
  family_generated_diploid <- hap2dip(data=family_generated, risk.variant.id=risk.variant.id, save.file=T)
  
  ##test based on new imputation framework
  result.trap <- family.test(data=family_generated, f=risk.variant.id)$p.value
  result.trafic_ext <- family.test.trafic.ext(data=family_generated, f=risk.variant.id)$p.value
  result.pedgene <- pedgene(ped=family_generated_diploid$ped, geno=family_generated_diploid$geno)
  result.pedgene.vc <- result.pedgene$pgdf$pval.kernel
  result.pedgene.burden <- result.pedgene$pgdf$pval.burden
  
  system(paste("/net/frodo/home/khlin/p3/FB-SKAT/FB-SKAT_v2.1 data_",seed,".ped variant_pass_",seed,".txt genes_",seed,".txt weight_",seed,".txt results_", seed, "_ mendelian_errors_",seed,"_ 10000 1 0.01 0 0", sep=""))
  system(paste("/net/frodo/home/khlin/p3/FB-SKAT/FB-SKAT_v2.1 data_",seed,".ped variant_pass_",seed,".txt genes_",seed,".txt weight_",seed,".txt results_", seed, "_ mendelian_errors_",seed,"_ 10000 1 0.01 1 0", sep=""), ignore.stdout = TRUE)
  x <- read.table(paste("results_", seed,"_1.000000_0.010000_0.txt", sep="")) #variance component
  result.fbskat.vc <- 1-pchisq(x[,4],x[,5])
  x <- read.table(paste("results_", seed,"_1.000000_0.010000_1.txt", sep="")) #burden
  result.fbskat.burden=1-pchisq(x[,4],x[,5])
  
  #case-control use burden test with equal weight
  data.cc <- gene_case_control_pe(n_case_control_pair=n_cc, p_dis=p_dis, Beta=Beta)
  data.cc.diploid <- hap2dip(data=list(data_family=data.cc), risk.variant.id=risk.variant.id, save.file=F)
  obj <- SKAT_Null_Model(data.cc.diploid$ped$trait ~ 1, out_type="D")
  result.cc <- SKAT(as.matrix(data.cc.diploid$geno[, -c(1:2)]), obj, weights.beta = c(1,1), r.corr = 1)$p.value
  
  #only report p.value
  c(result.trap, result.trafic_ext, 
    result.pedgene.vc, result.pedgene.burden,
    result.fbskat.vc, result.fbskat.burden, 
    result.cc)
})
##remove junk files produced by fbskat
system(paste("rm results_",seed,"*.txt", sep=""))
system(paste("rm mendelian_errors_",seed,"*.txt", sep=""))
system(paste("rm data_",seed,"*.ped", sep=""))
system(paste("rm genes_",seed,"*.ped", sep=""))
system(paste("rm weight_",seed,"*.ped", sep=""))
system(paste("rm variant_pass_",seed,"*.ped", sep=""))
## Write out your results to a csv file
result.df <- as.data.frame(t(sim_result))
colnames(result.df) <- c("result.trap", "result.trafic_ext", 
                         "result.pedgene.vc", "result.pedgene.burden",
                         "result.fbskat.vc", "result.fbskat.burden",
                         "result.cc")
result.df <- cbind(seed,f,r,p_dis,risk.variant.id=paste(c(risk.variant.id), collapse = "_"),risk.haplo.f,n_family,family_strct,result.df)
write.csv(result.df, paste("res_",seed,"_",r,"_",f,"_",n_family,".csv",sep=""), row.names=FALSE)
## R.miSniam
## End(Not run)
## Not run:

###########################################
##initialize
opts['seed']=1000 #starting seed number
opts['n_rep_total']=1000000 #total number of replications
opts['n_rep']=100 #number of replications in each parallele job
opts['n_ite']=opts['n_rep_total']/opts['n_rep'] #number of parallele jobs
opts['n_family']=1000 #number of family
opts['p_dis']=0.4 #prevalence

######################
#1.0. log the start time
######################
tgt = '{outputDir}/start.runmake.{jobName}.OK'.format(**opts)
dep = ''
cmd = ['[ ! -f {outputDir}/runmake_{jobName}_time.log ] && echo > {outputDir}/runmake_{jobName}_time.log; date | awk \'{{print "Simulations pipeline\\n\\nstart: "$$0}}\' >> {outputDir}/runmake_{jobName}_time.log'.format(**opts)]
makeJob('local', tgt, dep, cmd)

opts["param"] = "--time=1-12:0 --exclude=dl3601" #indicate this is a quick job
######################
#1.1. run simulations by calling mainSim.R
######################
inputFilesOK = []
opts['family_strct'] = '\"2g.3a.1u\"' #family structure
f = [0.01, 0.20] #rare and common
for h in f:
	for i in numpy.arange(1,1.05,0.1):
		for j in range(opts['n_ite']):
			opts['f'] = h
			opts['r'] = i
			tgt = 'callmainSim_{seed}.OK'.format(**opts)
			inputFilesOK.append(tgt)
			dep = ''
			cmd = ['R --vanilla --args seed {seed} n_rep {n_rep} r {r} f {f} n_family {n_family} p_dis {p_dis} family_strct {family_strct} < mainSim.R > mainSim{seed}_{n_rep}_{r}_{n_family}_{p_dis}_{f}_{family_strct}.Rout 2>&1'.format(**opts)]
			makeJob(opts['launchMethod'], tgt, dep, cmd)
			opts['seed'] += 1	

opts["param"] = "--mem-per-cpu=4096 --time=1-12:0 --exclude=dl3601" #indicate this is a quick job
opts['family_strct'] = '\"3g.3a.4u\"' #family structure
for h in f:
	for i in numpy.arange(1,1.05,0.1):
		for j in range(opts['n_ite']):
			opts['f'] = h
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
cmd = ['Rscript paste_mainSim_results.R']
makeJob('local', tgt, dep, cmd)

######################
#1.3. delete unnecesary files
######################
tgt = 'delTempFile.OK'
dep = 'pasteResults.OK'
cmd = ['rm -f result*.txt mendelian*.txt gene*.txt weight*.txt variant_pass*.txt data*.ped'.format(**opts)]
makeJob('local', tgt, dep, cmd)

######################
#1.4. copy the results to a folder with name $jobName
######################
tgt = 'backupResult.OK'
dep = 'delTempFile.OK'
cmd = ['cp * {jobName}/'.format(**opts)]
makeJob('local', tgt, dep, cmd)

######################
#1.5. log the end time
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
result <- read.csv("Table 1. type I error.csv", header=T)
result <- result %>% gather(key="method", value="p.value", -c(1:8))
result.plot <- result %>% group_by(f, risk.haplo.f, family_strct, r, method) %>% 
  summarise(n=n(), power=mean(p.value<2.5*10^-4, na.rm=T))
as.data.frame(result.plot)
# result.plot$risk.haplo.f <- factor(result.plot$risk.haplo.f, labels=c("f=0.0039","f=0.0178","f=0.202"))

#