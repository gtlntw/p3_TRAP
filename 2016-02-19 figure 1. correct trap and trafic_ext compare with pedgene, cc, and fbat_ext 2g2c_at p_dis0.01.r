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
  
  #simulation for two-generation families
  family_generated <<- gene_family_pe(family_strct=family_strct_ped, n_family=n_family, p_dis=p_dis, Beta=Beta, exact_affected = F) 
  #add pseudo null variant since FB-SKAT can only work with 2+ in a gene
  family_generated_diploid <- hap2dip(data=family_generated, risk.variant.id=risk.variant.id, save.file=T)
  
  #run tests
  result.trap <- family.test(data=family_generated, f=risk.variant.id)$p.value
  result.trafic.ext <- family.test.trafic.ext(data=family_generated, f=risk.variant.id)$p.value
  sim.fail <- tryCatch({  
      result.pedgene <- pedgene(ped=family_generated_diploid$ped, geno=family_generated_diploid$geno)
      print(result.pedgene)
      result.pedgene.vc <- result.pedgene$pgdf$pval.kernel
      result.pedgene.burden <- result.pedgene$pgdf$pval.burden
      F
    }, 
    error = function(e) return(T) 
  )
  
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

  #skip return when errors occured
  if(sim.fail==F) {
    #only report p.value
    sim_result[[i]] <- data.frame(result.trap, result.trafic.ext, 
                result.pedgene.vc, result.pedgene.burden,
                result.fbskat.vc, result.fbskat.burden,
                result.cc) 
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

###########################################
##initialize
opts['seed']=1000 #starting seed number
opts['n_rep_total']=1000 #total number of replications
opts['n_rep']=50 #number of replications in each parallele job
opts['n_ite']=opts['n_rep_total']/opts['n_rep'] #number of parallele jobs
opts['n_family']=1000 #number of family
opts['p_dis']=0.01 #prevalence

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
#family_strct=2g.3a.1u
opts['family_strct'] = '\"2g.2a.2u\"' #family structure
opts['f'] = 0.01 #rare
for i in numpy.arange(1,4.1,0.1):
  for j in range(opts['n_ite']):
    opts['r'] = i
    tgt = 'callmainSim_{seed}.OK'.format(**opts)
    inputFilesOK.append(tgt)
    dep = ''
    cmd = ['R --vanilla --args seed {seed} n_rep {n_rep} r {r} f {f} n_family {n_family} p_dis {p_dis} family_strct {family_strct} < mainSim.R > mainSim{seed}_{n_rep}_{r}_{n_family}_{p_dis}_{f}_{family_strct}.Rout 2>&1'.format(**opts)]
    makeJob(opts['launchMethod'], tgt, dep, cmd)
    opts['seed'] += 1	

opts['f'] = 0.05 #less rare
for i in numpy.arange(1,3.0,0.1):
  for j in range(opts['n_ite']):
    opts['r'] = i
    tgt = 'callmainSim_{seed}.OK'.format(**opts)
    inputFilesOK.append(tgt)
    dep = ''
    cmd = ['R --vanilla --args seed {seed} n_rep {n_rep} r {r} f {f} n_family {n_family} p_dis {p_dis} family_strct {family_strct} < mainSim.R > mainSim{seed}_{n_rep}_{r}_{n_family}_{p_dis}_{f}_{family_strct}.Rout 2>&1'.format(**opts)]
    makeJob(opts['launchMethod'], tgt, dep, cmd)
    opts['seed'] += 1	

opts['f'] = 0.20 #common
for i in numpy.arange(1,2.0,0.1):
  for j in range(opts['n_ite']):
    opts['r'] = i
    tgt = 'callmainSim_{seed}.OK'.format(**opts)
    inputFilesOK.append(tgt)
    dep = ''
    cmd = ['R --no-save --no-restore --args seed {seed} n_rep {n_rep} r {r} f {f} n_family {n_family} p_dis {p_dis} family_strct {family_strct} < mainSim.R > mainSim{seed}_{n_rep}_{r}_{n_family}_{p_dis}_{f}_{family_strct}.Rout 2>&1'.format(**opts)]
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
result <- read.csv("2016-02-19 figure 1. correct trap and trafic_ext compare with pedgene, cc, and fbat_ext 2g2c_at p_dis0.01.csv", header=T)
result <- result %>% gather(key="method", value="p.value", 9:15)
result.plot <- result %>% group_by(f, risk.variant.id, risk.haplo.f, r, method) %>% 
  summarise(n=n(), power=mean(p.value<2.5*10^-6, na.rm=T))
result.plot$risk.haplo.f <- factor(result.plot$risk.haplo.f, labels=c("f=0.01","f=0.05","f=0.20"))
result.plot$method <- factor(result.plot$method, levels=c("result.trap", "result.trafic.ext", "result.pedgene.vc", 
																													"result.pedgene.burden","result.fbskat.vc","result.fbskat.burden","result.cc"), 
														 labels = c("TRAP", "TRAFIC_EXT", "Pedgene.vc", "Pedgene", "FB-SKAT.vc", "FB-SKAT","CC"))

#figure 1
pd <- position_dodge(0.0)
filter(result.plot, !grepl("vc", method)) %>% ggplot(aes(x=r, y=power, ymax=max(power), group=method, col=method)) +
  #   geom_point(size=3, alpha=1) +
  facet_wrap(~risk.haplo.f, ncol=3, scale="free_x") +
  geom_line(size=1.2, alpha=0.7, position=pd) +
  geom_point(size=1.2, position=pd) +
#   ggtitle("f=0.202, consider effect size of risk haplotypes, TRAP") +
#   ggtitle("f=0.0178, consider effect size of risk haplotypes, TRAP") +
#   ggtitle("f=0.0039, consider effect size of risk haplotypes, TRAP") +
  labs(x="odds ratio r") +
  scale_y_continuous(limits=c(0,1)) +
  theme_gray(base_size = 20) +
  theme(legend.position="bottom", panel.background = element_rect(fill = 'grey85')) +
  scale_color_manual(values=cbbPalette)
