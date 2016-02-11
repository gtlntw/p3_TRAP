## ... your simualtion codek
#generate haplotype risk
haplotype.risk <- generate_haplorisk(r=r, het=T)

##gene drop simulation for two generations
if(family_strct=="2g.2a.2u") {
  family_strct <- data.frame(family=c(1,1,1,1), person=c(1,2,3,4), father=c(0,0,1,1), 
                             mother=c(0,0,2,2), sex=c(1,2,1,1), affect=c(1,1,2,2)) #1=male, 2=female, 1=unaffected, 2=affected
}

if(family_strct=="2g.3a.1u") {
  family_strct <- data.frame(family=c(1,1,1,1), person=c(1,2,3,4), father=c(0,0,1,1), 
                             mother=c(0,0,2,2), sex=c(1,2,1,1), affect=c(1,2,2,2)) #1=male, 2=female, 1=unaffected, 2=affected
}

if(family_strct=="2g.3a.2u") {
  family_strct <- data.frame(family=c(1,1,1,1,1), person=c(1,2,3,4,5), father=c(0,0,1,1,1), 
                             mother=c(0,0,2,2,2), sex=c(1,2,1,1,1), affect=c(1,1,2,2,2)) #1=male, 2=female, 1=unaffected, 2=affected
}

if(family_strct=="3g.3a.4u") {
  family_strct <- data.frame(family=c(1,1,1,1,1,1,1), person=c(1,2,3,4,5,6,7), father=c(0,0,1,1,0,4,4), 
                             mother=c(0,0,2,2,0,5,5), sex=c(1,2,1,1,2,1,1), affect=c(1,2,1,2,1,2,1)) #1=male, 2=female, 1=unaffected, 2=affected
}

n_cc <- length(family_strct$person)/2*n_family

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
  
  family_generated <<- gene_family(family_strct=family_strct, n=n_family, haplotype.risk=haplotype.risk) 
  
  #add pseudo null variant since FB-SKAT can only work with 2+ in a gene
  family_generated_diploid <- hap2dip(data=family_generated, risk.variant.id=risk.variant, save.file=T)
  
  ##test based on new imputation framework
  result.trap <- family.test(data=family_generated, f=risk.variant)$p.value
  result.trafic_ext <- family.test.trafic.ext(data=family_generated, f=risk.variant)$p.value
  result.pedgene <- pedgene(ped=family_generated_diploid$ped, geno=family_generated_diploid$geno)
  result.pedgene.vc <- result.pedgene$pgdf$pval.kernel
  result.pedgene.burden <- result.pedgene$pgdf$pval.burden
  
  system(paste("/net/frodo/home/khlin/p3/FB-SKAT/FB-SKAT_v2.1 data_",seed,".ped variant_pass.txt genes.txt weight.txt results_", seed, "_ mendelian_errors_",seed,"_ 10000 1 0.01 0 0", sep=""))
  system(paste("/net/frodo/home/khlin/p3/FB-SKAT/FB-SKAT_v2.1 data_",seed,".ped variant_pass.txt genes.txt weight.txt results_", seed, "_ mendelian_errors_",seed,"_ 10000 1 0.01 1 0", sep=""), ignore.stdout = TRUE)
  x <- read.table(paste("results_", seed,"_1.000000_0.010000_0.txt", sep="")) #variance component
  result.fbskat.vc <- 1-pchisq(x[,4],x[,5])
  x <- read.table(paste("results_", seed,"_1.000000_0.010000_1.txt", sep="")) #burden
  result.fbskat.burden=1-pchisq(x[,4],x[,5])
  
  #case-control
  data.cc <- gene_case_control(n_case_control_pair=n_cc)
  result.cc <- case_control.test(data=data.cc, f=risk.variant)[3]
  
  #only report p.value
  c(result.trap, result.trafic_ext, 
    result.pedgene.vc, result.pedgene.burden,
    result.fbskat.vc, result.fbskat.burden, 
    result.cc)
})
## Write out your results to a csv file
result.df <- as.data.frame(t(sim_result))
colnames(result.df) <- c("result.trap", "result.trafic_ext", 
                         "result.pedgene.vc", "result.pedgene.burden",
                         "result.fbskat.vc", "result.fbskat.burden",
                         "result.cc")
result.df <- cbind(seed,r,p_dis,risk.variant=paste(c(risk.variant), collapse = "_"),risk.haplo.f,n_family,result.df)
write.csv(result.df, paste("res_",r,"_",n_family,"_",seed,".csv",sep=""), row.names=FALSE)
## R.miSniam
## End(Not run)
## Not run:

###########################################
##initialize
opts['seed']=1000 #starting seed number
opts['n_rep_total']=1000 #total number of replications
opts['n_rep']=100 #number of replications in each parallele job
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

######################
#1.1. run simulations by calling mainSim.R
######################
inputFilesOK = []
#family_strct=2g.3a.1u
opts['riskVariant'] = '\"c(18,25,47)\"' #super rare
opts['family_strct'] = '\"2g.3a.1u\"' #family structure
for i in numpy.linspace(1,2.8,num=15):
	for j in range(opts['n_ite']):
		opts['r'] = i
		tgt = 'callmainSim_{seed}.OK'.format(**opts)
		inputFilesOK.append(tgt)
		dep = ''
		cmd = ['R --vanilla --args seed {seed} n_rep {n_rep} r {r} n_family {n_family} p_dis {p_dis} risk.variant {riskVariant} family_strct {family_strct}< mainSim.R > mainSim_{n_rep}_{r}_{n_family}_{p_dis}_{riskVariant}_{family_strct}.Rout{seed} 2>&1'.format(**opts)]
		makeJob(opts['launchMethod'], tgt, dep, cmd)
		opts['seed'] += 1	

opts['riskVariant'] = '\"c(4,15,32)\"' #rare
for i in numpy.linspace(1,2.3,num=15):
	for j in range(opts['n_ite']):
		opts['r'] = i
		tgt = 'callmainSim_{seed}.OK'.format(**opts)
		inputFilesOK.append(tgt)
		dep = ''
		cmd = ['R --vanilla --args seed {seed} n_rep {n_rep} r {r} n_family {n_family} p_dis {p_dis} risk.variant {riskVariant} family_strct {family_strct}< mainSim.R > mainSim_{n_rep}_{r}_{n_family}_{p_dis}_{riskVariant}_{family_strct}.Rout{seed} 2>&1'.format(**opts)]
		makeJob(opts['launchMethod'], tgt, dep, cmd)
		opts['seed'] += 1	

opts['riskVariant'] = '\"c(4,16,42)\"' #common
for i in numpy.linspace(1,1.8,num=15):
	for j in range(opts['n_ite']):
		opts['r'] = i
		tgt = 'callmainSim_{seed}.OK'.format(**opts)
		inputFilesOK.append(tgt)
		dep = ''
		cmd = ['R --vanilla --args seed {seed} n_rep {n_rep} r {r} n_family {n_family} p_dis {p_dis} risk.variant {riskVariant} family_strct {family_strct}< mainSim.R > mainSim_{n_rep}_{r}_{n_family}_{p_dis}_{riskVariant}_{family_strct}.Rout{seed} 2>&1'.format(**opts)]
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
result <- read.csv("figure 1. correct trap and trafic_ext compare with pedgene, cc, and fbat_ext 2g2c_at p_dis0.1.csv", header=T)
result <- result %>% gather(key="method", value="p.value", 7:13)
result.plot <- result %>% group_by(risk.variant, risk.haplo.f, r, method) %>% 
  summarise(n=n(), power=mean(p.value<2.5*10^-6, na.rm=T))
result.plot$risk.haplo.f <- factor(result.plot$risk.haplo.f, labels=c("f=0.0062","f=0.0116","f=0.2135"))
result.plot$method <- factor(result.plot$method, labels=c("TRAP","TRAFIC_EXT","Pedgene.vc", "Pedgene", "FB_SKAT.vc", "FB_SKAT", "Case-Control"))

#only super rare and common
pd <- position_dodge(0.0)
filter(result.plot, !grepl("vc", method), grepl("0.0062|0.2135", risk.haplo.f)) %>% ggplot(aes(x=r, y=power, ymax=max(power), group=method, col=method)) +
  #   geom_point(size=3, alpha=1) +
  facet_wrap(~risk.haplo.f, ncol=1, scale="free_x") +
  geom_line(size=3, alpha=0.7, position=pd) +
  # geom_point(size=2, position=pd) +
#   ggtitle("f=0.202, consider effect size of risk haplotypes, TRAP") +
#   ggtitle("f=0.0178, consider effect size of risk haplotypes, TRAP") +
#   ggtitle("f=0.0039, consider effect size of risk haplotypes, TRAP") +
  labs(x="relative risk") +
  scale_y_continuous(limits=c(0,1)) +
  theme_bw(base_size = 25) +
  theme(legend.position="bottom") +
  scale_color_manual(values=cbbPalette)

ggsave(file="figure 1.png", height=10, width=10)
