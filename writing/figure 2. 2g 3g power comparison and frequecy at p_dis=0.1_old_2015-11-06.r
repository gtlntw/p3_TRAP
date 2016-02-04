##load project related functions
source('function.r')

#need to take in parameter
## passed in from the command prompt.
parseCommandArgs()
## Will disply the default values on the first run,
print(seed)
print(n_rep)
print(r)
print(n_family) #number of family
print(p_dis) #prevelance
print(risk.variant)
print(family_strct)

## ... your simualtion codek
#generate haplotype risk
haplotype.risk <- generate_haplorisk(r=r, het=T)

##gene drop simulation for two generations
family_strct_ped <- select_family_strct(family_strct)

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

  #simulation for two-generation families
  family_generated <- gene_family(family_strct=family_strct_ped, n=n_family, haplotype.risk=haplotype.risk) 
  #add pseudo null variant since FB-SKAT can only work with 2+ in a gene
  family_generated_diploid <- hap2dip(data=family_generated, risk.variant.id=risk.variant, save.file=T)

  
  #run tests
  result.trap <- family.test(data=family_generated, f=risk.variant.id)$p.value
  result.trafic.ext <- family.test.trafic.ext(data=family_generated, f=risk.variant.id)$p.value
  result.pedgene <- pedgene(ped=family_generated_diploid$ped, geno=family_generated_diploid$geno)
  result.pedgene.vc <- result.pedgene$pgdf$pval.kernel
  result.pedgene.burden <- result.pedgene$pgdf$pval.burden
  
  system(paste("/net/frodo/home/khlin/p3/FB-SKAT/FB-SKAT_v2.1 data_",seed,".ped variant_pass.txt genes.txt weight.txt results_", seed, "_ mendelian_errors_",seed,"_ 10000 1 0.01 0 0", sep=""))
  system(paste("/net/frodo/home/khlin/p3/FB-SKAT/FB-SKAT_v2.1 data_",seed,".ped variant_pass.txt genes.txt weight.txt results_", seed, "_ mendelian_errors_",seed,"_ 10000 1 0.01 1 0", sep=""), ignore.stdout = TRUE)
  x <- read.table(paste("results_", seed,"_1.000000_0.010000_0.txt", sep="")) #variance component
  result.fbskat.vc <- 1-pchisq(x[,4],x[,5])
  x <- read.table(paste("results_", seed,"_1.000000_0.010000_1.txt", sep="")) #burden
  result.fbskat.burden=1-pchisq(x[,4],x[,5])

  #only report p.value
  c(result.trap, result.trafic.ext, 
    result.pedgene.vc, result.pedgene.burden,
    result.fbskat.vc, result.fbskat.burden)
})
## Write out your results to a csv file
result.df <- as.data.frame(t(sim_result))
colnames(result.df) <- c("result.trap", "result.trafic.ext", 
											    "result.pedgene.vc", "result.pedgene.burden",
											    "result.fbskat.vc", "result.fbskat.burden")
result.df <- cbind(seed,r,p_dis,risk.variant=paste(c(risk.variant), collapse = "_"),risk.haplo.f,n_family,family_strct,result.df)
write.csv(result.df, paste("res_",r,"_",n_family,"_",seed,".csv",sep=""), row.names=FALSE)
## R.miSniam
## End(Not run)
## Not run:

###############################################################################
####################commands to run jobs in parallel
##initialize
opts['seed']=1000 #starting seed number
opts['n_rep_total']=1000 #total number of replications
opts['n_rep']=100 #number of replications in each parallele job
opts['n_ite']=opts['n_rep_total']/opts['n_rep'] #number of parallele jobs
opts['n_family']=1000 #number of family
opts['p_dis']=0.1 #prevalence

#2g.3a.1u 3g.3a.4u 3g.2a.5u 3g.3a.5u
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
family_strct = ['2g.3a.1u','3g.3a.4u','3g.2a.5u','3g.3a.5u']
# opts['riskVariant'] = '\"c(18,25,47)\"' #super rare
# opts['family_strct'] = '\"2g.4a.1u\"' #family structure
# for i in numpy.linspace(1,2.0,num=15):
# 	for j in range(opts['n_ite']):
# 		opts['r'] = i
# 		tgt = 'callmainSim_{seed}.OK'.format(**opts)
# 		inputFilesOK.append(tgt)
# 		dep = ''
# 		cmd = ['R --vanilla --args seed {seed} n_rep {n_rep} r {r} n_family {n_family} p_dis {p_dis} risk.variant {riskVariant} family_strct {family_strct}< mainSim.R > mainSim_{n_rep}_{r}_{n_family}_{p_dis}_{riskVariant}_{family_strct}.Rout{seed} 2>&1'.format(**opts)]
# 		makeJob(opts['launchMethod'], tgt, dep, cmd)
# 		opts['seed'] += 1	

opts['riskVariant'] = '\"c(4,15,32)\"' #rare
for i in numpy.linspace(1,2.4,num=15):
	for j in range(opts['n_ite']):
		for k in family_strct:
			opts['family_strct'] = k
			opts['r'] = i
			tgt = 'callmainSim_{seed}.OK'.format(**opts)
			inputFilesOK.append(tgt)
			dep = ''
			cmd = ['--mem-per-cpu=4096  --time=2-0 R --vanilla --args seed {seed} n_rep {n_rep} r {r} n_family {n_family} p_dis {p_dis} risk.variant {riskVariant} family_strct {family_strct}< mainSim.R > mainSim_{n_rep}_{r}_{n_family}_{p_dis}_{riskVariant}_{family_strct}.Rout{seed} 2>&1'.format(**opts)]
			makeJob(opts['launchMethod'], tgt, dep, cmd)
			opts['seed'] += 1	

# opts['riskVariant'] = '\"c(4,16,42)\"' #common
# for i in numpy.linspace(1,1.4,num=15):
# 	for j in range(opts['n_ite']):
# 		opts['r'] = i
# 		tgt = 'callmainSim_{seed}.OK'.format(**opts)
# 		inputFilesOK.append(tgt)
# 		dep = ''
# 		cmd = ['R --vanilla --args seed {seed} n_rep {n_rep} r {r} n_family {n_family} p_dis {p_dis} risk.variant {riskVariant} family_strct {family_strct}< mainSim.R > mainSim_{n_rep}_{r}_{n_family}_{p_dis}_{riskVariant}_{family_strct}.Rout{seed} 2>&1'.format(**opts)]
# 		makeJob(opts['launchMethod'], tgt, dep, cmd)
# 		opts['seed'] += 1	

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


###############################################################
##power comparison between 
result <- read.csv("figure 2. 2g 3g power comparison and frequecy at p_dis=0.1_poster.csv", header=T)
result <- result %>% gather(key="method", value="p.value", 8:13)
result$family_strct <- factor(result$family_strct, levels=c("2g.3a.1u", "3g.3a.4u",
                                                            "3g.2a.5u", "3g.3a.5u"))
result$method <- factor(result$method, labels = c("TRAP", "TRAFIC_EXT", "Pedgene.vc", "Pedgene.burden", "FB-SKAT.vc", "FB-SKAT.burden"))
# result <- melt(id.vars=c("seed", "r", "p_dis","rish.haplo.f","n_family", ""), variable.name="", value.name="p.value")
result.plot <- result %>% group_by(r, family_strct, method) %>% 
  summarise(n=n(), power=mean(p.value < 2.5*10^-6))

#only TRAP and TRAFIC_EXT
filter(result.plot, grepl("TRAP|TRAFIC_EXT", method)) %>% ggplot(aes(x=r, y=power, col=family_strct)) +
  facet_wrap(~method) +
  #   geom_point(size=3, alpha=1) +
  geom_line(size=1.2, alpha=0.7) +
  ggtitle("f=0.0116") +
  labs(x="relative risk r") +
  theme_gray(base_size = 20) +
  theme(legend.position="bottom", panel.background = element_rect(fill="grey85"))

#all 6
ggplot(result.plot, aes(x=r, y=power, col=family_strct)) +
  facet_wrap(~method) +
  #   geom_point(size=3, alpha=1) +
  geom_line(size=1.2, alpha=0.7) +
  ggtitle("f=0.0116") +
  labs(x="relative risk r") +
  theme_gray(base_size = 20) +
  theme(legend.position="bottom", panel.background = element_rect(fill="grey85"))
