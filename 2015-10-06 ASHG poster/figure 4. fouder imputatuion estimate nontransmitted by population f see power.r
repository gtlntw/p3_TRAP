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

  family_generated <<- gene_family(family_strct=family_strct_ped, n=n_family, haplotype.risk=haplotype.risk) 
  
  #remove the founder from the data
  family_generated_c <- family_generated
  family_generated_founder <- list()
  temp.idx <- which(family_generated$data_family$person %in% 1:2) #take out founders
  family_generated_founder$data_family <- family_generated$data_family[temp.idx, ]
  family_generated_founder$tran_vec <- family_generated$tran_vec[temp.idx, ]
  family_generated_c$data_family <- family_generated$data_family[-temp.idx, ]
  family_generated_c$tran_vec <- family_generated$tran_vec[-temp.idx, ]
  
  ##test based on new imputation framework
  result.trap.3c.noimpute <- family.test(data=family_generated, f=risk.variant.id, nofounderphenotype=T)$p.value
  result.trap.3c.impute.founder.pop.f <- family.test.trap.impute.founder(data=family_generated_c, f=risk.variant.id, pop.f=T, pop.f.off=0)
  result.trap.3c.impute.founder.pop.f.off.5 <- family.test.trap.impute.founder(data=family_generated_c, f=risk.variant.id, pop.f=T, pop.f.off=5)
  result.trap.3c.impute.founder.pop.f.off.10 <- family.test.trap.impute.founder(data=family_generated_c, f=risk.variant.id, pop.f=T, pop.f.off=10)
  result.trap.3c.impute.founder.pop.f.off.20 <- family.test.trap.impute.founder(data=family_generated_c, f=risk.variant.id, pop.f=T, pop.f.off=20)
  result.trap.3c.impute.founder.pop.f.off.50 <- family.test.trap.impute.founder(data=family_generated_c, f=risk.variant.id, pop.f=T, pop.f.off=50)
  result.trap.3c.impute.founder.pop.f.off.n5 <- family.test.trap.impute.founder(data=family_generated_c, f=risk.variant.id, pop.f=T, pop.f.off=-5)
  result.trap.3c.impute.founder.pop.f.off.n10 <- family.test.trap.impute.founder(data=family_generated_c, f=risk.variant.id, pop.f=T, pop.f.off=-10)
  result.trap.3c.impute.founder.pop.f.off.n20 <- family.test.trap.impute.founder(data=family_generated_c, f=risk.variant.id, pop.f=T, pop.f.off=-20)
  result.trap.3c.impute.founder.pop.f.off.n50 <- family.test.trap.impute.founder(data=family_generated_c, f=risk.variant.id, pop.f=T, pop.f.off=-50)
  
  #only report p.value
  c(result.trap.3c.noimpute, result.trap.3c.impute.founder.pop.f,
    result.trap.3c.impute.founder.pop.f.off.5, result.trap.3c.impute.founder.pop.f.off.10,
    result.trap.3c.impute.founder.pop.f.off.20, result.trap.3c.impute.founder.pop.f.off.50,
    result.trap.3c.impute.founder.pop.f.off.n5, result.trap.3c.impute.founder.pop.f.off.n10,
    result.trap.3c.impute.founder.pop.f.off.n20, result.trap.3c.impute.founder.pop.f.off.n50)
})
## Write out your results to a csv file
result.df <- as.data.frame(t(sim_result))
colnames(result.df) <- c("result.trap.3c.noimpute", "result.trap.3c.impute.founder.pop.f",
                         "result.trap.3c.impute.founder.pop.f.off.5", "result.trap.3c.impute.founder.pop.f.off.10",
                         "result.trap.3c.impute.founder.pop.f.off.20", "result.trap.3c.impute.founder.pop.f.off.50",
											   "result.trap.3c.impute.founder.pop.f.off.n5", "result.trap.3c.impute.founder.pop.f.off.n10",
											   "result.trap.3c.impute.founder.pop.f.off.n20", "result.trap.3c.impute.founder.pop.f.off.n50")
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
opts['p_dis']=0.3 #prevalence

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
opts['riskVariant'] = '\"c(18,25,47)\"' #super rare
opts['family_strct'] = '\"2g.4a.1u\"' #family structure
for i in numpy.linspace(1,2.0,num=15):
	for j in range(opts['n_ite']):
		opts['r'] = i
		tgt = 'callmainSim_{seed}.OK'.format(**opts)
		inputFilesOK.append(tgt)
		dep = ''
		cmd = ['R --vanilla --args seed {seed} n_rep {n_rep} r {r} n_family {n_family} p_dis {p_dis} risk.variant {riskVariant} family_strct {family_strct}< mainSim.R > mainSim_{n_rep}_{r}_{n_family}_{p_dis}_{riskVariant}_{family_strct}.Rout{seed} 2>&1'.format(**opts)]
		makeJob(opts['launchMethod'], tgt, dep, cmd)
		opts['seed'] += 1	

opts['riskVariant'] = '\"c(4,15,32)\"' #rare
for i in numpy.linspace(1,1.6,num=15):
	for j in range(opts['n_ite']):
		opts['r'] = i
		tgt = 'callmainSim_{seed}.OK'.format(**opts)
		inputFilesOK.append(tgt)
		dep = ''
		cmd = ['R --vanilla --args seed {seed} n_rep {n_rep} r {r} n_family {n_family} p_dis {p_dis} risk.variant {riskVariant} family_strct {family_strct}< mainSim.R > mainSim_{n_rep}_{r}_{n_family}_{p_dis}_{riskVariant}_{family_strct}.Rout{seed} 2>&1'.format(**opts)]
		makeJob(opts['launchMethod'], tgt, dep, cmd)
		opts['seed'] += 1	

opts['riskVariant'] = '\"c(4,16,42)\"' #common
for i in numpy.linspace(1,1.4,num=15):
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
result <- read.csv("figure 4. fouder imputatuion estimate nontransmitted by population f see power.csv", header=T)
result <- result %>% gather(key="method", value="p.value", 8:17)
result.plot <- result %>% group_by(risk.variant, risk.haplo.f, r, method) %>% 
  summarise(n=n(), power=mean(p.value<0.05, na.rm=T))
result.plot$risk.haplo.f <- factor(result.plot$risk.haplo.f, labels=c("f=0.0062","f=0.0116","f=0.2135"))
result.plot$method <- factor(result.plot$method, labels=c("TRAP no_impute", "pop. allele freq. f", "off +5%", "off +10%", "off +20%", "off +50%", "off -5%", "off -10%", "off -20%", "off -50%"))

#exclude conditional and pop.off.10, pop.off.50 
pd <- position_dodge(0.0)
filter(result.plot, !grepl("conditional|10", method), grepl("0.0116", risk.haplo.f), !grepl("5%|50%", method)) %>% ggplot(aes(x=r, y=power, ymax=max(power), group=method, col=method)) +
  #   geom_point(size=3, alpha=1) +
  # facet_wrap(~risk.haplo.f, ncol=3, scale="free_x") +
  geom_line(size=3, alpha=0.7, position=pd) +
  #   ggtitle("f=0.202, consider effect size of risk haplotypes, TRAP") +
  #   ggtitle("f=0.0178, consider effect size of risk haplotypes, TRAP") +
  #   ggtitle("f=0.0039, consider effect size of risk haplotypes, TRAP") +
  labs(x="relative risk") +
	ggtitle("f=0.0116") + 
  coord_cartesian(xlim = c(1, 1.65), ylim = c(0, 1))+
  theme_bw(base_size = 27) +
  theme(legend.position="bottom") +
  scale_color_manual(values=cbbPalette)

ggsave(file="figure 4.png", height=7, width=10.5)
