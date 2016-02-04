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
family_strct.2g3c <- data.frame(family=c(1,1,1,1,1), person=c(1,2,3,4,5), father=c(0,0,1,1,1), 
                                mother=c(0,0,2,2,2), sex=c(1,2,1,1,1), affect=c(1,2,2,2,2)) #1=male, 2=female, 1=unaffected, 2=affected


rep.idx <<- 1
sim_result <- replicate(n_rep, {
  print(rep.idx)
  rep.idx <<- rep.idx + 1

  family_generated_2g3c <<- gene_family(family_strct=family_strct.2g3c, n=n_family, haplotype.risk=haplotype.risk) 
  
  #remove the founder from the data
  family_generated_3c <- family_generated_2g3c
  family_generated_founder <- list()
  temp.idx <- which(family_generated_2g3c$data_family$person %in% 1:2) #take out founders
  family_generated_founder$data_family <- family_generated_2g3c$data_family[temp.idx, ]
  family_generated_founder$tran_vec <- family_generated_2g3c$tran_vec[temp.idx, ]
  family_generated_3c$data_family <- family_generated_2g3c$data_family[-temp.idx, ]
  family_generated_3c$tran_vec <- family_generated_2g3c$tran_vec[-temp.idx, ]

  data <- family_generated_3c
  f <- risk.variant
  n_family <- max(data$data_family$family)
  n_family_member <- table(data$data_family$family)
  #check if founder's haplotype carries any variant's with f < 0.1
  if(length(f)==1 & f[1] <1) { #support allele frequency or list of snps
    snp2look.idx <-  which(snp$FREQ1 < f) # snp to look for
  } else(snp2look.idx <-  f)
  
  #disguitish those siblings with four founder's haplotype
  data.founder <- family_generated_2g3c
  n_family_member.w.founder <- n_family_member + 2
  #calculate the allele frequency
  founder.list <- list() #to store if the haplotype carrying a risk variant on multiple affected
  carrier.count.2 <- 0 
  chromosome.count.2 <- 0
  carrier.count.2.mis <- 0 
  chromosome.count.2.mis <- 0
  carrier.count.2.obs <- 0 
  chromosome.count.2.obs <- 0
  carrier.count.3 <- 0 
  chromosome.count.3 <- 0
  carrier.count.3.mis <- 0 
  chromosome.count.3.mis <- 0
  carrier.count.3.obs <- 0 
  chromosome.count.3.obs <- 0
  carrier.count.4 <- 0 
  chromosome.count.4 <- 0
  haplo.unique.list <- list()
  for(family.idx in 1:n_family) {
    current_row=sum(n_family_member.w.founder[1:family.idx]) - n_family_member.w.founder[family.idx]
    #assign a suitable transmission, not trivial to code up the algorithm assume is known for now
    tran_vec <- data.founder$tran_vec[(current_row+1):(current_row+n_family_member.w.founder[family.idx]),]
    ##tally the (un)ambiguous carrier haplotype
    h1 <- data.founder$data_family[(current_row+1):(current_row+n_family_member.w.founder[family.idx]),7:(6+n_snp)] #the first haplotype
    h2 <- data.founder$data_family[(current_row+1):(current_row+n_family_member.w.founder[family.idx]),-c(1:(6+n_snp))] #the second haplotype
    #observed allele count for each individual
    carrier.founder <- c(0,0,0,0)
    carrier.founder.mis <- c(0,0,0,0)
    carrier.founder.obs <- c(0,0,0,0)
    carrier.offspring <- c(0,0,0,0)
    for(a in 3:n_family_member.w.founder[family.idx]) { #offsprings
      idx.h1 <- match(tran_vec[a,"h1"], c("A", "B", "C", "D"))
      if(h1[a, snp2look.idx]==1) carrier.offspring[idx.h1] <- carrier.offspring[idx.h1] + 1
      idx.h2 <- match(tran_vec[a,"h2"], c("A", "B", "C", "D"))
      if(h2[a, snp2look.idx]==1) carrier.offspring[idx.h2] <- carrier.offspring[idx.h2] + 1
    }
    #the number of each unique haplotype occured on affected offsprings
    haplo.unique <- unique(as.vector(as.matrix(tran_vec[3:n_family_member.w.founder[family.idx] ,c("h1","h2")])))
    haplo.unique.count <- length(haplo.unique)
    haplo.unique.list[[family.idx]] <- length(haplo.unique)
    haplo.mis <- which(is.na(match(c("A", "B", "C", "D"), haplo.unique)))
    haplo.obs <- which(!is.na(match(c("A", "B", "C", "D"), haplo.unique)))
    
    for(a in 1:2) { #founder's 
      idx.h1 <- match(tran_vec[a,"h1"], c("A", "B", "C", "D"))
      if(h1[a, snp2look.idx]==1) carrier.founder[idx.h1] <- carrier.founder[idx.h1] + 1
      idx.h2 <- match(tran_vec[a,"h2"], c("A", "B", "C", "D"))
      if(h2[a, snp2look.idx]==1) carrier.founder[idx.h2] <- carrier.founder[idx.h2] + 1
      
      if(idx.h1 %in% haplo.mis & h1[a, snp2look.idx]==1) carrier.founder.mis[idx.h1] <- carrier.founder.mis[idx.h1] + 1
      if(idx.h2 %in% haplo.mis & h2[a, snp2look.idx]==1) carrier.founder.mis[idx.h2] <- carrier.founder.mis[idx.h2] + 1
      if(idx.h1 %in% haplo.obs & h1[a, snp2look.idx]==1) carrier.founder.obs[idx.h1] <- carrier.founder.obs[idx.h1] + 1
      if(idx.h2 %in% haplo.obs & h2[a, snp2look.idx]==1) carrier.founder.obs[idx.h2] <- carrier.founder.obs[idx.h2] + 1
    }
    
    #observed founder carrier based on offsprings
    founder.list[[family.idx]] <- c(observed=sum(carrier.offspring), carrier=sum(carrier.offspring>0), haplo.unique=haplo.unique)
    #steps to calculate true founder allele ferquency by the number of observed founder chromosome
    if(haplo.unique.count==2) {
      carrier.count.2 <- carrier.count.2 + sum(carrier.founder>=1)
      chromosome.count.2 <- chromosome.count.2 + 4
      carrier.count.2.mis <- carrier.count.2.mis + sum(carrier.founder.mis>=1)
      chromosome.count.2.mis <- chromosome.count.2.mis + 2
      carrier.count.2.obs <- carrier.count.2.obs + sum(carrier.founder.obs>=1)
      chromosome.count.2.obs <- chromosome.count.2.obs + 2
    }
    if(haplo.unique.count==3) {
      carrier.count.3 <- carrier.count.3 + sum(carrier.founder>=1)
      chromosome.count.3 <- chromosome.count.3 + 4
      carrier.count.3.mis <- carrier.count.3.mis + sum(carrier.founder.mis>=1)
      chromosome.count.3.mis <- chromosome.count.3.mis + 1
      carrier.count.3.obs <- carrier.count.3.obs + sum(carrier.founder.obs>=1)
      chromosome.count.3.obs <- chromosome.count.3.obs + 3      }
    if(haplo.unique.count==4) {
      carrier.count.4 <- carrier.count.4 + sum(carrier.founder>=1)
      chromosome.count.4 <- chromosome.count.4 + 4
    }
  }
  #the allele frequency
  (founder.freq.offspring.2 <- carrier.count.2/chromosome.count.2)
  (founder.freq.offspring.2.mis <- carrier.count.2.mis/chromosome.count.2.mis)
  (founder.freq.offspring.2.obs <- carrier.count.2.obs/chromosome.count.2.obs)
  (founder.freq.offspring.3 <- carrier.count.3/chromosome.count.3)
  (founder.freq.offspring.3.mis <- carrier.count.3.mis/chromosome.count.3.mis)
  (founder.freq.offspring.3.obs <- carrier.count.3.obs/chromosome.count.3.obs)
  (founder.freq.offspring.4 <- carrier.count.4/chromosome.count.4)

    #only report p.value
  c(founder.freq.offspring.2, founder.freq.offspring.2.mis, founder.freq.offspring.2.obs,
    founder.freq.offspring.3, founder.freq.offspring.3.mis, founder.freq.offspring.3.obs,
    founder.freq.offspring.4)
})
## Write out your results to a csv file
result.df <- as.data.frame(t(sim_result))
colnames(result.df) <- c("founder.freq.offspring.2", "founder.freq.offspring.2.mis", "founder.freq.offspring.2.obs",
                         "founder.freq.offspring.3", "founder.freq.offspring.3.mis", "founder.freq.offspring.3.obs",
                         "founder.freq.offspring.4")
result.df <- cbind(seed,r,p_dis,risk.variant,risk.haplo.f,n_family,result.df)
write.csv(result.df, paste("res_",r,"_",n_family,"_",seed,".csv",sep=""), row.names=FALSE)
## R.miSniam
## End(Not run)
## Not run:

##commands to run jobs in parallel
##passed in from the command prompt.
parseCommandArgs()
## Will disply the default values on the first run,
##initialize
seed=1000
#number of replications
n_rep=34
#number of family
n_family=1000
#prevalence
p_dis=0.3
# print(null) #null=FALSE
#family structure
# print(family_strct) #family_strct='family_strct.2g3c'

print(getwd())

##command line
parallel <- function(...) {
  names.args <- names(list(...))
  value.args <- as.vector(list(...))
  ##commands to run
  for(i in 1:3) {
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
for(i in seq(1,2, length.out = 5)) {
  parallel(n_rep=n_rep, r=i, n_family=n_family, p_dis=p_dis, risk.variant=2) #common
}
for(i in seq(1,2, length.out = 5)) {
  parallel(n_rep=n_rep, r=i, n_family=n_family, p_dis=p_dis, risk.variant=7) #rare
}
for(i in seq(1,2, length.out = 5)) {
  parallel(n_rep=n_rep, r=i, n_family=n_family, p_dis=p_dis, risk.variant=39) #super rare
}


###################################################
# The color-blind friednly palette
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
# The palette with black:
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
#check power using true, imputed founder carrier, minimum offspring carrier
#three versions -- command and rare and super rare #2 for common, 7 for rare, 39 for super rare
result <- read.csv("2015-06-30 imputation framework for TRAP to fix common variant see allele freqency on mis and obs.csv", header=T)
result <- result %>% melt(id.vars=c("seed", "r", "p_dis","risk.variant","risk.haplo.f","n_family"), variable.name="method", value.name="p.value")
result <- mutate(result, transmitted=ifelse(grepl("mis", method), "F", ifelse(grepl("obs", method), "T", "NA")))
result <- mutate(result, n_transmitted=ifelse(grepl("2", method), "2",ifelse(grepl("3", method), "3","4")))
result.plot <- result %>% group_by(risk.variant, risk.haplo.f, r, method, transmitted, n_transmitted) %>% 
  summarise(n=n(), maf=mean(p.value))

#full frequency by family who inherit different number of chromosomes
pd <- position_dodge(0.0)
filter(result.plot, risk.variant %in% c(2)) %>% 
  ggplot(aes(x=r, y=maf, ymax=max(maf), group=method, col=n_transmitted, lty=transmitted)) +
  facet_grid(.~risk.haplo.f,labeller = label_bquote(expr=f == .(x)), scale="free_y") + 
  #   geom_point(size=3, alpha=1) +
  geom_line(size=1.2, alpha=0.7, position=pd) +
  geom_point(size=1.2, position=pd) +
  #   ggtitle("f=0.202, consider effect size of risk haplotypes, TRAP") +
  #   ggtitle("f=0.0178, consider effect size of risk haplotypes, TRAP") +
  #   ggtitle("f=0.0039, consider effect size of risk haplotypes, TRAP") +
  labs(x="relative risk r", y="MAF") +
  coord_cartesian(ylim=c(0, 0.5))+
  theme_gray(base_size = 20) +
  theme(legend.position="right", panel.background = element_rect(fill = 'grey85')) +
  scale_color_manual(values=cbbPalette)


#only full frequency by family who inherit different number of chromosomes
pd <- position_dodge(0.0)
filter(result.plot, risk.variant %in% c(7), grepl("[0-9]$", method)) %>% 
  ggplot(aes(x=r, y=maf, ymax=max(maf), group=method, col=n_transmitted, lty=transmitted)) +
  facet_grid(.~risk.haplo.f,labeller = label_bquote(expr=f == .(x)), scale="free_y") + 
  #   geom_point(size=3, alpha=1) +
  geom_line(size=1.2, alpha=0.7, position=pd) +
  geom_point(size=1.2, position=pd) +
  #   ggtitle("f=0.202, consider effect size of risk haplotypes, TRAP") +
  #   ggtitle("f=0.0178, consider effect size of risk haplotypes, TRAP") +
  #   ggtitle("f=0.0039, consider effect size of risk haplotypes, TRAP") +
  labs(x="relative risk r", y="MAF") +
  coord_cartesian(ylim=c(0, 0.2))+
  theme_gray(base_size = 20)+
  theme(legend.position="right", panel.background = element_rect(fill = 'grey85')) +
  scale_color_manual(values=cbbPalette)
