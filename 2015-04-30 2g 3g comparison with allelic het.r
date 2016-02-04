source("function.r")
#need to take in parameter
## passed in from the command prompt.
parseCommandArgs()
## Will disply the default values on the first run,
print(seed)
print(n_rep)
print(r)
print(n_family) #number of family
print(p_dis) #prevelance

## ... your simualtion code
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
# risk.variant.id <- sort(sample(which(snp$FREQ1 < 0.02), 10))
risk.variant.id <- c(5,8,9,23,26,27,28,32,38,46)
risk.haplo.id <- which(apply(2-haplotype[, risk.variant.id+2], 1, sum)>0)
(rish.haplo.f <- mean(apply(2-haplotype[, risk.variant.id+2], 1, sum)>0)) #carrier haplotype frequency
print(risk.variant.id) #print risk variants
b <- rep(1, length=(ncol(haplotype)-2)) #initialize every haplotype is neutral
b[risk.variant.id] <- r #effect size is log10 of its allele frequency
#calculate the haplotype variants p(A|h)
haplotype.risk <- apply(2-haplotype[, -c(1:2)], 1, function(x) prod(b^x)*b0_sqrt)
#calculate the probability of drawing given the affected status of founder p(h|A)
# haplotype.draw_prob <- (haplotype.risk %o% haplotype.risk)/n_haplo^2/prevalence
(mu <<- mean(haplotype.risk[risk.haplo.id]/b0_sqrt)) #mean relative risk
(sigma2 <<- var(haplotype.risk[risk.haplo.id]/b0_sqrt)) #variance relative risk

##gene drop simulation for two generations
family_strct.2g.3a.1u <- data.frame(family=c(1,1,1,1), person=c(1,2,3,4), father=c(0,0,1,1), 
                                mother=c(0,0,2,2), sex=c(1,2,1,1), affect=c(1,2,2,2)) #1=male, 2=female, 1=unaffected, 2=affected
family_strct.3g.3a.4u <- data.frame(family=c(1,1,1,1,1,1,1), person=c(1,2,3,4,5,6,7), father=c(0,0,1,1,0,4,4), 
                                 mother=c(0,0,2,2,0,5,5), sex=c(1,2,1,1,2,1,1), affect=c(1,2,1,2,1,2,1)) #1=male, 2=female, 1=unaffected, 2=affected
family_strct.3g.2a.5u <- data.frame(family=c(1,1,1,1,1,1,1), person=c(1,2,3,4,5,6,7), father=c(0,0,1,1,0,4,4), 
                                 mother=c(0,0,2,2,0,5,5), sex=c(1,2,1,1,2,1,1), affect=c(1,2,1,1,1,2,1)) #1=male, 2=female, 1=unaffected, 2=affected
family_strct.3g.3a.5u <- data.frame(family=c(1,1,1,1,1,1,1,1), person=c(1,2,3,4,5,6,7,8), father=c(0,0,1,1,0,4,4,4), 
                                mother=c(0,0,2,2,0,5,5,5), sex=c(1,2,1,1,2,1,1,1), affect=c(1,2,1,2,1,2,1,1)) #1=male, 2=female, 1=unaffected, 2=affected

rep.idx <<- 1
sim_result <- replicate(n_rep, {
  print(rep.idx)
  rep.idx <<- rep.idx + 1
  
  #simulation for two-generation families
  family_generated_2g.3a.1u <- gene_family(family_strct=family_strct.2g.3a.1u, n=n_family, haplotype.risk=haplotype.risk) 
  result.trap.2g.3a.1u <- family.test(data=family_generated_2g.3a.1u, f=risk.variant.id)
  result.trafic.ext.2g.3a.1u <- family.test.trafic.ext(data=family_generated_2g.3a.1u, f=risk.variant.id)
    
  #simulation for three-generation families
  #three affected and four unaffected
  family_generated_3g.3a.4u <- gene_family(family_strct=family_strct.3g.3a.4u, n=n_family, haplotype.risk=haplotype.risk) 
  result.trap.3g.3a.4u <- family.test(data=family_generated_3g.3a.4u, f=risk.variant.id)
  result.trafic.ext.3g.3a.4u <- family.test.trafic.ext(data=family_generated_3g.3a.4u, f=risk.variant.id)
  #two affected and five unaffected
  family_generated_3g.2a.5u <- gene_family(family_strct=family_strct.3g.2a.5u, n=n_family, haplotype.risk=haplotype.risk) 
  result.trap.3g.2a.5u <- family.test(data=family_generated_3g.2a.5u, f=risk.variant.id)
  result.trafic.ext.3g.2a.5u <- family.test.trafic.ext(data=family_generated_3g.2a.5u, f=risk.variant.id)
  #three affected and five unaffected
  family_generated_3g.3a.5u <- gene_family(family_strct=family_strct.3g.3a.5u, n=n_family, haplotype.risk=haplotype.risk) 
  result.trap.3g.3a.5u <- family.test(data=family_generated_3g.3a.5u, f=risk.variant.id)
  result.trafic.ext.3g.3a.5u <- family.test.trafic.ext(data=family_generated_3g.3a.5u, f=risk.variant.id)
  #only report p.value
  c(result.trap.2g.3a.1u[7], result.trafic.ext.2g.3a.1u, 
    result.trap.3g.3a.4u[7], result.trafic.ext.3g.3a.4u,
    result.trap.3g.2a.5u[7], result.trafic.ext.3g.2a.5u,
    result.trap.3g.3a.5u[7], result.trafic.ext.3g.3a.5u)
})
## Write out your results to a csv file
result.df <- as.data.frame(t(sim_result))
colnames(result.df) <- c("result.trap.2g.3a.1u","result.trafic.ext.2g.3a.1u",
                         "result.trap.3g.3a.4u","result.trafic.ext.3g.3a.4u",
                         "result.trap.3g.2a.5u","result.trafic.ext.3g.2a.5u",
                         "result.trap.3g.3a.5u","result.trafic.ext.3g.3a.5u")
result.df <- cbind(seed,r,p_dis,n_family,result.df)
write.csv(result.df, paste("res_",r,"_",n_family,"_",seed,".csv",sep=""), row.names=FALSE)


#power curve without allelic heterogeneity
result <- read.csv("2015-04-20 2g 3g comparison.csv", header=T)
result <- select(result, -c(p_dis,n_family)) %>% melt(id.vars=c("seed", "r"), variable.name="method", value.name="p.value")
result.power <- result %>% group_by(r, method) %>% 
  summarise(n=n(), power=mean(p.value<0.05))
result.power <- mutate(result.power, type = grepl("trap", method))

p1 <- ggplot(result.power[grepl("result", result.power$method),], aes(x=r, y=power, col=method, lty=type)) +
  #   geom_point(size=3, alpha=1) +
  geom_line(size=1.2, alpha=0.7) +
  ggtitle("f=0.0102, consider effect size of risk haplotypes") +
  labs(x="relative risk r") +
  theme_gray(base_size = 20) +
  scale_color_discrete(name="Method", 
                       breaks=c("result.trap.3g.3a.4u","result.trap.2g.3a.1u",
                                "result.trap.3g.3a.5u","result.trap.3g.2a.5u",
                                "result.trafic.ext.2g.3a.1u","result.trafic.ext.3g.3a.4u",
                                "result.trafic.ext.3g.3a.5u","result.trafic.ext.3g.2a.5u")) +
  scale_linetype_discrete(name="Test", 
                          breaks=c("TRUE","FALSE"), labels=c("TRAP", "TRAFIC"))


#Load results with allelic heterogeneity
result <- read.csv("2015-04-30 2g 3g comparison with allelic het.csv", header=T)
result <- select(result, -c(p_dis,n_family)) %>% melt(id.vars=c("seed", "r"), variable.name="method", value.name="p.value")
result.power <- result %>% group_by(r, method) %>% 
  summarise(n=n(), power=mean(p.value<0.05))
result.power <- mutate(result.power, type = grepl("trap", method))

p2 <- ggplot(result.power[grepl("result", result.power$method),], aes(x=r, y=power, col=method, lty=type)) +
  #   geom_point(size=3, alpha=1) +
  geom_line(size=1.2, alpha=0.7) +
  ggtitle("f=0.0155, consider effect size of each risk variant. i.e. allelic heterogeneity") +
  labs(x="relative risk r") +
  theme_gray(base_size = 20) +
  scale_color_discrete(name="Method", 
                       breaks=c("result.trap.3g.3a.4u","result.trap.2g.3a.1u",
                                "result.trap.3g.3a.5u","result.trap.3g.2a.5u",
                                "result.trafic.ext.2g.3a.1u","result.trafic.ext.3g.3a.4u",
                                "result.trafic.ext.3g.3a.5u","result.trafic.ext.3g.2a.5u")) +
  scale_linetype_discrete(name="Test", 
                          breaks=c("TRUE","FALSE"), labels=c("TRAP", "TRAFIC"))

##make two plots in one
multiplot(p1,p2, cols=2)


##version 2 
result <- read.csv("2015-04-20 2g 3g comparison.csv", header=T)
result <- select(result, -c(p_dis,n_family)) %>% melt(id.vars=c("seed", "r"), variable.name="method", value.name="p.value")
result.power <- result %>% group_by(r, method) %>% 
  summarise(n=n(), power=mean(p.value<0.05))
result.power.p1 <- mutate(result.power, type = grepl("trap", method))
result.power.p1$het <- F

result <- read.csv("2015-04-30 2g 3g comparison with allelic het.csv", header=T)
result <- select(result, -c(p_dis,n_family)) %>% melt(id.vars=c("seed", "r"), variable.name="method", value.name="p.value")
result.power <- result %>% group_by(r, method) %>% 
  summarise(n=n(), power=mean(p.value<0.05))
result.power.p2 <- mutate(result.power, type = grepl("trap", method))
result.power.p2$het <- T

result.power <- rbind(result.power.p1, result.power.p2)

ggplot(result.power[grepl("result", result.power$method),], aes(x=r, y=power, col=method, lty=type)) +
  #   geom_point(size=3, alpha=1) +
  facet_grid(.~het) +
  geom_line(size=1.2, alpha=0.7) +
  ggtitle("f=0.0155, consider effect size of each risk variant. i.e. allelic heterogeneity") +
  labs(x="relative risk r") +
  theme_bw(base_size = 17) +
  scale_color_discrete(name="Method", 
                       breaks=c("result.trap.3g.3a.4u","result.trap.2g.3a.1u",
                                "result.trap.3g.3a.5u","result.trap.3g.2a.5u",
                                "result.trafic.ext.2g.3a.1u","result.trafic.ext.3g.3a.4u",
                                "result.trafic.ext.3g.3a.5u","result.trafic.ext.3g.2a.5u")) +
  scale_linetype_discrete(name="Test", 
                          breaks=c("TRUE","FALSE"), labels=c("TRAP", "TRAFIC"))

# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}