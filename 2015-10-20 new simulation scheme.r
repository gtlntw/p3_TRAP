source("function.r")
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

##2. gives structure and 
family_strct_ped <- select_family_strct("2g.3a.1u")

##3. simulate data using the new model
gene_family_pe <- function(family_strct=family_strct_ped, n_family= 1000, p_dis=0.1, Beta=Beta) {
	#get genotype matrix
	n_haplo <- nrow(haplotype)
	n_snp <- nrow(snp)
  n_family_member <- length(family_strct$person)
  data_family <- matrix(NA, nrow=n_family*n_family_member, ncol=(6+n_snp*2))
  tran_vec <- matrix(NA, nrow=n_family*n_family_member, ncol=3)
  affect_spec <- family_strct$affect-1 # afffect statuts
  haplotype.effect <- as.matrix(2- haplotype[, -(1:2)]) %*% Beta #precalculate to speed up
  #the strategy is generate the whole family and check affected status if fail then restart from the first individual
  data_family.idx <- 1 #index of the current family member being generated
  n_family.idx <- 1 #index of how many families have been generated
  while(n_family.idx <= n_family) {
    disease_vec <- matrix(NA, nrow=n_family_member, ncol=1) #store generated affected status
    family.haplo <- matrix(NA, nrow=n_family_member, ncol=2) #store haplotype info.
    ind.idx <- 1 # which individual is generating
    while(ind.idx <= n_family_member) { #until the generated matches the input
      #generate the current individual
      if(family_strct$father[ind.idx]==0 & family_strct$mother[ind.idx]==0) { #if founder
        family.haplo[ind.idx,] <- sample(1:n_haplo, 2, replace=T)
      } else{ #if not founder
        family.haplo[ind.idx,] <- c(sample(family.haplo[family_strct$father[ind.idx],], 1), sample(family.haplo[family_strct$mother[ind.idx],], 1))
      }
    	ind.idx <- ind.idx + 1
    }
		#save haplotype effect
    family.haplo.effect <- apply(family.haplo, 1, function(x) haplotype.effect[x[1]] + haplotype.effect[x[2]])
    
    #generate phenotype
    beta0 <- qlogis(p_dis) #baseline prevalence

		sigma_pe_sq <- 1.447215 #p*(1-p)/(exp(2*p)/(exp(p)+1)^4) to achieve 50% heritability
		K <- matrix(c(.5, 0, .25, .25,
									0, .5, .25, .25,
									.25, .25, .5, .25,
									.25, .25, .25, .5), nrow = 4, ncol = 4, byrow = T)
		PE <- 2*K*sigma_pe_sq #polygenic effect
    mu <- beta0 + family.haplo.effect
		
		p_star <- MASS::mvrnorm(1, mu=mu, Sigma=PE)
		(p <- plogis(p_star))
		(affect_status <- rbinom(rep(1,n_family_member), 1, p))
		
		if(all(affect_status==affect_spec)) {
	    #save the haplotype file
	    letter.idx <- 1 #indicator used in the transmission vector
	    for(i in 1:n_family_member) {
	      #store genope and transmission vector
	      data_family[data_family.idx, ] <- unlist(c(n_family.idx, family_strct[i,2:6], haplotype[family.haplo[i,1],-c(1:2)], haplotype[family.haplo[i,2],-c(1:2)]))
	      if(family_strct$father[i]==0 & family_strct$mother[i]==0) { #if founder
	        tran_vec[data_family.idx,] <- c(n_family.idx, LETTERS[letter.idx], LETTERS[letter.idx+1])
	        letter.idx <- letter.idx + 2
	      }
	      else{ #if not founder then compare with his/her parents
	        current_row <- (n_family.idx-1)*n_family_member
	        #print(current_row)
	        tran_vec[data_family.idx,] <- c(n_family.idx, ifelse(family.haplo[i,1] == family.haplo[family_strct$father[i],1], tran_vec[family_strct$father[i]+current_row, 2], tran_vec[family_strct$father[i]+current_row, 3]) 
	                                        ,ifelse(family.haplo[i,2] == family.haplo[family_strct$mother[i],1], tran_vec[family_strct$mother[i]+current_row, 2], tran_vec[family_strct$mother[i]+current_row, 3])) 
	      }
	      data_family.idx <- data_family.idx + 1
	    }
    #     print(n_family.idx)
    n_family.idx <- n_family.idx + 1	
    }
  }
  colnames(data_family) <- c("family","person","father","mother","sex","affect",rep(paste("SNP", 1:n_snp, sep=""),2))
  colnames(tran_vec) <- c("family","h1","h2")
  return(list(data_family=data.frame(data_family, stringsAsFactors=F), tran_vec=data.frame(tran_vec, stringsAsFactors=F)))
}


Beta <- rep(0, length.out=50)
Beta[7] <- log(1) #effect size of OR=2
gene_family_pe(family_strct=family_strct_ped, n_family= 100, p_dis=0.1, Beta=Beta)

#kinship matrix
#convert to probablility and generate phenotype
##4. 


# 
# #compare to old way
# ## ... your simualtion codek
# #generate haplotype risk
# haplotype.risk <- generate_haplorisk(r=r, het=T)
# 
# ##gene drop simulation for two generations
# family_strct_ped <- select_family_strct("2g.3a.1u")
# 
# 
# rep.idx <<- 1
# sim_result <- replicate(n_rep, {
#   print(rep.idx)
#   rep.idx <<- rep.idx + 1
# 
#   family_generated_2g2c2a <<- gene_family(family_strct=family_strct_ped, n=n_family, haplotype.risk=haplotype.risk) 
#   
#   ##trap w. and w.o. IBD information
#   result.trap.2g2c2a <- family.test(data=family_generated_2g2c2a, f=risk.variant)$p.value
#   result.trap.2g2c2a.noIBD <- family.test(data=family_generated_2g2c2a, f=risk.variant, noIBD=T)$p.value
#   
#   #only report p.value
#   c(result.trap.2g2c2a, result.trap.2g2c2a.noIBD)
# })
# ## Write out your results to a csv file
# result.df <- as.data.frame(t(sim_result))
# colnames(result.df) <- c("result.trap.2g2c2a", "result.trap.2g2c2a.noIBD")
# result.df <- cbind(seed,r,p_dis,risk.variant,risk.haplo.f,n_family,result.df)
# write.csv(result.df, paste("res_",r,"_",n_family,"_",seed,".csv",sep=""), row.names=FALSE)
# ## R.miSniam
# ## End(Not run)
# ## Not run:
