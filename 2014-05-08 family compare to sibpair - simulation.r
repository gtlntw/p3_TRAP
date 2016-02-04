source("function.r")

sim_result <- sim_family(m=1, var=0, f=0.01, SRR=5, n_family=10, rep=20)

## ... your simualtion code
sim_result <- sim_family(m=r, var=v, f=f, SRR=5, n_family=500, rep=n_rep)

## Write out your results to a csv file
write.csv(data.frame(seed=seed, f=f, r=r, v=v, true=sim_result),
          paste("res_",f,"_",r,"_",v,"_",seed,".csv",sep=""), row.names=FALSE)