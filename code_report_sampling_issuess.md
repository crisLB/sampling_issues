# Code report: Create incomplete sampling and taxonomic uncertainty in bipartite networks

Cristina Llopis-Belenguer

corresponding author: cristina.llopis.belenguer@gmail.com 


## 1. Utilities
### 1.1 Set a working directory
This is an optional step. Researchers could be interested in saving data and results in the same folder, the working directory. If data is stored in the working directory they will call data everytime with the name of the file and it will not be necessary to write the path to a file anymore.
``` r
setwd("path_to_the_directory")
```

### 1.2 Libraries and accessory functions
``` r
library(bipartite)
library(ade4)
library(parallel)
library(cluster)
library(plotrix)
library(ggplot2)
library(gridExtra)
library(viridis)
```
Load functions available in Fründ et al. 2016  https://doi.org/10.1111/oik.02256
``` r
source("get_skewabuns.R")
source("makeweb.R")
source("make_trueweb.R")
source("sampleweb.R")
```
Load functions to use in the taxonomic resolution analysis. After the allocation of parasite species to groups by the pairwise Euclidean similarities index (see analyses in section 2.3 below), these functions allow to aggregate interactions of parasites allocated to the same group by the Euclidean index.
``` r
source("taxonomic_resolution_functions.R")
```

### 1.3 Data availability
Researchers can create their own data. Alternatively, they may want to use data available in:

Llopis-Belenguer, Cristina, Balbuena, Juan Antonio, Blasco-Costa, Isabel, Karvonen, Anssi, Sarabeev, Volodimir, & Jokela, Jukka. (2022). Replication data for: Sensitivity of bipartite network analyses to incomplete sampling and taxonomic uncertainty [Data set]. Zenodo. https://doi.org/10.5281/zenodo.7386012



## 2. Create communities
### 2.1 Full communities
The full communities are weighted bipartite networks reflecting a chosen specialisation parameter. Follow this steps to simulate host and parasite species abundance distributions. Hosts are equally sampled. Parasite abundances follow a log-normal distribution.
``` r
nhost<-13
npara<-42
hostabun<-rep(1,nhost)
paraabun<-get_skewabuns(npara,5.89,1.45)
```
Simulate the mean number of interactions per parasite species.
``` r
ni<-86794
maxobs.rf<-round(ni/npara)
```
Create 10-replicate full communities.
``` r
myweb.p.list<-vector(mode = "list", length = 10)
myweb.p.list<-lapply(myweb.p.list, function (x) makeweb(specpar=55, npara, nhost, nicheshape="power"))
myweb.true.list<-lapply(myweb.p.list, function(x) make_trueweb(x, hostabun, paraabun))
myweb.full.list<-lapply(myweb.true.list, function(x) sampleweb(x, maxobs.rf, method="perweb"))
```
Alternatively, load the 10-replicate full communities in Llopis-Belenguer et al. (2022).
``` r
myweb.full.list<-readRDS("myweb.full.list.RDS")
```

### 2.2 Resampled communities with reduced effort in host sampling completeness
I simulated a decreasing number of sampled host individuals by means of reducing the mean number of interactions per parasite species in the 10-replicate full communities. I used myweb.true.list (see Create full communities) to build the resampled communities with a reduced effort in host sampling completeness. maxobs.rf was multiplied by the the mean percentage of interactions per parasite species that remained in all datasets when I reduced the number of host individuals in 10% steps starting from 90% and ending with 10% of the original host individuals.
``` r
# Mean percentage of interactions per parasite species that remained when we removed from 90% to 10% of host individuals
i.mean<-c(8.230,17.946,27.746,36.802,50.472,57.180,69.348,78.454,85.078,100.000) # mean number of interactions

# Resampled communities with 90% of host individuals (or 85.08% of interactions) with respect to the full communities
myweb.full.90.list<-lapply(myweb.true.list, function(x) sampleweb(x, as.numeric(round(maxobs.rf*i.mean[9]/100)), method='perweb'))

# Resampled communities with 80%-10% of host individuals with respect to the full communities
myweb.full.80.list<-lapply(myweb.true.list, function(x) sampleweb(x, as.numeric(round(maxobs.rf*i.mean[8]/100)), method='perweb'))
myweb.full.70.list<-lapply(myweb.true.list, function(x) sampleweb(x, as.numeric(round(maxobs.rf*i.mean[7]/100)), method='perweb'))
myweb.full.60.list<-lapply(myweb.true.list, function(x) sampleweb(x, as.numeric(round(maxobs.rf*i.mean[6]/100)), method='perweb'))
myweb.full.50.list<-lapply(myweb.true.list, function(x) sampleweb(x, as.numeric(round(maxobs.rf*i.mean[5]/100)), method='perweb'))
myweb.full.40.list<-lapply(myweb.true.list, function(x) sampleweb(x, as.numeric(round(maxobs.rf*i.mean[4]/100)), method='perweb'))
myweb.full.30.list<-lapply(myweb.true.list, function(x) sampleweb(x, as.numeric(round(maxobs.rf*i.mean[3]/100)), method='perweb'))
myweb.full.20.list<-lapply(myweb.true.list, function(x) sampleweb(x, as.numeric(round(maxobs.rf*i.mean[2]/100)), method='perweb'))
myweb.full.10.list<-lapply(myweb.true.list, function(x) sampleweb(x, as.numeric(round(maxobs.rf*i.mean[1]/100)), method='perweb'))
```
Alternatively, load the resampled communities with reduced effort in host sampling completeness in Llopis-Belenguer et al. (2022).
``` r
myweb.full.90.list<-readRDS("myweb.full.90.list.RDS")
myweb.full.80.list<-readRDS("myweb.full.80.list.RDS")
myweb.full.70.list<-readRDS("myweb.full.70.list.RDS")
myweb.full.60.list<-readRDS("myweb.full.60.list.RDS")
myweb.full.50.list<-readRDS("myweb.full.50.list.RDS")
myweb.full.40.list<-readRDS("myweb.full.40.list.RDS")
myweb.full.30.list<-readRDS("myweb.full.30.list.RDS")
myweb.full.20.list<-readRDS("myweb.full.20.list.RDS")
myweb.full.10.list<-readRDS("myweb.full.10.list.RDS")

```

### 2.3 Resampled communities with reduced effort in parasite taxonomic resolution
To simulate the effect of incomplete parasite taxonomic resolution, I first calculated pairwise Euclidean similarities between all parasite species based on their host-range overlap. I then used these pairwise distance estimates to distribute the parasite species to a requested number of groups and aggregate the interactions within groups. The number of groups ranged from 90% to 10% of the parasite species in each replicate of the full communities and resampled communities with reduced effort in host sampling completeness (the latter represent the crossed-biased communities), in 10% steps. These simulations gave a gradient of 90 resampled communities in a decreasing gradient of parasite taxonomic resolution. Additionally, I obtained 810 resampled communities, which were biased for both sampling issues to varying degrees.

``` r
# complete set of full and Resampled communities with reduced effort in host sampling complentess
# This a list of 10 lists of 10 matrices each
myweb.full.comp.list<-list(myweb.full.list, 
                           myweb.full.90.list,
                           myweb.full.80.list,
                           myweb.full.70.list,
                           myweb.full.60.list,
                           myweb.full.50.list,
                           myweb.full.40.list,
                           myweb.full.30.list,
                           myweb.full.20.list,
                           myweb.full.10.list)

# Give names to the original communities
hp.names<-myweb.full.comp.list
for (i in seq_along(hp.names)) {
  for (j in seq_along(hp.names[[i]])){
    row.names(hp.names[[i]][[j]])<-paste(rep("hos",nrow(hp.names[[i]][[j]])), 1:nrow(hp.names[[i]][[j]]), sep = ".")
    colnames(hp.names[[i]][[j]])<-paste(rep("par",ncol(hp.names[[i]][[j]])), 1:ncol(hp.names[[i]][[j]]), sep = ".")
  }
}


# 90% of Parasite taxonomic resolution
j.90<-vector('list', length(hp.names)) # j: jaccard
for (i in seq_along(hp.names)) {
  for (j in seq_along(hp.names[[i]])){
    j.90[[i]][[j]]<-pam(dist(x=t(hp.names[[i]][[j]]), method="euclidean"), k=ncol(hp.names[[i]][[j]])*0.90, diss=T)$clustering # cluster parasites
  }
}


cp.90 <- addgroupVar(hp.names, j.90)
cp.90 <- transposeAll(cp.90)
cp.90 <- nameGroup(cp.90)
cp.90 <- aggregateAll(cp.90)
for (i in seq_along(cp.90)) {
  for (j in seq_along(cp.90[[i]])) {
    cp.90[[i]][[j]]<-cp.90[[i]][[j]][,-1]
  }
}
cp.90 <- transposeAll(cp.90) # All these communities are biased for parasite taxonomic resolution at 90% level. 
                             # cp.90[[1]] are the full communities, but biased for parasite taxonomic resolution 
                             # at 90% level.
                             # From cp.90[[2]] to cp.90[[10]] are the Resampled communities with reduced effort 
                             # in host sampling completenes but biased for 
                             # parasite taxonomic resolution at 90% level.


# 80% of Parasite taxonomic resolution
j.80<-vector('list', length(hp.names))
for (i in seq_along(hp.names)) {
  for (j in seq_along(hp.names[[i]])){
    j.80[[i]][[j]]<-pam(dist(x=t(hp.names[[i]][[j]]), method="euclidean"), k=ncol(hp.names[[i]][[j]])*0.80, diss=T)$clustering # cluster parasites
  }
}

cp.80 <- addgroupVar(hp.names, j.80)
cp.80 <- transposeAll(cp.80)
cp.80 <- nameGroup(cp.80)
cp.80 <- aggregateAll(cp.80)
for (i in seq_along(cp.80)) {
  for (j in seq_along(cp.80[[i]])) {
    cp.80[[i]][[j]]<-cp.80[[i]][[j]][,-1]
  }
}
cp.80 <- transposeAll(cp.80)


# 70% of Parasite taxonomic resolution
j.70<-vector('list', length(hp.names))
for (i in seq_along(hp.names)) {
  for (j in seq_along(hp.names[[i]])){
    j.70[[i]][[j]]<-pam(dist(x=t(hp.names[[i]][[j]]), method="euclidean"), k=ncol(hp.names[[i]][[j]])*0.70, diss=T)$clustering # cluster parasites
  }
}

cp.70 <- addgroupVar(hp.names, j.70)
cp.70 <- transposeAll(cp.70)
cp.70 <- nameGroup(cp.70)
cp.70 <- aggregateAll(cp.70)
for (i in seq_along(cp.70)) {
  for (j in seq_along(cp.70[[i]])) {
    cp.70[[i]][[j]]<-cp.70[[i]][[j]][,-1]
  }
}
cp.70 <- transposeAll(cp.70)


# 60% of Parasite taxonomic resolution
j.60<-vector('list', length(hp.names))
for (i in seq_along(hp.names)) {
  for (j in seq_along(hp.names[[i]])){
    j.60[[i]][[j]]<-pam(dist(x=t(hp.names[[i]][[j]]), method="euclidean"), k=ncol(hp.names[[i]][[j]])*0.60, diss=T)$clustering # cluster parasites
  }
}

cp.60 <- addgroupVar(hp.names, j.60)
cp.60 <- transposeAll(cp.60)
cp.60 <- nameGroup(cp.60)
cp.60 <- aggregateAll(cp.60)
for (i in seq_along(cp.60)) {
  for (j in seq_along(cp.60[[i]])) {
    cp.60[[i]][[j]]<-cp.60[[i]][[j]][,-1]
  }
}
cp.60 <- transposeAll(cp.60)


# 50% of Parasite taxonomic resolution
j.50<-vector('list', length(hp.names))
for (i in seq_along(hp.names)) {
  for (j in seq_along(hp.names[[i]])){
    j.50[[i]][[j]]<-pam(dist(x=t(hp.names[[i]][[j]]), method="euclidean"), k=ncol(hp.names[[i]][[j]])*0.50, diss=T)$clustering # cluster parasites
  }
}

cp.50 <- addgroupVar(hp.names, j.50)
cp.50 <- transposeAll(cp.50)
cp.50 <- nameGroup(cp.50)
cp.50 <- aggregateAll(cp.50)
for (i in seq_along(cp.50)) {
  for (j in seq_along(cp.50[[i]])) {
    cp.50[[i]][[j]]<-cp.50[[i]][[j]][,-1]
  }
}
cp.50 <- transposeAll(cp.50)


# 40% of Parasite taxonomic resolution
j.40<-vector('list', length(hp.names))
for (i in seq_along(hp.names)) {
  for (j in seq_along(hp.names[[i]])){
    j.40[[i]][[j]]<-pam(dist(x=t(hp.names[[i]][[j]]), method="euclidean"), k=ncol(hp.names[[i]][[j]])*0.40, diss=T)$clustering # cluster parasites
  }
}

cp.40 <- addgroupVar(hp.names, j.40)
cp.40 <- transposeAll(cp.40)
cp.40 <- nameGroup(cp.40)
cp.40 <- aggregateAll(cp.40)
for (i in seq_along(cp.40)) {
  for (j in seq_along(cp.40[[i]])) {
    cp.40[[i]][[j]]<-cp.40[[i]][[j]][,-1]
  }
}
cp.40 <- transposeAll(cp.40)


# 30% of Parasite taxonomic resolution
j.30<-vector('list', length(hp.names))
for (i in seq_along(hp.names)) {
  for (j in seq_along(hp.names[[i]])){
    j.30[[i]][[j]]<-pam(dist(x=t(hp.names[[i]][[j]]), method="euclidean"), k=ncol(hp.names[[i]][[j]])*0.30, diss=T)$clustering # cluster parasites
  }
}

cp.30 <- addgroupVar(hp.names, j.30)
cp.30 <- transposeAll(cp.30)
cp.30 <- nameGroup(cp.30)
cp.30 <- aggregateAll(cp.30)
for (i in seq_along(cp.30)) {
  for (j in seq_along(cp.30[[i]])) {
    cp.30[[i]][[j]]<-cp.30[[i]][[j]][,-1]
  }
}
cp.30 <- transposeAll(cp.30)


# 20% of Parasite taxonomic resolution
j.20<-vector('list', length(hp.names))
for (i in seq_along(hp.names)) {
  for (j in seq_along(hp.names[[i]])){
    j.20[[i]][[j]]<-pam(dist(x=t(hp.names[[i]][[j]]), method="euclidean"), k=ncol(hp.names[[i]][[j]])*0.20, diss=T)$clustering # cluster parasites
  }
}

cp.20 <- addgroupVar(hp.names, j.20)
cp.20 <- transposeAll(cp.20)
cp.20 <- nameGroup(cp.20)
cp.20 <- aggregateAll(cp.20)
for (i in seq_along(cp.20)) {
  for (j in seq_along(cp.20[[i]])) {
    cp.20[[i]][[j]]<-cp.20[[i]][[j]][,-1]
  }
}
cp.20 <- transposeAll(cp.20)


# 10% of Parasite taxonomic resolution
j.10<-vector('list', length(hp.names))
for (i in seq_along(hp.names)) {
  for (j in seq_along(hp.names[[i]])){
    j.10[[i]][[j]]<-pam(dist(x=t(hp.names[[i]][[j]]), method="euclidean"), k=ncol(hp.names[[i]][[j]])*0.10, diss=T)$clustering # cluster parasites
  }
}

cp.10 <- addgroupVar(hp.names, j.10)
cp.10 <- transposeAll(cp.10)
cp.10 <- nameGroup(cp.10)
cp.10 <- aggregateAll(cp.10)
for (i in seq_along(cp.10)) {
  for (j in seq_along(cp.10[[i]])) {
    cp.10[[i]][[j]]<-cp.10[[i]][[j]][,-1]
  }
}
cp.10 <- transposeAll(cp.10)
```
Alternatively, load the resampled communities with reduced effort in  parasite taxonomic resolution in Llopis-Belenguer et al. (2022).
``` r
cp.90.list<-readRDS("cp.90.RDS")
cp.80.list<-readRDS("cp.80.RDS")
cp.70.list<-readRDS("cp.70.RDS")
cp.60.list<-readRDS("cp.60.RDS")
cp.50.list<-readRDS("cp.50.RDS")
cp.40.list<-readRDS("cp.40.RDS")
cp.30.list<-readRDS("cp.30.RDS")
cp.20.list<-readRDS("cp.20.RDS")
cp.10.list<-readRDS("cp.10.RDS")
```



## 3. Community and species level descriptors
### 3.1 Community descriptors
#### 3.1.1 Modularity and Nestedness

* Full communities

Note that I run R in parallel copies to minimise the computational time.

First, I created 1,000 null communities for each of the 10 full communities with the function swap.web.
Then, we measured modularity and nestedness of the full and null communities. Finally, we calculated the standardised modularity and nestedness of the full communities.
``` r
# Run R in parallel
expfuncs<-lsf.str(pos=search()[which(search()=="package:bipartite")])
cl <- makeCluster(7, type ="PSOCK")
clusterExport(cl,varlist=expfuncs)

# Null communities
nullwebs.list<-parLapply(cl, myweb.full.list, function(x) swap.web(1000,x))


# Modularity of full communities
webM.list<-parLapply(cl, myweb.full.list, function(x) computeModules(x))
webM.like.list<-parLapply(cl, webM.list, function(x) attributes(x)$likelihood)

# Modularity of null communities
nullM.list<-vector('list', length(nullwebs.list))
for(i in seq_along(nullwebs.list)){
  for(j in seq_along(nullwebs.list[[i]])){
    nullM.list[[i]][[j]] <- computeModules(nullwebs.list[[i]][[j]])@likelihood
  }
}
nullM.list<-parLapply(cl, nullM.list, function(x) unlist(x))

# Nestedness of full communities
webN.list<-parLapply(cl, myweb.full.list, function(x) networklevel(x,index="weighted NODF"))

# Nestedness of null communities
nullN.list<-vector("list", length(nullwebs.list))
for(i in seq_along(nullwebs.list)){
  for(j in seq_along(nullwebs.list[[i]])){
    nullN.list[[i]][[j]] <- networklevel(nullwebs.list[[i]][[j]],index="weighted NODF")
  }
}

# Stop running R in parallel
stopCluster(cl)

# Standardised modularity
webSM.list<-(unlist(webM.like.list)-unlist(lapply(nullM.list,mean)))/unlist(lapply(nullM.list,sd))

# Standardised nestedness
webSN.list<-(unlist(webN.list)-unlist(lapply(lapply(nullN.list, unlist),mean)))/unlist(lapply(lapply(nullN.list, unlist),sd))

```

* Resampled communities with reduced effort in host sampling completeness

I created 1,000 null communities for each of the 90 resampled communities with reduced effort in host sampling completeness. I measured modularity and nestedness of the resampled and null communities. Finally, I calculated standardised modularity and nestedness of the resampled communities with reduced effort in host sampling completeness.

``` r
# Run R in parallel
expfuncs<-lsf.str(pos=search()[which(search()=="package:bipartite")])
cl <- makeCluster(7, type ="PSOCK")
clusterExport(cl,varlist=expfuncs)

# Null communities
nullwebs90.list<-parLapply(cl, myweb.full.90.list, function(x) swap.web(1000,x))
nullwebs80.list<-parLapply(cl, myweb.full.80.list, function(x) swap.web(1000,x))
nullwebs70.list<-parLapply(cl, myweb.full.70.list, function(x) swap.web(1000,x))
nullwebs60.list<-parLapply(cl, myweb.full.60.list, function(x) swap.web(1000,x))
nullwebs50.list<-parLapply(cl, myweb.full.50.list, function(x) swap.web(1000,x))
nullwebs40.list<-parLapply(cl, myweb.full.40.list, function(x) swap.web(1000,x))
nullwebs30.list<-parLapply(cl, myweb.full.30.list, function(x) swap.web(1000,x))
nullwebs20.list<-parLapply(cl, myweb.full.20.list, function(x) swap.web(1000,x))
nullwebs10.list<-parLapply(cl, myweb.full.10.list, function(x) swap.web(1000,x))

# Modularity of resampled communities
# 90% of host sampling completeness
web90M.list<-parLapply(cl, myweb.full.90.list, function(x) computeModules(x))
web90M.like.list<-parLapply(cl, web90M.list, function(x) attributes(x)$likelihood)

# Modularity of null communities
# 90
null90M.list<-vector('list', length(nullwebs90.list))
for(i in seq_along(nullwebs90.list)){
  for(j in seq_along(nullwebs90.list[[i]])){
    null90M.list[[i]][[j]] <- computeModules(nullwebs90.list[[i]][[j]])@likelihood
  }
}
null90M.list<-parLapply(cl, null90M.list, function(x) unlist(x))

# Nestedness of resampled communities
# 90
web90N.list<-parLapply(cl, myweb.full.90.list, function(x) networklevel(x,index="weighted NODF"))

# Nestedness of resampled communities
null90N.list<-vector("list", length(nullwebs90.list))
for(i in seq_along(nullwebs90.list)){
  for(j in seq_along(nullwebs90.list[[i]])){
    null90N.list[[i]][[j]] <- networklevel(nullwebs90.list[[i]][[j]],index="weighted NODF")
  }
}

# Standardised modularity
web90SM.list<-(unlist(web90M.like.list)-unlist(lapply(null90M.list,mean)))/unlist(lapply(null90M.list,sd))

# Standardised nestedness
web90SN.list<-(unlist(web90N.list)-unlist(lapply(lapply(null90N.list,unlist),mean)))/unlist(lapply(lapply(null90N.list,unlist),sd))


# Modularity of resampled communities
# 80% of host sampling completeness
web80M.list<-parLapply(cl, myweb.full.80.list, function(x) computeModules(x))
web80M.like.list<-parLapply(cl, web80M.list, function(x) attributes(x)$likelihood)

# Modularity of null communities
# 80
null80M.list<-vector('list', length(nullwebs80.list))
for(i in seq_along(nullwebs80.list)){
  for(j in seq_along(nullwebs80.list[[i]])){
    null80M.list[[i]][[j]] <- computeModules(nullwebs80.list[[i]][[j]])@likelihood
  }
}
null80M.list<-parLapply(cl, null80M.list, function(x) unlist(x))

# Nestedness of resampled communities
# 80
web80N.list<-parLapply(cl, myweb.full.80.list, function(x) networklevel(x,index="weighted NODF"))

# Nestedness of resampled communities
null80N.list<-vector("list", length(nullwebs80.list))
for(i in seq_along(nullwebs80.list)){
  for(j in seq_along(nullwebs80.list[[i]])){
    null80N.list[[i]][[j]] <- networklevel(nullwebs80.list[[i]][[j]],index="weighted NODF")
  }
}

# Standardised modularity
web80SM.list<-(unlist(web80M.like.list)-unlist(lapply(null80M.list,mean)))/unlist(lapply(null80M.list,sd))

# Standardised nestedness
web80SN.list<-(unlist(web80N.list)-unlist(lapply(lapply(null80N.list,unlist),mean)))/unlist(lapply(lapply(null80N.list,unlist),sd))


# Modularity of resampled communities
# 70% of host sampling completeness
web70M.list<-parLapply(cl, myweb.full.70.list, function(x) computeModules(x))
web70M.like.list<-parLapply(cl, web70M.list, function(x) attributes(x)$likelihood)

# Modularity of null communities
# 70
null70M.list<-vector('list', length(nullwebs70.list))
for(i in seq_along(nullwebs70.list)){
  for(j in seq_along(nullwebs70.list[[i]])){
    null70M.list[[i]][[j]] <- computeModules(nullwebs70.list[[i]][[j]])@likelihood
  }
}
null70M.list<-parLapply(cl, null70M.list, function(x) unlist(x))

# Nestedness of resampled communities
# 70
web70N.list<-parLapply(cl, myweb.full.70.list, function(x) networklevel(x,index="weighted NODF"))

# Nestedness of resampled communities
null70N.list<-vector("list", length(nullwebs70.list))
for(i in seq_along(nullwebs70.list)){
  for(j in seq_along(nullwebs70.list[[i]])){
    null70N.list[[i]][[j]] <- networklevel(nullwebs70.list[[i]][[j]],index="weighted NODF")
  }
}

# Standardised modularity
web70SM.list<-(unlist(web70M.like.list)-unlist(lapply(null70M.list,mean)))/unlist(lapply(null70M.list,sd))

# Standardised nestedness
web70SN.list<-(unlist(web70N.list)-unlist(lapply(lapply(null70N.list,unlist),mean)))/unlist(lapply(lapply(null70N.list,unlist),sd))


# Modularity of resampled communities
# 60% of host sampling completeness
web60M.list<-parLapply(cl, myweb.full.60.list, function(x) computeModules(x))
web60M.like.list<-parLapply(cl, web60M.list, function(x) attributes(x)$likelihood)

# Modularity of null communities
# 60
null60M.list<-vector('list', length(nullwebs60.list))
for(i in seq_along(nullwebs60.list)){
  for(j in seq_along(nullwebs60.list[[i]])){
    null60M.list[[i]][[j]] <- computeModules(nullwebs60.list[[i]][[j]])@likelihood
  }
}
null60M.list<-parLapply(cl, null60M.list, function(x) unlist(x))

# Nestedness of resampled communities
# 60
web60N.list<-parLapply(cl, myweb.full.60.list, function(x) networklevel(x,index="weighted NODF"))

# Nestedness of resampled communities
null60N.list<-vector("list", length(nullwebs60.list))
for(i in seq_along(nullwebs60.list)){
  for(j in seq_along(nullwebs60.list[[i]])){
    null60N.list[[i]][[j]] <- networklevel(nullwebs60.list[[i]][[j]],index="weighted NODF")
  }
}

# Standardised modularity
web60SM.list<-(unlist(web60M.like.list)-unlist(lapply(null60M.list,mean)))/unlist(lapply(null60M.list,sd))

# Standardised nestedness
web60SN.list<-(unlist(web60N.list)-unlist(lapply(lapply(null60N.list,unlist),mean)))/unlist(lapply(lapply(null60N.list,unlist),sd))


# Modularity of resampled communities
# 50% of host sampling completeness
web50M.list<-parLapply(cl, myweb.full.50.list, function(x) computeModules(x))
web50M.like.list<-parLapply(cl, web50M.list, function(x) attributes(x)$likelihood)


# Modularity of null communities
# 50
null50M.list<-vector('list', length(nullwebs50.list))
for(i in seq_along(nullwebs50.list)){
  for(j in seq_along(nullwebs50.list[[i]])){
    null50M.list[[i]][[j]] <- computeModules(nullwebs50.list[[i]][[j]])@likelihood
  }
}
null50M.list<-parLapply(cl, null50M.list, function(x) unlist(x))

# Nestedness of resampled communities
# 50
web50N.list<-parLapply(cl, myweb.full.50.list, function(x) networklevel(x,index="weighted NODF"))

# Nestedness of resampled communities
null50N.list<-vector("list", length(nullwebs50.list))
for(i in seq_along(nullwebs50.list)){
  for(j in seq_along(nullwebs50.list[[i]])){
    null50N.list[[i]][[j]] <- networklevel(nullwebs50.list[[i]][[j]],index="weighted NODF")
  }
}

# Standardised modularity
web50SM.list<-(unlist(web50M.like.list)-unlist(lapply(null50M.list,mean)))/unlist(lapply(null50M.list,sd))

# Standardised nestedness
web50SN.list<-(unlist(web50N.list)-unlist(lapply(lapply(null50N.list,unlist),mean)))/unlist(lapply(lapply(null50N.list,unlist),sd))


# Modularity of resampled communities
# 40% of host sampling completeness
web40M.list<-parLapply(cl, myweb.full.40.list, function(x) computeModules(x))
web40M.like.list<-parLapply(cl, web40M.list, function(x) attributes(x)$likelihood)

# Modularity of null communities
# 40
null40M.list<-vector('list', length(nullwebs40.list))
for(i in seq_along(nullwebs40.list)){
  for(j in seq_along(nullwebs40.list[[i]])){
    null40M.list[[i]][[j]] <- computeModules(nullwebs40.list[[i]][[j]])@likelihood
  }
}
null40M.list<-parLapply(cl, null40M.list, function(x) unlist(x))

# Nestedness of resampled communities
# 40
web40N.list<-parLapply(cl, myweb.full.40.list, function(x) networklevel(x,index="weighted NODF"))

# Nestedness of resampled communities
null40N.list<-vector("list", length(nullwebs40.list))
for(i in seq_along(nullwebs40.list)){
  for(j in seq_along(nullwebs40.list[[i]])){
    null40N.list[[i]][[j]] <- networklevel(nullwebs40.list[[i]][[j]],index="weighted NODF")
  }
}

# Standardised modularity
web40SM.list<-(unlist(web40M.like.list)-unlist(lapply(null40M.list,mean)))/unlist(lapply(null40M.list,sd))

# Standardised nestedness
web40SN.list<-(unlist(web40N.list)-unlist(lapply(lapply(null40N.list,unlist),mean)))/unlist(lapply(lapply(null40N.list,unlist),sd))


# Modularity of resampled communities
# 30% of host sampling completeness
web30M.list<-parLapply(cl, myweb.full.30.list, function(x) computeModules(x))
web30M.like.list<-parLapply(cl, web30M.list, function(x) attributes(x)$likelihood)

# Modularity of null communities
# 30
null30M.list<-vector('list', length(nullwebs30.list))
for(i in seq_along(nullwebs30.list)){
  for(j in seq_along(nullwebs30.list[[i]])){
    null30M.list[[i]][[j]] <- computeModules(nullwebs30.list[[i]][[j]])@likelihood
  }
}
null30M.list<-parLapply(cl, null30M.list, function(x) unlist(x))

# Nestedness of resampled communities
# 30
web30N.list<-parLapply(cl, myweb.full.30.list, function(x) networklevel(x,index="weighted NODF"))

# Nestedness of resampled communities
null30N.list<-vector("list", length(nullwebs30.list))
for(i in seq_along(nullwebs30.list)){
  for(j in seq_along(nullwebs30.list[[i]])){
    null30N.list[[i]][[j]] <- networklevel(nullwebs30.list[[i]][[j]],index="weighted NODF")
  }
}

# Standardised modularity
web30SM.list<-(unlist(web30M.like.list)-unlist(lapply(null30M.list,mean)))/unlist(lapply(null30M.list,sd))

# Standardised nestedness
web30SN.list<-(unlist(web30N.list)-unlist(lapply(lapply(null30N.list,unlist),mean)))/unlist(lapply(lapply(null30N.list,unlist),sd))


# Modularity of resampled communities
# 20% of host sampling completeness
web20M.list<-parLapply(cl, myweb.full.20.list, function(x) computeModules(x))
web20M.like.list<-parLapply(cl, web20M.list, function(x) attributes(x)$likelihood)

# Modularity of null communities
# 20
null20M.list<-vector('list', length(nullwebs20.list))
for(i in seq_along(nullwebs20.list)){
  for(j in seq_along(nullwebs20.list[[i]])){
    null20M.list[[i]][[j]] <- computeModules(nullwebs20.list[[i]][[j]])@likelihood
  }
}
null20M.list<-parLapply(cl, null20M.list, function(x) unlist(x))

# Nestedness of resampled communities
# 20
web20N.list<-parLapply(cl, myweb.full.20.list, function(x) networklevel(x,index="weighted NODF"))

# Nestedness of resampled communities
null20N.list<-vector("list", length(nullwebs20.list))
for(i in seq_along(nullwebs20.list)){
  for(j in seq_along(nullwebs20.list[[i]])){
    null20N.list[[i]][[j]] <- networklevel(nullwebs20.list[[i]][[j]],index="weighted NODF")
  }
}

# Standardised modularity
web20SM.list<-(unlist(web20M.like.list)-unlist(lapply(null20M.list,mean)))/unlist(lapply(null20M.list,sd))

# Standardised nestedness
web20SN.list<-(unlist(web20N.list)-unlist(lapply(lapply(null20N.list,unlist),mean)))/unlist(lapply(lapply(null20N.list,unlist),sd))


# Modularity of resampled communities
# 10% of host sampling completeness
web10M.list<-parLapply(cl, myweb.full.10.list, function(x) computeModules(x))
web10M.like.list<-parLapply(cl, web10M.list, function(x) attributes(x)$likelihood)

# Modularity of null communities
# 10
null10M.list<-vector('list', length(nullwebs10.list))
for(i in seq_along(nullwebs10.list)){
  for(j in seq_along(nullwebs10.list[[i]])){
    null10M.list[[i]][[j]] <- computeModules(nullwebs10.list[[i]][[j]])@likelihood
  }
}
null10M.list<-parLapply(cl, null10M.list, function(x) unlist(x))

# Nestedness of resampled communities
# 10
web10N.list<-parLapply(cl, myweb.full.10.list, function(x) networklevel(x,index="weighted NODF"))

# Nestedness of resampled communities
null10N.list<-vector("list", length(nullwebs10.list))
for(i in seq_along(nullwebs10.list)){
  for(j in seq_along(nullwebs10.list[[i]])){
    null10N.list[[i]][[j]] <- networklevel(nullwebs10.list[[i]][[j]],index="weighted NODF")
  }
}

# Standardised modularity
web10SM.list<-(unlist(web10M.like.list)-unlist(lapply(null10M.list,mean)))/unlist(lapply(null10M.list,sd))

# Standardised nestedness
web10SN.list<-(unlist(web10N.list)-unlist(lapply(lapply(null10N.list,unlist),mean)))/unlist(lapply(lapply(null10N.list,unlist),sd))

# Stop running R in parallel
stopCluster(cl)
```

* Resampled communities with reduced effort in parasite taxonomic resolution

I created 1,000 null communities for each of the 90 resampled communities with reduced effort in parasite taxonomic resolution and the 810 resampled communities with reduced effort in both sampling issues (the crossed effect). I measured modularity and nestedness of the resampled and null communities. Finally, I calculated standardised modularity and nestedness of the resampled communities.

``` r
expfuncs<- lsf.str(pos=search()[which(search()=="package:bipartite")])
cl <- makeCluster(7, type ="PSOCK")
clusterExport(cl,varlist=expfuncs)

# 90% of parasite taxonomic resolution
# Null communities
nullwebs.list.cp.90<-vector('list', length(cp.90))
for(i in seq_along(cp.90)){
  for(j in seq_along(cp.90[[i]])){
    nullwebs.list.cp.90[[i]][[j]]<-swap.web(1000,cp.90[[i]][[j]])
  }
}

# saveRDS(nullwebs.list.cp.90,"nullwebs.list.cp.90.RDS")

# Modularity of null communities
nullM.list.cp.90<-vector('list', length(nullwebs.list.cp.90))
nullM.list.cp.90<-lapply(nullM.list.cp.90, function(x) vector('list', length(nullwebs.list.cp.90[[1]])))
for(i in seq_along(nullM.list.cp.90)){
  for(j in seq_along(nullM.list.cp.90[[i]])){
    nullM.list.cp.90[[i]][[j]]<-vector('list', length(nullwebs.list.cp.90[[1]][[1]]))
  }
}

for(i in seq_along(nullwebs.list.cp.90)){
  for(j in seq_along(nullwebs.list.cp.90[[i]])){
    for(k in seq_along(nullwebs.list.cp.90[[i]][[j]])){
      nullM.list.cp.90[[i]][[j]][[k]]<- computeModules(nullwebs.list.cp.90[[i]][[j]][[k]])@likelihood
    }
  }
}

for(i in seq_along(nullM.list.cp.90)){
  for(j in seq_along(nullM.list.cp.90[[i]])){
    nullM.list.cp.90[[i]][[j]]<-unlist(nullM.list.cp.90[[i]][[j]])
  }
}
# saveRDS(nullM.list.cp.90,"nullM.list.cp.90.RDS")

nullM.list.cp.90.2<-vector('list', length(nullM.list.cp.90))
nullM.list.cp.90.2.mean<-vector('list', length(nullM.list.cp.90))
nullM.list.cp.90.2.sd<-vector('list', length(nullM.list.cp.90))
for(i in seq_along(nullM.list.cp.90)){
  nullM.list.cp.90.2[[i]]<-parLapply(cl, nullM.list.cp.90[[i]], function(x) unlist(x))
  nullM.list.cp.90.2.mean[[i]]<-parLapply(cl, nullM.list.cp.90.2[[i]],function(k) mean(k))
  nullM.list.cp.90.2.mean[[i]]<-unlist(nullM.list.cp.90.2.mean[[i]])
  nullM.list.cp.90.2.sd[[i]]<-parLapply(cl, nullM.list.cp.90.2[[i]],function(k) sd(k))
  nullM.list.cp.90.2.sd[[i]]<-unlist(nullM.list.cp.90.2.sd[[i]])
}

# Modularity of resampled communities
webM.list.cp.90<-vector('list', length(cp.90))
for(i in seq_along(cp.90)){
  webM.list.cp.90[[i]]<-parLapply(cl, cp.90[[i]], function(k) computeModules(k))
}

webM.like.list.cp.90<-vector('list', length(webM.list.cp.90))
for(i in seq_along(webM.list.cp.90)){
  for(j in seq_along(webM.list.cp.90[[i]])){
    webM.like.list.cp.90[[i]]<-parLapply(cl, webM.list.cp.90[[i]], function(k) attributes(k)$likelihood)
  }
}
webM.like.list.cp.90<-lapply(webM.like.list.cp.90, unlist)
webM.like.list.cp.90<-lapply(webM.like.list.cp.90, as.numeric)

# Standardised modulariy of resampled communities
webSM.like.list.cp.90<-Map("/",Map("-", webM.like.list.cp.90,nullM.list.cp.90.2.mean),nullM.list.cp.90.2.sd)


# Nestedness 
# Nestedness of null communities
nullN.list.cp.90<-vector('list', length(nullwebs.list.cp.90))
nullN.list.cp.90<-lapply(nullN.list.cp.90, function(x) vector('list', length(nullwebs.list.cp.90[[1]])))
for(i in seq_along(nullN.list.cp.90)){
  for(j in seq_along(nullN.list.cp.90[[i]])){
    nullN.list.cp.90[[i]][[j]]<-vector('list', length(nullwebs.list.cp.90[[1]][[1]]))
  }
}

for(i in seq_along(nullwebs.list.cp.90)){
  for(j in seq_along(nullwebs.list.cp.90[[i]])){
    for(k in seq_along(nullwebs.list.cp.90[[i]][[j]])){
      nullN.list.cp.90[[i]][[j]][[k]] <- networklevel(nullwebs.list.cp.90[[i]][[j]][[k]],index="weighted NODF")
    }
  }
}

for(i in seq_along(nullN.list.cp.90)){
  for(j in seq_along(nullN.list.cp.90[[i]])){
    nullN.list.cp.90[[i]][[j]]<-as.numeric(unlist(nullN.list.cp.90[[i]][[j]]))
  }
}

# saveRDS(nullN.list.cp.90,"nullN.list.cp.90.RDS")

nullN.list.cp.90.2<-vector('list', length(nullN.list.cp.90))
nullN.list.cp.90.2.mean<-vector('list', length(nullN.list.cp.90))
nullN.list.cp.90.2.sd<-vector('list', length(nullN.list.cp.90))
for(i in seq_along(nullN.list.cp.90)){
  nullN.list.cp.90.2[[i]]<-parLapply(cl, nullN.list.cp.90[[i]], function(x) unlist(x))
  nullN.list.cp.90.2.mean[[i]]<-parLapply(cl, nullN.list.cp.90.2[[i]],function(k) mean(k))
  nullN.list.cp.90.2.mean[[i]]<-unlist(nullN.list.cp.90.2.mean[[i]])
  nullN.list.cp.90.2.sd[[i]]<-parLapply(cl, nullN.list.cp.90.2[[i]],function(k) sd(k))
  nullN.list.cp.90.2.sd[[i]]<-unlist(nullN.list.cp.90.2.sd[[i]])
}

# Nestedness of resampled communities
webN.list.cp.90<-vector('list', length(cp.90))
for(i in seq_along(cp.90)){
  webN.list.cp.90[[i]]<-parLapply(cl, cp.90[[i]], function(x) networklevel(x,index="weighted NODF"))
}

webN.list.cp.90<-lapply(webN.list.cp.90, unlist)

# Standardised nestedness of resampled communities
webSN.list.cp.90<-Map("/",Map("-", webN.list.cp.90,nullN.list.cp.90.2.mean),nullN.list.cp.90.2.sd)


# 80% of parasite taxonomic resolution
# Null communities
nullwebs.list.cp.80<-vector('list', length(cp.80))
for(i in seq_along(cp.80)){
  for(j in seq_along(cp.80[[i]])){
    nullwebs.list.cp.80[[i]][[j]]<-swap.web(1000,cp.80[[i]][[j]])
  }
}

# saveRDS(nullwebs.list.cp.80,"nullwebs.list.cp.80.RDS")

# Modularity of null communities
nullM.list.cp.80<-vector('list', length(nullwebs.list.cp.80))
nullM.list.cp.80<-lapply(nullM.list.cp.80, function(x) vector('list', length(nullwebs.list.cp.80[[1]])))
for(i in seq_along(nullM.list.cp.80)){
  for(j in seq_along(nullM.list.cp.80[[i]])){
    nullM.list.cp.80[[i]][[j]]<-vector('list', length(nullwebs.list.cp.80[[1]][[1]]))
  }
}

for(i in seq_along(nullwebs.list.cp.80)){
  for(j in seq_along(nullwebs.list.cp.80[[i]])){
    for(k in seq_along(nullwebs.list.cp.80[[i]][[j]])){
      nullM.list.cp.80[[i]][[j]][[k]]<- computeModules(nullwebs.list.cp.80[[i]][[j]][[k]])@likelihood
    }
  }
}

for(i in seq_along(nullM.list.cp.80)){
  for(j in seq_along(nullM.list.cp.80[[i]])){
    nullM.list.cp.80[[i]][[j]]<-unlist(nullM.list.cp.80[[i]][[j]])
  }
}
# saveRDS(nullM.list.cp.80,"nullM.list.cp.80.RDS")

nullM.list.cp.80.2<-vector('list', length(nullM.list.cp.80))
nullM.list.cp.80.2.mean<-vector('list', length(nullM.list.cp.80))
nullM.list.cp.80.2.sd<-vector('list', length(nullM.list.cp.80))
for(i in seq_along(nullM.list.cp.80)){
  nullM.list.cp.80.2[[i]]<-parLapply(cl, nullM.list.cp.80[[i]], function(x) unlist(x))
  nullM.list.cp.80.2.mean[[i]]<-parLapply(cl, nullM.list.cp.80.2[[i]],function(k) mean(k))
  nullM.list.cp.80.2.mean[[i]]<-unlist(nullM.list.cp.80.2.mean[[i]])
  nullM.list.cp.80.2.sd[[i]]<-parLapply(cl, nullM.list.cp.80.2[[i]],function(k) sd(k))
  nullM.list.cp.80.2.sd[[i]]<-unlist(nullM.list.cp.80.2.sd[[i]])
}

# Modularity of resampled communities
webM.list.cp.80<-vector('list', length(cp.80))
for(i in seq_along(cp.80)){
  webM.list.cp.80[[i]]<-parLapply(cl, cp.80[[i]], function(k) computeModules(k))
}

webM.like.list.cp.80<-vector('list', length(webM.list.cp.80))
for(i in seq_along(webM.list.cp.80)){
  for(j in seq_along(webM.list.cp.80[[i]])){
    webM.like.list.cp.80[[i]]<-parLapply(cl, webM.list.cp.80[[i]], function(k) attributes(k)$likelihood)
  }
}
webM.like.list.cp.80<-lapply(webM.like.list.cp.80, unlist)
webM.like.list.cp.80<-lapply(webM.like.list.cp.80, as.numeric)

# Standardised modulariy of resampled communities
webSM.like.list.cp.80<-Map("/",Map("-", webM.like.list.cp.80,nullM.list.cp.80.2.mean),nullM.list.cp.80.2.sd)


# Nestedness 
# Nestedness of null communities
nullN.list.cp.80<-vector('list', length(nullwebs.list.cp.80))
nullN.list.cp.80<-lapply(nullN.list.cp.80, function(x) vector('list', length(nullwebs.list.cp.80[[1]])))
for(i in seq_along(nullN.list.cp.80)){
  for(j in seq_along(nullN.list.cp.80[[i]])){
    nullN.list.cp.80[[i]][[j]]<-vector('list', length(nullwebs.list.cp.80[[1]][[1]]))
  }
}

for(i in seq_along(nullwebs.list.cp.80)){
  for(j in seq_along(nullwebs.list.cp.80[[i]])){
    for(k in seq_along(nullwebs.list.cp.80[[i]][[j]])){
      nullN.list.cp.80[[i]][[j]][[k]] <- networklevel(nullwebs.list.cp.80[[i]][[j]][[k]],index="weighted NODF")
    }
  }
}

for(i in seq_along(nullN.list.cp.80)){
  for(j in seq_along(nullN.list.cp.80[[i]])){
    nullN.list.cp.80[[i]][[j]]<-as.numeric(unlist(nullN.list.cp.80[[i]][[j]]))
  }
}

# saveRDS(nullN.list.cp.80,"nullN.list.cp.80.RDS")

nullN.list.cp.80.2<-vector('list', length(nullN.list.cp.80))
nullN.list.cp.80.2.mean<-vector('list', length(nullN.list.cp.80))
nullN.list.cp.80.2.sd<-vector('list', length(nullN.list.cp.80))
for(i in seq_along(nullN.list.cp.80)){
  nullN.list.cp.80.2[[i]]<-parLapply(cl, nullN.list.cp.80[[i]], function(x) unlist(x))
  nullN.list.cp.80.2.mean[[i]]<-parLapply(cl, nullN.list.cp.80.2[[i]],function(k) mean(k))
  nullN.list.cp.80.2.mean[[i]]<-unlist(nullN.list.cp.80.2.mean[[i]])
  nullN.list.cp.80.2.sd[[i]]<-parLapply(cl, nullN.list.cp.80.2[[i]],function(k) sd(k))
  nullN.list.cp.80.2.sd[[i]]<-unlist(nullN.list.cp.80.2.sd[[i]])
}

# Nestedness of resampled communities
webN.list.cp.80<-vector('list', length(cp.80))
for(i in seq_along(cp.80)){
  webN.list.cp.80[[i]]<-parLapply(cl, cp.80[[i]], function(x) networklevel(x,index="weighted NODF"))
}

webN.list.cp.80<-lapply(webN.list.cp.80, unlist)

# Standardised nestedness of resampled communities
webSN.list.cp.80<-Map("/",Map("-", webN.list.cp.80,nullN.list.cp.80.2.mean),nullN.list.cp.80.2.sd)


# 70% of parasite taxonomic resolution
# Null communities
nullwebs.list.cp.70<-vector('list', length(cp.70))
for(i in seq_along(cp.70)){
  for(j in seq_along(cp.70[[i]])){
    nullwebs.list.cp.70[[i]][[j]]<-swap.web(1000,cp.70[[i]][[j]])
  }
}

# saveRDS(nullwebs.list.cp.70,"nullwebs.list.cp.70.RDS")

# Modularity of null communities
nullM.list.cp.70<-vector('list', length(nullwebs.list.cp.70))
nullM.list.cp.70<-lapply(nullM.list.cp.70, function(x) vector('list', length(nullwebs.list.cp.70[[1]])))
for(i in seq_along(nullM.list.cp.70)){
  for(j in seq_along(nullM.list.cp.70[[i]])){
    nullM.list.cp.70[[i]][[j]]<-vector('list', length(nullwebs.list.cp.70[[1]][[1]]))
  }
}

for(i in seq_along(nullwebs.list.cp.70)){
  for(j in seq_along(nullwebs.list.cp.70[[i]])){
    for(k in seq_along(nullwebs.list.cp.70[[i]][[j]])){
      nullM.list.cp.70[[i]][[j]][[k]]<- computeModules(nullwebs.list.cp.70[[i]][[j]][[k]])@likelihood
    }
  }
}

for(i in seq_along(nullM.list.cp.70)){
  for(j in seq_along(nullM.list.cp.70[[i]])){
    nullM.list.cp.70[[i]][[j]]<-unlist(nullM.list.cp.70[[i]][[j]])
  }
}
# saveRDS(nullM.list.cp.70,"nullM.list.cp.70.RDS")

nullM.list.cp.70.2<-vector('list', length(nullM.list.cp.70))
nullM.list.cp.70.2.mean<-vector('list', length(nullM.list.cp.70))
nullM.list.cp.70.2.sd<-vector('list', length(nullM.list.cp.70))
for(i in seq_along(nullM.list.cp.70)){
  nullM.list.cp.70.2[[i]]<-parLapply(cl, nullM.list.cp.70[[i]], function(x) unlist(x))
  nullM.list.cp.70.2.mean[[i]]<-parLapply(cl, nullM.list.cp.70.2[[i]],function(k) mean(k))
  nullM.list.cp.70.2.mean[[i]]<-unlist(nullM.list.cp.70.2.mean[[i]])
  nullM.list.cp.70.2.sd[[i]]<-parLapply(cl, nullM.list.cp.70.2[[i]],function(k) sd(k))
  nullM.list.cp.70.2.sd[[i]]<-unlist(nullM.list.cp.70.2.sd[[i]])
}

# Modularity of resampled communities
webM.list.cp.70<-vector('list', length(cp.70))
for(i in seq_along(cp.70)){
  webM.list.cp.70[[i]]<-parLapply(cl, cp.70[[i]], function(k) computeModules(k))
}

webM.like.list.cp.70<-vector('list', length(webM.list.cp.70))
for(i in seq_along(webM.list.cp.70)){
  for(j in seq_along(webM.list.cp.70[[i]])){
    webM.like.list.cp.70[[i]]<-parLapply(cl, webM.list.cp.70[[i]], function(k) attributes(k)$likelihood)
  }
}
webM.like.list.cp.70<-lapply(webM.like.list.cp.70, unlist)
webM.like.list.cp.70<-lapply(webM.like.list.cp.70, as.numeric)

# Standardised modulariy of resampled communities
webSM.like.list.cp.70<-Map("/",Map("-", webM.like.list.cp.70,nullM.list.cp.70.2.mean),nullM.list.cp.70.2.sd)


# Nestedness 
# Nestedness of null communities
nullN.list.cp.70<-vector('list', length(nullwebs.list.cp.70))
nullN.list.cp.70<-lapply(nullN.list.cp.70, function(x) vector('list', length(nullwebs.list.cp.70[[1]])))
for(i in seq_along(nullN.list.cp.70)){
  for(j in seq_along(nullN.list.cp.70[[i]])){
    nullN.list.cp.70[[i]][[j]]<-vector('list', length(nullwebs.list.cp.70[[1]][[1]]))
  }
}

for(i in seq_along(nullwebs.list.cp.70)){
  for(j in seq_along(nullwebs.list.cp.70[[i]])){
    for(k in seq_along(nullwebs.list.cp.70[[i]][[j]])){
      nullN.list.cp.70[[i]][[j]][[k]] <- networklevel(nullwebs.list.cp.70[[i]][[j]][[k]],index="weighted NODF")
    }
  }
}

for(i in seq_along(nullN.list.cp.70)){
  for(j in seq_along(nullN.list.cp.70[[i]])){
    nullN.list.cp.70[[i]][[j]]<-as.numeric(unlist(nullN.list.cp.70[[i]][[j]]))
  }
}

# saveRDS(nullN.list.cp.70,"nullN.list.cp.70.RDS")

nullN.list.cp.70.2<-vector('list', length(nullN.list.cp.70))
nullN.list.cp.70.2.mean<-vector('list', length(nullN.list.cp.70))
nullN.list.cp.70.2.sd<-vector('list', length(nullN.list.cp.70))
for(i in seq_along(nullN.list.cp.70)){
  nullN.list.cp.70.2[[i]]<-parLapply(cl, nullN.list.cp.70[[i]], function(x) unlist(x))
  nullN.list.cp.70.2.mean[[i]]<-parLapply(cl, nullN.list.cp.70.2[[i]],function(k) mean(k))
  nullN.list.cp.70.2.mean[[i]]<-unlist(nullN.list.cp.70.2.mean[[i]])
  nullN.list.cp.70.2.sd[[i]]<-parLapply(cl, nullN.list.cp.70.2[[i]],function(k) sd(k))
  nullN.list.cp.70.2.sd[[i]]<-unlist(nullN.list.cp.70.2.sd[[i]])
}

# Nestedness of resampled communities
webN.list.cp.70<-vector('list', length(cp.70))
for(i in seq_along(cp.70)){
  webN.list.cp.70[[i]]<-parLapply(cl, cp.70[[i]], function(x) networklevel(x,index="weighted NODF"))
}

webN.list.cp.70<-lapply(webN.list.cp.70, unlist)

# Standardised nestedness of resampled communities
webSN.list.cp.70<-Map("/",Map("-", webN.list.cp.70,nullN.list.cp.70.2.mean),nullN.list.cp.70.2.sd)


# 60% of parasite taxonomic resolution
# Null communities
nullwebs.list.cp.60<-vector('list', length(cp.60))
for(i in seq_along(cp.60)){
  for(j in seq_along(cp.60[[i]])){
    nullwebs.list.cp.60[[i]][[j]]<-swap.web(1000,cp.60[[i]][[j]])
  }
}

# saveRDS(nullwebs.list.cp.60,"nullwebs.list.cp.60.RDS")

# Modularity of null communities
nullM.list.cp.60<-vector('list', length(nullwebs.list.cp.60))
nullM.list.cp.60<-lapply(nullM.list.cp.60, function(x) vector('list', length(nullwebs.list.cp.60[[1]])))
for(i in seq_along(nullM.list.cp.60)){
  for(j in seq_along(nullM.list.cp.60[[i]])){
    nullM.list.cp.60[[i]][[j]]<-vector('list', length(nullwebs.list.cp.60[[1]][[1]]))
  }
}

for(i in seq_along(nullwebs.list.cp.60)){
  for(j in seq_along(nullwebs.list.cp.60[[i]])){
    for(k in seq_along(nullwebs.list.cp.60[[i]][[j]])){
      nullM.list.cp.60[[i]][[j]][[k]]<- computeModules(nullwebs.list.cp.60[[i]][[j]][[k]])@likelihood
    }
  }
}

for(i in seq_along(nullM.list.cp.60)){
  for(j in seq_along(nullM.list.cp.60[[i]])){
    nullM.list.cp.60[[i]][[j]]<-unlist(nullM.list.cp.60[[i]][[j]])
  }
}
# saveRDS(nullM.list.cp.60,"nullM.list.cp.60.RDS")

nullM.list.cp.60.2<-vector('list', length(nullM.list.cp.60))
nullM.list.cp.60.2.mean<-vector('list', length(nullM.list.cp.60))
nullM.list.cp.60.2.sd<-vector('list', length(nullM.list.cp.60))
for(i in seq_along(nullM.list.cp.60)){
  nullM.list.cp.60.2[[i]]<-parLapply(cl, nullM.list.cp.60[[i]], function(x) unlist(x))
  nullM.list.cp.60.2.mean[[i]]<-parLapply(cl, nullM.list.cp.60.2[[i]],function(k) mean(k))
  nullM.list.cp.60.2.mean[[i]]<-unlist(nullM.list.cp.60.2.mean[[i]])
  nullM.list.cp.60.2.sd[[i]]<-parLapply(cl, nullM.list.cp.60.2[[i]],function(k) sd(k))
  nullM.list.cp.60.2.sd[[i]]<-unlist(nullM.list.cp.60.2.sd[[i]])
}

# Modularity of resampled communities
webM.list.cp.60<-vector('list', length(cp.60))
for(i in seq_along(cp.60)){
  webM.list.cp.60[[i]]<-parLapply(cl, cp.60[[i]], function(k) computeModules(k))
}

webM.like.list.cp.60<-vector('list', length(webM.list.cp.60))
for(i in seq_along(webM.list.cp.60)){
  for(j in seq_along(webM.list.cp.60[[i]])){
    webM.like.list.cp.60[[i]]<-parLapply(cl, webM.list.cp.60[[i]], function(k) attributes(k)$likelihood)
  }
}
webM.like.list.cp.60<-lapply(webM.like.list.cp.60, unlist)
webM.like.list.cp.60<-lapply(webM.like.list.cp.60, as.numeric)

# Standardised modulariy of resampled communities
webSM.like.list.cp.60<-Map("/",Map("-", webM.like.list.cp.60,nullM.list.cp.60.2.mean),nullM.list.cp.60.2.sd)


# Nestedness 
# Nestedness of null communities
nullN.list.cp.60<-vector('list', length(nullwebs.list.cp.60))
nullN.list.cp.60<-lapply(nullN.list.cp.60, function(x) vector('list', length(nullwebs.list.cp.60[[1]])))
for(i in seq_along(nullN.list.cp.60)){
  for(j in seq_along(nullN.list.cp.60[[i]])){
    nullN.list.cp.60[[i]][[j]]<-vector('list', length(nullwebs.list.cp.60[[1]][[1]]))
  }
}

for(i in seq_along(nullwebs.list.cp.60)){
  for(j in seq_along(nullwebs.list.cp.60[[i]])){
    for(k in seq_along(nullwebs.list.cp.60[[i]][[j]])){
      nullN.list.cp.60[[i]][[j]][[k]] <- networklevel(nullwebs.list.cp.60[[i]][[j]][[k]],index="weighted NODF")
    }
  }
}

for(i in seq_along(nullN.list.cp.60)){
  for(j in seq_along(nullN.list.cp.60[[i]])){
    nullN.list.cp.60[[i]][[j]]<-as.numeric(unlist(nullN.list.cp.60[[i]][[j]]))
  }
}

# saveRDS(nullN.list.cp.60,"nullN.list.cp.60.RDS")

nullN.list.cp.60.2<-vector('list', length(nullN.list.cp.60))
nullN.list.cp.60.2.mean<-vector('list', length(nullN.list.cp.60))
nullN.list.cp.60.2.sd<-vector('list', length(nullN.list.cp.60))
for(i in seq_along(nullN.list.cp.60)){
  nullN.list.cp.60.2[[i]]<-parLapply(cl, nullN.list.cp.60[[i]], function(x) unlist(x))
  nullN.list.cp.60.2.mean[[i]]<-parLapply(cl, nullN.list.cp.60.2[[i]],function(k) mean(k))
  nullN.list.cp.60.2.mean[[i]]<-unlist(nullN.list.cp.60.2.mean[[i]])
  nullN.list.cp.60.2.sd[[i]]<-parLapply(cl, nullN.list.cp.60.2[[i]],function(k) sd(k))
  nullN.list.cp.60.2.sd[[i]]<-unlist(nullN.list.cp.60.2.sd[[i]])
}

# Nestedness of resampled communities
webN.list.cp.60<-vector('list', length(cp.60))
for(i in seq_along(cp.60)){
  webN.list.cp.60[[i]]<-parLapply(cl, cp.60[[i]], function(x) networklevel(x,index="weighted NODF"))
}

webN.list.cp.60<-lapply(webN.list.cp.60, unlist)

# Standardised nestedness of resampled communities
webSN.list.cp.60<-Map("/",Map("-", webN.list.cp.60,nullN.list.cp.60.2.mean),nullN.list.cp.60.2.sd)


# 50% of parasite taxonomic resolution
# Null communities
nullwebs.list.cp.50<-vector('list', length(cp.50))
for(i in seq_along(cp.50)){
  for(j in seq_along(cp.50[[i]])){
    nullwebs.list.cp.50[[i]][[j]]<-swap.web(1000,cp.50[[i]][[j]])
  }
}

# saveRDS(nullwebs.list.cp.50,"nullwebs.list.cp.50.RDS")

# Modularity of null communities
nullM.list.cp.50<-vector('list', length(nullwebs.list.cp.50))
nullM.list.cp.50<-lapply(nullM.list.cp.50, function(x) vector('list', length(nullwebs.list.cp.50[[1]])))
for(i in seq_along(nullM.list.cp.50)){
  for(j in seq_along(nullM.list.cp.50[[i]])){
    nullM.list.cp.50[[i]][[j]]<-vector('list', length(nullwebs.list.cp.50[[1]][[1]]))
  }
}

for(i in seq_along(nullwebs.list.cp.50)){
  for(j in seq_along(nullwebs.list.cp.50[[i]])){
    for(k in seq_along(nullwebs.list.cp.50[[i]][[j]])){
      nullM.list.cp.50[[i]][[j]][[k]]<- computeModules(nullwebs.list.cp.50[[i]][[j]][[k]])@likelihood
    }
  }
}

for(i in seq_along(nullM.list.cp.50)){
  for(j in seq_along(nullM.list.cp.50[[i]])){
    nullM.list.cp.50[[i]][[j]]<-unlist(nullM.list.cp.50[[i]][[j]])
  }
}
# saveRDS(nullM.list.cp.50,"nullM.list.cp.50.RDS")

nullM.list.cp.50.2<-vector('list', length(nullM.list.cp.50))
nullM.list.cp.50.2.mean<-vector('list', length(nullM.list.cp.50))
nullM.list.cp.50.2.sd<-vector('list', length(nullM.list.cp.50))
for(i in seq_along(nullM.list.cp.50)){
  nullM.list.cp.50.2[[i]]<-parLapply(cl, nullM.list.cp.50[[i]], function(x) unlist(x))
  nullM.list.cp.50.2.mean[[i]]<-parLapply(cl, nullM.list.cp.50.2[[i]],function(k) mean(k))
  nullM.list.cp.50.2.mean[[i]]<-unlist(nullM.list.cp.50.2.mean[[i]])
  nullM.list.cp.50.2.sd[[i]]<-parLapply(cl, nullM.list.cp.50.2[[i]],function(k) sd(k))
  nullM.list.cp.50.2.sd[[i]]<-unlist(nullM.list.cp.50.2.sd[[i]])
}

# Modularity of resampled communities
webM.list.cp.50<-vector('list', length(cp.50))
for(i in seq_along(cp.50)){
  webM.list.cp.50[[i]]<-parLapply(cl, cp.50[[i]], function(k) computeModules(k))
}

webM.like.list.cp.50<-vector('list', length(webM.list.cp.50))
for(i in seq_along(webM.list.cp.50)){
  for(j in seq_along(webM.list.cp.50[[i]])){
    webM.like.list.cp.50[[i]]<-parLapply(cl, webM.list.cp.50[[i]], function(k) attributes(k)$likelihood)
  }
}
webM.like.list.cp.50<-lapply(webM.like.list.cp.50, unlist)
webM.like.list.cp.50<-lapply(webM.like.list.cp.50, as.numeric)

# Standardised modulariy of resampled communities
webSM.like.list.cp.50<-Map("/",Map("-", webM.like.list.cp.50,nullM.list.cp.50.2.mean),nullM.list.cp.50.2.sd)


# Nestedness 
# Nestedness of null communities
nullN.list.cp.50<-vector('list', length(nullwebs.list.cp.50))
nullN.list.cp.50<-lapply(nullN.list.cp.50, function(x) vector('list', length(nullwebs.list.cp.50[[1]])))
for(i in seq_along(nullN.list.cp.50)){
  for(j in seq_along(nullN.list.cp.50[[i]])){
    nullN.list.cp.50[[i]][[j]]<-vector('list', length(nullwebs.list.cp.50[[1]][[1]]))
  }
}

for(i in seq_along(nullwebs.list.cp.50)){
  for(j in seq_along(nullwebs.list.cp.50[[i]])){
    for(k in seq_along(nullwebs.list.cp.50[[i]][[j]])){
      nullN.list.cp.50[[i]][[j]][[k]] <- networklevel(nullwebs.list.cp.50[[i]][[j]][[k]],index="weighted NODF")
    }
  }
}

for(i in seq_along(nullN.list.cp.50)){
  for(j in seq_along(nullN.list.cp.50[[i]])){
    nullN.list.cp.50[[i]][[j]]<-as.numeric(unlist(nullN.list.cp.50[[i]][[j]]))
  }
}

# saveRDS(nullN.list.cp.50,"nullN.list.cp.50.RDS")

nullN.list.cp.50.2<-vector('list', length(nullN.list.cp.50))
nullN.list.cp.50.2.mean<-vector('list', length(nullN.list.cp.50))
nullN.list.cp.50.2.sd<-vector('list', length(nullN.list.cp.50))
for(i in seq_along(nullN.list.cp.50)){
  nullN.list.cp.50.2[[i]]<-parLapply(cl, nullN.list.cp.50[[i]], function(x) unlist(x))
  nullN.list.cp.50.2.mean[[i]]<-parLapply(cl, nullN.list.cp.50.2[[i]],function(k) mean(k))
  nullN.list.cp.50.2.mean[[i]]<-unlist(nullN.list.cp.50.2.mean[[i]])
  nullN.list.cp.50.2.sd[[i]]<-parLapply(cl, nullN.list.cp.50.2[[i]],function(k) sd(k))
  nullN.list.cp.50.2.sd[[i]]<-unlist(nullN.list.cp.50.2.sd[[i]])
}

# Nestedness of resampled communities
webN.list.cp.50<-vector('list', length(cp.50))
for(i in seq_along(cp.50)){
  webN.list.cp.50[[i]]<-parLapply(cl, cp.50[[i]], function(x) networklevel(x,index="weighted NODF"))
}

webN.list.cp.50<-lapply(webN.list.cp.50, unlist)

# Standardised nestedness of resampled communities
webSN.list.cp.50<-Map("/",Map("-", webN.list.cp.50,nullN.list.cp.50.2.mean),nullN.list.cp.50.2.sd)


# 40% of parasite taxonomic resolution
# Null communities
nullwebs.list.cp.40<-vector('list', length(cp.40))
for(i in seq_along(cp.40)){
  for(j in seq_along(cp.40[[i]])){
    nullwebs.list.cp.40[[i]][[j]]<-swap.web(1000,cp.40[[i]][[j]])
  }
}

# saveRDS(nullwebs.list.cp.40,"nullwebs.list.cp.40.RDS")

# Modularity of null communities
nullM.list.cp.40<-vector('list', length(nullwebs.list.cp.40))
nullM.list.cp.40<-lapply(nullM.list.cp.40, function(x) vector('list', length(nullwebs.list.cp.40[[1]])))
for(i in seq_along(nullM.list.cp.40)){
  for(j in seq_along(nullM.list.cp.40[[i]])){
    nullM.list.cp.40[[i]][[j]]<-vector('list', length(nullwebs.list.cp.40[[1]][[1]]))
  }
}

for(i in seq_along(nullwebs.list.cp.40)){
  for(j in seq_along(nullwebs.list.cp.40[[i]])){
    for(k in seq_along(nullwebs.list.cp.40[[i]][[j]])){
      nullM.list.cp.40[[i]][[j]][[k]]<- computeModules(nullwebs.list.cp.40[[i]][[j]][[k]])@likelihood
    }
  }
}

for(i in seq_along(nullM.list.cp.40)){
  for(j in seq_along(nullM.list.cp.40[[i]])){
    nullM.list.cp.40[[i]][[j]]<-unlist(nullM.list.cp.40[[i]][[j]])
  }
}
# saveRDS(nullM.list.cp.40,"nullM.list.cp.40.RDS")

nullM.list.cp.40.2<-vector('list', length(nullM.list.cp.40))
nullM.list.cp.40.2.mean<-vector('list', length(nullM.list.cp.40))
nullM.list.cp.40.2.sd<-vector('list', length(nullM.list.cp.40))
for(i in seq_along(nullM.list.cp.40)){
  nullM.list.cp.40.2[[i]]<-parLapply(cl, nullM.list.cp.40[[i]], function(x) unlist(x))
  nullM.list.cp.40.2.mean[[i]]<-parLapply(cl, nullM.list.cp.40.2[[i]],function(k) mean(k))
  nullM.list.cp.40.2.mean[[i]]<-unlist(nullM.list.cp.40.2.mean[[i]])
  nullM.list.cp.40.2.sd[[i]]<-parLapply(cl, nullM.list.cp.40.2[[i]],function(k) sd(k))
  nullM.list.cp.40.2.sd[[i]]<-unlist(nullM.list.cp.40.2.sd[[i]])
}

# Modularity of resampled communities
webM.list.cp.40<-vector('list', length(cp.40))
for(i in seq_along(cp.40)){
  webM.list.cp.40[[i]]<-parLapply(cl, cp.40[[i]], function(k) computeModules(k))
}

webM.like.list.cp.40<-vector('list', length(webM.list.cp.40))
for(i in seq_along(webM.list.cp.40)){
  for(j in seq_along(webM.list.cp.40[[i]])){
    webM.like.list.cp.40[[i]]<-parLapply(cl, webM.list.cp.40[[i]], function(k) attributes(k)$likelihood)
  }
}
webM.like.list.cp.40<-lapply(webM.like.list.cp.40, unlist)
webM.like.list.cp.40<-lapply(webM.like.list.cp.40, as.numeric)

# Standardised modulariy of resampled communities
webSM.like.list.cp.40<-Map("/",Map("-", webM.like.list.cp.40,nullM.list.cp.40.2.mean),nullM.list.cp.40.2.sd)


# Nestedness 
# Nestedness of null communities
nullN.list.cp.40<-vector('list', length(nullwebs.list.cp.40))
nullN.list.cp.40<-lapply(nullN.list.cp.40, function(x) vector('list', length(nullwebs.list.cp.40[[1]])))
for(i in seq_along(nullN.list.cp.40)){
  for(j in seq_along(nullN.list.cp.40[[i]])){
    nullN.list.cp.40[[i]][[j]]<-vector('list', length(nullwebs.list.cp.40[[1]][[1]]))
  }
}

for(i in seq_along(nullwebs.list.cp.40)){
  for(j in seq_along(nullwebs.list.cp.40[[i]])){
    for(k in seq_along(nullwebs.list.cp.40[[i]][[j]])){
      nullN.list.cp.40[[i]][[j]][[k]] <- networklevel(nullwebs.list.cp.40[[i]][[j]][[k]],index="weighted NODF")
    }
  }
}

for(i in seq_along(nullN.list.cp.40)){
  for(j in seq_along(nullN.list.cp.40[[i]])){
    nullN.list.cp.40[[i]][[j]]<-as.numeric(unlist(nullN.list.cp.40[[i]][[j]]))
  }
}

# saveRDS(nullN.list.cp.40,"nullN.list.cp.40.RDS")

nullN.list.cp.40.2<-vector('list', length(nullN.list.cp.40))
nullN.list.cp.40.2.mean<-vector('list', length(nullN.list.cp.40))
nullN.list.cp.40.2.sd<-vector('list', length(nullN.list.cp.40))
for(i in seq_along(nullN.list.cp.40)){
  nullN.list.cp.40.2[[i]]<-parLapply(cl, nullN.list.cp.40[[i]], function(x) unlist(x))
  nullN.list.cp.40.2.mean[[i]]<-parLapply(cl, nullN.list.cp.40.2[[i]],function(k) mean(k))
  nullN.list.cp.40.2.mean[[i]]<-unlist(nullN.list.cp.40.2.mean[[i]])
  nullN.list.cp.40.2.sd[[i]]<-parLapply(cl, nullN.list.cp.40.2[[i]],function(k) sd(k))
  nullN.list.cp.40.2.sd[[i]]<-unlist(nullN.list.cp.40.2.sd[[i]])
}

# Nestedness of resampled communities
webN.list.cp.40<-vector('list', length(cp.40))
for(i in seq_along(cp.40)){
  webN.list.cp.40[[i]]<-parLapply(cl, cp.40[[i]], function(x) networklevel(x,index="weighted NODF"))
}

webN.list.cp.40<-lapply(webN.list.cp.40, unlist)

# Standardised nestedness of resampled communities
webSN.list.cp.40<-Map("/",Map("-", webN.list.cp.40,nullN.list.cp.40.2.mean),nullN.list.cp.40.2.sd)


# 30% of parasite taxonomic resolution
# Null communities
nullwebs.list.cp.30<-vector('list', length(cp.30))
for(i in seq_along(cp.30)){
  for(j in seq_along(cp.30[[i]])){
    nullwebs.list.cp.30[[i]][[j]]<-swap.web(1000,cp.30[[i]][[j]])
  }
}

# saveRDS(nullwebs.list.cp.30,"nullwebs.list.cp.30.RDS")

# Modularity of null communities
nullM.list.cp.30<-vector('list', length(nullwebs.list.cp.30))
nullM.list.cp.30<-lapply(nullM.list.cp.30, function(x) vector('list', length(nullwebs.list.cp.30[[1]])))
for(i in seq_along(nullM.list.cp.30)){
  for(j in seq_along(nullM.list.cp.30[[i]])){
    nullM.list.cp.30[[i]][[j]]<-vector('list', length(nullwebs.list.cp.30[[1]][[1]]))
  }
}

for(i in seq_along(nullwebs.list.cp.30)){
  for(j in seq_along(nullwebs.list.cp.30[[i]])){
    for(k in seq_along(nullwebs.list.cp.30[[i]][[j]])){
      nullM.list.cp.30[[i]][[j]][[k]]<- computeModules(nullwebs.list.cp.30[[i]][[j]][[k]])@likelihood
    }
  }
}

for(i in seq_along(nullM.list.cp.30)){
  for(j in seq_along(nullM.list.cp.30[[i]])){
    nullM.list.cp.30[[i]][[j]]<-unlist(nullM.list.cp.30[[i]][[j]])
  }
}
# saveRDS(nullM.list.cp.30,"nullM.list.cp.30.RDS")

nullM.list.cp.30.2<-vector('list', length(nullM.list.cp.30))
nullM.list.cp.30.2.mean<-vector('list', length(nullM.list.cp.30))
nullM.list.cp.30.2.sd<-vector('list', length(nullM.list.cp.30))
for(i in seq_along(nullM.list.cp.30)){
  nullM.list.cp.30.2[[i]]<-parLapply(cl, nullM.list.cp.30[[i]], function(x) unlist(x))
  nullM.list.cp.30.2.mean[[i]]<-parLapply(cl, nullM.list.cp.30.2[[i]],function(k) mean(k))
  nullM.list.cp.30.2.mean[[i]]<-unlist(nullM.list.cp.30.2.mean[[i]])
  nullM.list.cp.30.2.sd[[i]]<-parLapply(cl, nullM.list.cp.30.2[[i]],function(k) sd(k))
  nullM.list.cp.30.2.sd[[i]]<-unlist(nullM.list.cp.30.2.sd[[i]])
}

# Modularity of resampled communities
webM.list.cp.30<-vector('list', length(cp.30))
for(i in seq_along(cp.30)){
  webM.list.cp.30[[i]]<-parLapply(cl, cp.30[[i]], function(k) computeModules(k))
}

webM.like.list.cp.30<-vector('list', length(webM.list.cp.30))
for(i in seq_along(webM.list.cp.30)){
  for(j in seq_along(webM.list.cp.30[[i]])){
    webM.like.list.cp.30[[i]]<-parLapply(cl, webM.list.cp.30[[i]], function(k) attributes(k)$likelihood)
  }
}
webM.like.list.cp.30<-lapply(webM.like.list.cp.30, unlist)
webM.like.list.cp.30<-lapply(webM.like.list.cp.30, as.numeric)

# Standardised modulariy of resampled communities
webSM.like.list.cp.30<-Map("/",Map("-", webM.like.list.cp.30,nullM.list.cp.30.2.mean),nullM.list.cp.30.2.sd)


# Nestedness 
# Nestedness of null communities
nullN.list.cp.30<-vector('list', length(nullwebs.list.cp.30))
nullN.list.cp.30<-lapply(nullN.list.cp.30, function(x) vector('list', length(nullwebs.list.cp.30[[1]])))
for(i in seq_along(nullN.list.cp.30)){
  for(j in seq_along(nullN.list.cp.30[[i]])){
    nullN.list.cp.30[[i]][[j]]<-vector('list', length(nullwebs.list.cp.30[[1]][[1]]))
  }
}

for(i in seq_along(nullwebs.list.cp.30)){
  for(j in seq_along(nullwebs.list.cp.30[[i]])){
    for(k in seq_along(nullwebs.list.cp.30[[i]][[j]])){
      nullN.list.cp.30[[i]][[j]][[k]] <- networklevel(nullwebs.list.cp.30[[i]][[j]][[k]],index="weighted NODF")
    }
  }
}

for(i in seq_along(nullN.list.cp.30)){
  for(j in seq_along(nullN.list.cp.30[[i]])){
    nullN.list.cp.30[[i]][[j]]<-as.numeric(unlist(nullN.list.cp.30[[i]][[j]]))
  }
}

# saveRDS(nullN.list.cp.30,"nullN.list.cp.30.RDS")

nullN.list.cp.30.2<-vector('list', length(nullN.list.cp.30))
nullN.list.cp.30.2.mean<-vector('list', length(nullN.list.cp.30))
nullN.list.cp.30.2.sd<-vector('list', length(nullN.list.cp.30))
for(i in seq_along(nullN.list.cp.30)){
  nullN.list.cp.30.2[[i]]<-parLapply(cl, nullN.list.cp.30[[i]], function(x) unlist(x))
  nullN.list.cp.30.2.mean[[i]]<-parLapply(cl, nullN.list.cp.30.2[[i]],function(k) mean(k))
  nullN.list.cp.30.2.mean[[i]]<-unlist(nullN.list.cp.30.2.mean[[i]])
  nullN.list.cp.30.2.sd[[i]]<-parLapply(cl, nullN.list.cp.30.2[[i]],function(k) sd(k))
  nullN.list.cp.30.2.sd[[i]]<-unlist(nullN.list.cp.30.2.sd[[i]])
}

# Nestedness of resampled communities
webN.list.cp.30<-vector('list', length(cp.30))
for(i in seq_along(cp.30)){
  webN.list.cp.30[[i]]<-parLapply(cl, cp.30[[i]], function(x) networklevel(x,index="weighted NODF"))
}

webN.list.cp.30<-lapply(webN.list.cp.30, unlist)

# Standardised nestedness of resampled communities
webSN.list.cp.30<-Map("/",Map("-", webN.list.cp.30,nullN.list.cp.30.2.mean),nullN.list.cp.30.2.sd)


# 20% of parasite taxonomic resolution
# Null communities
nullwebs.list.cp.20<-vector('list', length(cp.20))
for(i in seq_along(cp.20)){
  for(j in seq_along(cp.20[[i]])){
    nullwebs.list.cp.20[[i]][[j]]<-swap.web(1000,cp.20[[i]][[j]])
  }
}

# saveRDS(nullwebs.list.cp.20,"nullwebs.list.cp.20.RDS")

# Modularity of null communities
nullM.list.cp.20<-vector('list', length(nullwebs.list.cp.20))
nullM.list.cp.20<-lapply(nullM.list.cp.20, function(x) vector('list', length(nullwebs.list.cp.20[[1]])))
for(i in seq_along(nullM.list.cp.20)){
  for(j in seq_along(nullM.list.cp.20[[i]])){
    nullM.list.cp.20[[i]][[j]]<-vector('list', length(nullwebs.list.cp.20[[1]][[1]]))
  }
}

for(i in seq_along(nullwebs.list.cp.20)){
  for(j in seq_along(nullwebs.list.cp.20[[i]])){
    for(k in seq_along(nullwebs.list.cp.20[[i]][[j]])){
      nullM.list.cp.20[[i]][[j]][[k]]<- computeModules(nullwebs.list.cp.20[[i]][[j]][[k]])@likelihood
    }
  }
}

for(i in seq_along(nullM.list.cp.20)){
  for(j in seq_along(nullM.list.cp.20[[i]])){
    nullM.list.cp.20[[i]][[j]]<-unlist(nullM.list.cp.20[[i]][[j]])
  }
}
# saveRDS(nullM.list.cp.20,"nullM.list.cp.20.RDS")

nullM.list.cp.20.2<-vector('list', length(nullM.list.cp.20))
nullM.list.cp.20.2.mean<-vector('list', length(nullM.list.cp.20))
nullM.list.cp.20.2.sd<-vector('list', length(nullM.list.cp.20))
for(i in seq_along(nullM.list.cp.20)){
  nullM.list.cp.20.2[[i]]<-parLapply(cl, nullM.list.cp.20[[i]], function(x) unlist(x))
  nullM.list.cp.20.2.mean[[i]]<-parLapply(cl, nullM.list.cp.20.2[[i]],function(k) mean(k))
  nullM.list.cp.20.2.mean[[i]]<-unlist(nullM.list.cp.20.2.mean[[i]])
  nullM.list.cp.20.2.sd[[i]]<-parLapply(cl, nullM.list.cp.20.2[[i]],function(k) sd(k))
  nullM.list.cp.20.2.sd[[i]]<-unlist(nullM.list.cp.20.2.sd[[i]])
}

# Modularity of resampled communities
webM.list.cp.20<-vector('list', length(cp.20))
for(i in seq_along(cp.20)){
  webM.list.cp.20[[i]]<-parLapply(cl, cp.20[[i]], function(k) computeModules(k))
}

webM.like.list.cp.20<-vector('list', length(webM.list.cp.20))
for(i in seq_along(webM.list.cp.20)){
  for(j in seq_along(webM.list.cp.20[[i]])){
    webM.like.list.cp.20[[i]]<-parLapply(cl, webM.list.cp.20[[i]], function(k) attributes(k)$likelihood)
  }
}
webM.like.list.cp.20<-lapply(webM.like.list.cp.20, unlist)
webM.like.list.cp.20<-lapply(webM.like.list.cp.20, as.numeric)

# Standardised modulariy of resampled communities
webSM.like.list.cp.20<-Map("/",Map("-", webM.like.list.cp.20,nullM.list.cp.20.2.mean),nullM.list.cp.20.2.sd)


# Nestedness 
# Nestedness of null communities
nullN.list.cp.20<-vector('list', length(nullwebs.list.cp.20))
nullN.list.cp.20<-lapply(nullN.list.cp.20, function(x) vector('list', length(nullwebs.list.cp.20[[1]])))
for(i in seq_along(nullN.list.cp.20)){
  for(j in seq_along(nullN.list.cp.20[[i]])){
    nullN.list.cp.20[[i]][[j]]<-vector('list', length(nullwebs.list.cp.20[[1]][[1]]))
  }
}

for(i in seq_along(nullwebs.list.cp.20)){
  for(j in seq_along(nullwebs.list.cp.20[[i]])){
    for(k in seq_along(nullwebs.list.cp.20[[i]][[j]])){
      nullN.list.cp.20[[i]][[j]][[k]] <- networklevel(nullwebs.list.cp.20[[i]][[j]][[k]],index="weighted NODF")
    }
  }
}

for(i in seq_along(nullN.list.cp.20)){
  for(j in seq_along(nullN.list.cp.20[[i]])){
    nullN.list.cp.20[[i]][[j]]<-as.numeric(unlist(nullN.list.cp.20[[i]][[j]]))
  }
}

# saveRDS(nullN.list.cp.20,"nullN.list.cp.20.RDS")

nullN.list.cp.20.2<-vector('list', length(nullN.list.cp.20))
nullN.list.cp.20.2.mean<-vector('list', length(nullN.list.cp.20))
nullN.list.cp.20.2.sd<-vector('list', length(nullN.list.cp.20))
for(i in seq_along(nullN.list.cp.20)){
  nullN.list.cp.20.2[[i]]<-parLapply(cl, nullN.list.cp.20[[i]], function(x) unlist(x))
  nullN.list.cp.20.2.mean[[i]]<-parLapply(cl, nullN.list.cp.20.2[[i]],function(k) mean(k))
  nullN.list.cp.20.2.mean[[i]]<-unlist(nullN.list.cp.20.2.mean[[i]])
  nullN.list.cp.20.2.sd[[i]]<-parLapply(cl, nullN.list.cp.20.2[[i]],function(k) sd(k))
  nullN.list.cp.20.2.sd[[i]]<-unlist(nullN.list.cp.20.2.sd[[i]])
}

# Nestedness of resampled communities
webN.list.cp.20<-vector('list', length(cp.20))
for(i in seq_along(cp.20)){
  webN.list.cp.20[[i]]<-parLapply(cl, cp.20[[i]], function(x) networklevel(x,index="weighted NODF"))
}

webN.list.cp.20<-lapply(webN.list.cp.20, unlist)

# Standardised nestedness of resampled communities
webSN.list.cp.20<-Map("/",Map("-", webN.list.cp.20,nullN.list.cp.20.2.mean),nullN.list.cp.20.2.sd)


# 10% of parasite taxonomic resolution
# Null communities
nullwebs.list.cp.10<-vector('list', length(cp.10))
for(i in seq_along(cp.10)){
  for(j in seq_along(cp.10[[i]])){
    nullwebs.list.cp.10[[i]][[j]]<-swap.web(1000,cp.10[[i]][[j]])
  }
}

# saveRDS(nullwebs.list.cp.10,"nullwebs.list.cp.10.RDS")

# Modularity of null communities
nullM.list.cp.10<-vector('list', length(nullwebs.list.cp.10))
nullM.list.cp.10<-lapply(nullM.list.cp.10, function(x) vector('list', length(nullwebs.list.cp.10[[1]])))
for(i in seq_along(nullM.list.cp.10)){
  for(j in seq_along(nullM.list.cp.10[[i]])){
    nullM.list.cp.10[[i]][[j]]<-vector('list', length(nullwebs.list.cp.10[[1]][[1]]))
  }
}

for(i in seq_along(nullwebs.list.cp.10)){
  for(j in seq_along(nullwebs.list.cp.10[[i]])){
    for(k in seq_along(nullwebs.list.cp.10[[i]][[j]])){
      nullM.list.cp.10[[i]][[j]][[k]]<- computeModules(nullwebs.list.cp.10[[i]][[j]][[k]])@likelihood
    }
  }
}

for(i in seq_along(nullM.list.cp.10)){
  for(j in seq_along(nullM.list.cp.10[[i]])){
    nullM.list.cp.10[[i]][[j]]<-unlist(nullM.list.cp.10[[i]][[j]])
  }
}
# saveRDS(nullM.list.cp.10,"nullM.list.cp.10.RDS")

nullM.list.cp.10.2<-vector('list', length(nullM.list.cp.10))
nullM.list.cp.10.2.mean<-vector('list', length(nullM.list.cp.10))
nullM.list.cp.10.2.sd<-vector('list', length(nullM.list.cp.10))
for(i in seq_along(nullM.list.cp.10)){
  nullM.list.cp.10.2[[i]]<-parLapply(cl, nullM.list.cp.10[[i]], function(x) unlist(x))
  nullM.list.cp.10.2.mean[[i]]<-parLapply(cl, nullM.list.cp.10.2[[i]],function(k) mean(k))
  nullM.list.cp.10.2.mean[[i]]<-unlist(nullM.list.cp.10.2.mean[[i]])
  nullM.list.cp.10.2.sd[[i]]<-parLapply(cl, nullM.list.cp.10.2[[i]],function(k) sd(k))
  nullM.list.cp.10.2.sd[[i]]<-unlist(nullM.list.cp.10.2.sd[[i]])
}

# Modularity of resampled communities
webM.list.cp.10<-vector('list', length(cp.10))
for(i in seq_along(cp.10)){
  webM.list.cp.10[[i]]<-parLapply(cl, cp.10[[i]], function(k) computeModules(k))
}

webM.like.list.cp.10<-vector('list', length(webM.list.cp.10))
for(i in seq_along(webM.list.cp.10)){
  for(j in seq_along(webM.list.cp.10[[i]])){
    webM.like.list.cp.10[[i]]<-parLapply(cl, webM.list.cp.10[[i]], function(k) attributes(k)$likelihood)
  }
}
webM.like.list.cp.10<-lapply(webM.like.list.cp.10, unlist)
webM.like.list.cp.10<-lapply(webM.like.list.cp.10, as.numeric)

# Standardised modulariy of resampled communities
webSM.like.list.cp.10<-Map("/",Map("-", webM.like.list.cp.10,nullM.list.cp.10.2.mean),nullM.list.cp.10.2.sd)


# Nestedness 
# Nestedness of null communities
nullN.list.cp.10<-vector('list', length(nullwebs.list.cp.10))
nullN.list.cp.10<-lapply(nullN.list.cp.10, function(x) vector('list', length(nullwebs.list.cp.10[[1]])))
for(i in seq_along(nullN.list.cp.10)){
  for(j in seq_along(nullN.list.cp.10[[i]])){
    nullN.list.cp.10[[i]][[j]]<-vector('list', length(nullwebs.list.cp.10[[1]][[1]]))
  }
}

for(i in seq_along(nullwebs.list.cp.10)){
  for(j in seq_along(nullwebs.list.cp.10[[i]])){
    for(k in seq_along(nullwebs.list.cp.10[[i]][[j]])){
      nullN.list.cp.10[[i]][[j]][[k]] <- networklevel(nullwebs.list.cp.10[[i]][[j]][[k]],index="weighted NODF")
    }
  }
}

for(i in seq_along(nullN.list.cp.10)){
  for(j in seq_along(nullN.list.cp.10[[i]])){
    nullN.list.cp.10[[i]][[j]]<-as.numeric(unlist(nullN.list.cp.10[[i]][[j]]))
  }
}

# saveRDS(nullN.list.cp.10,"nullN.list.cp.10.RDS")

nullN.list.cp.10.2<-vector('list', length(nullN.list.cp.10))
nullN.list.cp.10.2.mean<-vector('list', length(nullN.list.cp.10))
nullN.list.cp.10.2.sd<-vector('list', length(nullN.list.cp.10))
for(i in seq_along(nullN.list.cp.10)){
  nullN.list.cp.10.2[[i]]<-parLapply(cl, nullN.list.cp.10[[i]], function(x) unlist(x))
  nullN.list.cp.10.2.mean[[i]]<-parLapply(cl, nullN.list.cp.10.2[[i]],function(k) mean(k))
  nullN.list.cp.10.2.mean[[i]]<-unlist(nullN.list.cp.10.2.mean[[i]])
  nullN.list.cp.10.2.sd[[i]]<-parLapply(cl, nullN.list.cp.10.2[[i]],function(k) sd(k))
  nullN.list.cp.10.2.sd[[i]]<-unlist(nullN.list.cp.10.2.sd[[i]])
}

# Nestedness of resampled communities
webN.list.cp.10<-vector('list', length(cp.10))
for(i in seq_along(cp.10)){
  webN.list.cp.10[[i]]<-parLapply(cl, cp.10[[i]], function(x) networklevel(x,index="weighted NODF"))
}

webN.list.cp.10<-lapply(webN.list.cp.10, unlist)

# Standardised nestedness of resampled communities
webSN.list.cp.10<-Map("/",Map("-", webN.list.cp.10,nullN.list.cp.10.2.mean),nullN.list.cp.10.2.sd)

# Stop running R in parallel
stopCluster(cl)

```


#### 3.1.2 Connectance and Specialisation (H~2~^'^)

* Full and Resampled communities with reduced effort in host sampling completeness

I measured connectance and specialisation (H~2~^'^) of the 10 full and the 90 resampled communities with reduced effort in host sampling completeness.

``` r
# Connectance
webC.list<-vector('list', length(myweb.full.comp.list))
for(i in seq_along(myweb.full.comp.list)){
  for(j in seq_along(myweb.full.comp.list[[i]])){
    webC.list[[i]][[j]] <- networklevel(myweb.full.comp.list[[i]][[j]],
                                        c("weighted connectance"))
    }
}


# Specialisation
webH.list<-vector('list', length(myweb.full.comp.list))
for(i in seq_along(myweb.full.comp.list)){
  for(j in seq_along(myweb.full.comp.list[[i]])){
    webH.list[[i]][[j]] <- networklevel(myweb.full.comp.list[[i]][[j]], c("H2"))
    }
}

```

* Resampled communities with reduced effort in parasite taxonomic resolution

I measured connectance and specialisation (H~2~^'^) of the 90 resampled communities with reduced effort in parasite taxonomic resolution and the 810 resampled communities with reduced effort in both sampling issues (the crossed effect).

``` r
# 90% of parasite taxonomic resolution
# Connectance
webC.list.cp.90<-vector('list', length(cp.90))
for(i in seq_along(cp.90)){
  for(j in seq_along(cp.90[[i]])){
    webC.list.cp.90[[i]][[j]] <- networklevel(cp.90[[i]][[j]],
                                              c("weighted connectance"))
    webC.list.cp.90[[i]]<-unlist(webC.list.cp.90[[i]])
  }
}


# Specialisation
webH.list.cp.90<-vector('list', length(cp.90))
for(i in seq_along(cp.90)){
  for(j in seq_along(cp.90[[i]])){
    webH.list.cp.90[[i]][[j]] <- networklevel(cp.90[[i]][[j]],
                                        c("H2"))
    webH.list.cp.90[[i]]<-unlist(webH.list.cp.90[[i]])
  }
}


# 80% of parasite taxonomic resolution
# Connectance
webC.list.cp.80<-vector('list', length(cp.80))
for(i in seq_along(cp.80)){
  for(j in seq_along(cp.80[[i]])){
    webC.list.cp.80[[i]][[j]] <- networklevel(cp.80[[i]][[j]],
                                              c("weighted connectance"))
    webC.list.cp.80[[i]]<-unlist(webC.list.cp.80[[i]])
  }
}


# Specialisation
webH.list.cp.80<-vector('list', length(cp.80))
for(i in seq_along(cp.80)){
  for(j in seq_along(cp.80[[i]])){
    webH.list.cp.80[[i]][[j]] <- networklevel(cp.80[[i]][[j]],
                                              c("H2"))
    webH.list.cp.80[[i]]<-unlist(webH.list.cp.80[[i]])
  }
}


# 70% of parasite taxonomic resolution
# Connectance
webC.list.cp.70<-vector('list', length(cp.70))
for(i in seq_along(cp.70)){
  for(j in seq_along(cp.70[[i]])){
    webC.list.cp.70[[i]][[j]] <- networklevel(cp.70[[i]][[j]],
                                              c("weighted connectance"))
    webC.list.cp.70[[i]]<-unlist(webC.list.cp.70[[i]])
  }
}


# Specialisation
webH.list.cp.70<-vector('list', length(cp.70))
for(i in seq_along(cp.70)){
  for(j in seq_along(cp.70[[i]])){
    webH.list.cp.70[[i]][[j]] <- networklevel(cp.70[[i]][[j]],
                                              c("H2"))
    webH.list.cp.70[[i]]<-unlist(webH.list.cp.70[[i]])
  }
}


# 60% of parasite taxonomic resolution
# Connectance
webC.list.cp.60<-vector('list', length(cp.60))
for(i in seq_along(cp.60)){
  for(j in seq_along(cp.60[[i]])){
    webC.list.cp.60[[i]][[j]] <- networklevel(cp.60[[i]][[j]],
                                              c("weighted connectance"))
    webC.list.cp.60[[i]]<-unlist(webC.list.cp.60[[i]])
  }
}


# Specialisation
webH.list.cp.60<-vector('list', length(cp.60))
for(i in seq_along(cp.60)){
  for(j in seq_along(cp.60[[i]])){
    webH.list.cp.60[[i]][[j]] <- networklevel(cp.60[[i]][[j]],
                                              c("H2"))
    webH.list.cp.60[[i]]<-unlist(webH.list.cp.60[[i]])
  }
}


# 50% of parasite taxonomic resolution
# Connectance
webC.list.cp.50<-vector('list', length(cp.50))
for(i in seq_along(cp.50)){
  for(j in seq_along(cp.50[[i]])){
    webC.list.cp.50[[i]][[j]] <- networklevel(cp.50[[i]][[j]],
                                              c("weighted connectance"))
    webC.list.cp.50[[i]]<-unlist(webC.list.cp.50[[i]])
  }
}


# Specialisation
webH.list.cp.50<-vector('list', length(cp.50))
for(i in seq_along(cp.50)){
  for(j in seq_along(cp.50[[i]])){
    webH.list.cp.50[[i]][[j]] <- networklevel(cp.50[[i]][[j]],
                                              c("H2"))
    webH.list.cp.50[[i]]<-unlist(webH.list.cp.50[[i]])
  }
}


# 40% of parasite taxonomic resolution
# Connectance
webC.list.cp.40<-vector('list', length(cp.40))
for(i in seq_along(cp.40)){
  for(j in seq_along(cp.40[[i]])){
    webC.list.cp.40[[i]][[j]] <- networklevel(cp.40[[i]][[j]],
                                              c("weighted connectance"))
    webC.list.cp.40[[i]]<-unlist(webC.list.cp.40[[i]])
  }
}


# Specialisation
webH.list.cp.40<-vector('list', length(cp.40))
for(i in seq_along(cp.40)){
  for(j in seq_along(cp.40[[i]])){
    webH.list.cp.40[[i]][[j]] <- networklevel(cp.40[[i]][[j]],
                                              c("H2"))
    webH.list.cp.40[[i]]<-unlist(webH.list.cp.40[[i]])
  }
}


# 30% of parasite taxonomic resolution
# Connectance
webC.list.cp.30<-vector('list', length(cp.30))
for(i in seq_along(cp.30)){
  for(j in seq_along(cp.30[[i]])){
    webC.list.cp.30[[i]][[j]] <- networklevel(cp.30[[i]][[j]],
                                              c("weighted connectance"))
    webC.list.cp.30[[i]]<-unlist(webC.list.cp.30[[i]])
  }
}


# Specialisation
webH.list.cp.30<-vector('list', length(cp.30))
for(i in seq_along(cp.30)){
  for(j in seq_along(cp.30[[i]])){
    webH.list.cp.30[[i]][[j]] <- networklevel(cp.30[[i]][[j]],
                                              c("H2"))
    webH.list.cp.30[[i]]<-unlist(webH.list.cp.30[[i]])
  }
}


# 20% of parasite taxonomic resolution
# Connectance
webC.list.cp.20<-vector('list', length(cp.20))
for(i in seq_along(cp.20)){
  for(j in seq_along(cp.20[[i]])){
    webC.list.cp.20[[i]][[j]] <- networklevel(cp.20[[i]][[j]],
                                              c("weighted connectance"))
    webC.list.cp.20[[i]]<-unlist(webC.list.cp.20[[i]])
  }
}


# Specialisation
webH.list.cp.20<-vector('list', length(cp.20))
for(i in seq_along(cp.20)){
  for(j in seq_along(cp.20[[i]])){
    webH.list.cp.20[[i]][[j]] <- networklevel(cp.20[[i]][[j]],
                                              c("H2"))
    webH.list.cp.20[[i]]<-unlist(webH.list.cp.20[[i]])
  }
}


# 10% of parasite taxonomic resolution
# Connectance
webC.list.cp.10<-vector('list', length(cp.10))
for(i in seq_along(cp.10)){
  for(j in seq_along(cp.10[[i]])){
    webC.list.cp.10[[i]][[j]] <- networklevel(cp.10[[i]][[j]],
                                              c("weighted connectance"))
    webC.list.cp.10[[i]]<-unlist(webC.list.cp.10[[i]])
  }
}


# Specialisation
webH.list.cp.10<-vector('list', length(cp.10))
for(i in seq_along(cp.10)){
  for(j in seq_along(cp.10[[i]])){
    webH.list.cp.10[[i]][[j]] <- networklevel(cp.10[[i]][[j]],
                                              c("H2"))
    webH.list.cp.10[[i]]<-unlist(webH.list.cp.10[[i]])
  }
}
```


### 3.2 Species descriptors

* Full communities

I assessed the pattern of interactions of the full communities at the species level. I measured alpha diversities of the parasite communities of each host species, and two species level network descriptors of centrality for both host and parasite species independently: weighted betweenness and weighted closeness.

``` r
# Alpha diversity
ad<-lapply(myweb.full.list, data.frame) # ad: alpha diversity
dpcoa.ad<-lapply(ad, function(x) dpcoa(x, scannf = F, full = T))
alpha.ad<-lapply(dpcoa.ad, function (x) 1/(1-x$RaoDiv))


# Betweenness
webB.list<-vector('list', length(myweb.full.comp.list))
for(i in seq_along(myweb.full.comp.list)){
        for(j in seq_along(myweb.full.comp.list[[i]])){
                webB.list[[i]][[j]] <- specieslevel(myweb.full.comp.list[[i]][[j]],
                index = "betweenness")
        }
}

# betweenness parasites
webB.list.p<-vector('list', length(webB.list))
# betweenness hosts
webB.list.h<-vector('list', length(webB.list))

for(i in seq_along(webB.list)){
  for(j in seq_along(webB.list[[i]])){
    webB.list.p[[i]][[j]] <- webB.list[[i]][[j]]$`higher level`[2][,1]
    webB.list.h[[i]][[j]] <- webB.list[[i]][[j]]$`lower level`[2][,1]
  }
}


# Closeness
webE.list<-vector('list', length(myweb.full.comp.list)) # E: cEntrality
for(i in seq_along(myweb.full.comp.list)){
    for(j in seq_along(myweb.full.comp.list[[i]])){
      webE.list[[i]][[j]] <- specieslevel(myweb.full.comp.list[[i]][[j]], index = "closeness")
        }
}

# closeness parasites
webE.list.p<-vector('list', length(webE.list)) 
# closeness hosts
webE.list.h<-vector('list', length(webE.list)) 

for(i in seq_along(webE.list)){
  for(j in seq_along(webE.list[[i]])){
    webE.list.p[[i]][[j]] <- webE.list[[i]][[j]]$`higher level`[2][,1]
    webE.list.h[[i]][[j]] <- webE.list[[i]][[j]]$`lower level`[2][,1]
    }
}
```


## 4. Figure

I first created a data frame with the results to plot the effects of decreasing host sampling completeness and parasite taxonomic resolution on community descriptors: modularity, nestedness, connectance and specialisation.

``` r
results<-data.frame(
  sam=rev(as.character(rep(seq(from=10,to=100, by=10), 10))), # sam: host sampling completeness
  res=as.character(sort(rep(seq(from=10,to=100, by=10), 10), T)), # res: Parasite.Taxonomic.Resolution
  
  M=  # M: Standardised Modularity
    # Means of the standardised modularities of the full and resampled communities with reduced effort in 
    # host sampling complteness
    c(unlist(lapply(list(webSM.list,web90SM.list,web80SM.list,web70SM.list,web60SM.list,
                         web50SM.list,web40SM.list,web30SM.list,web20SM.list,web10SM.list),mean)),
      # Means of the standardised modularities of the resampled communities with reduced effort in 
      # parasite taxonomic resolution
      unlist(lapply(webSM.like.list.cp.90, mean)),
      unlist(lapply(webSM.like.list.cp.80, mean)),unlist(lapply(webSM.like.list.cp.70, mean)),
      unlist(lapply(webSM.like.list.cp.60, mean)),unlist(lapply(webSM.like.list.cp.50, mean)),
      unlist(lapply(webSM.like.list.cp.40, mean)),unlist(lapply(webSM.like.list.cp.30, mean)),
      unlist(lapply(webSM.like.list.cp.20, mean)),unlist(lapply(webSM.like.list.cp.10, mean))
    ),
  
  Mse= # Mse: standardised error of the mean of modularity
    # Standard errors of the means of the standardised modularities of the full and resampled communities with reduced effort in 
    # host sampling complteness
    c(unlist(lapply(list(webSM.list,web90SM.list,web80SM.list,web70SM.list,web60SM.list,web50SM.list,web40SM.list,
                         web30SM.list,web20SM.list,web10SM.list), std.error)),
      # Standard errors of the means  of the standardised modularities of the resampled communities with reduced effort in 
      # parasite taxonomic resolution
      unlist(lapply(webSM.like.list.cp.90, std.error)),
      unlist(lapply(webSM.like.list.cp.80, std.error)),unlist(lapply(webSM.like.list.cp.70, std.error)),
      unlist(lapply(webSM.like.list.cp.60, std.error)),unlist(lapply(webSM.like.list.cp.50, std.error)),
      unlist(lapply(webSM.like.list.cp.40, std.error)),unlist(lapply(webSM.like.list.cp.30, std.error)),
      unlist(lapply(webSM.like.list.cp.20, std.error)),unlist(lapply(webSM.like.list.cp.10, std.error))
    ),
  
  N=  # N: Standardised Nestedness
    # Means of the standardised nestedness of the full and resampled communities with reduced effort in 
    # host sampling complteness
    c(unlist(lapply(list(webSN.list,web90SN.list,web80SN.list,web70SN.list,web60SN.list,
                         web50SN.list,web40SN.list,web30SN.list,web20SN.list,web10SN.list),mean)),
      # Means of the standardised nestedness of the resampled communities with reduced effort in 
      # parasite taxonomic resolution
      unlist(lapply(webSN.list.cp.90, mean)),
      unlist(lapply(webSN.list.cp.80, mean)),unlist(lapply(webSN.list.cp.70, mean)),
      unlist(lapply(webSN.list.cp.60, mean)),unlist(lapply(webSN.list.cp.50, mean)),
      unlist(lapply(webSN.list.cp.40, mean)),unlist(lapply(webSN.list.cp.30, mean)),
      unlist(lapply(webSN.list.cp.20, mean)),unlist(lapply(webSN.list.cp.10, mean))
    ),
  
  Nse= # Nse: standardised error of the mean of nestedness
    # Standard errors of the means of the standardised nestedness of the full and resampled communities with reduced effort in 
    # host sampling complteness  
    c(unlist(lapply(list(webSN.list,web90SN.list,web80SN.list,web70SN.list,web60SN.list,web50SN.list,web40SN.list,
                         web30SN.list,web20SN.list,web10SN.list), std.error)),
      # Standard errors of the means of the standardised nestedness of the resampled communities with reduced effort in 
      # parasite taxonomic resolution
      unlist(lapply(webSN.list.cp.90, std.error)),
      unlist(lapply(webSN.list.cp.80, std.error)),unlist(lapply(webSN.list.cp.70, std.error)),
      unlist(lapply(webSN.list.cp.60, std.error)),unlist(lapply(webSN.list.cp.50, std.error)),
      unlist(lapply(webSN.list.cp.40, std.error)),unlist(lapply(webSN.list.cp.30, std.error)),
      unlist(lapply(webSN.list.cp.20, std.error)),unlist(lapply(webSN.list.cp.10, std.error))
    ),
  
  C= # C: connectance
    # Means of connectance of the full and resampled communities with reduced effort in 
    # host sampling complteness
    c(unlist(lapply(lapply(webC.list, unlist),mean)),
      # Means of the connectances of the resampled communities with reduced effort in 
      # parasite taxonomic resolution
      unlist(lapply(webC.list.cp.90, mean)),
      unlist(lapply(webC.list.cp.80, mean)),unlist(lapply(webC.list.cp.70, mean)),
      unlist(lapply(webC.list.cp.60, mean)),unlist(lapply(webC.list.cp.50, mean)),
      unlist(lapply(webC.list.cp.40, mean)),unlist(lapply(webC.list.cp.30, mean)),
      unlist(lapply(webC.list.cp.20, mean)),unlist(lapply(webC.list.cp.10, mean))
    ),
  
  Cse=  # Cse: standardised error of the mean of connectance
    # Standard errors of the means of the connectances of the full and resampled communities with reduced effort in 
    # host sampling complteness  
    c(unlist(lapply(lapply(webC.list, unlist),std.error)),
      # Standard errors of the means of the connectances of the resampled communities with reduced effort in 
      # parasite taxonomic resolution
      unlist(lapply(webC.list.cp.90, std.error)),
      unlist(lapply(webC.list.cp.80, std.error)),unlist(lapply(webC.list.cp.70, std.error)),
      unlist(lapply(webC.list.cp.60, std.error)),unlist(lapply(webC.list.cp.50, std.error)),
      unlist(lapply(webC.list.cp.40, std.error)),unlist(lapply(webC.list.cp.30, std.error)),
      unlist(lapply(webC.list.cp.20, std.error)),unlist(lapply(webC.list.cp.10, std.error))
    ),
  
  H= # H: specialisation
    # Means of specialisation of the full and resampled communities with reduced effort in 
    # host sampling complteness
    c(unlist(lapply(lapply(webH.list, unlist),mean)),
      # Means of the specialisations of the resampled communities with reduced effort in 
      # parasite taxonomic resolution
      unlist(lapply(webH.list.cp.90, mean)),
      unlist(lapply(webH.list.cp.80, mean)),unlist(lapply(webH.list.cp.70, mean)),
      unlist(lapply(webH.list.cp.60, mean)),unlist(lapply(webH.list.cp.50, mean)),
      unlist(lapply(webH.list.cp.40, mean)),unlist(lapply(webH.list.cp.30, mean)),
      unlist(lapply(webH.list.cp.20, mean)),unlist(lapply(webH.list.cp.10, mean))
    ),
  
  Hse=  # Hse: standardised error of the mean of specialisation
    # Standard errors of the means of the specialisations of the full and resampled communities with reduced effort in 
    # host sampling complteness  
    c(unlist(lapply(lapply(webH.list, unlist),std.error)),
      # Standard errors of the means of the specialisations of the resampled communities with reduced effort in 
      # parasite taxonomic resolution
      unlist(lapply(webH.list.cp.90, std.error)),
      unlist(lapply(webH.list.cp.80, std.error)),unlist(lapply(webH.list.cp.70, std.error)),
      unlist(lapply(webH.list.cp.60, std.error)),unlist(lapply(webH.list.cp.50, std.error)),
      unlist(lapply(webH.list.cp.40, std.error)),unlist(lapply(webH.list.cp.30, std.error)),
      unlist(lapply(webH.list.cp.20, std.error)),unlist(lapply(webH.list.cp.10, std.error))
    ))
```
Finally, we created a figure with the results.

``` r
Mx<-ggplot(results,
           aes(x=sam, y=M, group=res, color=res)) +
  geom_point(size=3) +
  geom_line() +
  geom_errorbar(aes(ymin=M-Mse, ymax=M+Mse), width=.2,
                position=position_dodge(0.05)) +
  ggtitle("(a) Standardised Modularity") +
  xlab("Host Sampling Completeness (%)") +
  ylab(expression(paste("Standardised Modularity", phantom(.)%+-%phantom(.), "S.E.M."))) +
  theme_bw() +
  scale_color_viridis(discrete = TRUE, direction=-1, "Parasite Taxonomic Resolution (%)") +
  scale_fill_viridis(discrete = TRUE, direction=-1)

Nx<-ggplot(results,
           aes(x=sam, y=N, group=res, color=res)) +
  geom_point(size=3) +
  geom_line() +
  geom_errorbar(aes(ymin=N-Nse, ymax=N+Nse), width=.2,
                position=position_dodge(0.05)) +
  ggtitle("(b) Standardised Nestedness") +
  xlab("Host Sampling Completeness (%)") +
  ylab(expression(paste("Standardised Nestedness", phantom(.)%+-%phantom(.), "S.E.M."))) +
  theme_bw() +
  scale_color_viridis(discrete = TRUE, direction=-1, "Parasite Taxonomic Resolution (%)") +
  scale_fill_viridis(discrete = TRUE, direction=-1)

Cx<-ggplot(results,
           aes(x=sam, y=C, group=res, color=res)) +
  geom_point(size=3) +
  geom_line() +
  geom_errorbar(aes(ymin=C-Cse, ymax=C+Cse), width=.2,
                position=position_dodge(0.05)) +
  ggtitle("(c) Connectance") +
  xlab("Host Sampling Completeness (%)") +
  ylab(expression(paste("Connectance", phantom(.)%+-%phantom(.), "S.E.M."))) +
  theme_bw() +
  scale_color_viridis(discrete = TRUE, direction=-1, "Parasite Taxonomic Resolution (%)") +
  scale_fill_viridis(discrete = TRUE, direction=-1)

Hx<-ggplot(results,
           aes(x=sam, y=H, group=res, color=res)) +
  geom_point(size=3) +
  geom_line() +
  geom_errorbar(aes(ymin=H-Hse, ymax=H+Hse), width=.2,
                position=position_dodge(0.05)) +
  ggtitle("(d) Specialisation") +
  xlab("Host Sampling Completeness (%)") +
  ylab(expression(paste("H2", phantom(.)%+-%phantom(.), "S.E.M."))) +
  theme_bw() +
  scale_color_viridis(discrete = TRUE, direction=-1, "Parasite Taxonomic Resolution (%)") +
  scale_fill_viridis(discrete = TRUE, direction=-1)
  
grid.arrange(Mx,Nx,Cx,Hx)
```

## 5. Two-way ANOVA
I analysed community descriptors in four two-way ANOVAs with host sampling completeness and parasite taxonomic resolution as fixed factors.

``` r
a2<-data.frame( # a2: two-way anova data
  sampling=rep(sort(rep(c("wfull","web90","web80","web70","web60","web50","web40","web30","web20","web10"),10),T),10),
  taxonomic=sort(rep(c("wfull","web90","web80","web70","web60","web50","web40","web30","web20","web10"),100),T),
  
  M=c(c(webSM.list,web90SM.list,web80SM.list,
        web70SM.list,web60SM.list,web50SM.list,
        web40SM.list,web30SM.list,web20SM.list,
        web10SM.list),
      unlist(webSM.like.list.cp.90),
      unlist(webSM.like.list.cp.80),unlist(webSM.like.list.cp.70),
      unlist(webSM.like.list.cp.60),unlist(webSM.like.list.cp.50),
      unlist(webSM.like.list.cp.40),unlist(webSM.like.list.cp.30),
      unlist(webSM.like.list.cp.20),unlist(webSM.like.list.cp.10)),
  
  N=c(c(webSN.list,web90SN.list,web80SN.list,
      web70SN.list,web60SN.list,web50SN.list,
      web40SN.list,web30SN.list,web20SN.list,
      web10SN.list),
    unlist(webSN.list.cp.90),
    unlist(webSN.list.cp.80),unlist(webSN.list.cp.70),
    unlist(webSN.list.cp.60),unlist(webSN.list.cp.50),
    unlist(webSN.list.cp.40),unlist(webSN.list.cp.30),
    unlist(webSN.list.cp.20),unlist(webSN.list.cp.10)),
  
  C=c(unlist(webC.list),
      unlist(webC.list.cp.90),
      unlist(webC.list.cp.80),unlist(webC.list.cp.70),
      unlist(webC.list.cp.60),unlist(webC.list.cp.50),
      unlist(webC.list.cp.40),unlist(webC.list.cp.30),
      unlist(webC.list.cp.20),unlist(webC.list.cp.10)),
  
  H=c(unlist(webH.list),
      unlist(webH.list.cp.90),
      unlist(webH.list.cp.80),unlist(webH.list.cp.70),
      unlist(webH.list.cp.60),unlist(webH.list.cp.50),
      unlist(webH.list.cp.40),unlist(webH.list.cp.30),
      unlist(webH.list.cp.20),unlist(webH.list.cp.10))
  )


# Modularity
anova(lm(M ~ sampling + taxonomic + sampling * taxonomic, data=a2))
plot(lm(M ~ sampling + taxonomic + sampling * taxonomic, data=a2))

# Nestedness
anova(lm(N ~ sampling + taxonomic + sampling * taxonomic, data=a2))
plot(lm(N ~ sampling + taxonomic + sampling * taxonomic, data=a2))

# Connectance
anova(lm(C ~ sampling + taxonomic + sampling * taxonomic, data=a2))
plot(lm(C ~ sampling + taxonomic + sampling * taxonomic, data=a2))

# H2
anova(lm(H ~ sampling + taxonomic + sampling * taxonomic, data=a2))
plot(lm(H ~ sampling + taxonomic + sampling * taxonomic, data=a2))

````
