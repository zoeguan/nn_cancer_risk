### simulate AJ families based on structures sampled from CGN

library(BayesMendel)
library(dplyr)
source("fam.gen.R")
source("geno.gen.R")
source("pheno.gen.R")
load("cgn.rel.counts.RData")

script_num = 1
# get script_num from command line
args = commandArgs(trailingOnly=T)
if (length(args)>=1) {
  script_num = as.numeric(args[1])
}

# number of families to generate
N = 20000

# allow up to 5 relatives of a given type
max.count = 5
cgn.rel.counts[, 2:11] = apply(cgn.rel.counts[, 2:11], 2, function(c) pmin(c, max.count))
cgn.rel.counts[, c("HALF_SIS", "HALF_BRO")] = NULL
cgn.rel.counts = na.omit(cgn.rel.counts)

# allow proband ages between 20 and 84
cgn.ages = cgn.ages[which(cgn.ages>=20 & cgn.ages<=84)]

# sample CGN families
set.seed(script_num)
ind.sample = sample(1:nrow(cgn.rel.counts), N, replace=T)
ind.sample.ages = sample(1:length(cgn.ages), N, replace=T)
probands = data.frame(FamID=1:N, agecur=cgn.ages[ind.sample.ages])
count.name = c("FULL_BRO", "FULL_SIS", "TOT_SONS", "TOT_DAUGHTERS", 
               "TOT_PAT_UNCLES", "TOT_PAT_AUNTS", "TOT_MAT_UNCLES", "TOT_MAT_AUNTS")
probands = cbind(probands, cgn.rel.counts[ind.sample, count.name])

# simulate N families with the maximal structure allowed
# use default BRCAPRO AJ prevalences (v2.1-6)
# mean.age.diff=27, sd.age.diff=6 based on distribution of age differences between probands and their mothers in CGN
# mean.age.cens=80, sd.age.cens=15 based on U.S. life expectancy statistics
set.seed(script_num)
fam = fam.gen(n=N,
              n.unc.pat=max.count, n.aunt.pat=max.count, n.unc.mat=max.count, n.aunt.mat=max.count, 
              n.bro=max.count, n.sis=max.count,
              n.son=max.count, n.daught=max.count,
              n.neph.bro=0, n.nie.bro=0, n.neph.sis=0, n.nie.sis=0,
              ethnic="AJ",
              age.pro=cgn.ages[ind.sample.ages], 
              mean.age.diff=27, sd.age.diff=6, min.age.diff=14, max.age.diff=40,
              max.age=94,
              mean.age.cens=80, sd.age.cens=15, bl.yr=2001)
fam = geno.gen(fam, pbrca1=0.01366243, pbrca2=0.01168798)
fam = pheno.gen(fam, penetrance = penet.brca.net)


# prune the simulated families so that they have the same structure as the sampled CGN families
rel.df = data.frame(firstID=c(60, 70, 80, 90, 10, 20, 30, 40), 
           relationship=c("Brother", "Sister", "Son", "Daughter", 
                          "Paternal Uncle", "Paternal Aunt", "Maternal Uncle", "Maternal Aunt"), 
           count.name=count.name)


fam = left_join(fam, probands[, c("FamID", count.name)], by="FamID")

fam$remove = 0
for (i in 1:length(count.name)) {
  fam$remove[which(fam$relationship %in% rel.df$relationship[i] & fam$ID >= rel.df$firstID[i] + fam[, count.name[i]])] = 1
}


fam = fam[which(fam$remove==0 & fam$agecur>0),]
fam[ , c("remove", count.name, "birth.yr.m", "birth.yr.f")] = NULL

fam$FamID = fam$FamID + (script_num-1)*N
probands$FamID = probands$FamID + (script_num-1)*N

save(probands, fam, file=paste0("../families/simFamAJ_", script_num, ".RData"))
