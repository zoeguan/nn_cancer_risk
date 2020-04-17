script_num=1
args = commandArgs(trailingOnly=T)
if (length(args)>=1) {
  script_num = as.numeric(args[1])
}

load(paste0("../families/simFamAJ_", script_num, ".RData"))

# relative types
first.deg.f = c("Mother", "Sister", "Daughter")
second.deg.f = c("Paternal Grandmother", "Maternal Grandmother", 
                 "Paternal Aunt", "Maternal Aunt", "Niece")
first.deg.m = c("Father", "Brother", "Son")
second.deg.m = c("Paternal Grandfather", "Maternal Grandfather", 
                 "Paternal Uncle", "Maternal Uncle", "Nephew")

### function to generate false negatives and false positives for cancer diagnoses
# ind.c: indices of true positives
# ind.nc: indices of true negatives
# cancer.name: cancer site ("Breast" or "Ovary")
# fnr: false negative rate
# fpr: false positive rate
# true.ages: vector of diagnosis ages of true positives
perturb.canc = function(fam, ind.c, ind.nc, cancer.name="Breast", fnr, fpr, true.ages) {
  fam[ind.c, paste0("Affected", cancer.name)] = rbinom(length(ind.c), 1, 1-fnr)
  fam[ind.nc, paste0("Affected", cancer.name)] = rbinom(length(ind.nc), 1, fpr)
  # for false positives, sample diagnosis ages from diagnosis ages of true positives
  ind.update.age = intersect(ind.nc, which(fam[ind.nc, paste0("Affected", cancer.name)] %in% 1))
  fam[ind.update.age, paste0("Age", cancer.name)] = pmin(sample(true.ages, length(ind.update.age), replace=T), fam$agecur[ind.update.age] )
  return(fam)
}


### function to perturb ages
# ind: indices of individuals whose ages are to be perturbed
# age.name: name of age variable to be perturbed ("AgeBreast", "AgeOvary", or "agecur")
# mean.diff: mean difference between true and reported ages
# sd.diff: standard deviation of difference between true and reported ages
## default mean.diff is from Braun et al. 2017
perturb.age = function(fam, ind, age.name="AgeBreast", mean.diff=4, sd.diff=3, min.age=2, max.age=94) {
  if (!age.name %in% "agecur") {
    max.age = fam$agecur[ind]
  }
  age.diff = round(abs(rnorm(length(ind), mean.diff, sd.diff)))
  
  fam[ind, age.name] = pmin(max.age, 
                            pmax(min.age, 
                                 fam[ind, age.name] + sample(c(-1, 1), length(ind), replace=T)*age.diff ))
  return(fam)
}


### function for perturbing families (simulate misreported diagnoses and ages/diagnosis ages)
# BC.error.rates: dataframe where row 1 contains the false positive and false negative rates for reported breast cancer diagnoses among first-degree relatives and row 2 the rates for reported breast cancer diagnoses among second-degree relatives
# OC.error.rates: dataframe where row 1 contains the false positive and false negative rates for reported ovarian cancer diagnoses among first-degree relatives and row 2 the rates for reported ovarian cancer diagnoses among second-degree relatives
## default error rates are from Ziogas and Anton-Culver 2003 
add.noise = function(fam, 
                     BC.error.rates=data.frame(degree=c(1,2), FNR=c(5, 18)/100, FPR=c(3, 3)/100), 
                     OC.error.rates=data.frame(degree=c(1,2), FNR=c(17, 56)/100, FPR=c(1, 2)/100),
                     seed=1) {
  set.seed(seed)
  
  ##### misreported diagnoses
  
  ### 1st-degree BC
  ind.BC.1 = which(fam$relationship %in% c(first.deg.f, first.deg.m) & 
                     fam$AffectedBreast==1)
  ind.noBC.1 = setdiff(which(fam$relationship %in% c(first.deg.f, first.deg.m)), ind.BC.1)
  BC.ages.1 = fam$AgeBreast[ind.BC.1]
  fam = perturb.canc(fam, ind.BC.1, ind.noBC.1, "Breast", BC.error.rates$FNR[1], BC.error.rates$FPR[1], BC.ages.1)
  
  ### 2nd-degree BC
  ind.BC.2 = which(fam$relationship %in% c(second.deg.f, second.deg.m) & 
                     fam$AffectedBreast==1)
  ind.noBC.2 = setdiff(which(fam$relationship %in% c(second.deg.f, second.deg.m)), ind.BC.2)
  BC.ages.2 = fam$AgeBreast[ind.BC.2]
  fam = perturb.canc(fam, ind.BC.2, ind.noBC.2, "Breast", BC.error.rates$FNR[2], BC.error.rates$FPR[2], BC.ages.2)
  
  ### 1st-degree OC
  ind.OC.1 = which(fam$relationship %in% first.deg.f & fam$AffectedOvary==1)
  ind.noOC.1 = setdiff(which(fam$relationship %in% first.deg.f), ind.OC.1)
  OC.ages.1 = fam$AgeOvary[ind.OC.1]
  fam = perturb.canc(fam, ind.OC.1, ind.noOC.1, "Ovary", OC.error.rates$FNR[1], OC.error.rates$FPR[1], OC.ages.1)
  
  ### 2nd-degree OC
  ind.OC.2 = which(fam$relationship %in% second.deg.f & fam$AffectedOvary==1)
  ind.noOC.2 = setdiff(which(fam$relationship %in% second.deg.f), ind.OC.2)
  OC.ages.2 = fam$AgeOvary[ind.OC.2]
  fam = perturb.canc(fam, ind.OC.2, ind.noOC.2, "Ovary", OC.error.rates$FNR[2], OC.error.rates$FPR[2], OC.ages.2)
  
  ##### misreported ages
  
  ### current age 
  # rates of misreporting are based on accuracy estimates from the CGN
  fam$agecur.orig = fam$agecur
  ind.age.1 = sample(which(fam$relationship %in% c(first.deg.f, first.deg.m)), round(0.02*length(which(fam$relationship %in% c(first.deg.f, first.deg.m)))), replace=F)
  fam = perturb.age(fam, ind.age.1, "agecur")
  ind.age.2 = sample(which(fam$relationship %in% c(second.deg.f, second.deg.m)), round(0.2*length(which(fam$relationship %in% c(second.deg.f, second.deg.m)))), replace=F)
  fam = perturb.age(fam, ind.age.2, "agecur")
  fam$agecur[which(fam$agecur<fam$AgeBreast)] = fam$AgeBreast[which(fam$agecur<fam$AgeBreast)]
  fam$agecur[which(fam$agecur<fam$AgeOvary)] = fam$AgeOvary[which(fam$agecur<fam$AgeOvary)]
  
  ### BC age 
  # misreporting rate from Braun et al. 2017
  ind.BC.age.1 = sample(which(fam$AffectedBreast==1 & fam$relationship %in% c(first.deg.f, first.deg.m)), round(0.03*length(which(fam$AffectedBreast==1 & fam$relationship %in% c(first.deg.f, first.deg.m)))), replace=F)
  fam = perturb.age(fam, ind.BC.age.1, "AgeBreast")
  
  ind.BC.age.2 = sample(which(fam$AffectedBreast==1 & fam$relationship %in% c(second.deg.f, second.deg.m)), round(0.03*length(which(fam$AffectedBreast==1 & fam$relationship %in% c(second.deg.f, second.deg.m)))), replace=F)
  fam = perturb.age(fam, ind.BC.age.2, "AgeBreast")
  
  ### OC age
  # misreporting rate from Braun et al. 2017
  ind.OC.age.1 = sample(which(fam$AffectedOvary==1 & fam$relationship %in% c(first.deg.f)), round(0.04*length(which(fam$AffectedOvary==1 & fam$relationship %in% c(first.deg.f)))), replace=F)
  fam = perturb.age(fam, ind.OC.age.1, "AgeOvary")
  
  ind.OC.age.2 = sample(which(fam$AffectedOvary==1 & fam$relationship %in% c(second.deg.f)), round(0.04*length(which(fam$AffectedOvary==1 & fam$relationship %in% c(second.deg.f)))), replace=F)
  fam = perturb.age(fam, ind.OC.age.2, "AgeOvary")
  
  ### for those reported as unaffected, set AgeBreast/AgeOvary to current age
  fam$AgeBreast[which(fam$AffectedBreast==0)] = fam$agecur[which(fam$AffectedBreast==0)]
  fam$AgeOvary[which(fam$AffectedOvary==0 & fam$Gender==0)] = fam$agecur[which(fam$AffectedOvary==0 & fam$Gender==0)]
  
  
  return(fam)
}

fam.misrep = add.noise(fam)

save(fam.misrep, file=paste0("../families_misrep/simFamNoise_", script_num, ".RData"))


# ## check proportions
# 
# ind = which(fam$relationship %in% c(first.deg.f, first.deg.m) & fam$AffectedBreast==1)
# mean(fam$AffectedBreast[ind] > fam.misrep$AffectedBreast[ind])
# BC.error.rates$FNR[1]
# 
# ind = which(fam$relationship %in% c(first.deg.f, first.deg.m) & fam$AffectedBreast==0)
# mean(fam$AffectedBreast[ind] < fam.misrep$AffectedBreast[ind])
# BC.error.rates$FPR[1]
# 
# ind = which(fam$relationship %in% c(second.deg.f, second.deg.m) & fam$AffectedBreast==1)
# mean(fam$AffectedBreast[ind] > fam.misrep$AffectedBreast[ind])
# BC.error.rates$FNR[2]
# 
# ind = which(fam$relationship %in% c(second.deg.f, second.deg.m) & fam$AffectedBreast==0)
# mean(fam$AffectedBreast[ind] < fam.misrep$AffectedBreast[ind])
# BC.error.rates$FPR[2]
# 
# ind = which(fam$relationship %in% c(first.deg.f, first.deg.m) & fam$AffectedOvary==1)
# mean(fam$AffectedOvary[ind] > fam.misrep$AffectedOvary[ind])
# OC.error.rates$FNR[1]
# 
# ind = which(fam$relationship %in% c(first.deg.f, first.deg.m) & fam$AffectedOvary==0)
# mean(fam$AffectedOvary[ind] < fam.misrep$AffectedOvary[ind])
# OC.error.rates$FPR[1]
# 
# ind = which(fam$relationship %in% c(second.deg.f, second.deg.m) & fam$AffectedOvary==1)
# mean(fam$AffectedOvary[ind] > fam.misrep$AffectedOvary[ind])
# OC.error.rates$FNR[2]
# 
# ind = which(fam$relationship %in% c(second.deg.f, second.deg.m) & fam$AffectedOvary==0)
# mean(fam$AffectedOvary[ind] < fam.misrep$AffectedOvary[ind])
# OC.error.rates$FPR[2]
