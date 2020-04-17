
library(BayesMendel)

cancer.age.gen = function(m, gender.cancer="fFX", genotype="B00", penetrance) {
  penet = penetrance[[gender.cancer]]
  w = append(penet[, genotype], 1-sum(penet[, genotype]))
  cancer.ages = rmultinom(m, 1, prob=w)
  cancer.ages.vec = rep(NA, m)
  for (i in 1:length(w)) {
    cancer.ages.vec[which(cancer.ages[i, ]==1)] = i
  }
  return(cancer.ages.vec)
  #return(apply(cancer.ages, 2, function(c) which(c==1)))
}



### generate breast cancer and ovarian cancer phenotypes
pheno.gen = function(fam, penetrance = penet.brca.net, max.age=94) {
  
  #########################################
  ### breast cancer ages
  #########################################
  
  fam$AgeBreast.lt = NA
  
  # males
  ind12.male = which(fam$BRCA1==1 & fam$BRCA2==1 & fam$Gender==1)
  fam$AgeBreast.lt[ind12.male] = cancer.age.gen(length(ind12.male), "fMX", "B11", penetrance)
  
  ind1.male = which(fam$BRCA1==1 & fam$BRCA2==0 & fam$Gender==1)
  fam$AgeBreast.lt[ind1.male] = cancer.age.gen(length(ind1.male), "fMX", "B10", penetrance)
  
  ind2.male = which(fam$BRCA1==0 & fam$BRCA2==1 & fam$Gender==1)
  fam$AgeBreast.lt[ind2.male] = cancer.age.gen(length(ind2.male), "fMX", "B01", penetrance)
  
  ind0.male = which(fam$BRCA1==0 & fam$BRCA2==0 & fam$Gender==1)
  fam$AgeBreast.lt[ind0.male] = cancer.age.gen(length(ind0.male), "fMX", "B00", penetrance)
  
  # females
  ind12.female = which(fam$BRCA1==1 & fam$BRCA2==1 & fam$Gender==0)
  fam$AgeBreast.lt[ind12.female] = cancer.age.gen(length(ind12.female), "fFX", "B11", penetrance)
  
  ind1.female = which(fam$BRCA1==1 & fam$BRCA2==0 & fam$Gender==0)
  fam$AgeBreast.lt[ind1.female] = cancer.age.gen(length(ind1.female), "fFX", "B10", penetrance)
  
  ind2.female = which(fam$BRCA1==0 & fam$BRCA2==1 & fam$Gender==0)
  fam$AgeBreast.lt[ind2.female] = cancer.age.gen(length(ind2.female), "fFX", "B01", penetrance)
  
  ind0.female = which(fam$BRCA1==0 & fam$BRCA2==0 & fam$Gender==0)
  fam$AgeBreast.lt[ind0.female] = cancer.age.gen(length(ind0.female), "fFX", "B00", penetrance)
  
  #########################################
  ### ovarian cancer ages
  #########################################
  
  fam$AgeOvary.lt = NA
  fam$AgeOvary.lt[ind12.female] = cancer.age.gen(length(ind12.female), "fFY", "B11", penetrance)
  fam$AgeOvary.lt[ind1.female] = cancer.age.gen(length(ind1.female), "fFY", "B10", penetrance)
  fam$AgeOvary.lt[ind2.female] = cancer.age.gen(length(ind2.female), "fFY", "B01", penetrance)
  fam$AgeOvary.lt[ind0.female] = cancer.age.gen(length(ind0.female), "fFY", "B00", penetrance)
  
  #########################################
  ### AffectedBreast, AffectedOvary
  #########################################
  
  fam$AffectedBreast.lt = 0
  fam$AffectedOvary.lt = 0
  
  ind.BC = which(fam$AgeBreast.lt<=length(penetrance$fFX[, 1]))
  ind.OC = which(fam$AgeOvary.lt<=length(penetrance$fFX[, 1]))
  fam$AffectedBreast.lt[ind.BC] = 1
  fam$AffectedOvary.lt[ind.OC] = 1
  
  
  fam$AffectedBreast = as.numeric(fam$AgeBreast.lt <= fam$agecur)
  fam$AffectedOvary = as.numeric(fam$AgeOvary.lt <= fam$agecur)
  fam$AgeBreast = fam$agecur
  fam$AgeOvary = fam$agecur
  
  fam$AgeBreastContralateral = 0
  fam$AgeBreast[which(fam$AffectedBreast==1)] = fam$AgeBreast.lt[which(fam$AffectedBreast==1)]
  fam$AgeOvary[which(fam$AffectedOvary==1)] = fam$AgeOvary.lt[which(fam$AffectedOvary==1)]
  
  return(fam)
  
}