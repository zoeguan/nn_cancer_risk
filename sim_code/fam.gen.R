
library(dplyr)

### generate n families with the same structure
# n: number of families to generate
# n.unc.pat: number of paternal uncles 
# n.aunt.pat: number of paternal aunts
# n.unc.mat: number of maternal uncles
# n.aunt.mat: number of maternal aunts
# n.bro: number of brothers
# n.sis: number of sisters
# n.son: number of sons
# n.daught: number of daughters
# n.neph.bro: number of sons for each brother or vector of length n.bro containing number of sons for each brother
# n.nie.bro: number of daughters for each brother or vector of length n.bro containing number of daughters for each brother
# n.neph.sis: number of sons for each sister or vector of length n.sis containing numbers of sons for each sister
# n.nie.sis: number of daughters for each sister or vector of length n.sis containing numbers of daughters for each sister
# ethnic: ethnicity (see BayesMendel formatting)
# age.pro: proband's current age (scalar or vector)
# mean.age.diff, sd.age.diff: ages of non-probands will be generated assuming the age difference between parent and offpsring follows N(mean.age.diff, sd.age.diff), bounded by min.age.diff and max.age.diff
# min.age.diff: minimum age difference between parent and offspring
# max.age.diff: maximum age difference between parent and offspring
# max.age: maximum possible age
# mean.age.cens, sd.age.cens: censoring ages will be generated from N(mean.age.cens, sd.age.cens), bounded by max.age.diff and max.age
# bl.yr: baseline year for risk assessment
fam.gen = function(n=1,
                   n.unc.pat=0, n.aunt.pat=0, n.unc.mat=0, n.aunt.mat=0, 
                   n.bro=0, n.sis=0,
                   n.son=0, n.daught=0,
                   n.neph.bro=0, n.nie.bro=0, n.neph.sis=0, n.nie.sis=0,
                   ethnic="AJ",
                   age.pro=50, 
                   mean.age.diff=27, sd.age.diff=6, min.age.diff=14, max.age.diff=40,
                   max.age=94,
                   mean.age.cens=NULL, sd.age.cens=NULL, bl.yr=NULL) {
  
  if (length(n.neph.bro)==1) {
    n.neph.bro = rep(n.neph.bro, n.bro)
  }
  if (length(n.nie.bro)==1) {
    n.nie.bro = rep(n.nie.bro, n.bro)
  }
  if (length(n.neph.sis)==1) {
    n.neph.sis = rep(n.neph.sis, n.sis)
  }
  if (length(n.nie.sis)==1) {
    n.nie.sis = rep(n.nie.sis, n.sis)
  }
  
  # don't allow more than 10 of a given relative type
  n.unc.pat = pmin(10, n.unc.pat)
  n.aunt.pat = pmin(10, n.aunt.pat)
  n.unc.mat = pmin(10, n.unc.mat)
  n.aunt.mat = pmin(10, n.aunt.mat)
  
  # family size
  fam.size = 8 + n.unc.pat + n.aunt.pat + n.unc.mat + n.aunt.mat + n.bro + n.sis + n.son + n.daught + sum(c(n.neph.bro, n.nie.bro, n.neph.sis, n.nie.sis))
  
  ### set up family dataframe
  fam = data.frame(FamID=rep(1:n, each=fam.size), 
                   ID=rep(c(1:5,
                            rep(10:(9+n.unc.pat), n.unc.pat>0),
                            rep(20:(19+n.aunt.pat), n.aunt.pat>0),
                            6,
                            rep(30:(29+n.unc.mat), n.unc.mat>0),
                            rep(40:(39+n.aunt.mat), n.aunt.mat>0),
                            50,
                            rep(60:(59+n.bro), n.bro>0),
                            rep(70:(69+n.sis), n.sis>0),
                            51,
                            rep(80:(79+n.son), n.son>0),
                            rep(90:(89+n.daught), n.daught>0),
                            rep(100:(99+sum(n.neph.bro)), sum(n.neph.bro)>0),
                            rep(200:(199+sum(n.nie.bro)), sum(n.nie.bro)>0),
                            rep(300:(299+sum(n.neph.sis)), sum(n.neph.sis)>0),
                            rep(400:(399+sum(n.nie.sis)), sum(n.nie.sis)>0)),
                          times=n),
                   relationship=rep(c("Paternal Grandfather", "Paternal Grandmother", "Maternal Grandfather", "Maternal Grandmother", "Father",
                                      rep("Paternal Uncle", n.unc.pat),
                                      rep("Paternal Aunt", n.aunt.pat),
                                      "Mother",
                                      rep("Maternal Uncle", n.unc.mat),
                                      rep("Maternal Aunt", n.aunt.mat),
                                      "Self",
                                      rep("Brother", n.bro),
                                      rep("Sister", n.sis),
                                      "Spouse",
                                      rep("Son", n.son),
                                      rep("Daughter", n.daught),
                                      rep("Nephew", sum(n.neph.bro)),
                                      rep("Niece", sum(n.nie.bro)),
                                      rep("Nephew", sum(n.neph.sis)),
                                      rep("Niece", sum(n.nie.sis))), times=n))
  ### Gender
  fam$Gender = 0
  fam$Gender[which(fam$relationship %in% c("Paternal Grandfather", "Maternal Grandfather", "Father", "Paternal Uncle", "Maternal Uncle", "Brother", "Spouse", "Son", "Nephew"))] = 1
  ### FatherID
  fam$FatherID = 0
  fam$FatherID[which(fam$relationship %in% c("Father", "Paternal Uncle", "Paternal Aunt"))] = 1
  fam$FatherID[which(fam$relationship %in% c("Mother", "Maternal Uncle", "Maternal Aunt"))] = 3
  fam$FatherID[which(fam$relationship %in% c("Self", "Brother", "Sister"))] = 5
  fam$FatherID[which(fam$relationship %in% c("Son", "Daughter"))] = 51
  fam$FatherID[which(fam$relationship %in% c("Nephew"))] = rep(c(rep(rep(60:(59+n.bro), n.bro>0), n.neph.bro), 0*rep(rep(70:(69+n.sis), n.sis>0), n.neph.sis)), n)
  fam$FatherID[which(fam$relationship %in% c("Niece"))] = rep(c(rep(rep(60:(59+n.bro), n.bro>0), n.nie.bro), 0*rep(rep(70:(69+n.sis), n.sis>0), n.nie.sis)), n)
  ### MotherID
  fam$MotherID = 0
  fam$MotherID[which(fam$relationship %in% c("Father", "Paternal Uncle", "Paternal Aunt"))] = 2
  fam$MotherID[which(fam$relationship %in% c("Mother", "Maternal Uncle", "Maternal Aunt"))] = 4
  fam$MotherID[which(fam$relationship %in% c("Self", "Brother", "Sister"))] = 6
  fam$MotherID[which(fam$relationship %in% c("Son", "Daughter"))] = 50
  fam$MotherID[which(fam$relationship %in% c("Nephew"))] = rep(c(0*rep(rep(60:(59+n.bro), n.bro>0), n.neph.bro), rep(rep(70:(69+n.sis), n.sis>0), n.neph.sis)), n)
  fam$MotherID[which(fam$relationship %in% c("Niece"))] = rep(c(0*rep(rep(60:(59+n.bro), n.bro>0), n.nie.bro), rep(rep(70:(69+n.sis), n.sis>0), n.nie.sis)), n)
  
  ### ethnicity
  fam$ethnic = ethnic
  
  ### twins
  fam$Twins = 0
  
  
  
  
  
  ### current age 
  fam$agecur = NA
  fam$birth.yr = NA
  fam$bl.yr = NA
  ## proband
  fam$agecur[which(fam$relationship=="Self")] = age.pro
  # set bl.yr and birth.yr
  fam$bl.yr = as.numeric(format(Sys.Date(), "%Y")) - 10
  if (!is.null(bl.yr)) {
    fam$bl.yr = bl.yr
  }
  fam$birth.yr[which(fam$relationship=="Self")] = fam$bl.yr[which(fam$relationship=="Self")] - fam$agecur[which(fam$relationship=="Self")]
  # age difference with mother
  fam$age.diff.m = pmax(pmin(round(rnorm(nrow(fam), mean.age.diff, sd.age.diff)), max.age.diff), min.age.diff)
  fam$age.diff.f = pmax(pmin(round(rnorm(nrow(fam), mean.age.diff, sd.age.diff)), max.age.diff), min.age.diff)
  ## mother
  fam$birth.yr[which(fam$relationship %in% c("Mother"))] = fam$birth.yr[which(fam$relationship%in% c("Self"))] - fam$age.diff.m[which(fam$relationship%in% c("Self"))]
  # maternal grandparents
  fam$birth.yr[which(fam$relationship %in% c("Maternal Grandmother"))] = fam$birth.yr[which(fam$relationship %in% c("Mother"))] - fam$age.diff.m[which(fam$relationship%in% c("Mother"))]
  fam$birth.yr[which(fam$relationship %in% c("Maternal Grandfather"))] = fam$birth.yr[which(fam$relationship %in% c("Mother"))] - fam$age.diff.f[which(fam$relationship%in% c("Mother"))]
  # mother's siblings
  fam$birth.yr[which(fam$relationship %in% c("Maternal Aunt", "Maternal Uncle"))] = pmax(rep(fam$birth.yr[which(fam$relationship  %in% c("Maternal Grandmother"))], each=n.unc.mat+n.aunt.mat) + fam$age.diff.m[which(fam$relationship %in% c("Maternal Aunt", "Maternal Uncle"))], 
  rep(fam$birth.yr[which(fam$relationship  %in% c("Maternal Grandfather"))], each=n.unc.mat+n.aunt.mat) + fam$age.diff.f[which(fam$relationship %in% c("Maternal Aunt", "Maternal Uncle"))] )
  # father
  fam$birth.yr[which(fam$relationship %in% c("Father"))] = fam$birth.yr[which(fam$relationship=="Self")] - fam$age.diff.f[which(fam$relationship%in% c("Self"))]
  # paternal grandparents
  fam$birth.yr[which(fam$relationship %in% c("Paternal Grandmother"))] = fam$birth.yr[which(fam$relationship %in% c("Father"))] - fam$age.diff.m[which(fam$relationship%in% c("Father"))]
  fam$birth.yr[which(fam$relationship %in% c("Paternal Grandfather"))] = fam$birth.yr[which(fam$relationship %in% c("Father"))] - fam$age.diff.f[which(fam$relationship%in% c("Father"))]
  # father's siblings
  fam$birth.yr[which(fam$relationship %in% c("Paternal Aunt", "Paternal Uncle"))] = pmax(rep(fam$birth.yr[which(fam$relationship  %in% c("Paternal Grandmother"))], each=n.unc.pat+n.aunt.pat) + fam$age.diff.m[which(fam$relationship %in% c("Paternal Aunt", "Paternal Uncle"))], 
rep(fam$birth.yr[which(fam$relationship  %in% c("Paternal Grandfather"))], each=n.unc.pat+n.aunt.pat) + fam$age.diff.f[which(fam$relationship %in% c("Paternal Aunt", "Paternal Uncle"))])
  # siblings
  fam$birth.yr[which(fam$relationship %in% c("Sister", "Brother"))] = pmax(rep(fam$birth.yr[which(fam$relationship  %in% c("Mother"))], each=n.sis+n.bro) + fam$age.diff.m[which(fam$relationship %in% c("Sister", "Brother"))], 
                                                                           rep(fam$birth.yr[which(fam$relationship  %in% c("Father"))], each=n.sis+n.bro) + fam$age.diff.f[which(fam$relationship %in% c("Sister", "Brother"))])
  # spouse
  fam$birth.yr[which(fam$relationship %in% c("Spouse"))] = fam$birth.yr[which(fam$relationship %in% c("Self"))] + pmin(5, pmax(-5, round(rnorm(length(which(fam$relationship=="Spouse"))))))
  # children, nephews, nieces
  fam$unique.id = paste0(fam$FamID, "_", fam$ID)
  fam$unique.mid = paste0(fam$FamID, "_", fam$MotherID)
  fam$unique.mid[which(fam$MotherID==0)] = NA
  fam$unique.fid = paste0(fam$FamID, "_", fam$FatherID)
  fam$unique.fid[which(fam$FatherID==0)] = NA
  fam = left_join(fam, fam[, c("unique.id", "birth.yr")], by=c("unique.fid" = "unique.id"), suffix=c("", ".f"))
  fam = left_join(fam, fam[, c("unique.id", "birth.yr")], by=c("unique.mid" = "unique.id"), suffix=c("", ".m"))
  fam$birth.yr[which(fam$FatherID %in% 60:69)] = fam$birth.yr.f[which(fam$FatherID %in% 60:69)] + fam$age.diff.f[which(fam$FatherID %in% 60:69)]
  fam$birth.yr[which(fam$MotherID %in% c(50, 70:79))] = pmax(fam$birth.yr.m[which(fam$MotherID %in% c(50, 70:79))] + fam$age.diff.m[which(fam$MotherID %in% c(50, 70:79))],
                                                             fam$birth.yr.f[which(fam$MotherID %in% c(50, 70:79))] + fam$age.diff.f[which(fam$MotherID %in% c(50, 70:79))])
  
  fam$age.diff.f = NULL
    
  # birth.yr should be <= bl.yr
  fam$birth.yr[which(fam$birth.yr > fam$bl.yr)] = fam$bl.yr[which(fam$birth.yr > fam$bl.yr)]
  fam$agecur = fam$bl.yr - fam$birth.yr
  
  # age 1 treated as unknown in BayesMendel, so change ages that are 1 to 2
  ind.1 = which(fam$agecur==1)
  fam$agecur[ind.1] = fam$agecur[ind.1]+1
  fam$birth.yr[ind.1] = fam$birth.yr[ind.1]-1
  
  fam$Death = NA
  fam$AgeDeath = NA
  
  
  # deal with ages over max.age
  if (is.null(mean.age.cens)) {
    fam$agecur[which(fam$agecur>94)] = 94
  } else {
    fam$AgeDeath = pmax(max.age.diff, pmin(94, round(rnorm(nrow(fam), mean.age.cens, sd.age.cens))))
    fam$Death[which(fam$AgeDeath <= fam$agecur)] = 1
    fam$agecur[which(fam$Death==1)] = fam$AgeDeath[which(fam$Death==1)]
  }
  
  fam$relationship = as.character(fam$relationship)
  
  return(fam)

}
