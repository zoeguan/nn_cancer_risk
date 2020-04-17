library(BayesMendel)
packageVersion("BayesMendel")
library(foreach)
library(doMC)
library(parallel)
registerDoMC(cores=16)

script_num=1
args = commandArgs(trailingOnly=T)
if (length(args)>=1) {
  script_num = as.numeric(args[1])
}

load(paste0("../families_misrep/simFamNoise_", script_num, ".RData"))


ids = unique(fam.misrep$FamID)

# set proband's baseline breast cancer status to 0 so BRCAPRO doesn't calculate contralateral BC risk (will eventually exclude the probands affected at baseline)
fam.misrep$AffectedBreast[which(fam.misrep$ID==50)] = 0
fam.misrep$AgeBreast[which(fam.misrep$ID==50)] = fam.misrep$agecur[which(fam.misrep$ID==50)]
fam.misrep$AgeOvary[which(fam.misrep$AffectedOvary==1 & fam.misrep$AgeOvary==1)] = 2

bparams = brcaparams(age.by=5, age.to=94)
bparams$comprisk = 0*bparams$comprisk + 1

run.brcapro = function(family, impute=F, seed=1) {
  set.seed(seed)
  out = tryCatch(brcapro(family[, c("ID", "Gender", "MotherID", "FatherID", "ethnic", "Twins", "AffectedBreast", "AgeBreast", "AffectedOvary", "AgeOvary", "AgeBreastContralateral")], counselee.id=50, params=bparams, print=F, 
                         imputeAges=impute, imputeRelatives=impute, net=T),
                 error = function(e) NA)
  return(tryCatch(unlist(c(out@probs[1, c(1, 2, 3, 4)], out@predictions[1:2,2])), error = function(e) rep(NA, 6)))
}


ind.start = 1
ind.end = length(ids)
results = foreach (i = ind.start:ind.end, .combine=rbind) %dopar% {
  c(ids[i], run.brcapro(fam.misrep[which(fam.misrep$FamID==ids[i]),]))
}

brcapro.results.misrep = data.frame(results)
names(brcapro.results.misrep) = c("FamID", "Prob.Carrier", "Prob.BRCA1", "Prob.BRCA2", "Prob.BRCA12", "brcapro5", "brcapro10")

save(brcapro.results.misrep, file=paste0("../brcapro_results_misrep/brcapro_noise_", script_num, ".RData"))


