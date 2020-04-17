library(BayesMendel)
packageVersion("BayesMendel")
library(foreach)
library(doMC)
library(parallel)

registerDoMC(cores=16)

script_num=1

# get script_num from command line
args = commandArgs(trailingOnly=T)
if (length(args)>=1) {
  script_num = as.numeric(args[1])
}

load(paste0("../families/simFamAJ_", script_num, ".RData"))

ids = unique(fam$FamID)

# set proband's baseline breast cancer status to 0 so BRCAPRO doesn't calculate contralateral BC risk (will eventually exclude the probands affected at baseline)
fam$AffectedBreast[which(fam$ID==50)] = 0
fam$AgeBreast[which(fam$ID==50)] = fam$agecur[which(fam$ID==50)]
# change diagnosis age 1 to 2 so BRCAPRO doesn't treat it as unknown
fam$AgeOvary[which(fam$AffectedOvary==1 & fam$AgeOvary==1)] = 2

bparams = brcaparams(age.by=5, age.to=94)
bparams$comprisk = 0*bparams$comprisk + 1

# get carrier probabilities, 5-year risk and 10-year risk
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
  print(i)
  c(ids[i], run.brcapro(fam[which(fam$FamID==ids[i]),]))
}

brcapro.results = data.frame(results)
names(brcapro.results) = c("FamID", "Prob.Carrier", "Prob.BRCA1", "Prob.BRCA2", "Prob.BRCA12", "brcapro5", "brcapro10")


save(brcapro.results, file=paste0("../brcapro_results/brcapro_AJ_", script_num, ".RData"))
