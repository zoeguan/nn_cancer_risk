### standardize perturbed pedigrees for NN

source('standardize_ped.R')

pro.list = vector("list", length = 50)
fam.list = vector("list", length = 50)

script_num = 1
load(paste0("../families_misrep/simFamNoise_", script_num, ".RData"))
load(paste0("../brcapro_results_misrep/brcapro_noise_", script_num, ".RData"))

pro = fam.misrep[which(fam.misrep$ID==50),]
pro = pro[which(pro$AffectedBreast==0),]
pro = left_join(pro, brcapro.results.misrep)
fam.std = standardize_ped(fam.misrep)

pro.list[[script_num]] = pro
fam.list[[script_num]] = fam.std

for (script_num in 2:25) {
  print(script_num)
  load(paste0("../families_misrep/simFamNoise_", script_num, ".RData"))
  load(paste0("../brcapro_results_misrep/brcapro_noise_", script_num, ".RData"))
  
  pro = fam.misrep[which(fam.misrep$ID==50),]
  pro = pro[which(pro$AffectedBreast==0),]
  pro = left_join(pro, brcapro.results.misrep)
  fam.std = standardize_ped(fam.misrep)
  
  pro.list[[script_num]] = pro
  fam.list[[script_num]] = fam.std
}

for (script_num in 26:50) {
  print(script_num)
  load(paste0("../families_misrep/simFamNoise_", script_num, ".RData"))
  
  pro = fam.misrep[which(fam.misrep$ID==50),]
  pro = pro[which(pro$AffectedBreast==0),]
  pro[, setdiff(names(brcapro.results.misrep), "FamID")] = 0
  fam.std = standardize_ped(fam.misrep)
  
  pro.list[[script_num]] = pro
  fam.list[[script_num]] = fam.std
}


library(data.table)
pro.all = rbindlist(pro.list)
nn.input = rbindlist(fam.list)
pro.all = pro.all[which(pro.all$Death %in% NA),]
nn.input = nn.input[which(nn.input$FamID %in% pro.all$FamID),]


pro.all$BC.5 = pro.all$AffectedBreast.lt==1 & pro.all$AgeBreast.lt <= pro.all$agecur+5
pro.all$BC.10 = pro.all$AffectedBreast.lt==1 & pro.all$AgeBreast.lt <= pro.all$agecur+10

nn.input.noise = data.frame(nn.input)

load("../nn_input/nn_input.RData")

nn.input.noise[, c("FDR", "SDR")] = nn.input[, c("FDR", "SDR")]

nn.input.noise = left_join(nn.input.noise, pro.all[, c("FamID", "brcapro5", "brcapro10", "BC.5", "BC.10")], by="FamID")


nn.input.noise$BC.5 = as.numeric(nn.input.noise$BC.5)
nn.input.noise$BC.10 = as.numeric(nn.input.noise$BC.10)

nn.input.noise[, which(grepl("dau3", names(nn.input.noise)) | grepl("son3", names(nn.input.noise)))] = NULL

save(nn.input.noise, file="../nn_input/nn_input_noise.RData")


write.csv(nn.input.noise, file="../nn_input/nn_input_noise.csv")



cnn.input.noise = nn.input.noise[, which(!grepl("FDR|SDR|brcapro|BC.", names(nn.input.noise)))]
cnn.input.noise[, c("Gender_null", "AffB_null", "AgeB_null", "AffO_null", "AgeO_null", "Age_null", "Miss_null")] = 0
cnn.input.noise$Miss_null = 1
cnn.input.noise[, c("FDR", "SDR", "brcapro5", "brcapro10", "BC.5", "BC.10")] = nn.input.noise[, c("FDR", "SDR", "brcapro5", "brcapro10", "BC.5", "BC.10")]

write.csv(cnn.input.noise, file="../nn_input/cnn_input_noise.csv")
