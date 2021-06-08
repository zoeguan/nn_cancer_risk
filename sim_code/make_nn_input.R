### standardize pedigrees for NN

library(data.table)

source('standardize_ped.R')


make_nn_input = function(rel.types.df, suffix="med") {
  pro.list = vector("list", length = 50)
  fam.list = vector("list", length = 50)
  
  script_num = 1
  load(paste0("../families/simFamAJ_", script_num, ".RData"))
  load(paste0("../brcapro_results/brcapro_AJ_", script_num, ".RData"))
  
  pro = fam[which(fam$ID==50),]
  pro = pro[which(pro$AffectedBreast==0),]
  pro = left_join(pro, brcapro.results)
  fam.std = standardize_ped(fam, rel.types=rel.types.df)
  
  pro.list[[script_num]] = pro
  fam.list[[script_num]] = fam.std
  
  for (script_num in 2:50) {
    print(script_num)
    load(paste0("../families/simFamAJ_", script_num, ".RData"))
    load(paste0("../brcapro_results/brcapro_AJ_", script_num, ".RData"))
    
    pro = fam[which(fam$ID==50),]
    pro = pro[which(pro$AffectedBreast==0),]
    pro = left_join(pro, brcapro.results)
    fam.std = standardize_ped(fam, rel.types=rel.types.df)
    
    pro.list[[script_num]] = pro
    fam.list[[script_num]] = fam.std
  }
  
  
  pro.all = rbindlist(pro.list)
  nn.input = rbindlist(fam.list)
  pro.all = pro.all[which(pro.all$Death %in% NA),]
  nn.input = nn.input[which(nn.input$FamID %in% pro.all$FamID),]
  
  
  pro.all$BC.5 = pro.all$AffectedBreast.lt==1 & pro.all$AgeBreast.lt <= pro.all$agecur+5
  pro.all$BC.10 = pro.all$AffectedBreast.lt==1 & pro.all$AgeBreast.lt <= pro.all$agecur+10
  
  nn.input = data.frame(nn.input)
  nn.input$FDR = rowSums(nn.input[, which(grepl("AffB_mom|AffB_dad|AffB_sis|AffB_bro|AffB_dau|AffB_son", names(nn.input)))])
  nn.input$SDR = rowSums(nn.input[, which(grepl("AffB", names(nn.input)))]) - nn.input$FDR
  nn.input = left_join(nn.input, pro.all[, c("FamID", "brcapro5", "brcapro10", "BC.5", "BC.10")], by="FamID")
  
  
  nn.input$BC.5 = as.numeric(nn.input$BC.5)
  nn.input$BC.10 = as.numeric(nn.input$BC.10)
  
  #nn.input[, which(grepl("dau3", names(nn.input)) | grepl("son3", names(nn.input)))] = NULL
  
  save(nn.input, file=paste0("../nn_input/nn_input_", suffix, ".RData"))
  
  write.csv(nn.input, file=paste0("../nn_input/nn_input_", suffix, ".csv"))
  
  cnn.input = nn.input[, which(!grepl("FDR|SDR|brcapro|BC.", names(nn.input)))]
  cnn.input[, c("Gender_null", "AffB_null", "AgeB_null", "AffO_null", "AgeO_null", "Age_null", "Miss_null")] = 0
  cnn.input$Miss_null = 1
  cnn.input[, c("FDR", "SDR", "brcapro5", "brcapro10", "BC.5", "BC.10")] = nn.input[, c("FDR", "SDR", "brcapro5", "brcapro10", "BC.5", "BC.10")]
  write.csv(cnn.input, file=paste0("../nn_input/cnn_input_", suffix, ".csv"))
  
}





